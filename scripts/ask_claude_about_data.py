"""Ask Claude questions about LARGE local datasets without uploading them.

The data never leaves your computer. Claude gets a few small, read-only tools
(list files, read a matrix header, summarise a table column). It calls them,
your computer runs them, and only the short text results go back to Claude.

Setup:   pip install anthropic pandas
         set ANTHROPIC_API_KEY (or run `ant auth login`)
Run:     DATA_ROOT=/path/to/data python scripts/ask_claude_about_data.py "How many cells are in each library?"
"""
import gzip
import os
import sys
from pathlib import Path

import anthropic
import pandas as pd
from anthropic import beta_tool

DATA_ROOT = Path(os.environ.get("DATA_ROOT", ".")).resolve()   # Claude can only look inside this folder
MAX_CHARS = 4000                                               # cap on what any one tool returns


def _inside_root(path: str) -> Path:
    p = (DATA_ROOT / path).resolve()
    if p != DATA_ROOT and DATA_ROOT not in p.parents:
        raise ValueError("path is outside the data folder")
    return p


def _short(text: str) -> str:
    return text if len(text) <= MAX_CHARS else text[:MAX_CHARS] + "\n...[truncated]"


def list_files(folder: str = ".") -> str:
    """List files and their sizes (MB) in a folder inside the data folder.

    Args:
        folder: Folder relative to the data folder.
    """
    rows = [f"{f.name}\t{f.stat().st_size / 1e6:.1f} MB" for f in sorted(_inside_root(folder).iterdir()) if f.is_file()]
    return _short("\n".join(rows) or "no files")


def mtx_dimensions(path: str) -> str:
    """Read only the header of a 10x matrix.mtx or matrix.mtx.gz file.

    Args:
        path: File relative to the data folder.
    """
    p = _inside_root(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt") as handle:
        for line in handle:
            if not line.startswith("%"):
                genes, cells, nonzero = line.split()
                return f"genes={genes} cells={cells} nonzero_entries={nonzero}"
    return "no header found"


def table_overview(path: str) -> str:
    """Show the columns and first 5 rows of a CSV or TSV table (reads only 5 rows).

    Args:
        path: File relative to the data folder.
    """
    p = _inside_root(path)
    sep = "," if ".csv" in p.name else "\t"
    head = pd.read_csv(p, sep=sep, nrows=5)
    return _short(f"{len(head.columns)} columns: {list(head.columns)}\n\n{head.to_string()}")


def column_summary(path: str, column: str) -> str:
    """Summarise one column of a CSV or TSV table (reads only that column).

    Args:
        path: File relative to the data folder.
        column: Column name.
    """
    p = _inside_root(path)
    sep = "," if ".csv" in p.name else "\t"
    values = pd.read_csv(p, sep=sep, usecols=[column])[column]
    if pd.api.types.is_numeric_dtype(values):
        return _short(values.describe().to_string())
    return _short(f"missing={values.isna().sum()}  distinct={values.nunique()}\n{values.value_counts().head(15).to_string()}")


TOOLS = [beta_tool(f) for f in (list_files, mtx_dimensions, table_overview, column_summary)]

SYSTEM = (
    "You help a bioinformatician explore large local datasets. You cannot see the data: "
    "use the tools, and report only what they return. If the tools cannot answer the question, "
    "say so and suggest what to run on the computer instead. Never ask for the data to be uploaded."
)


def main() -> None:
    question = " ".join(sys.argv[1:]) or "List the files and describe what is in them."
    client = anthropic.Anthropic()
    runner = client.beta.messages.tool_runner(
        model="claude-opus-5-5",
        max_tokens=16000,
        thinking={"type": "adaptive"},
        system=SYSTEM,
        tools=TOOLS,
        messages=[{"role": "user", "content": question}],
    )
    final = None
    for message in runner:       # the runner calls your tools and loops until Claude is done
        final = message
    print("".join(block.text for block in final.content if block.type == "text"))


if __name__ == "__main__":
    main()
