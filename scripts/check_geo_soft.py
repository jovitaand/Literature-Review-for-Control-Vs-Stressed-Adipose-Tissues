#!/usr/bin/env python3
"""Check a GEO family SOFT file (e.g. GSE128889_family.soft or .gz) against the problem statement using the sample records only.
Usage: python3 scripts/check_geo_soft.py FILE [--out report.csv]
Reports: species, human samples, depots (subcutaneous / visceral), lean/obese/diabetes words, single-cell vs bulk, and a verdict table.
It cannot count cells (SOFT files contain no counts)."""
import re, sys, gzip, csv, collections
def read(p): return (gzip.open(p, 'rt', errors='replace') if p.endswith('.gz') else open(p, errors='replace')).read()
def main(p):
    t = read(p); title = (re.findall(r'!Series_title = (.*)', t) or [''])[0]; acc = (re.findall(r'!Series_geo_accession = (.*)', t) or [''])[0]
    S = []
    for b in re.split(r'(?m)^\^SAMPLE = ', t)[1:]:
        g = lambda k: ' ; '.join(re.findall(r'!Sample_' + k + r'\w* = (.*)', b))
        S.append(dict(gsm=b.split('\n')[0].strip(), title=g('title'), organism=g('organism_ch1'), text=' '.join([g('title'), g('source_name_ch1'), g('characteristics_ch1'), g('description')]).lower(), desc=g('description'), strategy=g('library_strategy')))
    hs = [s for s in S if 'sapiens' in s['organism'].lower()]
    has = lambda L, pat: [s['gsm'] for s in L if re.search(pat, s['text'])]
    sc = [s['gsm'] for s in hs if re.search(r'single.?cell|scrna|10x|single.?nucle|snrna|drop', s['text'] + s['title'].lower())]
    V = [("Human adipose tissue", f"{len(hs)} of {len(S)} samples are human" + ("" if hs else " (NOT MET)")),
         ("Single-cell / single-nucleus samples (human)", f"{len(sc)} sample(s): {', '.join(sc) or 'none'}"),
         ("Subcutaneous", f"{len(has(hs, r'subcut|abdominal|inguinal'))} human sample(s)"), ("Visceral", f"{len(has(hs, r'visc|omen|mesent|perigon'))} human sample(s)" + ("" if has(hs, r'visc|omen') else " (NOT MET)")),
         ("Lean / obese / BMI labels (human)", f"lean {len(has(hs, r'lean'))}, obese {len(has(hs, r'obes'))}, BMI {len(has(hs, r'bmi'))}"),
         ("Type 2 diabetes labels (human)", f"{len(has(hs, r't2d|type 2|diabet'))} sample(s)" + ("" if has(hs, r't2d|type 2|diabet') else " (NOT MET)")),
         ("Cell counts / endothelial counts", "not available in a SOFT file; needs the count files (RAW.tar / matrices)")]
    print(f"{acc}: {title}\n"); [print(f"{a:48s} {b}") for a, b in V]
    return V
if __name__ == '__main__':
    V = main(sys.argv[1])
    if '--out' in sys.argv: csv.writer(open(sys.argv[sys.argv.index('--out') + 1], 'w', newline='')).writerows([("requirement", "result")] + V)
