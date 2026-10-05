# How to download the datasets

These commands are untested: GEO, EBI and SRA were blocked from the environment where this was compiled. Accessions are the
ones read in the papers (GSE86473, GSE154126, SRP075970). The FTP path pattern and tool names are from memory. If a path fails,
copy the exact link from the GEO page.

## Start here: Lawlor 2017, processed (GEO GSE86473)

### Option A: browser (simplest)
Open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86473. Scroll to **Supplementary file** at the bottom and use the (ftp) or (http) links. Also open **Series Matrix File(s)**: it carries the sample labels (control vs T2D, donor ID), which you must confirm.

### Option B: command line
GEO series folders follow a fixed pattern: the last three digits of the number become `nnn`.

```
mkdir -p data/GSE86473 && cd data/GSE86473
wget -r -np -nH --cut-dirs=4 -R "index.html*" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86473/suppl/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86473/matrix/
```

### Option C: Python
```
import GEOparse            # pip install GEOparse
gse = GEOparse.get_GEO("GSE86473", destdir="data/")
meta = gse.phenotype_data   # sample table: look for diagnosis and donor columns
gse.download_supplementary_files(directory="data/GSE86473_suppl")
```

### Option D: R
`GEOquery::getGEO("GSE86473")` and `getGEOSuppFiles("GSE86473")`.

Most supplementary files of this kind are plain-text count tables. Load with `pandas.read_csv(..., sep="\t", index_col=0)`, then build an AnnData object with `anndata.AnnData(df.T)`. Check first whether rows are genes or cells.

## Raw reads (FASTQ): only if you need to realign
Raw data for Lawlor is in SRA under SRP075970 (BioProject PRJNA323853).

```
# install sra-tools: conda install -c bioconda sra-tools
# get the run list: SRA Run Selector -> search SRP075970 -> download Accession List
prefetch --option-file SRR_Acc_List.txt
fasterq-dump --split-files SRRxxxxxxx
```
Raw FASTQ files are large: check disk space first. For a first project, use the processed matrix.

## Avrahami 2020 (GEO GSE154126)
Same steps as above. The FTP folder is `.../geo/series/GSE154nnn/GSE154126/`.

## Segerstolpe 2016
- Easiest: the authors' portal at http://sandberg.cmb.ki.se/pancreas.
- ArrayExpress E-MTAB-5061 (accession unverified): open the study page, go to the Files tab, download the processed counts file and the metadata (SDRF) file. The SDRF is the sample table: look in it for the disease label and the donor.

## HPAP (main project dataset)
HPAP is not a one-click download.
1. Go to PANC-DB (https://hpap.pmacs.upenn.edu, URL from memory) and register.
2. Accept the data-use terms. This may need an institutional signature or PI approval.
3. Download donor metadata and the scRNA-seq files. The Elgamal paper says raw FASTQ files come from there.
4. Ask your PI now: approval can take time.

Processing 60-plus donors from FASTQ also needs a cluster. Check whether Khalifa University HPC is available to you.

## Before you analyze
- Metadata first: build a donor table (donor ID, diagnosis, sex, age, BMI, batch). Count donors per group, not cells.
- Record every file's source URL and download date in `data/README.md`.
- Keep raw downloads read-only and never edit them.
