#!/usr/bin/env python3
"""Merge the GSE129363 discovery-cohort files into one cell table, one gene table and one sample table, and write the schema.
Usage: python3 scripts/merge_gse129363_schema.py ANNOTATION.txt.gz BARCODES.tsv.gz GENES.tsv.gz OUT_DIR
Inputs: CellAnnotation (CellID, SampleName, Condition, Tissue), barcodes (matrix column order), genes (Ensembl ID, symbol; matrix row order).
Also uses data/metadata/GSE129363_sample_map_and_counts.csv (barcode suffix -> sample -> GEO record)."""
import sys, pandas as pd, pathlib
an_f, bc_f, ge_f, out = sys.argv[1:5]; out = pathlib.Path(out); out.mkdir(parents=True, exist_ok=True)
root = pathlib.Path(__file__).resolve().parent.parent
an = pd.read_csv(an_f, sep='\t'); bc = pd.read_csv(bc_f, header=None, names=['barcode']); ge = pd.read_csv(ge_f, sep='\t', header=None, names=['ensembl_id', 'symbol'])
smap = pd.read_csv(root / 'data/metadata/GSE129363_sample_map_and_counts.csv')
# ---- cells (one row per barcode = one matrix column)
bc['matrix_col'] = range(1, len(bc) + 1); bc['barcode_suffix'] = bc.barcode.str.rsplit('-', n=1).str[1].astype(int)
cells = bc.merge(an.rename(columns={'CellID': 'barcode', 'SampleName': 'sample_name', 'Condition': 'condition_annotation', 'Tissue': 'tissue'}), on='barcode', how='left')
cells['in_annotation'] = cells.sample_name.notna()
s = smap.rename(columns={'sample_name': 'sample_name_map', 'condition_in_annotation': 'condition_annotation_map', 'condition_in_GEO': 'condition_geo', 'tissue': 'tissue_map', 'GEO_sample': 'geo_sample_title', 'GSM': 'gsm'})
s['person'] = s.geo_sample_title.str.split('-').str[1]; s['condition_geo'] = s.condition_geo.str.replace(' ', '').str.replace('NonDiabetic', 'NonDiabetic')
s['label_conflict'] = s.condition_geo != s.condition_annotation_map
cells = cells.merge(s[['barcode_suffix', 'gsm', 'geo_sample_title', 'person', 'condition_geo', 'label_conflict']], on='barcode_suffix', how='left')
cells['group_annotation'] = cells.condition_annotation.map({'NonDiabetic': 'control (obese, no T2D)', 'Diabetic': 'T2D (obese + T2D)'})
cells['group_geo'] = cells.condition_geo.map({'NonDiabetic': 'control (obese, no T2D)', 'Diabetic': 'T2D (obese + T2D)'})
cells = cells[['matrix_col', 'barcode', 'barcode_suffix', 'in_annotation', 'sample_name', 'gsm', 'geo_sample_title', 'person', 'tissue', 'condition_annotation', 'condition_geo', 'label_conflict', 'group_annotation', 'group_geo']]
# ---- genes (one row per gene = one matrix row)
ge['matrix_row'] = range(1, len(ge) + 1); ge['symbol_is_duplicated'] = ge.symbol.duplicated(keep=False)
ge['symbol_unique'] = ge.symbol.where(~ge.symbol.duplicated(), ge.symbol + '.' + ge.groupby('symbol').cumcount().astype(str))
GL = {'ERG': 'ERG', 'PECAM1': 'CD31', 'CDH5': 'VE-cadherin', 'NOS3': 'eNOS', 'VCAM1': 'VCAM1', 'ICAM1': 'ICAM1', 'HIF1A': 'HIF-1alpha', 'ANGPT2': 'ANGPT2'}
ge['gene_list_name'] = ge.symbol.map(GL); ge['in_gene_list'] = ge.gene_list_name.notna()
# ---- samples (one row per barcode suffix)
cnt = cells.groupby('barcode_suffix').agg(barcodes=('barcode', 'size'), cells_annotated=('in_annotation', 'sum')).reset_index()
samples = s[['barcode_suffix', 'sample_name_map', 'gsm', 'geo_sample_title', 'person', 'tissue_map', 'condition_annotation_map', 'condition_geo', 'label_conflict']].rename(columns={'sample_name_map': 'sample_name', 'tissue_map': 'tissue', 'condition_annotation_map': 'condition_annotation'}).merge(cnt, on='barcode_suffix')
cells.to_csv(out / 'GSE129363_cells_merged.csv', index=False); ge.to_csv(out / 'GSE129363_genes_merged.csv', index=False); samples.to_csv(out / 'GSE129363_samples_merged.csv', index=False)
print('cells', cells.shape, '| annotated', int(cells.in_annotation.sum()), '| genes', ge.shape, '| duplicated symbols', int(ge.symbol_is_duplicated.sum()), '| samples', samples.shape)
print('gene-list rows:'); print(ge[ge.in_gene_list][['matrix_row', 'ensembl_id', 'symbol', 'gene_list_name']].to_string(index=False))
chk = cells[cells.in_annotation].merge(smap[['barcode_suffix', 'tissue']].rename(columns={'tissue': 'tissue_map'}), on='barcode_suffix'); print('tissue consistent with sample map:', bool((chk.tissue == chk.tissue_map).all()))
print('label conflicts: samples', int(samples.label_conflict.sum()), '| annotated cells affected', int(cells[cells.in_annotation].label_conflict.sum()))
print(cells.groupby(['tissue', 'group_annotation']).size().unstack().to_string()); print(cells.groupby(['tissue', 'group_geo']).size().unstack().to_string())
