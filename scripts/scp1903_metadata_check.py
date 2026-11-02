#!/usr/bin/env python3
"""SCP1903 (Single Cell Portal, human adipose single-nucleus RNA-seq) criteria check from the portal metadata files.
Usage: python3 scripts/scp1903_metadata_check.py METADATA.tsv RAW_GENES.tsv.gz OUT.xlsx
Inputs: metadata.tsv (cell-level: cluster, donor, depot, wtstatus ...), raw_genes.tsv.gz (gene list). No count matrix is needed for the counts below."""
import sys, gzip, pandas as pd
meta, genes, out = sys.argv[1:4]
m = pd.read_csv(meta, sep='\t', skiprows=[1]); g = [l.strip() for l in gzip.open(genes, 'rt')]
GL = ['PECAM1', 'CDH5', 'NOS3', 'HIF1A', 'ERG', 'VCAM1', 'ICAM1', 'ANGPT2']
blood = m[m.seurat_clusters.isin(['EC', 'EC-P'])]; lec = m[m.seurat_clusters == 'LEC']
by_cluster = pd.crosstab(m.seurat_clusters, [m.depot, m.wtstatus], margins=True)
ec_tab = m[m.seurat_clusters.isin(['EC', 'EC-P', 'LEC'])].groupby(['depot', 'wtstatus', 'seurat_clusters']).size().unstack(fill_value=0)
ec_tab['blood_EC'] = ec_tab.get('EC', 0) + ec_tab.get('EC-P', 0)
per_sample = m.groupby('biosample_id').agg(donor=('donor_id', 'first'), depot=('depot', 'first'), wtstatus=('wtstatus', 'first'), sex=('sex', 'first'), nuclei=('NAME', 'size')).join(blood.groupby('biosample_id').size().rename('blood_EC')).fillna(0)
per_sample['pct_blood_EC'] = (100 * per_sample.blood_EC / per_sample.nuclei).round(1)
group = m.groupby(['depot', 'wtstatus']).agg(nuclei=('NAME', 'size'), donors=('donor_id', 'nunique')).join(blood.groupby(['depot', 'wtstatus']).agg(blood_EC=('NAME', 'size')))
genes_tab = pd.DataFrame({'gene': GL, 'in_raw_gene_list': [x in g for x in GL]})
disease = m.groupby(['disease__ontology_label', 'wtstatus']).size().rename('nuclei').reset_index()
verdict = pd.DataFrame([
 ('Control group (lean or obese without diabetes)', f"MET (lean): {m[m.wtstatus=='lean'].donor_id.nunique()} lean donors, {m[m.wtstatus=='obese'].donor_id.nunique()} obese donors; no diabetes status recorded"),
 ('Type 2 diabetes group', "NOT MET: no diabetes label; the 'disease' field is 'obesity disorder' for every cell"),
 ('Endothelial cells present', f"MET: clusters EC, EC-P (blood, {len(blood)} nuclei) and LEC (lymphatic, {len(lec)})"),
 ('More than 10,000 endothelial cells', f"NOT MET: {len(blood)} blood endothelial nuclei ({len(blood)+len(lec)} with lymphatic) of {len(m)} nuclei"),
 ('Human adipose tissue', 'MET: 6 donors, 10x 3 prime v2 single-nucleus RNA-seq'),
 ('Both subcutaneous and visceral', f"MET: visceral {int((m.depot=='vat').sum())} and subcutaneous {int((m.depot=='sat').sum())} nuclei; 5 donors with both depots"),
 ('Gene list', f"{int(genes_tab.in_raw_gene_list.sum())} of 8 genes in the gene list (PECAM1 absent); expression needs the count matrix")], columns=['requirement', 'result'])
with pd.ExcelWriter(out) as w:
    verdict.to_excel(w, sheet_name='Verdict', index=False); group.to_excel(w, sheet_name='Nuclei_by_depot_group'); ec_tab.to_excel(w, sheet_name='Endothelial_by_depot_group'); per_sample.to_excel(w, sheet_name='Per_sample')
    by_cluster.to_excel(w, sheet_name='Clusters'); genes_tab.to_excel(w, sheet_name='Gene_list_presence', index=False); disease.to_excel(w, sheet_name='Disease_field', index=False)
print(verdict.to_string(index=False))
