# GSE176171 (Emont et al., Nature 2022): criteria check (2 Nov 2026)

Drive file shared: `GSE176171.gz`, 45,122,572 bytes. **Its contents could not be read** (the Drive connector passes files under 10 MB; Drive and GEO file hosts are blocked). The check below uses the GEO series and 24 human sample records, the GEO supplementary file list, and the Single Cell Portal page (SCP1376).

The Drive file does not match any single file in the GEO supplementary directory (cell_metadata.tsv.gz is 27 MB; Hs10X counts.mtx.gz is 738 MB), so it is probably a re-saved or merged file. Run `notebooks/GSE176171_criteria_check.ipynb` on it to find out what it is.

| Requirement | Result | Evidence |
|---|---|---|
| Human adipose tissue | Met | 24 human samples (11 visceral/omental, 13 subcutaneous) from 13 people; 166,149 human nuclei (SCP1376) |
| Both subcutaneous and visceral | Met | both depots in the sample records |
| Control group | Met, by BMI only | lean BMI < 25: 3 people (22.2, 24.0, 24.3); overweight 1; obese 9 (BMI 30 to 49) |
| Type 2 diabetes group | **Not met** | no diabetes field in any sample record. UPMC participants (bariatric patients and lean controls) excluded anyone with type 1 or 2 diabetes. |
| Endothelial cells | Present | vascular cells reported (22,734 human vascular cells on SCP1376, which also includes pericytes, smooth muscle and lymphatic cells) |
| More than 10,000 endothelial cells | **Not established** | the endothelial-only human count was not found; the vascular total is 22,734 |
| Gene list | Not checked | needs the count matrix; the human features file is available at GEO |

**Verdict:** does not meet the problem statement, because there is no type 2 diabetes group. It would work as a lean vs obese comparison with both depots.
Raw human data are in dbGaP; processed counts are at GEO (Hs10X files, 738 MB matrix).
