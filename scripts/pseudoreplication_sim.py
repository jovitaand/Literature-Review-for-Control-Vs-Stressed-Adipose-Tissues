# Toy simulation: false-positive rate of cell-level vs pseudobulk tests when there is NO condition effect.
# Donor-level random effect on the mean (log-scale SD), negative-binomial counts. Requires numpy and scipy.
# Run: python3 scripts/pseudoreplication_sim.py
import numpy as np
from scipy import stats
rng=np.random.default_rng(1)
def run(donors=10,cells=200,sd=0.3,genes=2000,mu=1.0,disp=0.5):
    fp_cell=fp_w=fp_pb=0
    for g in range(genes):
        eff=rng.normal(0,sd,2*donors)           # donor-level random effect, no condition effect
        lam=mu*np.exp(eff)
        X=[]
        for d in range(2*donors):
            r=1/disp; p=r/(r+lam[d])
            X.append(rng.negative_binomial(r,p,cells))
        X=np.array(X)
        grp=np.repeat([0,1],donors)
        a=np.log1p(X[grp==0].ravel()); b=np.log1p(X[grp==1].ravel())
        fp_cell+= stats.ttest_ind(a,b,equal_var=False).pvalue<0.05
        fp_w+= stats.mannwhitneyu(a,b).pvalue<0.05
        pb=np.log1p(X.sum(1)/cells*1000)
        fp_pb+= stats.ttest_ind(pb[grp==0],pb[grp==1],equal_var=False).pvalue<0.05
    return fp_cell/genes,fp_w/genes,fp_pb/genes
for d,c,sd in [(10,200,0.3),(10,1000,0.3),(10,200,0.1),(3,200,0.3)]:
    print(d,c,sd,run(d,c,sd,genes=1500))
