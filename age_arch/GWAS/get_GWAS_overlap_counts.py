import glob
import os, sys
from pybedtools import BedTool
import subprocess

GWAS = "/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/gwasCatalog_2019-09-24_hg19_unique_cleaned_LDEx_p5e-8.bed"
PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/architecture_coordinates/all_fantom/"
FS = glob.glob(f"{PATH}*.bed")

for F in FS:
    obs_sum = 0

    A = BedTool(F)
    B = BedTool(GWAS)
    obs_intersect = A.intersect(B, wo = True)
    for line in obs_intersect:
        obs_sum += int(line[-1]) # sum snps in obs regions
    print(F, obs_sum)
    
