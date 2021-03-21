import matplotlib.pyplot as plt
import os, sys
import pandas as pd
import pybedtools as pb
from scipy import stats
import seaborn as sns
import subprocess


PATH = "/data/hodges_lab/ATAC-STARR_V2/data/ATAC-STARR_cts_matricies/"
FILE = "GM12878inGM12878_counts.tsv"
OUTPATH = "/dors/capra_lab/users/fongsl/tyler/data/"

F = os.path.join(PATH, FILE)
F


#%%

# cut the coordinates


def cut_bed(f, outpath):

    outf = f'{outpath}GM12878inGM12878_counts.bed'

    if os.path.exists(outf) == False:

        # cut only bed coordinates and geneID
        cut_cmd = '''awk -F'\t' -v OFS="\t" '{ print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10}'  %s > %s''' % (f, outf)

        subprocess.call(cut_cmd, shell = True)

        os.chdir(outpath)

        # remove the header from the cut file
        tail_cmd = f"tail -n+3 {outf}> t && mv t {outf}"
        subprocess.call(tail_cmd, shell = True)

    return outf

def merge_bed(outf, outpath, f):

    merged_f = f"{outpath}merged.bed"
    merged_peak = f"{outpath}merged_peaks.bed"

    if os.path.exists(outf) == False:
        # merge the coordinates
        cmd = f"bedtools merge -i {outf} > {merged_f}"
        subprocess.call(cmd, shell = True)

        # intersect peaks and coordinates
        int_cmd = f"bedtools intersect -a {f} -b {merged_f} -wa -wb > {merged_peak}"
        subprocess.call(int_cmd, shell = True)

    return merged_peak

outf = cut_bed(F, OUTPATH)
merged_f = merge_bed(outf, OUTPATH, F)

#%%


cols = ["peakid", "chr", "start",
"end","strand", "len",
"dna_rep1", "dna_rep2", "rna_rep1", "rna_rep2"]

peaks = pd.read_csv(F, sep = '\t', skiprows =1)
peaks.columns = cols

merged_cols = ["chr", "start", "end"]

merged = pd.read_csv(merged_f, sep = '\t', names = merged_cols)
merged.head()

# 1.14e7 rows

#%%

peaks.iloc[:, 5:].describe() # look at just the length and replicate data
"""
        len	dna_rep1	dna_rep2	rna_rep1	rna_rep2
count	1.140877e+07	1.140877e+07	1.140877e+07	1.140877e+07	1.140877e+07
mean	2.044606e+02	6.855084e+02	5.931834e+02	1.142335e+03	1.386095e+03
std	9.589663e+01	6.152083e+02	5.363281e+02	1.060707e+03	1.293021e+03
min	2.900000e+01	2.000000e+00	2.000000e+00	0.000000e+00	0.000000e+00
25%	1.270000e+02	2.320000e+02	1.980000e+02	3.810000e+02	4.590000e+02
50%	1.820000e+02	4.940000e+02	4.240000e+02	7.980000e+02	9.660000e+02
75%	2.730000e+02	9.500000e+02	8.250000e+02	1.565000e+03	1.902000e+03
max	7.190000e+02	5.989000e+03	5.240000e+03	1.031800e+04	1.340000e+04

"""

#%%
peaks.head(10)


#%%


x = peaks["len"]

sns.distplot(x = x)
