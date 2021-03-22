import matplotlib.pyplot as plt
import os, sys
import pandas as pd
import pybedtools as pb
from scipy import stats
import seaborn as sns
import subprocess


PATH = "/dors/capra_lab/users/fongsl/tyler/data/"
FILE = "GG-LL_species-specific_OCRs_rank.bed"

F = os.path.join(PATH, FILE)


HARPATH = "/dors/capra_lab/data/human_accelerated_regions/"
HARFILE = "nchaes_merged_hg38.bed"
HARFILE = "Doan2016_HARs.liftOver.to.hg38.bed"
HARF = os.path.join(HARPATH, HARFILE)


#%%

def bed_intersect(f, harf, path):

    sid = (f.split("/")[-1]).split(".bed")[0]
    harid = (harf.split("/")[-1]).split(".")[0]
    outf = f"{path}{sid}_x_{harid}.bed"

    cmd = f"bedtools intersect -a {f} -b {harf} -wao > {outf}" # write all the overlapping and non-overlapping hars

    if os.path.exists(outf) == False:

        subprocess.call(cmd, shell = True)

    return outf



def open_df(outf):
    cols = ["chr", "start", "end", "species", "chr_h", "start_h", "end_h",
    "region", "gene_assoc", "something", "HAR_id", "overlap"]

    df = pd.read_csv(outf, sep = '\t', header = None)
    df.columns = cols

    df["overlap_bin"] = 0
    df.loc[df["overlap"] > 0, "overlap_bin"] = 1

    return df
#%%

# run intersection
outf = bed_intersect(F, HARF, PATH)


#%% open file
df = open_df(outf)
df.head()
df.shape
df.groupby("species")["overlap_bin"].sum()
"""
overlap_bin
GM12878_specific    3
LCL8664_specific    4
Name: overlap_bin, dtype: int64
"""

obs = [[3,4999], [4, 5002]]
OR, P = stats.fisher_exact(obs)
print(OR, P)
# OR = 0.75 P = 1.0
#%%
df.loc[df.overlap_bin>0]

"""
chr	start	end	species	chr_h	start_h	end_h	region	gene_assoc	something	HAR_id	overlap	overlap_bin
957	chr12	12582215	12583614	GM12878_specific	chr12	12582992	12583305	intergenic	DUSP16,CREBL2	B	HAR_Merge50-00400	313	1
4002	chr6	107668961	107669765	GM12878_specific	chr6	107668779	107669052	intergenic	SOBP,SCML4	P	HAR_Merge50-02190	91	1
4929	chrX	112667277	112667504	GM12878_specific	chrX	112667447	112668119	intronic	LHFPL1	P	HAR_Merge50-02701	57	1
5782	chr11	109532762	109533534	LCL8664_specific	chr11	109533149	109533429	intergenic	C11orf87,ZC3H12C	B	HAR_Merge50-00350	280	1
7942	chr3	87601002	87601287	LCL8664_specific	chr3	87600953	87601135	intergenic	POU1F1,HTR1F	B	HAR_Merge50-01634	133	1
8440	chr5	113266907	113267338	LCL8664_specific	chr5	113266925	113267105	intronic	MCC	P	HAR_Merge50-01998	180	1
9322	chr7	114345397	114345885	LCL8664_specific	chr7	114345587	114345780	intergenic	NONE,FOXP2	P	HAR_Merge50-02337	193	1
"""
