import os, sys
import pandas as pd
from scipy import stats
import subprocess


ENCODE = "/dors/capra_lab/projects/enhancer_ages/encode/data/ELS_combined_HepG2/breaks/truncated.bed"
BUILD = "hg38"
CL = "ELS_combined_HepG2"


RE = f"/dors/capra_lab/projects/enhancer_ages/encode/results/intronic_{CL}/"

if os.path.exists(RE) == False:
    os.mkdir(RE)

#%% Functions

def get_ensembl(build):
    path = "/dors/capra_lab/users/fongsl/data/ensembl/"
    edict = {
    "hg38":f"{path}ensGene_v32_hg38_coding_exons_autosome_merged.bed",
    "hg19":f"{path}ensGene_v36lift37_hg19_coding_exons_autosome_merged.bed"
    }

    ensembl = edict[build]
    return ensembl


def bedtools_intersect(encode, ensembl, cl, re, build):
    outf = f"{re}{cl}_x_ensembl_{build}.tsv"

    if os.path.exists(outf) == False: # if you haven't done this intersection already...
        cmd = f"bedtools intersect -a {encode} -b {ensembl} -wao > {outf}"
        subprocess.call(cmd, shell = True)
        print(cmd)

    else:
        print("you've done this intersection before")

    return outf

#%% get ensembl coding merged file, Run intersection

ENSEMBL = get_ensembl(BUILD)
outf = bedtools_intersect(ENCODE, ENSEMBL, CL, RE, BUILD)
outf
#%%
df = pd.read_csv(outf, sep = '\t', header = None)
cols = ["#chr_enh", "start_enh", "end_enh", "enh_id",
"core_remodeling","arch", "mrca", "#chr_ex", "start_ex",
"end_ex", "ex_len"]
df.columns = cols
df.head()
#%%
df['ex_bin'] = 0
df.loc[df.ex_len >0, 'ex_bin'] = 1

obsdf = df.groupby(["core_remodeling", "ex_bin"])["enh_id"].count().reset_index()

obsdf
complex_ex = obsdf["enh_id"].iloc[3]
complex_noex = obsdf["enh_id"].iloc[2]
simple_ex = obsdf["enh_id"].iloc[1]
simple_noex = obsdf["enh_id"].iloc[0]

obs = [[complex_ex, complex_noex], [simple_ex, simple_noex]]

a#%% are complex enhancers enriched for exon overlap?
stats.fisher_exact(obs) #(1.0391750501101307, 0.0723880988809173)

# complex enhancers may be slightly enriched for exon overlap
