import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import subprocess

#%% dictionary of values

vals = {"enh_path": "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/CL_0000182_hepatocyte_expressed_enhancers/non-genic/no-exon_CL_0000182/ages",
"arens_path":"/dors/capra_lab/data/mpra/arensbergen18/",
"outdata_path": "/dors/capra_lab/projects/enhancer_ages/arensbergen2019/data/",
"sid": "CL_0000182_hepatocyte_expressed_enhancers"}

vals["enh"] = "%s/syn_breaks_no-exon_CL_0000182_ages.bed"% vals["enh_path"]
vals["arens"] = "%sSuRE_SNP_table_181029_hg19.bed" % vals["arens_path"]
vals["outfile"] = "%sarens_%s.bed" % (vals["outdata_path"], vals["sid"])


# bedtools intersection arensbergen19 x enhancers


cmd = "bedtools intersect -a %s -b %s -wao > %s" % (vals["arens"],
vals["enh"], vals["outfile"])

subprocess.call(cmd, shell = True)


#%% make the dataframe


columns = ["chr","start_snp", "end_snp",
           "SNP_ID","ref.element.count", "alt.element.count",
           "k562.ref.mean","k562.alt.mean", "k562.wilcoxp","k562.wilcoxp_random",
           "hepg2.ref.mean", "hepg2.alt.mean","hepg2.wilcoxp","hepg2.wilcoxp_random",
           "chr_syn", "start_syn", "end_syn", "enh_id",
           "chr_enh", "start_enh", "end_enh",
           "seg_index", "core_remodeling", "core",
           "mrca", "overlap"]

df = pd.read_csv(vals["outfile"], sep = '\t', header = None)

# rename columns
df.columns = columns

# syn_id
df["syn_id"]  = df.chr_syn + ":" + df.start_syn.map(str) + "-" + df.end_syn.map(str)


# mark all the significant variants
df["sig_hepg2_fdr05"] = 0

# 5% FDR from authors.
df.loc[df["hepg2.wilcoxp"]<0.006192715, "sig_hepg2_fdr05"] =1

df.head()


#%%
# 325 enhancers overlap arensbergen of 668 total enhancers)

len(df.enh_id.unique())

#%% 257 variants overlap simple enhancers, 149 enhancers.
#%% 310 variants overlap complex enhancers, 175 enhancers total


df.groupby(["core_remodeling"])["chr"].count()


#%%

len(df.loc[df.sig_hepg2_fdr05>0, "SNP_ID"].unique())

#%% Focus on only the enhancers that overlap SNPs

overlap = df.loc[df.overlap >0 , ["enh_id" ,"SNP_ID", "core_remodeling", "sig_hepg2_fdr05", "overlap", "core", "syn_id"]]
overlap.head(20)


#%% count the number of enhancers that overlap SNPs


n_simple_enhancer = len(overlap.loc[overlap.core_remodeling.astype(int) ==0]["enh_id"].unique())
n_complex_enhancer = len(overlap.loc[overlap.core_remodeling.astype(int) ==1]["enh_id"].unique())


n_complex_enhancer_core = len(overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.core.astype(int) ==1)].groupby("enh_id")["overlap"].count())

n_complex_enhancer_der = len(overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.core.astype(int) ==0)].groupby("enh_id")["overlap"].count())

print(n_simple_enhancer, n_complex_enhancer, n_complex_enhancer_core, n_complex_enhancer_der)

#%% count the number of enhancers that overlap significant SNPs

# sig SNPs in simple enhancers
n_sig_simple_enhancer = overlap.loc[(overlap.core_remodeling.astype(int) ==0)\
 & (overlap.sig_hepg2_fdr05.astype(int) > 0 )]["enh_id"].count()

 # sig SNPs in complex enhancers
n_sig_complex_enhancer = overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.sig_hepg2_fdr05.astype(int) > 0 )]["enh_id"].count()

 # sig SNPs in complex cores
n_sig_complex_core = overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.sig_hepg2_fdr05.astype(int) > 0 ) &(overlap.core.astype(int) ==1)]["syn_id"].count()

 # sig SNPs in complex cores
n_sig_complex_der = overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.sig_hepg2_fdr05.astype(int) > 0 ) &(overlap.core.astype(int) ==0)]["syn_id"].count()

print(n_sig_complex_enhancer, n_sig_complex_core,
 n_sig_complex_der, n_sig_simple_enhancer)


#%% simple enhancers have a higher fraction of significant variants than complex
# % of simple

print("% simple", n_sig_simple_enhancer/n_simple_enhancer)
9/149

# % of complex
print("% complex enhancer", n_sig_complex_enhancer/n_complex_enhancer)
10/175

# % of complex cores
print("% complex core", n_sig_complex_core/n_complex_enhancer_core)
6/100

# % of complex derived
print("% complex derived", n_sig_complex_der/n_complex_enhancer_der)
4/104

#%% all complex enhancers w/ sig SNP overlap.
# There are no cases where sig SNPS overlap both core and derived regions of same enhancer

overlap.loc[(overlap.core_remodeling.astype(int) ==1) & (overlap.sig_hepg2_fdr05.astype(int) > 0 )]



#%% get number of SNPS in simple, complex core, complex derived segments


overlap.groupby(["core_remodeling", "core"])["syn_id"].count()


#%% background % significant variants NOT overlapping enhancers.


no_overlap = len(df.loc[df.overlap.astype(int) == 0])
sig_no_overlap = len(df.loc[(df.overlap.astype(int) == 0) & (df.sig_hepg2_fdr05.astype(float) == 1)])
print(no_overlap,sig_no_overlap, sig_no_overlap/no_overlap)
