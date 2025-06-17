import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import subprocess

#%% dictionary of values

vals = {"enh_path": "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh/ages",
"arens_path":"/dors/capra_lab/data/mpra/arensbergen18/",
"outdata_path": "/dors/capra_lab/projects/enhancer_ages/arensbergen2019/data/",
"sid": "all_fantom_enh"}

vals["enh"] = "%s/syn_breaks_all_fantom_enh_ages.bed"% vals["enh_path"]
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


# mark all the significant variants w/ binary
# 5% FDR from authors.
if "hepato" in vals["sid"]:
    df["sig_fdr05"] = 0
    df.loc[df["hepg2.wilcoxp"]<0.006192715, "sig_fdr05"] =1
    print("HepG2")

elif "granulocyte" in vals["sid"]:
    df["sig_fdr05"] = 0
    df.loc[df["k562.wilcoxp"]<0.00173121, "sig_fdr05"] =1
    print("K562")

else:
    df["sig_fdr05_k562"] = 0
    df.loc[df["k562.wilcoxp"]<0.00173121, "sig_fdr05_k562"] =1
    df["sig_fdr05_hepg2"] = 0
    df.loc[df["hepg2.wilcoxp"]<0.006192715, "sig_fdr05_hepg2"] =1


#%%
# 13265 enhancers overlap arensbergen of 30474 total enhancers)

len(df.enh_id.unique())

#%% 257 variants overlap simple enhancers, 149 enhancers.
#%% 310 variants overlap complex enhancers, 175 enhancers total


df.groupby(["core_remodeling"])["chr"].count()


#%%

len(df.loc[df.sig_fdr05_k562>0, "SNP_ID"].unique())
len(df.loc[df.sig_fdr05_hepg2>0, "SNP_ID"].unique())

#%% Focus on only the enhancers that overlap SNPs

overlap = df.loc[df.overlap >0 , ["enh_id" ,"SNP_ID", "core_remodeling",

overlap.head()


#%% count the number of enhancers that overlap SNPs


n_simple_enhancer = len(overlap.loc[overlap.core_remodeling.astype(int) ==0]["enh_id"].unique())
n_complex_enhancer = len(overlap.loc[overlap.core_remodeling.astype(int) ==1]["enh_id"].unique())


n_complex_enhancer_core = len(overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.core.astype(int) ==1)].groupby("enh_id")["overlap"].count())

n_complex_enhancer_der = len(overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.core.astype(int) ==0)].groupby("enh_id")["overlap"].count())

print(n_simple_enhancer, n_complex_enhancer, n_complex_enhancer_core, n_complex_enhancer_der)



#%% count all SNPS overlapping simple, complex,

n_simple_enhancer = len(overlap.loc[overlap.core_remodeling.astype(int) ==0])
n_complex_enhancer = len(overlap.loc[overlap.core_remodeling.astype(int) ==1])


n_complex_enhancer_core = len(overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.core.astype(int) ==1)])

n_complex_enhancer_der = len(overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.core.astype(int) ==0)])

print(n_simple_enhancer, n_complex_enhancer, n_complex_enhancer_core, n_complex_enhancer_der)

###
# HEPG2
###

#%% count the number of significant SNPs

# sig SNPs in simple enhancers
n_sig_simple_enhancer = len(overlap.loc[(overlap.core_remodeling.astype(int) ==0)\
 & (overlap.sig_fdr05_hepg2.astype(int) > 0 )])

 # sig SNPs in complex enhancers
n_sig_complex_enhancer = len(overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.sig_fdr05_hepg2.astype(int) > 0 )])

 # sig SNPs in complex cores
n_sig_complex_core = overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.sig_fdr05_hepg2.astype(int) > 0 ) &(overlap.core.astype(int) ==1)]["syn_id"].count()

 # sig SNPs in complex cores
n_sig_complex_der = overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.sig_fdr05_hepg2.astype(int) > 0 ) &(overlap.core.astype(int) ==0)]["syn_id"].count()

print(n_sig_simple_enhancer, n_sig_complex_enhancer, n_sig_complex_core,
 n_sig_complex_der)

#%% simple enhancers have a higher fraction of significant variants than complex
# % of simple

print("% simple", n_sig_simple_enhancer/n_simple_enhancer)


# % of complex
print("% complex enhancer", n_sig_complex_enhancer/n_complex_enhancer)


# % of complex cores
print("% complex core", n_sig_complex_core/n_complex_enhancer_core)


# % of complex derived
print("% complex derived", n_sig_complex_der/n_complex_enhancer_der)


#%% all complex enhancers w/ sig SNP overlap.
# There are no cases where sig SNPS overlap both core and derived regions of same enhancer

overlap.loc[(overlap.core_remodeling.astype(int) ==1) & (overlap.sig_fdr05_k562.astype(int) > 0 )]



#%% get number of SNPS in simple, complex core, complex derived segments


overlap.groupby(["core_remodeling", "core"])["syn_id"].count()


#%% background % significant variants NOT overlapping enhancers.


no_overlap = len(df.loc[df.overlap.astype(int) == 0])
sig_no_overlap = len(df.loc[(df.overlap.astype(int) == 0) & (df.sig_fdr05_hepg2.astype(float) == 1)])
print(no_overlap,sig_no_overlap, sig_no_overlap/no_overlap)

###
# k562
###

#%% count the number of significant SNPs

# sig SNPs in simple enhancers
n_sig_simple_enhancer = len(overlap.loc[(overlap.core_remodeling.astype(int) ==0)\
 & (overlap.sig_fdr05_k562.astype(int) > 0 )])

 # sig SNPs in complex enhancers
n_sig_complex_enhancer = len(overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.sig_fdr05_k562.astype(int) > 0 )])

 # sig SNPs in complex cores
n_sig_complex_core = overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.sig_fdr05_k562.astype(int) > 0 ) &(overlap.core.astype(int) ==1)]["syn_id"].count()

 # sig SNPs in complex cores
n_sig_complex_der = overlap.loc[(overlap.core_remodeling.astype(int) ==1)\
 & (overlap.sig_fdr05_k562.astype(int) > 0 ) &(overlap.core.astype(int) ==0)]["syn_id"].count()

print(n_sig_simple_enhancer, n_sig_complex_enhancer, n_sig_complex_core,
 n_sig_complex_der)

#%% simple enhancers have a higher fraction of significant variants than complex
# % of simple

print("% simple", n_sig_simple_enhancer/n_simple_enhancer)


# % of complex
print("% complex enhancer", n_sig_complex_enhancer/n_complex_enhancer)


# % of complex cores
print("% complex core", n_sig_complex_core/n_complex_enhancer_core)


# % of complex derived
print("% complex derived", n_sig_complex_der/n_complex_enhancer_der)



#%% get number of SNPS in simple, complex core, complex derived segments


overlap.groupby(["core_remodeling", "core"])["syn_id"].count()


#%% background % significant variants NOT overlapping enhancers.


no_overlap = len(df.loc[df.overlap.astype(int) == 0])
sig_no_overlap = len(df.loc[(df.overlap.astype(int) == 0) & (df.sig_fdr05_k562.astype(float) == 1)])
print(no_overlap,sig_no_overlap, sig_no_overlap/no_overlap)
