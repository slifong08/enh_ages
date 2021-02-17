import argparse
from collections import Counter
import glob
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy import stats
import seaborn as sns



FANTOMPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh/ages/"
FANTOMFILE = "syn_breaks_all_fantom_enh_ages.bed"
FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

SAMPLE_ID = "all_fantom_enh_hg19"

OUTPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh/"
RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/gc_content/"

#%%

def getfasta(test_enh, sample_id, test_path, results_path):
    os.chdir(test_path)

    chr_cmd = "awk '{print >$1\"_%s_temp.bed\"}' %s" % (sample_id, test_enh) # split test into chrN.bed files
    os.system(chr_cmd)

    enh_chr_list = glob.glob("%s/chr*_%s_temp.bed" % (test_path, sample_id)) # glob chromosomes

    for enh_chr in enh_chr_list: # intersect test chrN.bed w/ syntenic block

        chr_num = (enh_chr.split("/")[-1]).split("_")[0]

        dna_path = "/dors/capra_lab/data/dna/human/hg19/"

        fasta_in = ("%s%s.fa" % (dna_path, chr_num))
        fasta_out = "%s%s_%s.fa" % (results_path, sample_id, chr_num)

        cmd = "bedtools getfasta -fi %s -bed %s -fo %s" %(fasta_in, enh_chr, fasta_out)

        os.system(cmd)

        cleanup_cmd = "rm %s" % enh_chr # clean up chr temp file
        os.system(cleanup_cmd)


    cat_fasta = "%s%s.fa" % (results_path, sample_id)
    cat_cmd = "cat %s%s_chr*.fa > %s" % (results_path, sample_id, cat_fasta)
    os.system(cat_cmd)

    cleanup_cmd = "rm %s%s_chr*.fa" % (results_path, sample_id)
    os.system(cleanup_cmd)

    return cat_fasta



def countgc(sequence):
    gc = []
    total = []
    letters = ["G", "C", "g", "c"]
    counts = Counter(sequence)

    for letter in letters:
        gc.append(int(counts[letter]))

    gc_sum = sum(gc)
    return gc_sum


#%%


outfasta = getfasta(FANTOM, SAMPLE_ID, FANTOMPATH, OUTPATH)


#%%


cols = ["seq", "syn_id"]
df = pd.read_csv(outfasta, sep = ">",  header =None, names =cols )

df["syn_id"] = df["syn_id"].shift(1) # need to shift col 1 by 1.
df = df.dropna() # drop na
df

#%%

df["GC"] = df["seq"].apply(lambda x: countgc(x)) # get GC content
df["motif_len"] = df["seq"].apply(lambda x: len(x)) # get sequence len
df["GC_density"] = df["GC"].divide(df["motif_len"] )



df
#%% layer architecture info


syn_cols = ["chr_syn", "start_syn", "end_syn",
"enh_id",
"chr", "start", "end",
"seg_index", "core_remodeling", "core",
"mrca",]

syn = pd.read_csv(FANTOM, sep ='\t', header = None, names = syn_cols)
syn["syn_id"] = syn.chr_syn + ":" + syn.start_syn.map(str) + "-" + syn.end_syn.map(str)
syn.head()

#%%

# age and taxon file
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]] # whittle down the df
syn["mrca"] = syn["mrca"].round(3) # round the ages
syn = pd.merge(syn, syn_gen_bkgd, how = "left", on = "mrca")

#%% merge GC content and architecture

merged = pd.merge(df, syn, how = "left", on = "syn_id")

merged
#%%


merged["arch"] = "simple"
merged.loc[(merged.core_remodeling ==1) & (merged.core ==1), "arch"] = "complex_core"
merged.loc[(merged.core_remodeling ==1) & (merged.core ==0), "arch"] = "complex_derived"


#%%
fig, ax = plt.subplots(figsize = (6,6))

hue_order = ["simple", "complex_core", "complex_derived"]
x = "arch"
y = "GC_density"
outfile = os.path.join(RE, "GC_density_all_FANTOM_arch.pdf")
print(outfile)

sns.barplot(x, y, data = merged, order = hue_order)
ax.set(ylim = (0.45, 0.5))
plt.savefig(outfile, bbox_inches = 'tight')

#%%


core = merged.loc[merged.arch == "complex_core", "GC_density"]
derived = merged.loc[merged.arch == "complex_derived", "GC_density"]
simple = merged.loc[merged.arch == "simple", "GC_density"]

cstat, cp = stats.mannwhitneyu(core, derived)
print("core v. derived", cstat, cp )

stat, p = stats.mannwhitneyu(simple, core)
print("complex core v simple", stat, p )

merged.groupby("arch")["GC_density"].mean()

#%%


fig, ax = plt.subplots(figsize = (15,10))

x = "mrca_2"
y = "GC_density"
hue = "arch"
outfile = os.path.join(RE, "GC_density_all_FANTOM_arch_mrca.pdf")

sns.barplot(x, y, data = merged, hue = hue, hue_order = hue_order)
ax.set(ylim = (0.35, 0.6))
plt.savefig(outfile, bbox_inches = 'tight')

#%%

count_agearch = merged.groupby(["arch", "mrca_2"])["syn_id"].count().reset_index()

fig, ax = plt.subplots(figsize = (10,10))
x = "mrca_2"
y = "syn_id"
hue = "arch"
outfile = os.path.join(RE, "count_all_FANTOM_arch.pdf")

sns.barplot(x, y, data = count_agearch, hue = hue, hue_order = hue_order)
plt.savefig(outfile, bbox_inches = 'tight')
