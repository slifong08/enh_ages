import argparse
from collections import Counter
import glob
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy import stats
import seaborn as sns


SPECIES = 'hg19'


FANTOMPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"

FANTOMFILE = "syn_breaks_all_fantom_enh_ages.bed"
FANTOM = os.path.join(FANTOMPATH, "all_fantom_enh/ages", FANTOMFILE)

SAMPLE_ID = "all_fantom_enh_%s" % SPECIES

OUTPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh/"
RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/gc_content/"

MERGE_DISTANCES = 0

if MERGE_DISTANCES == 0:
    PROMOTERFILE = "promoter_data_hg19.bed"
    PROMOTERSAMPLE_ID = "all_fantom_pr_merged_%s" % SPECIES
else:
    PROMOTERFILE = "merged_%s_bp_promoters.bed" % MERGE_DISTANCES
    PROMOTERSAMPLE_ID = "all_fantom_pr_merged%s_%s" % (MERGE_DISTANCES, SPECIES)

PROMOTER = os.path.join(FANTOMPATH, PROMOTERFILE)



#%%


def getfasta(test_enh, sample_id, test_path, results_path, species):


    cat_fasta = "%s%s.fa" % (results_path, sample_id)

    if os.path.exists(cat_fasta) == False:
        os.chdir(test_path)

        chr_cmd = "awk '{print >$1\"_%s_temp.bed\"}' %s" % (sample_id, test_enh) # split test into chrN.bed files
        os.system(chr_cmd)

        enh_chr_list = glob.glob("%s/chr*_%s_temp.bed" % (test_path, sample_id)) # glob chromosomes

        for enh_chr in enh_chr_list: # intersect test chrN.bed w/ syntenic block

            chr_num = (enh_chr.split("/")[-1]).split("_")[0]

            dna_path = "/dors/capra_lab/data/dna/human/%s/" % species

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

    else:
        print("already done!")

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


def format_fa(fasta_file):

    cols = ["seq", "syn_id"]
    df = pd.read_csv(fasta_file, sep = ">",  header =None, names =cols )

    df["syn_id"] = df["syn_id"].shift(1) # need to shift col 1 by 1.
    df = df.dropna() # drop na


    df["GC"] = df["seq"].apply(lambda x: countgc(x)) # get GC content
    df["motif_len"] = df["seq"].apply(lambda x: len(x)) # get sequence len
    longer_than5 = df.loc[df.motif_len>=6].copy()
    longer_than5["GC_density"] = longer_than5["GC"].divide(longer_than5["motif_len"])

    return longer_than5


def format_syndf(enh_age_file):

    syn_cols = ["chr_syn", "start_syn", "end_syn",
    "enh_id",
    "chr", "start", "end",
    "seg_index", "core_remodeling", "core",
    "mrca",]

    syn = pd.read_csv(FANTOM, sep ='\t', header = None, names = syn_cols)

    syn["syn_id"] = syn.chr_syn + ":" + syn.start_syn.map(str) + "-" + syn.end_syn.map(str)
    syn["enh_len"] = syn.end - syn.start


    # age and taxon file
    syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

    syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]] # whittle down the df
    syn["mrca"] = syn["mrca"].round(3) # round the ages

    syn = pd.merge(syn, syn_gen_bkgd, how = "left", on = "mrca")

    syn["arch"] = "simple"

    syn.loc[(syn.core_remodeling ==1) & (syn.core ==1), "arch"] = "complex_core"
    syn.loc[(syn.core_remodeling ==1) & (syn.core ==0), "arch"] = "complex_derived"

    return syn


def merge_fa_ages(fa, syn):

    merged = pd.merge(fa, syn, how = "left", on = "syn_id")

    return merged


def make_pdf(file_name, RE):
    outfile = file_name +".pdf"
    outf = os.path.join(RE, outfile)

    return outf


def custom_round(x, base=10):
    return int(base * round(float(x)/base))


def prep_match_len_df(df, id_col, len_col):

    df = df[[id_col, len_col]].drop_duplicates() # get columns

    columns_names = ["matching_ids", "matching_len"]
    df.columns = columns_names # format column names

    return df


def match_len(df1, df2, base_len, balanced_ids_per_bin):

    # round_lengths

    df1.matching_len = df1.matching_len.astype(float).apply(lambda x: custom_round(x, base=base_len)) # round to the nearest 100bp

    df2.matching_len = df2.matching_len.astype(float).apply(lambda x: custom_round(x, base=base_len))

    # set of lengths in both datasets
    lens = set(list(df1.matching_len.unique()) + list(df2.matching_len.unique())) 

    match_dict = {}

    # per length bin
    for length in lens:

        # count the min number ids with that length bin in both datasets
        df2n = len(df2.loc[df2.matching_len == length])
        df1n = len(df1.loc[df1.matching_len == length])

        # choose the number of ids to match on
        # by scaling to the smaller number () of ids.
        min_ids = min(df2n, df1n)

        if length > 0 and min_ids > 0:

            # get ids per length-matched bins
            if balanced_ids_per_bin == 1: # make a balanced set with the same number of ids per bin
                df2_ids = df2.loc[df2.matching_len == length].sample(n = min_ids, replace = False) # sample w/ replacement
                df1_ids = df1.loc[df1.matching_len == length].sample(n = min_ids, replace = False) # sample w/ replacement

            else: # get all ids in length bins

                df2_ids = df2.loc[df2.matching_len == length].sample(n = df2n, replace = False) # sample w/ replacement
                df1_ids = df1.loc[df1.matching_len == length].sample(n = df1n, replace = False) # sample w/ replacement

            balanced = pd.concat([df1_ids, df2_ids])
            match_dict[length] = balanced

    final_matched_id = pd.concat(match_dict.values())

    return final_matched_id.matching_ids.unique()

def compare_gc_arch(df):

    core = df.loc[df.arch == "complex_core", "GC_density"]
    derived = df.loc[df.arch == "complex_derived", "GC_density"]
    simple = df.loc[df.arch == "simple", "GC_density"]

    cstat, cp = stats.mannwhitneyu(core, derived)
    print("core v. derived", cstat, cp )

    stat, p = stats.mannwhitneyu(simple, core)
    print("complex core v simple", stat, p )

    print(df.groupby("arch")["GC_density"].mean(), "mean")
#%% make fasta files


enhfasta = getfasta(FANTOM, SAMPLE_ID, FANTOMPATH, OUTPATH, SPECIES)
prfasta = getfasta(PROMOTER, PROMOTERSAMPLE_ID, FANTOMPATH, OUTPATH, SPECIES)


#%% format fasta files and calculate GC density

fa = format_fa(enhfasta)

prfa = format_fa(prfasta)

prfa["arch"] = "promoter"
prfa = prfa.loc[prfa.motif_len>=50]


syn = format_syndf(FANTOM)


#%% merge GC content and architecture

df = merge_fa_ages(fa, syn)


#%% evaluate length distribution
fig, ax =plt.subplots()
sns.kdeplot(prfa.motif_len, label = "promoter")
sns.kdeplot(fa.motif_len, label = "enhancer")
ax.legend()

#%% select only the records that match lengths

enh_to_match = prep_match_len_df(df, "enh_id", "enh_len")


pr_to_match = prep_match_len_df(prfa, "syn_id", "motif_len")

# match promoter lengths to enhancer lengths
base_len = 10
balanced_ids_per_bin = 1
matching_len_ids = match_len(enh_to_match, pr_to_match, base_len, balanced_ids_per_bin)

len(matching_len_ids) #29780 matching IDs


#%% plot matching lengths

match_prfa = prfa.loc[prfa["syn_id"].isin(matching_len_ids)]

match_enh = df.loc[df["enh_id"].isin(matching_len_ids)]

fig, ax =plt.subplots()
sns.kdeplot(match_prfa["motif_len"], label = "promoter")
sns.kdeplot(match_enh["enh_len"], label = "enhancer")
ax.legend()
#%%


all_matched = pd.concat([match_enh, match_prfa])

all_matched.describe()
#%%

fig, ax = plt.subplots(figsize = (6,6))

hue_order = ["simple", "complex_core", "complex_derived", "promoter"]
x = "arch"
y = "GC_density"
data = all_matched
file_name = "GC_density_all_FANTOM_arch"
outfile = make_pdf(file_name, RE)

sns.set("poster")
xlabs = ["simple", "complex\ncore", "complex\nderived", "promoter"]
sns.barplot(x, y, data = data, order = hue_order)
ax.set(ylim = (0.35, 0.7), xticklabels = xlabs)
plt.savefig(outfile, bbox_inches = 'tight')

#%% all enhancers

compare_gc_arch(df)

compare_gc_arch(match_enh)

#%%


fig, ax = plt.subplots(figsize = (15,10))

x = "mrca_2"
y = "GC_density"
hue = "arch"
file_name = "GC_density_all_FANTOM_arch_mrca"
outfile = make_pdf(file_name, RE)

sns.barplot(x, y, data = data, hue = hue, hue_order = hue_order)
#ax.set(ylim = (0.35, 0.6))
plt.savefig(outfile, bbox_inches = 'tight')

#%%

count_agearch = all_matched.groupby(["arch", "mrca_2"])["syn_id"].count().reset_index()

fig, ax = plt.subplots(figsize = (10,10))
x = "mrca_2"
y = "syn_id"
hue = "arch"
data = count_agearch
file_name = "count_all_FANTOM_arch"
outfile = make_pdf(file_name, RE)

sns.barplot(x, y, data = count_agearch, hue = hue, hue_order = hue_order)
plt.savefig(outfile, bbox_inches = 'tight')
#%%

prfa.describe()
