from itertools import combinations, combinations_with_replacement, product
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm
import subprocess


from joblib import Parallel, delayed
import multiprocessing


ENHBASE = "/dors/capra_lab/projects/enhancer_ages/encode/hepg2/data/"

ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/"

RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/cCRE_x_tfbs_encode3/HepG2"



colors = [ "amber", "dusty purple", "windows blue"]
PAL = sns.xkcd_palette(colors)
sns.palplot(PAL)

colors = [ "windows blue"]
DERPAL = sns.xkcd_palette(colors)
sns.palplot(DERPAL)

all_simple_jaccard
#%% Functions


def get_cell_lines():
    sample_dict = {
    "HepG2": "no-exon_dELS_combined",
    }

    return sample_dict


def get_paths(cell_line, file_tag, fantombase, encodepath):

    FANTOMPATH = os.path.join(fantombase, file_tag, "ages")
    FANTOMFILE = "syn_breaks_%s_ages.bed" % file_tag
    FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

    if "CL" in cell_line:
        ENCODEFILE = "cells/" + cell_line.split("_CL")[0] + ".bed.gz"
    elif cell_line == "all_fantom_enh":
        ENCODEFILE = "trimmed_encRegTfbsClusteredWithCells.bed"
    else:
        ENCODEFILE = "cells/" + cell_line + ".bed.gz"

    ENCODE = os.path.join(encodepath, ENCODEFILE)


    INTERSECTIONPATH = os.path.join(fantombase, file_tag)
    INTERSECTIONFILE = "%s_x_ENCODE.bed" % file_tag
    INTERSECTION = os.path.join(INTERSECTIONPATH, INTERSECTIONFILE)

    return FANTOM, ENCODE, INTERSECTION


def bed_intersect(fantom, encode, intersection):

    if os.path.exists(intersection) == False:
        cmd = "bedtools intersect -a %s -b %s -wao > %s" % (fantom, encode, intersection)

        subprocess.call(cmd, shell = True)

        print(cmd)
    else:
        print("previously done enh x encode intersection")


def format_df(intersection_file):

    cols = ["chr_syn", "start_syn", "end_syn",
    "enh_id","chr", "start", "end",
    "seg_index", "core_remodeling", "core",
    "mrca",
    "chr_tf", "start_tf", "end_tf",
    "tf_id", "peak_len",
    "tf", "cell_line" , "overlap"
    ]

    df = pd.read_csv(intersection_file,
    sep = '\t',
    header = None).drop_duplicates()

    df.columns = cols # add column names

    df["tf"] = df["tf_id"].apply(lambda x: x.split("_")[0])

    # add architecture label - core, derived, simple
    df["arch"] = "simple"
    df.loc[(df.core_remodeling ==1) & (df.core ==1), "arch"] = "complex_core"
    df.loc[(df.core_remodeling ==1) & (df.core ==0), "arch"] = "complex_derived"

    # add architecture label - complex, simple
    df["overallarch"] = "simple"
    df.loc[df.core_remodeling ==1, "overallarch"] = "complex"

    # add syn identifier
    df["syn_id"] = df.chr_syn + ":" + df.start_syn.map(str) + "-" + df.end_syn.map(str)

    #calculate enhancer and syntenic block length
    df["enh_len"] = df.end - df.start
    df["syn_len"] = df.end_syn - df.start_syn
    df.loc[df.syn_len <6, "syn_len"] = 0


    # binary for TF overlap, any TF that overlaps less than 6bp is not counted.
    df["tfoverlap_bin"] = 1
    df.loc[df.tf == ".", "tfoverlap_bin"] = 0
    df.loc[df.overlap <6 , "tfoverlap_bin"] = 0

    return df


def get_df(cell_line, val, fantombase, encodepath):

    print(cell_line, val)
    fantom, encode, intersection = get_paths(cell_line, val, fantombase, encodepath)

    #Bed command
    bed_intersect(fantom, encode, intersection)

    #dataframe
    df = format_df(intersection)

    return df


def get_mrca2(df):

    SYN_GROUP = "/dors/capra_lab/projects/enhancer_ages/hg38_syn_taxon.bed"
    syn = pd.read_csv(SYN_GROUP, sep = '\t')

    # round all values
    syn[["mrca", "mrca_2"]] = syn[["mrca", "mrca_2"]].round(3)
    df.mrca = df.mrca.round(3)

    df = pd.merge(df, syn, how = "left")


    cores = df.groupby("enh_id")["mrca_2", 'taxon2'].max().reset_index()
    cores.columns = ["enh_id", "coremrca_2", "coretaxon2"]

    df = pd.merge(df, cores, how = "left")

    return df


def describe_tf_overlaps(df):

    #%% count how many tfs have how many overlaps

    tf_count = df.groupby("tf")["syn_id"].count().reset_index() # count all the TFs in simple
    tf_count.columns = ["tf", "tf_count"]

    len(tf_count) # there are 117 total TFs in simple eutherian enhancers

    median = tf_count.tf_count.median()


    fig, ax =plt.subplots(figsize = (6,6))
    sns.histplot(tf_count.tf_count)

    ax.set( xlabel = f"peaks per dataset\nmedian peaks/dataset = {median}",
    ylabel = "TF datasets",
     title = "all HepG2 TF counts")

    outfile = f'{RE}HepG2_ENCODE_TFBS_tf_count.pdf'

    plt.savefig(outfile, bbox_inches = 'tight')


# get TFs with min number of instances.
def get_tf_combos(taxon2, arch, min_instances, df ):

    #constrain df to age and architecture of interest


    if "complex" in arch: # get TF counts based on core age

        age_arch = df.loc[(df.coretaxon2 == taxon2) & (df.arch.str.contains(arch))]

    elif taxon2 == "all":
        age_arch = df.loc[(df.arch.str.contains(arch))]

    else:
        age_arch = df.loc[(df.taxon2 == taxon2) & (df.arch.str.contains(arch))]

    # count how many enhancers are in this category
    enh_count = len(age_arch.enh_id.unique())

    print(enh_count, "enhancers in", taxon2, arch)

    # count the number of peaks per TF in architecture
    tf_count = age_arch.groupby("tf")["syn_id"].count().reset_index()
    tf_count.columns = ["tf", "tf_count"] # rename columns

    # get all possible combos of TF w/ peak count greater than min number of instances of peak.
    tf_min_list = tf_count.loc[(tf_count.tf_count >= MIN_INSTANCES) & (tf_count.tf !="."), "tf"]


    # create all combinations of TFs w/ counts >MIN_INSTANCES
    #tf_combos = list(combinations(tf_min_list, 2)) # AB AC AD BC BD CD
    #tf_combos = list(combinations_with_replacement(tf_min_list, 2)) # AA AB AC AD BB BC BD CC CD DD

    tf_combos = list(product(tf_min_list, repeat = 2)) #AA AB AC AD BA BB BC BD CA CB CC CD DA DB DC DD


    print(len(tf_min_list), "TFs over min_instance thresh", len(tf_combos), "TF combinations")

    return tf_combos


def run_parallel_jaccard(taxon2, arch, df, arch_comparison, enhbase, tf_combo_list):

    ### RUN PARALLEL Jaccards  ASSEMBLY ###

    num_cores = multiprocessing.cpu_count()
    print("number of cores", num_cores)

    # run parallel jobs
    outdir = os.path.join(enhbase, arch_comparison)

    if os.path.exists(outdir) ==False: # make a directory for results
        os.mkdir(outdir)

    jaccardfile = f"summary_jaccard_{arch_comparison}.txt"
    jaccardf = os.path.join(outdir, jaccardfile)

    if os.path.exists(jaccardf) == False:

        Parallel(n_jobs=num_cores, verbose=100, prefer="threads")(delayed(run_jaccard)(taxon2, arch, df, arch_comparison, outdir, tf_pair) for tf_pair in tf_combo_list)

    return outdir


def run_jaccard(taxon2, arch, df, arch_comparison, outdir, tf):

    tf1 = tf[0]
    tf2 = tf[1]

    if taxon2 != "all": # filter dataframe by age and architecture?
        test = df.loc[(df.coretaxon2 == taxon2) & (df.arch.str.contains(arch))]

    else:
        test = df.loc[(df.arch.str.contains(arch))] # or filter just on architecture.


    if "core_v_der" in arch_comparison: # if you're comparing core v. derived
        df1 = test.loc[(test.core ==1)]
        df2 = test.loc[(test.core !=1)]

    else:
        df1 = test
        df2 = test

    # make a comparison name
    comparison_name = tf1 + "/" + tf2 + "-" + taxon2


    results = jaccard_function(tf1, tf2, df1, df2)

    results["arch_comparison"] = arch_comparison

    outfile = f"{arch_comparison}_{taxon2}_{tf1}_{tf2}.tsv"
    outf = os.path.join(outdir, outfile)

    results.to_csv(outf, sep ='\t', header = False, index = False)

    return age_arch_results


def jaccard_function(tf1, tf2, df1, df2):

    tf1Set = set(df1.loc[(df1.tf == tf1), "enh_id"])
    tf2Set = set(df2.loc[(df2.tf == tf2), "enh_id"])


    intersection = len(tf1Set.intersection(tf2Set))
    union = (len(tf1Set) + len(tf2Set)) - intersection

    jaccard_index = float(intersection) / union

    jaccard_df = pd.DataFrame({"tf1": [tf1],
    "tf2": [tf2],
    "tf1Set_len": [len(tf1Set)],
    "tf2Set_len": [len(tf2Set)],
    "intersection": [intersection],
    "union": [union],
    "jaccard_index":[jaccard_index]
    })

    return jaccard_df


def open_df(outdir, arch_comparison):

    jaccardfile = f"summary_jaccard_{arch_comparison}.txt"
    jaccardf = os.path.join(outdir, jaccardfile)

    if os.path.exists(jaccardf) == False:
        cmd = f"cat {outdir}/*.tsv > {jaccardf}"
        subprocess.call(cmd, shell = True)

    print(jaccardf)
    cols = ["tf1", "tf2", "tf1Set_len", "tf2Set_len",
    "intersection", "union", "jaccard_index", "arch_comparison"]

    jaccardDF = pd.read_csv(jaccardf, sep = '\t', header = None, names =cols)

    if len(jaccardDF) >5: # clean up the other files.
        cmd = f'rm {outdir}/*.tsv'
        subprocess.call(cmd, shell = True)

    return jaccardDF_complex


def plot_jaccard(jaccardDF, arch_comparison, re):

    # HISTOGRAM

    sns.set_context("talk")
    fig, ax = plt.subplots(figsize = (6,6))

    sns.histplot(jaccardDF.jaccard_index)
    median = jaccardDF.jaccard_index.median()
    ax.set( xlabel = f"Jaccard index\nmedian Jaccard = {median}",
    ylabel = "count",
    title = f"{arch_comparison}")

    outfile = f'{re}HepG2_ENCODE_TFBS_jaccard_hist_{arch_comparison}.pdf'

    plt.savefig(outfile, bbox_inches = 'tight')


    jaccard_table = pd.pivot(jaccardDF.sort_values(by="tf1"),
                            values = "jaccard_index",
                            index = "tf2",
                            columns = "tf1")


    sns.set_context("paper")

    mask = jaccard_table < 0 # mask any zero values.


    # CLUSTER MAP
    g = sns.clustermap(jaccard_table.fillna(-0.001),
    figsize = (12,12),
    #row_cluster = False,
    mask = mask)

    g.fig.suptitle(f'{arch_comparison} TF Co-binding\n Jaccard Index')

    outfile = f'{RE}{arch_comparison}_HepG2_ENCODE_TFBS_co-binding_heatmap_clustermap.pdf'
    plt.savefig(outfile, bbox_inches = 'tight')


#%%


sample_dict = get_cell_lines()

ALPHA = 0.1
MIN_INSTANCES = 1000

cell_line = "HepG2"
val = sample_dict[cell_line]


df = get_df(cell_line, val, ENHBASE, ENCODEPATH)

print("this many TFS overlap dataset", len(df.tf.unique()))

#%%

df = get_mrca2(df)
describe_tf_overlaps(df)

#%% run all simple


TAXON2 = "all"
ARCH = "simple"
ARCH_COMPARISON = "simple"
tf_combos  = get_tf_combos(TAXON2, ARCH, MIN_INSTANCES, df)


outdir = run_parallel_jaccard(TAXON2, ARCH, df, ARCH_COMPARISON, ENHBASE, tf_combos)
outdir
new = open_df(outdir, ARCH_COMPARISON)

# Jaccard hist
plot_jaccard(new, ARCH_COMPARISON, RE)

del ARCH, ARCH_COMPARISON, tf_combos, outdir, new

#%% Run all complex


TAXON2 = "all"
ARCH = "complex"
ARCH_COMPARISON = "complex"
tf_combos  = get_tf_combos(TAXON2, ARCH, MIN_INSTANCES, df)


outdir = run_parallel_jaccard(TAXON2, ARCH, df, ARCH_COMPARISON, ENHBASE, tf_combos)
outdir
new = open_df(outdir, ARCH_COMPARISON)

# Jaccard hist
plot_jaccard(new, ARCH_COMPARISON, RE)

del ARCH, ARCH_COMPARISON, tf_combos, outdir, new


#%%


TAXON2 = "all"
ARCH = "complex"
ARCH_COMPARISON = "core_v_der"

tf_combos  = get_tf_combos(TAXON2, ARCH, MIN_INSTANCES, df)


outdir = run_parallel_jaccard(TAXON2, ARCH, df, ARCH_COMPARISON, ENHBASE, tf_combos)

new = open_df(outdir, ARCH_COMPARISON)

# Jaccard hist
plot_jaccard(new, ARCH_COMPARISON, RE)

del ARCH, ARCH_COMPARISON, tf_combos, outdir, new

#%%


TAXON2 = "Eutheria"
MIN_INSTANCES = 500
ARCH = "simple"
ARCH_COMPARISON = "simple_eutherian"

tf_combos  = get_tf_combos(TAXON2, ARCH, MIN_INSTANCES, df)


outdir = run_parallel_jaccard(TAXON2, ARCH, df, ARCH_COMPARISON, ENHBASE, tf_combos)

new = open_df(outdir, ARCH_COMPARISON)

# Jaccard hist
plot_jaccard(new, ARCH_COMPARISON, RE)

del ARCH, ARCH_COMPARISON, tf_combos, outdir, new

#%%

TAXON2 = "Eutheria"
ARCH = "complex"
ARCH_COMPARISON = "complex_eutherian"

tf_combos  = get_tf_combos(TAXON2, ARCH, MIN_INSTANCES, df)


outdir = run_parallel_jaccard(TAXON2, ARCH, df, ARCH_COMPARISON, ENHBASE, tf_combos)

new = open_df(outdir, ARCH_COMPARISON)

# Jaccard hist
plot_jaccard(new, ARCH_COMPARISON, RE)

del ARCH, ARCH_COMPARISON, tf_combos, outdir, new

#%%

TAXON2 = "Eutheria"
ARCH = "complex"
ARCH_COMPARISON = "core_v_der_eutherian"

tf_combos  = get_tf_combos(TAXON2, ARCH, MIN_INSTANCES, df)

outdir = run_parallel_jaccard(TAXON2, ARCH, df, ARCH_COMPARISON, ENHBASE, tf_combos)

new = open_df(outdir, ARCH_COMPARISON)

# Jaccard hist
plot_jaccard(new, ARCH_COMPARISON, RE)

del ARCH, ARCH_COMPARISON, tf_combos, outdir, new
