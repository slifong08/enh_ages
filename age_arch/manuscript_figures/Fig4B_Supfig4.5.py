from collections import Counter
import glob
import pandas as pd
from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
import numpy as np
import os, sys
import pybedtools as pb
from scipy import stats
import seaborn as sns
import statsmodels.api as sm

# colors
faded_green = "#7bb274"
amber = "#feb308"


# sns colors
arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
arch_palette = sns.xkcd_palette(arch_colors)

colors = [ "amber", "faded green", "dusty purple", "windows blue","greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(arch_palette)



#%% contrls


RUN_BED = 0
RUN_SHUF=0
RUN_COUNTS = 0
RUN_OR = 0


#%% paths


sid_dict = {"UBERON_0002107": "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/UBERON_0002107_liver_expressed_enhancers/UBERON_0002107_liver_expressed_enhancers_enh_ages_enh_age_arch_summary_matrix.bed",
            "CL_0000094": "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/CL_0000094_granulocyte_expressed_enhancers/CL_0000094_granulocyte_expressed_enhancers_enh_ages_enh_age_arch_summary_matrix.bed",
            "E123": "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E123/trimmed/no-exon_trimmed-310-E123.bed",
            "E118": "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E118/trimmed/no-exon_trimmed-310-E118.bed",
            "all_fantom" : "/dors/capra_lab/projects/enhancer_ages/fantom/data/trimmed_all_unique_fantom_enh_112_tissue.bed"
           }

arens_path = "/dors/capra_lab/data/mpra/arensbergen18/"
outdata_path = "/dors/capra_lab/projects/enhancer_ages/arensbergen2019/data/"
outdata_shuffle_path = "/dors/capra_lab/projects/enhancer_ages/arensbergen2019/data/shuffle/"
RE = "/dors/capra_lab/projects/enhancer_ages/arensbergen2019/results/"

#%% intersect enhancers, shuffle, w/ arensbergen
# only do this once. Takes a long time


if RUN_BED ==1:
    arens = pb.BedTool("%sSuRE_SNP_table_181029_hg19.bed" % (arens_path)).sort()

    for sid, f in sid_dict.items():

        print(sid)

        outfile = "%sarens_%s.bed" % (outdata_path, sid)

        if RUN_BED ==1:
            enh = pb.BedTool(f).sort()
            arens.intersect(enh, wao = True).saveas(outfile)


if RUN_SHUF ==1:

    shuf_path_dict = {"all_fantom_enh":"/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/first_round_breaks/"}
    shuf_filename_dict = {"all_fantom_enh":"shuf-all_fantom_enh_age_breaks_summary_matrix_noexon-"}

    arens = pb.BedTool("%sSuRE_SNP_table_181029_hg19.bed" % (arens_path)).sort()

    for sid, shufpath in shuf_path_dict.items():


        for i in np.arange(1,100):
            shuffile = shuf_filename_dict[sid] + str(i) + ".bed"

            SHUFF = os.path.join(shufpath, shuffile)

            outfile = "%sarens_shuf-%s-%s.bed" % (outdata_shuffle_path, sid, i)

            if RUN_SHUF ==1:
                print(i)
                enh = pb.BedTool(SHUFF).sort()
                arens.intersect(enh, wao = True).saveas(outfile)

# ran up to shuffle iteration 32, then quit because there is no shuffle file iteration 33.


#%% FUNCTIONS


def assign_columns(sid, df):

    # different files have different column lengths. This is to address those differences.

    fantom_sids = ["UBERON_0002107", "CL_0000094"]
    roadmap_sids = ["E123", "E118"]

    if sid in fantom_sids:
        df = df[[0, 1,2,3,4,5,6,7,8,9,10,11,12,13, 17, 19, 20, 21, 22, 23, 32]]

        columns = ["chr","start_snp", "end_snp",
                   "SNP_ID","ref.element.count", "alt.element.count",
                   "k562.ref.mean","k562.alt.mean", "k562.wilcoxp","k562.wilcoxp_random",
                   "hepg2.ref.mean", "hepg2.alt.mean","hepg2.wilcoxp","hepg2.wilcoxp_random",
               "enh_id", "sid", "core_remodeling","arch", "seg_index", "mrca","overlap"]
        df.columns = columns
        df["core_remodeling"] =0

    elif sid in roadmap_sids:
        df = df[[0, 1,2,3,4,5,6,7,8,9,10,11,12,13, 17, 19, 20, 21, 30]]

        columns = ["chr","start_snp", "end_snp", "SNP_ID","ref.element.count",
              "alt.element.count","k562.ref.mean","k562.alt.mean",
              "k562.wilcoxp","k562.wilcoxp_random","hepg2.ref.mean",
              "hepg2.alt.mean","hepg2.wilcoxp","hepg2.wilcoxp_random",
                "enh_id", "arch",  "seg_index","mrca",
               "overlap"]

        df.columns = columns
        df["core_remodeling"] =0
        df.loc[df.seg_index ==1, "core_remodeling"] = 1


    elif sid == "all_fantom":
        columns = ["chr","start_snp", "end_snp",
           "SNP_ID","ref.element.count","alt.element.count",
           "k562.ref.mean","k562.alt.mean", "k562.wilcoxp","k562.wilcoxp_random",
           "hepg2.ref.mean", "hepg2.alt.mean","hepg2.wilcoxp","hepg2.wilcoxp_random",
           "chr_enh", "start_enh", "end_enh",
             "enh_id", "core_remodeling", "mrca", "arch", "dataset", "dataset2", "overlap"]
        df.columns = columns
        df["seg_index"] = 1

    df = df.drop_duplicates()

    return df


def format_df(sid, df):

    fantom_sids = ["UBERON_0002107", "CL_0000094"]
    roadmap_sids = ["E123", "E118"]
    # name the columns. This is annoying bc I made different file formats.

    df = assign_columns(sid, df)

    # replace all the missing data with a 1
    df[["k562.ref.mean","k562.alt.mean",
          "k562.wilcoxp","k562.wilcoxp_random","hepg2.ref.mean",
          "hepg2.alt.mean","hepg2.wilcoxp","hepg2.wilcoxp_random"]] = df[["k562.ref.mean","k562.alt.mean",
          "k562.wilcoxp","k562.wilcoxp_random","hepg2.ref.mean",
          "hepg2.alt.mean","hepg2.wilcoxp","hepg2.wilcoxp_random"]].replace(".",1)

    # float columns
    df[["k562.ref.mean","k562.alt.mean", "k562.wilcoxp","hepg2.ref.mean",
                  "hepg2.alt.mean", "hepg2.wilcoxp"]]\
    = df[["k562.ref.mean","k562.alt.mean", "k562.wilcoxp","hepg2.ref.mean",
                  "hepg2.alt.mean", "hepg2.wilcoxp"]].astype(float)

    # reduce dataframe
    df2 = df[["enh_id",
               "k562.ref.mean", "k562.alt.mean","k562.wilcoxp", "core_remodeling",
               "hepg2.ref.mean", "hepg2.alt.mean", "hepg2.wilcoxp", "overlap", "seg_index"]]\
    .drop_duplicates()


    # replace infinite values with 0
    df2 = df2.replace([np.inf, -np.inf], 0)

    # created significant activity 0|1 binaries
    # based on nominal p-value = 0.05 or 5% FDR p-value cut off determined in arensbergen 2019 publication
    df2["sig_k562_05"],  df2["sig_hepg2_05"], df2["sig_k562_fdr05"], df2["sig_hepg2_fdr05"] = 0, 0, 0, 0

    df2.loc[df2["k562.wilcoxp"]<0.05, "sig_k562_05"] =1
    df2.loc[df2["hepg2.wilcoxp"]<0.05, "sig_hepg2_05"] =1

    df2.loc[df2["k562.wilcoxp"]<0.00173121, "sig_k562_fdr05"] =1
    df2.loc[df2["hepg2.wilcoxp"]<0.006192715, "sig_hepg2_fdr05"] =1

    return df2


def assign_core_remodeling(sid, df):

    # reassign enhancer architectures (aka core_remodeling). If not enhancer, arch = -1

    if sid != "all_fantom":
        df.loc[df.seg_index !="1", "core_remodeling"] = 1
        df.loc[df.seg_index =="1", "core_remodeling"] = 0
        df.loc[df.seg_index ==".", "core_remodeling"] = -1
    print(df.core_remodeling.unique())

    return df


def bootstrapCI(list): # bootstrap C.I.s

    xbar = np.mean(list) # get the real median fold-change from 500 shuffles
    n = len(list)
    nboot = 10000 # resample 10000 times
    val = 0
    bs_means = []

    while val < nboot:

        bs_dist = np.random.choice(list, replace = True, size = n)
        bsmean = np.mean(bs_dist)
        bs_means.append(bsmean)
        val +=1

    bs = pd.DataFrame(data = bs_means, index = np.arange(nboot), columns = ["bs_means"])

    bs["deltas"] = bs.bs_means - xbar

    bs = bs.sort_values(by = "deltas", ascending= False)

    low = bs.deltas.quantile(0.025)
    high = bs.deltas.quantile(0.975)
    print(low, high)

    ci = xbar - [high, low]
    stdev = np.std(bs_means)
    return ci, stdev


def get_arch_snp_counts(df, sig):

    # Simple
    simple_list = df.loc[(df.core_remodeling ==0) & df.overlap ==1, sig].to_list()

    simple_total = len(simple_list)
    simple_nonsig = Counter(simple_list)[0] # count non
    simple_sig = Counter(simple_list)[1] # count sig overlap

    # Complex
    complex_list = df.loc[(df.core_remodeling ==1) & df.overlap ==1, sig].to_list()

    complex_total = len(complex_list)
    complex_nonsig =  Counter(complex_list)[0]
    complex_sig =  Counter(complex_list)[1]

    # bkgd (w/o enhancer overlap)
    bkgd_list =  df.loc[(df.core_remodeling ==-1) & df.overlap ==0, sig].to_list() # no core_remodeling, no enh overlap

    bkdg_sig = Counter(bkgd_list)[1]
    bkdg_nonsig = Counter(bkgd_list)[0]

    count_list = [simple_sig, simple_nonsig, simple_total, complex_sig, complex_nonsig, complex_total, bkdg_sig, bkdg_nonsig]

    return count_list


def make_count_df(cell_model, sig, count_list, sid):

    count_id_list = ["simple_sig", "simple_nonsig", "simple_total",
                    "complex_sig","complex_nonsig", "complex_total",
                    "bkdg_sig", "bkdg_nonsig"]

    newdf = pd.DataFrame({
         "count_id": count_id_list,
         "counts": count_list,
          })

    newdf["cell_model"], newdf["sig_metric"], newdf["sid"] = cell_model, sig, sid
    newdf["count_list"] = [count_list]

    return newdf


def estimate_OR(count_list, cell_model):

    simple_sig, simple_nonsig, simple_total = count_list[0], count_list[1], count_list[2]
    complex_sig,complex_nonsig, complex_total = count_list[3], count_list[4], count_list[5]
    bkdg_sig, bkdg_nonsig = count_list[6], count_list[7]


    # Test 1 - sig simple v. all background
    #[[simple sig, simple not sig], [allbkgd sig(minus enh), allbkgd not sig(minus enh)]]


    COMPARISON1 = "sig simple v. bkgd"
    obs1 = np.array([[simple_sig, simple_nonsig], [bkdg_sig, bkdg_nonsig]]) # make a contigency table

    simple_v_bkgd = get_OR(obs1, COMPARISON1, cell_model) # calculate OR

    key = cell_model + str(nom_p) + COMPARISON1
    f_results[key] = simple_v_bkgd


    # Test 2 - sig complex v. all background
    #[[complex sig, complex not sig], [allbkgd sig(minus enh), allbkgd not sig(minus enh)]]


    COMPARISON2 = "sig complex v. bkgd"
    obs2 = np.array([[complex_sig, complex_nonsig], [bkdg_sig, bkdg_nonsig]]) # make a contigency table

    complex_v_bkgd = get_OR(obs2, COMPARISON2, cell_model)# calculate OR

    key = cell_model + str(nom_p) + COMPARISON2
    f_results[key] = complex_v_bkgd


    # Test 3 - sig simple v. Signsigificant complex
    #[[complex sig, complex not sig], [allbkgd sig(minus enh), allbkgd not sig(minus enh)]]


    COMPARISON3 = "sig simple v. sig complex"
    obs3 = np.array([[simple_sig, simple_nonsig], [complex_sig, complex_nonsig]]) # make a contigency table

    simple_v_complex = get_OR(obs3, COMPARISON3, cell_model)# calculate OR

    key = cell_model + str(nom_p) + COMPARISON3
    f_results[key] = simple_v_complex


    concat_results = pd.concat(f_results.values()) # concat all results for this cell line

    return concat_results


def get_OR(obs, comparison, cell_model):
    print(obs, comparison, cell_model)
    odds, p = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()


    results = pd.DataFrame({"sid": [sid],
            "cell_model": [cell_model],
            "nominal_p":[nom_p],
            "a": [obs[0][0]],
            "b": [obs[0][1]],
            "c": [obs[1][0]],
            "d": [obs[1][1]],
            "OR": [odds],
            "p": [p],
            "test":[comparison],
            "ci_lower" :[odds_ci[0]],
            "ci_upper" :[odds_ci[1]]})

    return results


def format_OR_dataframe(df):
    # match cell-specific enhancers with cell model background tested in MPRA

    df["cell_line_match"] = 0 # create dummy column

    hepg2_sids = ["UBERON_0002107", "E118", "all_fantom"] # hepg2/liver facet ids

    k562_sids = ["CL_0000094", "E123", "all_fantom"] # k562/granulocyte facet ids

    df.loc[(df.cell_model =="HEPG2") & (df.sid.isin(hepg2_sids)), "cell_line_match"] = 1
    df.loc[(df.cell_model =="K562") & (df.sid.isin(k562_sids)), "cell_line_match"] = 1

    df["-log10p"] = -1*np.log10(df.p) # -log10p

    df["log2or"] = np.log2(df.OR) # log2(OR)

    # dataset names
    df['dataset'] = "Fantom_all"
    df.loc[df.sid.str.contains("E1"), "dataset"] = "Roadmap-Matched Cell line"
    df.loc[df.sid.str.contains("000"), "dataset"] = "Fantom-Matched Tissue"
    df["test2"] = df["cell_model"] + "-" + df["test"] # label for plotting Figure 4a

    return df


def count_OR_df(count):


    count["SUMab"] = count.a + count.b # total alleles in simple enhancers

    count["SUMcd"] = count.c + count.d # total alleles in complex enhancers
    count["total"] = count.a + count.b + count.c + count.d # total alleles in either enhancer

    count["FRACa"] = count.a.divide(count.SUMab) # fraction of sig alleles in simple enhancers
    count["FRACab"] = count.SUMab.divide(count.total) # fraction of simple alleles in total allele dataset
    count["FRACc"] = count.c.divide(count.SUMcd) # fraction of sig alleles in complex enhancers
    count["FRACcd"] = count.SUMcd.divide(count.total) # fraction of complex alleles in total allele dataset

    return count


def sim_and_bootstrap_freq_overlap(df):


    sig, non_sig = df["a"].iloc[0], df["b"].iloc[0] # get counts of sig, non sig enhancer overlap

    freq_list = []

    # make a list of 0s and 1s
    sample = [1]*sig + [0]*non_sig
    #this is a vector of all enhancer architectures overlapping sig/non-sig
    #MPRA variants in arensbergen

    total = sig + non_sig

    # create 10k randomly sampled vectors to estimate frequency of overlaps

    # sample w/ replacement!
    # calculate frequency of 1s per randomly sampled dataset
    # append frequency to frequency_list

    for n in np.arange(1000): # sample 1000x

        sample_list = np.random.choice(sample, replace = True, size =total)
        # size = total # of enhancers
        # with replacement
        # sample 1s and 0s
        random_sig = Counter(sample_list)[1] # count 1s
        freq = random_sig/total # get frequency of sig. MPRA variants randomly sampled
        freq_list.append(freq) # append frequency to list.

    sns.countplot(freq_list)
    ci, std = bootstrapCI(freq_list) # calculate bootstrapped cis

    df["ci_lower"], df["ci_upper"] =  ci[0], ci[1]

    df["stdev"] = std # stdev on the bootstrappedCI distribution

    return df


def plot_freq(test, cell_model, enh_dataset, outf):


    simple = test.loc[test.test == "sig simple v. bkgd"][["FRACa", "yerr", "FRACc"]]
    simVal, simErr = simple.FRACa.iloc[0], simple.yerr.iloc[0]

    complexenh = test.loc[test.test == "sig complex v. bkgd"][["FRACa", "yerr"]]
    comVal, comErr = complexenh.FRACa.iloc[0], complexenh.yerr.iloc[0]

    bkgd = simple.FRACc.iloc[0] # % bkgd sig variants (not overlapping enh)

    sids = cell_model


    ind = np.arange(1)  # the x locations for the groups
    width = 0.2  # the width of the bars
    barWidth = 0.25


    fig, ax = plt.subplots(figsize = (6,6))

    # Set position of bar on X axis
    r1 = 1
    r2 = 1 + barWidth


    # Make the plot
    sns.set("poster")
    sns.set_style("white")


    plt.bar(r1, simVal, color=amber, width=barWidth, edgecolor='white',
    label='simple', yerr =simErr, error_kw=dict(lw=5))
    plt.bar(r2, comVal, color=faded_green, width=barWidth, edgecolor='white',
    label='complexenh', yerr =comErr, error_kw=dict(lw=5))


    ax.axhline(bkgd, color = "k", ls = "--", label = "bkgd")

    ax.set(ylabel = "% sig. MRPA variants\nin architecture",
            xlabel = cell_model,
            title = enh_dataset,
            xticklabels = ["",  "simple", "complex"])
    if enh_dataset == "UBERON_0002107":
        ax.set_ylim(0,0.08)
    else:
        ax.set_ylim(0, 0.05)
    plt.savefig(outf, bbox_inches = "tight")

#%% gather all the files


fs = glob.glob("%sarens_*.bed" % (outdata_path))


#%% filter for the files you want to analyze?


select_sids =["UBERON_0002107", "E123", "CL_0000094", "E118", "all_fantom"]

fsnew = {}

for f in fs:

    f_sid = "_".join(((f.split("/")[-1]).split(".")[0]).split("_")[1:])
    if f_sid in select_sids:

        fsnew[f_sid] = f
fsnew


#%% prepare to count overlaps


MODELS = ["K562", "HEPG2"]
nom_p = "FDR05"
sig_dict= {'K562': "sig_k562_fdr05", 'HEPG2':"sig_hepg2_fdr05"}

ci = {}
all_OR_results = {} # collect all results for both K562 and HepG2 cell lines
all_counts_results = {}


#%% count overlaps


RUN_COUNTS = 0


if RUN_COUNTS ==1:

    for sid, f in fsnew.items():

        print(sid, f)

        key = str(sid+"-"+nom_p)

        if key not in all_counts_results.keys():

            df = pd.read_csv(f, sep = '\t', header = None) # open the file

            df = format_df(sid, df) # format dataframe

            df = assign_core_remodeling(sid, df) # assign core_remodeling


            for CELL_MODEL in MODELS:

                sig = sig_dict[CELL_MODEL]

                ### Get numbers ###
                count_list = get_arch_snp_counts(df, sig)
                all_counts_results[key] = make_count_df(CELL_MODEL, sig, count_list, key)



    counts = pd.concat(all_counts_results.values())
    counts.to_csv("%sall_trimmed_counts.tsv"% (outdata_path), sep = '\t', header= True, index = True)

else:

    counts = pd.read_csv("%sall_trimmed_counts.tsv"% (outdata_path), sep = '\t')


#%% calculate OR


if RUN_OR ==1:


    for sid in df.sid.unique():
        count_list = df.loc[df.sid == sid, "counts"].to_list()
        all_OR_results[key] = estimate_OR(count_list, sid)


    df = pd.concat(all_OR_results.values())
    df.sort_values(by ="sid")
    df.to_csv("%sall_trimmed_or.tsv"% (outdata_path), sep = '\t', header= True, index = False)

elif RUN_OR ==0:
    df = pd.read_csv("%sall_trimmed_or.tsv"% (outdata_path), sep = '\t')

ordf = format_OR_dataframe(df)
ordf = ordf.loc[ordf.cell_line_match ==1] # get only the matching enhnacer cell lines and cell models


#%% simulate and estimate confidence intervals for frequency overlap


zipped = zip(ordf.test2, ordf.sid)

ci_dict = {}
#%%
val = 0
for test2, sid in zipped:


    if "bkgd" in test2:

        test_df = ordf.loc[ (ordf.test2 == test2)
         & (ordf.sid == sid)]

        print(test2, sid, len(test_df))

        ci_results = sim_and_bootstrap_freq_overlap(test_df)

        key = val

        #ci_dict[key] = ci_results

        val +=1
        break




#%%
ci_dict.keys()
len(new_ordf)
new_ordf = pd.concat(ci_dict.values()).drop_duplicates()

new_ordf["yerr"] = new_ordf.ci_upper- new_ordf.ci_lower

count = count_OR_df(new_ordf) # focus on allele overlap between simple and complex only.

count.dataset.unique()


#%%


plot_zip = zip(count.cell_model, count.sid)

for cell_model, enh_dataset in plot_zip:

    test = count.loc[(count.cell_model == cell_model) &(count.sid == enh_dataset)]
    outf = "%sfig4b-%s_arens_FRAC_nomp-%s.pdf" % (RE, enh_dataset, cell_model)
    plot_freq(test, cell_model, enh_dataset, outf)

#%%
count.sid.unique()
fantom_all = count.loc[count.dataset == "Fantom_all"]
fantom_all
E123 = count.loc[count.sid == "E123"]
E123
E118 = count.loc[count.sid == "E118"]
CL = count.loc[count.sid.str.contains("CL_0000094")]
CL
UB = count.loc[count.sid.str.contains("UBERON_0002107")]
UB
def ttest_perm(df):

    for cell_model in df.cell_model.unique():
        test = df.loc[(df.cell_model == cell_model)]
        simple = test.loc[test.test== "sig simple v. bkgd"]
        complexenh = test.loc[test.test == "sig complex v. bkgd"]

        simple_mean, simple_std = simple["FRACa"].iloc[0], simple["stdev"].iloc[0]
        complex_mean, complex_std = complexenh["FRACa"].iloc[0], complexenh["stdev"].iloc[0]
        result, p = stats.ttest_ind_from_stats(mean1 = simple_mean,
                                           std1 = simple_std,
                                           nobs1 = 1000,
                                          mean2 = complex_mean,
                                           std2 = complex_std,
                                           nobs2 = 1000,
                                           equal_var = False)

        print(cell_model, result, p)


ttest_perm(fantom_all)
ttest_perm(E123)
ttest_perm(E118)
ttest_perm(CL)
ttest_perm(UB)

'''

two-tailed permutation tests


all fantom enh
K562 2113.424347967078 0.0
HEPG2 1122.5971469691642 0.0

ttest_perm(E123) - K562 1892.273538458676 0.0

ttest_perm(E118) - HEPG2 1756.6936782193204 0.0

ttest_perm(CL) - K562 1413.6372582057693 0.0

ttest_perm(UB) - HEPG2 1951.1964120469022 0.0

'''
