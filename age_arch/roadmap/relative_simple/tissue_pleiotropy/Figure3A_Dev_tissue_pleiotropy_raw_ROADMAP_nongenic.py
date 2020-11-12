#!/usr/bin/env python
# coding: utf-8

#%%


import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
import statsmodels


colors = [ "amber", "faded green"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)
plt.rcParams.update({'font.size': 15})
sns.set_style("white")

import datetime
LAST_RUN = datetime.datetime.now()
TODAY = (datetime.date.today())
RE = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/pleiotropy/"

print("last run", LAST_RUN)


# # get enhancer file tissue/cell line descriptions


desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
desc_df= pd.read_csv(desc_file,header = None)
desc_df.columns = ["sid", "desc"]


# get fetal tissues

devsid_list = desc_df.loc[desc_df.desc.str.contains("Fetal"), "sid"].to_list()

#%%


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t')
syn_gen_bkgd["mrca"] = syn_gen_bkgd["mrca"].round(3)
syn_gen_bkgd["mrca_2"] = syn_gen_bkgd["mrca_2"].round(3)
syn_gen_bkgd["mya"] = syn_gen_bkgd["taxon2"].apply(lambda x: x.split(" ")[-1])

syn_gen_bkgd.head()

def get_relative_simple(sid, pdf):
    simple_break_percentile = 0.5
    if sid in pdf.sid2.unique():
        test = pdf.loc[(pdf.sid2 == sid)]

        relative_simple = pdf.loc[(pdf.sid2 == sid) &\
         (pdf.pct < simple_break_percentile), "seg_index"].max().astype(int)

    else:
        return "no_relative_simple",


#%%
for sid in devsid_list:

    multipath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/non-genic/"
    multifile = "%smultiintersect/trim-%s_multiintersect_0.5_count.bed"% (multipath,  sid) # all 98 tissues
    #multifile = "%smultiintersect/trim-%s_multiintersect_0.5_count_dev_only.bed"% (multipath,  sid) # 8 dev tissues

    # open the enhancer age file (architecture == raw architectures)
    agefile = "%sbreaks/no-exon_%s_parallel_breaks_enh_age_arch_summary_matrix.bed"% (multipath,  sid)
    if os.path.exists(agefile)==True:
        enh = pd.read_csv(agefile, sep = '\t', header = None, usecols=[3, 5, 7, 8, 10, 13, 14])


        enh.columns =[ "old_enh_id", "dataset", "arch", "seg_index", "old_len", "taxon2", "mrca_2"]

        # get relative percentiles and reassign "arch"
        enh_percentiles = "%sall_noexon_ROADMAP_breaks_percentiles.bed" % multipath
        pdf = pd.read_csv(enh_percentiles, sep = '\t')



        # RELATIVE SIMPLE DEF

        relative_simple = get_relative_simple(sid, pdf)
        print(sid, relative_simple)
        enh["relative_arch"] = "rel_simple"
        enh.loc[enh.seg_index.astype(float) >=relative_simple, "relative_arch"] = "rel_complex"
        enh.loc[enh.relative_arch == "rel_simple", "arch"] = "simple"


        # import pleiotropy multi-tissue overlap
        multi = pd.read_csv(multifile, sep ='\t', header = None)



        multi.columns = ["chr_enh", "start_enh", "end_enh", \
        "old_enh_id", "old_len","mean_len", "count_overlap"] # rename columns
        multi = pd.merge(enh, multi[["old_enh_id","old_len","count_overlap"]], how = "left", on = ["old_enh_id", "old_len"])
        multi = multi.loc[~multi.dataset.str.contains("shuf")]

        #multi = multi.loc[multi.chr_enh != "chrX"]


        # count some


        multi = multi.loc[multi.arch != ""]

        """
        sns.distplot(multi["count_overlap"].loc[multi["arch"].str.contains("simple")], label = "simple")
        sns.distplot(multi["count_overlap"].loc[multi["arch"].str.contains("complex")], label = "complex")
        plt.legend()
        """




        t, pval = stats.mannwhitneyu(multi["count_overlap"].loc[multi["arch"].str.contains("complex")],
                                           multi["count_overlap"].loc[multi["arch"].str.contains("simple")])
        print(t, pval)

        simple_med = multi["count_overlap"].loc[multi["arch"].str.contains("simple")].median()
        complex_med = multi["count_overlap"].loc[multi["arch"].str.contains("complex")].median()
        print( "simple, complex medians", simple_med, complex_med)




        multi.groupby("arch")["old_enh_id"].count()




        multi.groupby(["mrca_2","arch"])["count_overlap"].count()




        from matplotlib import gridspec
        from matplotlib.ticker import MultipleLocator

        order = ["simple", "complexenh"]
        order1 = ["simple", "complexenh"]

        fig = plt.figure(figsize = (12, 8))
        gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
        ax0 = plt.subplot(gs[0])

        sns.barplot(x = "arch", y = "count_overlap", data = multi,
                    palette = palette, order = order,
                    #showfliers = False, notch = True,
                    ax = ax0)
        #ax0.set_xlabel("")
        simplen = multi.loc[multi["arch"] == "simple"]["count_overlap"].count()
        simplemean = multi.loc[multi["arch"] == "simple"]["count_overlap"].mean().round(0)

        complexn = multi.loc[multi["arch"] != "simple"]["count_overlap"].count()
        complexmean = multi.loc[multi["arch"] != "simple"]["count_overlap"].mean().round(0)

        ax0.set(xticklabels =["simple n=%s\nmean = %s" %(simplen, simplemean),
         "complex n=%s\nmean = %s"% (complexn,complexmean)],
        ylim = (-1, 20), ylabel="Number of Tissue/Cell Datasets")

        ax0.set_xticklabels(ax0.get_xticklabels(), rotation = 90)

        ax0.yaxis.set_major_locator(MultipleLocator(4))
        sns.set("poster")

        ax2 = plt.subplot(gs[1])
        sns.barplot(x = "taxon2", y = "count_overlap", hue = "arch",
                      data = multi.sort_values(by = "mrca_2"),
                        palette = palette,
                    ax = ax2)

        ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90)
        desc = desc_df.loc[desc_df.sid == sid, "desc"].iloc[0] + "-" + sid
        ax2.set (xlabel = "", ylim = (-1, 20), ylabel = "",
        yticklabels = "", title = "%s" %desc )

        ax2.legend().remove()


        ax2.yaxis.set_major_locator(MultipleLocator(4))
        #ax2.set_xticklabels("")


        sns.set("poster")
        plt.savefig("%sFig3a-JOINT_barplot_fantom_sample_overlap_x_mrca_2_%s.pdf" % (RE,desc), bbox_inches = "tight" )


##%%
relative_simple
