#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pybedtools as pb
from scipy import stats
import seaborn as sns
order = ["simple", "complex"]
get_ipython().run_line_magic('matplotlib', 'inline')


# sns colors
arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
arch_palette = sns.xkcd_palette(arch_colors)


colors = [ "amber", "faded green", "dusty purple", "windows blue","greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(arch_palette)


# sns graphing preferences
sns.set(color_codes=True)
sns.set(font_scale=1.5)
sns.set_style("white")
sns.despine(bottom=True, left=True)


import datetime
last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())

RE = "/dors/capra_lab/projects/enhancer_ages/villar15/results/"


# # load TE, syngenbkgd

# In[2]:


te_path = "/dors/capra_lab/users/abraha1/projects/transposable_elements/data/"
tefam = "%shg19_TE_counts_wlin.txt" % te_path
dfam_file="%sdfam/hg38_dfam.nrph.hits_genomicCoordinates_tidy.bed" % te_path
repeatmasker = "%sfiltered_formated_hg19fromhg38-TE_coords.tsv" % te_path


# In[3]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df


# In[4]:


def custom_round(x, base=10):
    return int(base * round(float(x)/base))


def match_len(simple_df, complex_df, base_len):

    columns = ["enh_id", "enh_len"]
    columns_names = ["matching_ids", "matching_len"]

    simple = simple_df[columns].drop_duplicates()

    simple.columns = columns_names
    simple.matching_len = simple.matching_len.astype(float).apply(lambda x: custom_round(x, base=base_len)) # round to the nearest 100bp

    complex = complex_df[columns].drop_duplicates()
    complex.columns = columns_names

    complex.matching_len = complex.matching_len.astype(float).apply(lambda x: custom_round(x, base=base_len))

    lens = set(list(simple.matching_len.unique()) + list(complex.matching_len.unique()))

    match_dict = {}

    for length in lens:
        complexn = len(complex.loc[complex.matching_len == length])
        simplen = len(simple.loc[simple.matching_len == length])

        sn = min(complexn, simplen)

        if length > 0 and sn > 0:


            # find 100 matched enhancer samples

            complex_ids = complex.loc[complex.matching_len == length].sample(n = complexn, replace = False) # sample w/ replacement
            simple_ids = simple.loc[simple.matching_len == length].sample(n = simplen, replace = False) # sample w/ replacement
            balanced = pd.concat([simple_ids, complex_ids])
            match_dict[sn] = balanced

    final_matched_id = pd.concat(match_dict.values())

    return final_matched_id.matching_ids.unique()


def get_counts(df, strat):

    if strat == 1:
        counts = df.groupby(["mrca_2", "core_remodeling"])["enh_id"].count().reset_index()

        # add empty dataframe for complex human enhancers (do not exist by our def)
        empty = pd.DataFrame({"mrca_2": [0.000],
        "core_remodeling": [1],
        "enh_id": [0]
        })

        counts = pd.concat([counts, empty]) # concat the dataframe

        counts =counts.sort_values(by = ["core_remodeling", 'mrca_2']).reset_index()

    else:
        counts = df.groupby("arch")["enh_id"].count().reset_index()
        counts =counts.sort_values(by = "arch", ascending = False).reset_index()

    counts["enh_id"] = counts["enh_id"].astype(int)
    counts = counts.drop(["index"], axis = 1)

    return counts


# In[5]:


# 'VSPEC' = Villar species alignments. These fragments were mapped to human reference genome
villar_f = "/dors/capra_lab/projects/enhancer_uniqueness/data2/villar15/conservation/Hsap_Enhancers_conservationInOthers.bed"
vspec = pd.read_table("/dors/capra_lab/projects/enhancer_uniqueness/data2/villar15/conservation/Hsap_Enhancers_conservationInOthers.bed", sep = '\t')

hq_species =['Mmul', 'Cjac', 'Mmus', 'Rnor', 'Ocun', 'Btau', 'Sscr', 'Cfam', 'Fcat']

species = vspec.columns[4:22]

#Step 1- sum activity across species
# "-" = alignable, not active

vspec["all_map"] = (vspec[species] == '-').sum(axis =1)
vspec["hq_map"] = (vspec[hq_species] == '-').sum(axis =1)

vspec[species] = vspec[species].replace({'H3': 1}, regex = True) # replace all the active enhancers with '1' to sum
vspec["all_act"] = (vspec[species] == 1).sum(axis=1)
vspec["hq_act"] = (vspec[hq_species] ==1).sum(axis = 1)

#Step 2- sum alignable across species

vspec["all_map_act_sum"] = vspec["all_map"] + vspec["all_act"]
vspec["hq_map_act_sum"]= vspec["hq_map"] + vspec["hq_act"]

print(len(vspec), "villar enhancers")
vspec.head()


# In[6]:


# rename columns
vspec.columns=['chr_enh', 'start_enh', 'end_enh', 'IDs', 'Fcat',
 'Csab', 'Mdom', 'Ddel', 'Mbid', 'Cpor', 'Bbor', 'Rnor', 'Btau', 'Cjac',
 'Cfam', 'Mmus', 'Sscr', 'Ocun', 'Tbel', 'Shar', 'Mmul', 'Mfur', 'all_map',
 'hq_map', 'all_act', 'hq_act', 'all_map_act_sum', 'hq_map_act_sum']
vspec["enh_id"] = vspec["chr_enh"] +":"+ vspec["start_enh"].map(str) + "-" + vspec["end_enh"].map(str)
vspec["enh_len"] = vspec["end_enh"] - vspec["start_enh"]
vspec.head()


# In[7]:


# Fantom x TE intersection
BED_INTERSECTION = 0

sid = (villar_f.split("/")[-1]).split(".")[0]
print(sid)
outpath = "/home/fongsl/enhancer_ages/villar/data/"
outf = "%s%s_x_te.bed" % (outpath, sid)


# bedtools intersect
if BED_INTERSECTION == 1:

    f = pb.BedTool(villar_f).sort()
    r = pb.BedTool(repeatmasker).sort()

    f.intersect(r, wao = True).saveas(outf)
    cmd = '''tr -s '\t' <%s > %stemp.bed && mv %stemp.bed %s''' % (outf, outpath, outpath, outf) # remove double tabs
    os.system(cmd)


# # evaluate tissue overlap v. species overlap

# In[8]:


# roadmap tissue overlap
f_multi = "/dors/capra_lab/projects/enhancer_ages/villar15/data/multiintersect/VILLAR_multiintersect_0.5_count.bed"
multi = pd.read_csv(f_multi, sep = '\t', header = None)
multi.columns = ["chr_enh", "start_enh", "end_enh", "mrca", "enh_id", "core_remodeling","count_tissue_overlap"]

multi.head()


# In[9]:


# smaller villar df
vdf = vspec[['IDs', 'all_map',
 'hq_map', 'all_act', 'hq_act', 'all_map_act_sum', 'hq_map_act_sum', "enh_id", "enh_len"]].drop_duplicates()

len(vdf)


# In[10]:


merged = pd.merge(vdf,multi)
#merged = merged.loc[~merged.count_tissue_overlap.isna()]
merged.head()


# # match on length

# In[74]:


# add in max seg index information
breaks = "/dors/capra_lab/projects/enhancer_ages/villar15/data/breaks/HSap_age_seg_count.bed"
b = pd.read_csv(breaks, header = None)

b["seg_index"] = b[0].apply(lambda x: x.split("chr")[0])
b["enh_id"] = "chr" + b[0].apply(lambda x: x.split("chr")[1])
b.head()


# In[75]:


# set relative simple definition
merged = pd.merge(merged, b[["enh_id", "seg_index"]], how = "left", on = "enh_id")

# redefine simple architecture according to the median age breaks.
merged.loc[merged.seg_index.astype(int) < merged.seg_index.astype(int).median(), "core_remodeling"] = 0

merged.groupby("core_remodeling")["enh_id"].count()


# In[76]:


simpledf = merged.loc[merged.core_remodeling ==0]

complex_df = merged.loc[(merged.core_remodeling ==1)]

matched_id = match_len(simpledf, complex_df, 100)


matched = merged.loc[merged.enh_id.isin(matched_id)]

matched["arch"] = "simple"
matched.loc[matched.core_remodeling ==1, "arch"] = "complexenh"
matched.mrca = matched.mrca.round(3)
matched = pd.merge(matched, syn_gen_bkgd[["mrca", "mrca_2", "taxon2"]], how = "left", on ="mrca")
matched = matched.loc[matched.chr_enh != "chrX"]
matched.head()




matched.groupby("core_remodeling")["enh_id"].count()


# In[79]:



mwu, p = stats.mannwhitneyu(matched.loc[matched.core_remodeling == 1, "hq_act"],\
matched.loc[matched.core_remodeling != 1, "hq_act"])
print(mwu, p)


#%%
matched.groupby("core_remodeling")["hq_act"].mean()


order = ["simple", "complexenh"]

matched.head()

#%%

from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator


fig = plt.figure(figsize = (12, 8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
ax0 = plt.subplot(gs[0])

splot = sns.barplot(x = "arch", y = "hq_act", data = matched,
            palette = palette, order = order,
            ax = ax0)
STRAT = 0
counts = get_counts(matched, STRAT)
for n, p in enumerate(splot.patches):
    value = counts.iloc[n]["enh_id"]
    splot.annotate(value,
                   (p.get_x() + p.get_width() / 2.,0.10),
                   ha = 'center', va = 'baseline',
                   size=15,
                   rotation = 90,
                   color = "white",
                   xytext = (0, 1),
                   textcoords = 'offset points'
                   )

ax0.set_xticklabels(["simple", "complex"], rotation = 90)
ax0.set(xlabel="", ylabel ="Number of Active Species", ylim=(0,3))
ax0.yaxis.set_major_locator(MultipleLocator(1))

sns.set("poster")

ax2 = plt.subplot(gs[1])
mplot = sns.barplot(x = "taxon2", y = "hq_act", hue = "core_remodeling",
              data = matched.sort_values(by = "mrca_2"),
                palette = palette,
            ax = ax2)
STRAT = 1
counts = get_counts(matched, STRAT)
counts
for n, p in enumerate(mplot.patches):
    value = counts.iloc[n]["enh_id"].astype(int)
    mplot.annotate(value,
                   (p.get_x() + p.get_width() / 2.,0.05),
                   ha = 'center', va = 'baseline',
                   size=15,
                   rotation = 90,
                   color = "white",
                   xytext = (0, 1),
                   textcoords = 'offset points'
                   )
xlab = ['Homo', 'Prim', 'Euar', 'Bore', 'Euth', 'Ther', 'Mam',
'Amni', 'Tetr', 'Vert']
ax2.set_xticklabels(xlab, rotation = 90)
sns.set("poster")
ax2.yaxis.set_major_locator(MultipleLocator(1))
ax2.set(ylabel="",  ylim=(0,3), yticklabels = "")
ax2.legend().remove()
plt.savefig("%sFig3b-JOINT_barplot_villar_cross_species_overlap_x_mrca_2.pdf" % RE, bbox_inches = "tight" )

#%%


# In[ ]:


te = pd.read_csv(outf, sep = "\t", header = None)

te = te [[0,1,2,3, 22,23,24,25,26]] # keep only important info
te.head()


# In[ ]:


te.columns = ['chr_enh', 'start_enh', 'end_enh', 'IDs', "chr_t", "start_t", "end_t", "te", "len_te"]
te.head()


# In[ ]:


matchedcomplex_te = pd.merge(matchedcomplex, te[["IDs", "chr_t", "start_t", "end_t", "te", "len_te"]], how = "left", on="IDs")
matchedcomplex_te["te_bin"] = 0
matchedcomplex_te.loc[matchedcomplex_te.start_t>0, "te_bin"] = 1
matchedcomplex_te = matchedcomplex_te.drop_duplicates()

simpledf_te = pd.merge(simpledf, te[["IDs", "chr_t", "start_t", "end_t", "te", "len_te"]], how = "left", on="IDs").drop_duplicates()
simpledf_te["te_bin"] = 0
simpledf_te.loc[simpledf_te.start_t>0, "te_bin"] = 1
simpledf_te = simpledf_te.drop_duplicates()


# In[ ]:


matchedcomplex_te = matchedcomplex_te.groupby(["enh_id", "IDs", "enh_len"])["te_bin", "mrca", "hq_act"].max().reset_index()
simpledf_te = simpledf_te.groupby(["enh_id", "IDs", "enh_len"])["te_bin", "mrca", "hq_act"].max().reset_index()


# In[ ]:


print(simpledf_te.enh_len.median())
print(matchedsimple.enh_len.median())

matchedcomplex_te.enh_len.median()


# In[ ]:


sns.distplot(simpledf.enh_len, kde = False, norm_hist = True, label = "simple")
sns.distplot(matchedcomplex.enh_len, kde = False, norm_hist = True, label = "complex")
plt.legend()


# In[ ]:


pleiotropy = sns.jointplot( y= "count_tissue_overlap", x ="hq_map_act_sum", data = merged)
pleiotropy.annotate(stats.pearsonr, fontsize = 10)


# In[ ]:


fig = plt.figure(figsize = (16,6))
fig.subplots_adjust(hspace=0.4, wspace=0.4)
ax1 = fig.add_subplot(1, 2, 1)
s = sns.jointplot( y= "count_tissue_overlap", x ="hq_map_act_sum",
                  data = merged.loc[merged["core_remodeling"] == 1], color = "g", kind = "reg", ax = ax1)
ax2 = fig.add_subplot(1, 2, 2)
c = sns.jointplot( y= "count_tissue_overlap", x ="hq_map_act_sum",
                  data = merged.loc[merged["core_remodeling"] == 0], color = "y", kind = "reg", ax = ax2)
s.annotate(stats.pearsonr, fontsize = 10)
c.annotate(stats.pearsonr, fontsize = 10)


# # Active enhancers across species

# In[ ]:


# calculate the frequency of enhancers mapping to species

# simple

simple_enh_freq = (simpledf.groupby("hq_act")["enh_id"].count().reset_index())
simple_enh_freq.columns = ["hq_act", "count"]
simple_enh_freq["freq"] = simple_enh_freq["count"]/(simple_enh_freq["count"].sum())
simple_enh_freq["code"] = "simple"

# complex
complex_enh_freq = (matchedcomplex.groupby("hq_act")["enh_id"].count().reset_index())
complex_enh_freq.columns = ["hq_act", "count"]
complex_enh_freq["freq"] = complex_enh_freq["count"]/(complex_enh_freq["count"].sum())
complex_enh_freq["code"] = "complexenh"

# merge the dataframes
enh_freq = pd.concat([simple_enh_freq,complex_enh_freq])

colors = [ "amber", "faded green"]
palette = sns.xkcd_palette(colors)
sns.palplot(arch_palette)
fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (16,6))


hue_order = ["simple", "complexenh"]
sns.barplot(x = "hq_act", y = "freq", data = enh_freq, hue = "code",  ax = ax1,
            hue_order = hue_order, palette = palette)
ax1.set_ylabel("Freq")
ax1.set_xlabel("# overlapping species")
ax1.set_title("Frequency of architectures ACTIVE in other species")

sns.barplot(x = "hq_act", y = "count", data = enh_freq, hue = "code",  ax = ax2,
            hue_order = hue_order, palette = palette)
ax2.set_ylabel("Count")
ax2.set_xlabel("# overlapping species")
ax2.set_title("Count of architectures ACTIVE in other species")
m, p = stats.mannwhitneyu(simpledf["hq_act"],
                   matchedcomplex["hq_act"])
print(m,p)
plt.tight_layout()

#plt.savefig("%sVILLAR_enharch_freq_species_activity_bar_plot.pdf" %RE, bbox_inches = "tight")


# In[ ]:


mrcas = pd.concat([matchedsimple, matchedcomplex])
mrcas.mrca = mrcas.mrca.round(3)
mrcas["code"] = "simple"
mrcas["code"].loc[mrcas.core_remodeling ==1] = "complexenh"
mrcas.head()


# In[ ]:


mrcas = pd.merge(mrcas, syn_gen_bkgd, how = "left" )


# In[ ]:


mrcas.groupby(["core_remodeling"])["hq_act"].median()


# In[ ]:


mrcas.groupby(["core_remodeling"])["IDs"].count()


# In[ ]:


counts = mrcas.groupby(["mrca_2", "core_remodeling", "taxon2"])["IDs"].count().reset_index()

counts["IDs"] = counts["IDs"].astype(str) # turn counts into a str

countLabel = counts.groupby(["taxon2", "mrca_2"])["IDs"].apply(','.join).reset_index()
countLabel["label"] = countLabel.taxon2 +"\n(" + countLabel.IDs +")"

labels = countLabel.sort_values(by = "mrca_2")["label"].tolist()
labels


# In[ ]:


mrcas.head()


# In[ ]:


print(len(mrcas))
#mrcas = mrcas.loc[mrcas.mrca_2>0.131]
print(len(mrcas))


# In[ ]:


mrcas.code.unique()


# In[ ]:


mrcas.groupby("core_remodeling")["mrca"].count()


# In[ ]:


from matplotlib import gridspec

order = ["simple", "complex"]
order1 = ["simple", "complex"]

fig = plt.figure(figsize = (12, 8))
sns.set("poster")
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
ax0 = plt.subplot(gs[0])

hue_order = ["simple", "complexenh"]
sns.barplot(x = "core_remodeling", y = "hq_act", data = mrcas.sort_values(by = "mrca_2"),
             ax = ax0,
            hue_order = hue_order,
            palette = palette,)
ax0.set_xlabel("")
ax0.set_xticklabels("")
ax0.set_ylabel("Number of Active Species")
ax0.set_ylim(-0.5,5.5)

sns.set("poster")

ax2 = plt.subplot(gs[1])
sns.barplot(x = "taxon2", y = "hq_act", hue = "core_remodeling",
              data = mrcas.sort_values(by = "mrca_2"),
                palette = palette, ax = ax2)

ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 270)
ax2.set_xlabel("")

ax2.legend().remove()
sns.set("poster")
ax2.set_ylim(-0.5, 5.5)

#ax2.set_xticklabels("")
ax2.set_ylabel("")

#plt.savefig("%sfig3b-JOINT VILLAR_species_activity_arch_mrca_olderthaneuar.pdf" % RE, bbox_inches = "tight")


# In[ ]:


stats.mannwhitneyu(mrcas.loc[mrcas.code == "simple", "hq_act"],
                  mrcas.loc[mrcas.code != "simple", "hq_act"])


# In[ ]:


mrcas.groupby("core_remodeling")["hq_act"].mean()


# In[ ]:


mrcas.groupby(["mrca_2","core_remodeling"])["enh_id"].count()


# In[ ]:


arch_te = [ "amber", "lemon", "faded green", "pale green", ]
arch_te_pal = sns.xkcd_palette(arch_te)
sns.palplot(arch_te_pal)


# In[ ]:


# calculate the frequency of enhancers mapping to species

# simple

simple_enh_freq = (simpledf_te.groupby(["te_bin", "hq_act", ])["enh_id"].count().reset_index())
simple_enh_freq.columns = ["te_bin", "hq_act", "count"]
simple_enh_freq["freq"] = simple_enh_freq["count"]/(simple_enh_freq["count"].sum())
simple_enh_freq["code"] = "simple"
simple_enh_freq["code2"] = simple_enh_freq["code"] +"-TE_overlap_"+ simple_enh_freq["te_bin"].map(str)

# complex
complex_enh_freq = (matchedcomplex_te.groupby(["te_bin", "hq_act", ])["enh_id"].count().reset_index())
complex_enh_freq.columns = ["te_bin", "hq_act", "count"]
complex_enh_freq["freq"] = complex_enh_freq["count"]/(complex_enh_freq["count"].sum())
complex_enh_freq["code"] = "complexenh"
complex_enh_freq["code2"] = complex_enh_freq["code"] +"-TE_overlap_"+ complex_enh_freq["te_bin"].map(str)
# merge the dataframes

enh_freq = pd.concat([simple_enh_freq,complex_enh_freq])

fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (16,6))


#hue_order = ["simple", "complexenh"]
sns.barplot(x = "hq_act", y = "freq", data = enh_freq, hue = "code2",  ax = ax1, palette = arch_te_pal)
#            hue_order = hue_order, palette = palette)
ax1.set_ylabel("Freq")
ax1.set_xlabel("# overlapping species")
ax1.set_title("Frequency of architectures ACTIVE in other species")
ax1.legend(loc = "upper right")
sns.barplot(x = "hq_act", y = "count", data = enh_freq, hue = "code2",  ax = ax2, palette = arch_te_pal)
#            hue_order = hue_order, palette = palette)
ax2.set_ylabel("Count")
ax2.set_xlabel("# overlapping species")
ax2.set_title("Count of architectures ACTIVE in other species")
ax2.legend(loc = "upper right")
m, p = stats.mannwhitneyu(simpledf["hq_act"],
                   matchedcomplex["hq_act"])
print(m,p)
plt.tight_layout()

plt.savefig("/home/fongsl/enhancer_ages/villar/VILLAR_enharch_freq_species_activity_te.pdf", bbox_inches = "tight")


# In[ ]:


within_te = enh_freq.groupby(["code", "te_bin"])["freq", "count"].sum().reset_index()
within_te.columns = ["code", "te_bin", "withinarch_te_freq", "withinarch_te_count"]

within_te.head()


# In[ ]:


enh_freq = pd.merge(enh_freq, within_te, how = "left")
enh_freq["withinarch_te_freq2"] = enh_freq["count"].divide(enh_freq.withinarch_te_count)
enh_freq.head()


# In[ ]:


fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (16,6))

sns.barplot(x = "hq_act", y = "withinarch_te_freq2", data = enh_freq, hue = "code2",  ax = ax1, palette = arch_te_pal)

ax1.set_ylabel("Freq")
ax1.set_xlabel("# overlapping species")
ax1.set_title("Frequency of architectures ACTIVE in other species")
ax1.legend(loc = "upper right")
sns.barplot(x = "hq_act", y = "count", data = enh_freq, hue = "code2",  ax = ax2, palette = arch_te_pal)

ax2.set_ylabel("Count")
ax2.set_xlabel("# overlapping species")
ax2.set_title("Count of architectures ACTIVE in other species")
ax2.legend(loc = "upper right")
m, p = stats.kruskal(enh_freq.loc[enh_freq["code2"] == "simple-TE_overlap_0", "withinarch_te_freq2"],
                   enh_freq.loc[enh_freq["code2"] == "simple-TE_overlap_1", "withinarch_te_freq2"],
                    enh_freq.loc[enh_freq["code2"] == "complexenh-TE_overlap_0", "withinarch_te_freq2"],
                    enh_freq.loc[enh_freq["code2"] == "complexenh-TE_overlap_1", "withinarch_te_freq2"])
print(m,p)
plt.tight_layout()

plt.savefig("/home/fongsl/enhancer_ages/villar/VILLAR_enharch_freq2_species_activity_te.pdf",
            bbox_inches = "tight")


# #
