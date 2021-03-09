import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/age_breaks/"

es_colors = [ "slate grey","greyish"]
es_pal = sns.xkcd_palette(es_colors)
sns.palplot(es_pal)

colors = [ "amber", "faded green", "dusty purple", "windows blue","greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)

#%% Files
path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"

enh = "%sFANTOM_enh_age_arch_full_matrix.tsv" % path
summaryEnh = "%sFANTOM_enh_age_arch_summary_matrix.tsv" % path

shuf = "%sSHUFFLED_FANTOM_enh_age_arch_full_matrix.tsv" % path
summaryShuf = "%sSHUFFLE_FANTOM_enh_age_arch_summary_matrix.tsv" % path

#%% other summary files

# age and taxon file
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df
syn_gen_bkgd

# tissue/cell line dataset descriptions
desc_file = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/sample_id_descriptions.txt"
desc_df= pd.read_csv(desc_file, sep = '\t', header = None)

#%% LOAD Files
enh = "%sFANTOM_enh_age_arch_full_matrix.tsv" % path
summaryEnh = "%sFANTOM_enh_age_arch_summary_matrix.tsv" % path

shuf = "%sSHUFFLED_FANTOM_enh_age_arch_full_matrix.tsv" % path
summaryShuf = "%sSHUFFLE_FANTOM_enh_age_arch_summary_matrix.tsv" % path

shuffle = pd.read_csv(shuf, sep = '\t')
final_merge = pd.read_csv(enh, sep = '\t')


#%%


enh_lens = final_merge.groupby(["enh_id", "core_remodeling","enh_len", "arch"])["mrca_2"].max().reset_index()
enh_lens.mrca_2 = enh_lens.mrca_2.round(3)

enh_lens = pd.merge(enh_lens, syn_gen_bkgd[[ "mrca_2","taxon2", "mya2"]], how = "left", on = "mrca_2")
enh_lens = enh_lens.drop_duplicates()
print(enh_lens.shape)

enh_lens["datatype"]="FANTOM"
enh_lens.head()
#%%

shuf_len = shuffle.groupby(["enh_id", "core_remodeling",\
 "enh_len", "shuf_id", "arch"])["mrca_2"].max().reset_index()

print(shuf_len.shape)

shuf_len["datatype"]="SHUFFLE"
shuf_len.mrca_2 = shuf_len.mrca_2.round(3)
shuf_len = pd.merge(shuf_len, syn_gen_bkgd[["mrca_2", "taxon2", "mya2"]],\
 how = "left", on = "mrca_2")

shuf_len = shuf_len.drop_duplicates()
print(shuf_len.shape)


#%% sample 1/10th of shuffle dataset and combine w/ enhancers for comparison


sample_shuf = shuf_len.sample(frac = 0.1)

lens = pd.concat([enh_lens, sample_shuf])
lens.head()

#%% simple v complex

fig, ax = plt.subplots(figsize = (8, 8))
sns.boxplot(x = "core_remodeling", y = "enh_len",
            data = enh_lens, palette = palette,
           showfliers = False)
ax.set(xticklabels = ["Simple", "Complex"], xlabel = "Architecture",
ylabel= "Enhancer Length (bp)")
ax.legend(bbox_to_anchor=(1.0, 1.0))


# Get medians


lens_data = lens.groupby(["core_remodeling", "datatype"])["enh_len"].median().reset_index()
lens_data_mean = lens.groupby(["core_remodeling", "datatype"])["enh_len"].mean().round().reset_index()

print("medians", lens_data,)
print("means", lens_data_mean,)
ms, msp = stats.mannwhitneyu(enh_lens.loc[enh_lens.core_remodeling==1, "enh_len"],
                            enh_lens.loc[enh_lens.core_remodeling==0, "enh_len"])
print("simple", ms, msp)

plt.savefig("%sfantom_enh_arch_lens.pdf" %RE, bbox_inches = "tight")


""" RESULTS MWU SIMPLE V. COMPLEX ENHANCER LENGTHS
simple 71739682.5 0.0
"""


""" RESULTS
No handles with labels found to put in legend.
medians    core_remodeling datatype  enh_len
0              0.0   FANTOM      259
1              0.0  SHUFFLE      255
2              1.0   FANTOM      347
3              1.0  SHUFFLE      339
means    core_remodeling datatype  enh_len
0              0.0   FANTOM    277.0
1              0.0  SHUFFLE    272.0
2              1.0   FANTOM    370.0
3              1.0  SHUFFLE    363.0
simple 71739682.5 0.0
"""

#%% Plot FANTOM v. SHUFFLE simple and complex lengths


fig, ax = plt.subplots(figsize = (8, 8))
sns.boxplot(x = "core_remodeling", y = "enh_len", hue = "datatype",
            data = lens, palette = es_pal,
           showfliers = False, notch = True)
ax.set(xticklabels = ["Simple", "Complex"], xlabel = "Architecture",
ylabel= "Enhancer Length (bp)")
ax.legend(bbox_to_anchor=(1.0, 1.0))

# Get medians

lens_data = lens.groupby(["core_remodeling", "datatype"])["enh_len"].median().reset_index()
lens_data_mean = lens.groupby(["core_remodeling", "datatype"])["enh_len"].mean().reset_index()
print("medians", lens_data,)
print("means", lens_data_mean,)

ax.get_legend().remove()
plt.savefig("%sfantom_enh_lens.pdf" %RE, bbox_inches = "tight")


""" RESULTS
medians    core_remodeling datatype  enh_len
0              0.0   FANTOM      259
1              0.0  SHUFFLE      255
2              1.0   FANTOM      347
3              1.0  SHUFFLE      339
means    core_remodeling datatype     enh_len
0              0.0   FANTOM  276.759779
1              0.0  SHUFFLE  271.579385
2              1.0   FANTOM  369.700593
3              1.0  SHUFFLE  363.292826
"""

#%% Old (>placental) v. young (<= placental) fantom median lengths


print(stats.mannwhitneyu(lens.loc[(lens.mrca_2>0.175) & (lens.datatype == "FANTOM"),"enh_len"],
                  lens.loc[(lens.mrca_2<=0.175) & (lens.datatype == "FANTOM"),"enh_len"],))

print("\nolder than mammalian lengths\n")
print(lens.loc[lens.mrca_2>0.175].groupby("datatype")["enh_len"].median())
print("\nmammalian and younger lengths\n")
lens.loc[lens.mrca_2<=0.175].groupby("datatype")["enh_len"].median()


""" RESULTS
MannwhitneyuResult(statistic=88089072.0, pvalue=8.184782752230459e-122)

older than mammalian lengths

datatype
FANTOM     321
SHUFFLE    311
Name: enh_len, dtype: int64

mammalian and younger lengths

datatype
FANTOM     277
SHUFFLE    286
Name: enh_len, dtype: int64
"""

#%% Young fantom v. shuffle lengths


print("\nYoung fantom v. shuffle lengths\n")
print(stats.mannwhitneyu(lens.loc[(lens.mrca_2<=0.175) & (lens.datatype == "FANTOM"),"enh_len"],
                  lens.loc[(lens.mrca_2<=0.175) & (lens.datatype == "SHUFFLE"),"enh_len"],))

""" RESULTS Young fantom v. shuffle lengths

MannwhitneyuResult(statistic=2021614693.5, pvalue=2.7548239106496657e-16)

"""

#%% FANTOM enhancer v shuffled lengths. 


hue_order = ["FANTOM", "SHUFFLE"]
fig, ax = plt.subplots(figsize = (8, 8))
sns.barplot(x = "taxon2", y = "enh_len", hue = "datatype",
            data = lens.sort_values(by = 'mrca_2'),\
            palette = es_pal, hue_order = hue_order)

ax.set(ylim = (225,375), xlabel= "Architecture", ylabel= "Enhancer Length")
ax.set_xticklabels(ax.get_xticklabels(),\
 rotation = 90, horizontalalignment = "left")
ax.legend(bbox_to_anchor=(1.0, 1.0))
# Get medians

lens_data = lens.groupby(["core_remodeling", "datatype"])\
["enh_len"].median().reset_index()
lens_data_mean = lens.groupby(["core_remodeling", "datatype"])\
["enh_len"].mean().reset_index()

print("medians", lens_data,)
print("means", lens_data_mean,)

ax.get_legend().remove()
plt.savefig("%sfig1.2b-fantom_enh_mrca_lens.pdf" %RE, bbox_inches = "tight")


""" RESULTS MEDIAN AND MEAN FANTOM/SHUFFLE LENGTHS
medians    core_remodeling datatype  enh_len
0              0.0   FANTOM      259
1              0.0  SHUFFLE      255
2              1.0   FANTOM      347
3              1.0  SHUFFLE      339
means    core_remodeling datatype     enh_len
0              0.0   FANTOM  276.759779
1              0.0  SHUFFLE  271.579385
2              1.0   FANTOM  369.700593
3              1.0  SHUFFLE  363.292826
"""


#%% MWU FANTOM v. SHUFFLE simple enhancer lengths


stats.mannwhitneyu(lens.loc[(lens.arch == "simple") & (lens.datatype == "FANTOM"), "enh_len"],
                  lens.loc[(lens.arch == "simple") & (lens.datatype == "SHUFFLE"), "enh_len"])
""" RESULTS MWU FANTOM v. SHUFFLE simple enhancer lengths
MannwhitneyuResult(statistic=1526143108.5, pvalue=8.179141620886248e-05)
"""
#%% MRCA linear regression

xs = enh_lens.mrca_2.loc[enh_lens.arch.str.contains("simple")]
ys = enh_lens.enh_len.loc[enh_lens.arch.str.contains("simple")]

slope_s, intercept_s, r_value_s, p_value_s, std_err_s = stats.linregress(xs,ys)
f
xc = enh_lens.mrca_2.loc[enh_lens.arch.str.contains("complex")]
yc = enh_lens.enh_len.loc[enh_lens.arch.str.contains("complex")]

slope_c, intercept_c, r_value_c, p_value_c, std_err_c  = stats.linregress(xc,yc)
slope_c, intercept_c, r_value_c, p_value_c, std_err_c

xsshuf = shuf_len.mrca_2.loc[shuf_len.arch.str.contains("simple")]
ysshuf = shuf_len.enh_len.loc[shuf_len.arch.str.contains("simple")]

slope_ss, intercept_ss, r_value_ss, p_value_ss, std_err_ss = stats.linregress(xsshuf,ysshuf)

xcs = shuf_len.mrca_2.loc[shuf_len.arch.str.contains("complex")]
ycs = shuf_len.enh_len.loc[shuf_len.arch.str.contains("complex")]

slope_cs, intercept_cs, r_value_cs, p_value_cs, std_err_cs = stats.linregress(xcs,ycs)

print("Simple enhancer slope = %s,\
 r2 = %s, pval = %s\nComplex enhancer slope = %s, r2 = %s, pval = %s"\
              % (round(slope_s, 0), round(r_value_s,2), p_value_s,\
               round(slope_c,0),  round(r_value_c,2), p_value_c))
print("SHUFFLE Simple enhancer slope = %s, r2 = %s, \
pval = %s\nSHUFFLE Complex enhancer slope = %s, r2 = %s, pval = %s"\
              % (round(slope_ss, 0), round(r_value_ss,2), p_value_ss,\
               round(slope_cs,0),  round(r_value_cs,2), p_value_cs))

""" RESULTS MRCA LINEAR REGRESSION
Simple enhancer slope = -12.0, r2 = -0.01, pval = 0.0965096502617304
Complex enhancer slope = 81.0, r2 = 0.1, pval = 5.238360778796791e-24
SHUFFLE Simple enhancer slope = -24.0, r2 = -0.03, pval = 6.836831404608524e-296
SHUFFLE Complex enhancer slope = 39.0, r2 = 0.05, pval = 0.0
"""

#%% MYA LINEAR REGRESSION
xs = enh_lens.mya2.loc[enh_lens.arch.str.contains("simple")]
ys = enh_lens.enh_len.loc[enh_lens.arch.str.contains("simple")]

slope_s, intercept_s, r_value_s, p_value_s, std_err_s = stats.linregress(xs,ys)

xc = enh_lens.mya2.loc[enh_lens.arch.str.contains("complex")]
yc = enh_lens.enh_len.loc[enh_lens.arch.str.contains("complex")]

slope_c, intercept_c, r_value_c, p_value_c, std_err_c  = stats.linregress(xc,yc)
slope_c, intercept_c, r_value_c, p_value_c, std_err_c

xsshuf = shuf_len.mya2.loc[shuf_len.arch.str.contains("simple")]
ysshuf = shuf_len.enh_len.loc[shuf_len.arch.str.contains("simple")]

slope_ss, intercept_ss, r_value_ss, p_value_ss, std_err_ss = stats.linregress(xsshuf,ysshuf)

xcs = shuf_len.mya2.loc[shuf_len.arch.str.contains("complex")]
ycs = shuf_len.enh_len.loc[shuf_len.arch.str.contains("complex")]

slope_cs, intercept_cs, r_value_cs, p_value_cs, std_err_cs = stats.linregress(xcs,ycs)

print("Simple enhancer slope = %s, r2 = %s, pval = %s\
\nComplex enhancer slope = %s, r2 = %s, pval = %s"\
              % (round(slope_s, 3), round(r_value_s,2), p_value_s,\
               round(slope_c,3),  round(r_value_c,2), p_value_c))
print("SHUFFLE Simple enhancer slope = %s, r2 = %s, pval = %s\
\nSHUFFLE Complex enhancer slope = %s, r2 = %s, pval = %s"\
              % (round(slope_ss, 3), round(r_value_ss,2), p_value_ss, \
              round(slope_cs,3),  round(r_value_cs,2), p_value_cs))
""" RESULTS MYA LINEAR REGRESSION
Simple enhancer slope = -0.008, r2 = -0.0, pval = 0.4899561844561535
Complex enhancer slope = 0.107, r2 = 0.08, pval = 6.303958714096115e-18
SHUFFLE Simple enhancer slope = -0.031, r2 = -0.02, pval = 7.439211772479368e-206
SHUFFLE Complex enhancer slope = 0.059, r2 = 0.05, pval = 0.0
"""
