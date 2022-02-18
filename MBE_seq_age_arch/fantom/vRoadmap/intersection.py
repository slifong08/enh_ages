import os, sys
import pandas as pd
import subprocess

fPath = "/dors/capra_lab/projects/enhancer_ages/fantom/data/breaks/"
fF= "%sall_unique_fantom_erna_112_tissue.bed"% fPath

rPath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E072/cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_E072/breaks/"
rF = "%sno-exon_E072.bed" % rPath

oP = "/home/fongsl/fantom_roadmap/"
oF = "%sfantom_x_E072.bed" % oP

cmd = "bedtools intersect -a %s -b %s -wa -wb > %s" % (fF, rF, oF)
subprocess.call(cmd, shell = True)

print(cmd)

cutF = "%scut.bed" % oP
cut = "cut -f 4,5,6,7,8,22,23,24,25,26 %s > %s" %(oF, cutF)
subprocess.call(cut, shell = True)

#%%
df = pd.read_csv(cutF, sep ='\t', header = None)
df.head()
cols = ["fantom_enh_id", "id2", "fantom_id", "arch_binF", "archF", "E003_enh_id", "roadmap_id", "roadmap_seg", "arch_binR", "archR"]
df.columns = cols
df["consistent_arch"] = 0
df.loc[df.arch_binF == df.arch_binR, "consistent_arch"] =1

df.groupby(df["consistent_arch"])["id2"].count()
