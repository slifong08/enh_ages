# sarahfong
# 20210204

# obj: make a bed file out of the kircher19 saturating mutatgenesis MPRA assay

import pandas as pd
import subprocess

vals = {"Kpath":"/dors/capra_lab/data/mpra/kircher19/",
"sid": "GRCh37_ALL",
"AGEpath" : "/dors/capra_lab/projects/enhancer_ages/kircher19/"}

vals["Kfile"] = "%sGRCh37_ALL.tsv" % vals["Kpath"]

df = pd.read_csv(vals["Kfile"], sep ='\t')

# make a bed file!

# make new chr column

df["chr"] = "chr" + df.Chromosome

# make start column at position -1
df["start"] = df.Position - 1

# make end column at position
df["end"] = df.Position

outdf = df[["chr", "start", "end",
"Position", "Ref", "Alt", "Tags",
"DNA", "RNA", "Value", "P-Value",
"Element"]]

# write the bed file
outf = "%s%s.bed" %(vals["Kpath"], vals["sid"])

outdf.to_csv(outf, sep = '\t', header = True, index = False)

# move to project file

cmd = "cp %s %s" % (outf, vals["AGEpath"])

subprocess.call(cmd, shell = True)


#%% get enhancer coordinates
enh_end = outdf.groupby("Element")[["chr","end"]].max().reset_index()
enh_start= outdf.groupby("Element")[["chr","start"]].min().reset_index()
enh = pd.merge(enh_start, enh_end)
enh = enh[["chr", "start", "end", "Element"]]
enh["len"] = enh.end - enh.start

#%% write the bed file
enhf = "%s%s_elements.bed" %(vals["Kpath"], vals["sid"])

enh.to_csv(enhf, sep = '\t', header = False, index = False)

cmd = "cp %s %s" % (enhf, vals["AGEpath"])

subprocess.call(cmd, shell = True)
