#!/usr/bin/env python
# coding: utf-8

#%% In[1]:


#%% Intro
# created 2019-08-05

# sarahfong

# Goal - multiintersect all roadmap histone enhancers with a count for the number of overlapping coordinates between datasets

# Method - multiintersect

## step 1 - make bed file with all non-redundant enhancers, named "all_fantom_enhancers". This will be the -a arguement
## step 2 - bedtools intersect -a [all_roadmap_enhancers] -b [E005_core_ages.bed, E004_core_ages.bed ... etc.]

## -c (count number of overlapping features)

import glob
import os, sys
import pandas as pd

fraction_overlap = sys.argv[1]
fraction_overlap = 0.5


#%% In[2]:


path = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/"
outpath = "/dors/capra_lab/projects/enhancer_ages/fantom/data/multiintersect/"
datapath = "/dors/capra_lab/projects/enhancer_ages/fantom/data/architecture_coordinates/"

# -b base

# Entire Enhancer #
# Trim enhancers to mean lengths #
infiles = glob.glob("%s*_complexenh.bed" % datapath)
#inpath = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"
#infiles =["%sall_fantom_enh_112_tissues.bed" % inpath]


# -a
all_uniq = "%sall_unique_fantom_erna_112_tissue.bed" % outpath # 31089 unique enhancer IDs

# -a columns = chr start end enh_id maxmrca maxbreaks arch_type

# -b
#path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/multiintersect/trimmed/"
all_beds = glob.glob("%s**/*.bed" % path, recursive = True)
#all_beds = glob.glob("%s*.bed" % path)
all_beds_str = ' '.join(str(bed) for bed in all_beds)



all_beds_str


#%% In[26]: TRIM the enhancers to mean length per enhancer dataset.


for infile in infiles:
    sid = "".join((infile.split("/")[-1]).split("_")[:1])

    if sid != "FANTOM":

        df = pd.read_csv(infile, sep ='\t', header = -1)

        df.columns = ["chr", "start", "end", "mrca", "enh_id",  "core_remodeling"]

        df["len"]= df.end - df.start # calculate enhancer length

        df["midpoint"] = (df.start + (df.len)/2).astype(int) #%% Identify the midpoint of each enhancer

        mean = int(df["len"].mean()) # mean enhancer length in dataset
        print(sid, mean)
        df["start_new"] =((df.midpoint - (mean/2)).round(0)).astype(int) # calculate new start as the midpoint - (mean length/2)

        df["end_new"] = ((df.midpoint + (mean/2)).round(0)).astype(int)

        df["desc"] = sid

        trimmed = df[["chr", "start_new", "end_new", "len", "core_remodeling", "mrca", "desc"]].drop_duplicates()
         # -a
        test_file = "%strimmed/trimmed_%s_FANTOM_complexenh.bed" % (outpath, sid) # trimmed file for multiintersect

        trimmed.to_csv(test_file, header = False, index = False, sep = '\t') # save the trimmed file.

        # -b
        all_beds = glob.glob("%s**/*.bed" % path, recursive = True)
        all_beds_str = ' '.join(str(bed) for bed in all_beds)

        # recursive one tissue v. all tissue multi intersect
        multi_cmd = "bedtools intersect -a %s -b %s -f %s -c > %strimmed_%s_multiintersect_%s_count.bed" % (test_file, all_beds_str, fraction_overlap,  outpath, sid, fraction_overlap)
        print(multi_cmd)
        os.system(multi_cmd)


# # Multi intersect all_fantom_enh_112_tissues.bed

#%% In[6]:


infile = "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh_112_tissues.bed"
sid = (infile.split("/")[-1]).split(".")[0]


df = pd.read_csv(infile, sep ='\t', header = -1)


df.columns = ["chr", "start", "end", "seg_index",  "core_remodeling",  "mrca",]

df["len"]= df.end - df.start # calculate enhancer length

df["midpoint"] = (df.start + (df.len)/2).astype(int) #%% Identify the midpoint of each enhancer

mean = int(df["len"].mean()) # mean enhancer length in dataset
print(sid, mean)
df["start_new"] =((df.midpoint - (mean/2)).round(0)).astype(int) # calculate new start as the midpoint - (mean length/2)

df["end_new"] = ((df.midpoint + (mean/2)).round(0)).astype(int)

df["desc"] = sid

trimmed = df[["chr", "start_new", "end_new", "len", "core_remodeling", "mrca", "desc"]].drop_duplicates()
 # -a
test_file = "%strimmed/trimmed_%s_FANTOM_complexenh.bed" % (outpath, sid) # trimmed file for multiintersect

trimmed.to_csv(test_file, header = False, index = False, sep = '\t') # save the trimmed file.

# -b
all_beds = glob.glob("%s**/*.bed" % path, recursive = True)
all_beds_str = ' '.join(str(bed) for bed in all_beds)

# multi intersect
multi_cmd = "bedtools intersect -a %s -b %s -f %s -c >\
 %strimmed_%s_multiintersect_%s_count.bed" % (test_file, all_beds_str, \
 fraction_overlap,  outpath, sid, fraction_overlap)
print(multi_cmd)
os.system(multi_cmd)
