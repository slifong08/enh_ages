#!/usr/bin/env python
# coding: utf-8

#%% In[1]:


#%% Intro

# created 2019-08-05
# updated 2020-10-08
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


#%% set paths, files


path = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/"
outpath = "/dors/capra_lab/projects/enhancer_ages/fantom/data/multiintersect/"
datapath = "/dors/capra_lab/projects/enhancer_ages/fantom/data/architecture_coordinates/"

 # -a = infile

infile = "/dors/capra_lab/projects/enhancer_ages/fantom/data/FANTOM_enh_age_arch_summary_matrix.bed"

#%% # -b
# make -b base from raw enhancer tissue datasets


all_beds = glob.glob("%s*_complexenh.bed" % datapath)

all_beds_str = ' '.join(str(bed) for bed in all_beds)

#all_beds_str # create a string of all the raw bed files.


#%% multi intersect raw length FANTOM_enh_age_arch_summary_matrix.bed


sid = (infile.split("/")[-1]).split(".")[0]
print(sid)


#%% # multi intersect


multi_cmd = "bedtools intersect -a %s -b %s -f %s -c >\
 %sraw_%s_multiintersect_%s_count.bed" % (infile, all_beds_str, \
 fraction_overlap,  outpath, sid, fraction_overlap)
print(multi_cmd)
os.system(multi_cmd)
