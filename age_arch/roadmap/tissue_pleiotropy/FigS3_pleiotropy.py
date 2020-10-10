# intro

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


fraction_overlap = 0.5

basepath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/"
fpath = "%sarchitecture_coordinates/" %basepath
path = "%s/download/h3k27ac_plus_h3k4me3_minus_peaks/"%basepath

# -a
test_list =["E050", "E029", "E034", "E069", "E072", "E073", "E118", "E123", "E116"]

#%%
# -a columns = chr start end enh_id maxmrca maxbreaks arch_type

# -b base
all_beds = glob.glob("%sHsap_H3K27ac_plus_H3K4me3_minus_*.bed.gz" % path, recursive = True)

#all_beds = glob.glob("%sE*_ROADMAP_complexenh.bed" % fpath, recursive = True)
all_beds
#%%
corepath = "%sarchitecture_coordinates/" % basepath
outpath = "%smultiintersect/" % basepath
for test in test_list:

    print(test)
    # Entire Enhancer #
    infile = "%s%s_ROADMAP_complexenh.bed" % (fpath, test)

    all_beds_minus_one = []

    for all_bed_f in all_beds: # remove CORE file from intersection
        if test not in all_bed_f:
            all_beds_minus_one.append(all_bed_f)

    all_beds_str = ' '.join(str(bed) for bed in all_beds_minus_one)

    # Enhancer + ROADMAP Bedtools intersection
    multi_cmd = "bedtools intersect -a %s -b %s -f %s -c > %s%s_multiintersect_%s_count.bed" %\
    (infile, all_beds_str, fraction_overlap,  outpath, test, fraction_overlap)

    print(multi_cmd)
    os.system(multi_cmd)
