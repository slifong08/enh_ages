# intro

# created 2019-08-05
# updated 2020-10-09 for mean trimmed enhancer length per dataset.
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

#%% path to trimmed means
fpath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/"


# -a
trimmed_means = glob.glob("%snon-genic/trimmed/*.bed" %fpath)
test_list = []
for f in trimmed_means:
    sid = ((f.split("/")[-1]).split('-')[-1]).split(".")[0]
    if sid not in test_list:
        test_list.append(sid)
print(test_list)
#%%

# -a columns = chr start end enh_id maxmrca maxbreaks arch_type

# -b base
all_beds = glob.glob("%s/*/Hsap_H3K27ac_plus_H3K4me3_minus_E*.bed" % fpath)

outpath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/multiintersect/"

for test in test_list:

    print(test)
    # Entire Enhancer #
    # -a

    infile = "%snon-genic/trimmed/trimmed-mean-%s.bed" % (fpath, test)

    # ALL OTHER ROADMAP TISSUES #
    # -b

    all_beds_minus_one = []
    mask_file = "%sHsap_H3K27ac_plus_H3K4me3_minus_%s.bed" % (fpath, test)
    for all_bed_f in all_beds: # remove CORE file from intersection
        if all_bed_f != mask_file:
            all_beds_minus_one.append(all_bed_f)

    all_beds_str = ' '.join(str(bed) for bed in all_beds_minus_one)

    # Enhancer + ROADMAP Bedtools intersection
    multi_cmd = "bedtools intersect -a %s -b %s -f %s -c > %strim-%s_multiintersect_%s_count.bed" %\
    (infile, all_beds_str, fraction_overlap,  outpath, test, fraction_overlap)

    os.system(multi_cmd)
#%%
