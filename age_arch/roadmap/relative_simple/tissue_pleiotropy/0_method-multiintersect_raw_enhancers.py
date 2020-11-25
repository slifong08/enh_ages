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

#%%
fpath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/non-genic/"

desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
desc_df= pandas.read_csv(desc_file, header = None)
desc_df.columns = ["sid", "desc"]
desc_df.head()

dev = desc_df.loc[desc_df.desc.str.contains("Fetal"), "sid"]

#%%
# -a
test_list =["E050", "E029", "E034", "E069", "E072", "E073", "E118", "E123", "E116"]


# -a columns = chr start end enh_id maxmrca maxbreaks arch_type

# -b base
all_beds = glob.glob("%sno-exon_E*.bed" % fpath)

outpath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/multiintersect/"
outpath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/non-genic/multiintersect/"

for test in dev:

    print(test)
    # Entire Enhancer #
    # -a

    infile = "%sno-exon_%s.bed" % (fpath, test)

    # ALL OTHER ROADMAP TISSUES #
    # -b

    all_beds_minus_one = []
    mask_file = "%ssno-exon_%s.bed" % (fpath, test)
    for all_bed_f in all_beds: # remove CORE file from intersection
        if all_bed_f != mask_file:
            all_beds_minus_one.append(all_bed_f)

    all_beds_str = ' '.join(str(bed) for bed in all_beds_minus_one)

    # Enhancer + ROADMAP Bedtools intersection
    multi_cmd = "bedtools intersect -a %s -b %s -f %s -c > %sraw-%s_multiintersect_%s_count.bed" %\
    (infile, all_beds_str, fraction_overlap,  outpath, test, fraction_overlap)

    os.system(multi_cmd)
#%%

dev.tolist()