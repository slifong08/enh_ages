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
NONGENIC = 1

#%%
if NONGENIC == 1:
    fpath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/non-genic/"
    datapath = fpath
elif NONGENIC ==0:
    fpath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/"
    datapath = "/dors/capra_lab/data/roadmap_epigenomics/release_v9/consolidated/histone/h3k27ac_plus_h3k4me3_minus_peaks/"

# -a
test_list =["E050", "E029", "E034", "E069", "E072", "E073", "E118", "E123", "E116"]

desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
desc_df= pd.read_csv(desc_file,header = None)
desc_df.columns = ["sid", "desc"]


#test_list = desc_df.loc[desc_df.desc.str.contains("Fetal"), "sid"].to_list()

#%%
desc_df.loc[desc_df.desc.str.contains("Fetal")]
#%%
# -a columns = chr start end enh_id maxmrca maxbreaks arch_type

# -b base
all_beds_dev = []
if NONGENIC ==1:
    all_beds = glob.glob("%sbreaks/no-exon_E*.bed" % fpath)
elif NONGENIC ==0:
    all_beds = glob.glob("%strimmed/trimmed-mean-*.bed" %fpath)
print(len(all_beds))
#%%
all_beds
#%%
#test_list = desc_df.sid.to_list()
for test in test_list:
    for bed in all_beds:

        if NONGENIC ==1:
            bed_id = ((bed.split("/")[-1]).split("_")[1])
            print(bed_id)

            if test == bed_id:
                all_beds_dev.append(bed)

        elif NONGENIC ==0:
            bed_id = ((bed.split("/")[-1]).split("-")[-1]).split(".")[0]


            if test == bed_id:
                print(test)
                all_beds_dev.append(bed)

#%%
test_list
#%%
len(all_beds_dev)
#%%
outpath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/breaks/multiintersect/"

for test in test_list:


    # Entire Enhancer #
    # -a

    infile = "%strimmed/trimmed-mean-%s.bed" % (fpath, test)

    if os.path.exists(infile)==True:
        # ALL OTHER ROADMAP TISSUES #
        # -b
        print(test)
        all_beds_minus_one = []

        if NONGENIC ==1:
            mask_file = "%sno-exon_%s.bed" % (fpath, test)
        elif NONGENIC ==0:
            mask_file = "%sHsap_H3K27ac_plus_H3K4me3_minus_%s.bed" % (fpath, test)

        #for all_bed_f in all_beds_dev:

        for all_bed_f in all_beds: # remove CORE file from intersection #
            if all_bed_f != mask_file:
                all_beds_minus_one.append(all_bed_f)

        all_beds_str = ' '.join(str(bed) for bed in all_beds_minus_one)

        # Enhancer + ROADMAP Bedtools intersection
        multi_cmd = "bedtools intersect -a %s -b %s -f %s -c > %strim-%s_multiintersect_%s_count.bed" %\
        (infile, all_beds_str, fraction_overlap,  outpath, test, fraction_overlap)

        os.system(multi_cmd)
#%%
multi_cmd
