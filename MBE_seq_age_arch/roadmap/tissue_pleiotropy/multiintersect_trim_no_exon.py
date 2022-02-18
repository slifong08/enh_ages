# intro

# created 2019-08-05
# updated 2020-10-09 for mean trimmed enhancer length per dataset.
# sarahfong

# Goal - multiintersect all roadmap histone enhancers with a count for the number of overlapping coordinates between datasets

# Method - multiintersect

## bedtools intersect -a [all_roadmap_enhancers] -b [E005_core_ages.bed, E004_core_ages.bed ... etc.]

## -c (count number of overlapping features)

import glob
import os, sys
import pandas as pd
import subprocess

FRAC_OVERLAP = 0.5

#%% path to trimmed means
PATH = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/"

#%% Function


def mask_file(a, b_list):
    print("masking", a)
    b_list_new = []

    for b in b_list: # remove CORE file from intersection
        if b != a:
            b_list_new.append(b)

    return b_list_new


def multi_intersect(a, b_list, fraction_overlap, sid):

    outpath = "/".join(a.split("/")[:-1]) + "/multiintersect/"

    if os.path.exists(outpath) == False:
        os.mkdir(outpath)

    if "Hsap_H3K27ac_plus_H3K4me3_minus_E" in a:
        b_list_masked = mask_file(a, b_list)
        all_beds_str = ' '.join(str(bed) for bed in b_list_masked)

    else:
        # make a string of the multi_intersect files
        all_beds_str = ' '.join(str(bed) for bed in b_list)

    # make afile to write to
    outf = f"{outpath}{sid}_multiintersect_{fraction_overlap}_count.bed"

    # do the intersection
    multi_cmd = f"bedtools intersect -a {a} -b {all_beds_str} -f {fraction_overlap} -c > {outf}"
    #print(multi_cmd)
    subprocess.call(multi_cmd, shell = True)

    print("made file", outf)

    return outf


#%% multi intersect for all roadmap enhancers

trims = ["310", "1934"]

for trim in trims:
    # -a
    F = f"{PATH}all_roadmap_enh/trimmed/trim_{trim}_no-exon_all_roadmap_enh.bed"
    sid = (F.split("/")[-1]).split(".bed")[0] + "_mean"

    # -b base
    all_beds = glob.glob(f"{PATH}/*/Hsap_H3K27ac_plus_H3K4me3_minus_E*.bed")

    # multi intersect
    outf = multi_intersect(F, all_beds, FRAC_OVERLAP, sid)

#%% for each roadmap enhancer

pathpath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E*/non-genic/trimmed/trim_*_E*.bed"

FS = glob.glob(pathpath)

print(len(FS))

all_beds_trim_mean = []
all_beds_trim_310 = []

# sort out the means from the 310 trimmed files
for f in FS:
    if "310" in f:
        all_beds_trim_310.append(f)
    else:
        all_beds_trim_mean.append(f)

print(len(all_beds_trim_310), len(all_beds_trim_mean))

#%% run multi intersect on means
for F in all_beds_trim_mean:
    sid = (F.split("_")[-1]).split(".")[0] + "_mean"
    outf = multi_intersect(F, all_beds_trim_mean, FRAC_OVERLAP, sid)

#%% run multi intersect on 310
for F in all_beds_trim_310:
    sid = (F.split("_")[-1]).split(".")[0]+ "_310"

    outf = multi_intersect(F, all_beds_trim_mean, FRAC_OVERLAP, sid)
