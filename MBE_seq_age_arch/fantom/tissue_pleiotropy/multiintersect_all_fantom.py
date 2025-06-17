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
PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"

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

    if "000" in a:
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

trims = ["310"]

for trim in trims:
    # -a
    F = f"{PATH}non-genic/trimmed/trim_310_all_fantom_enh.bed"
    sid = (F.split("/")[-1]).split(".bed")[0]

    # -b base
    all_beds = glob.glob(f"{PATH}download/*/non-genic/no-exon_*.bed")

    # multi intersect
    outf = multi_intersect(F, all_beds, FRAC_OVERLAP, sid)
