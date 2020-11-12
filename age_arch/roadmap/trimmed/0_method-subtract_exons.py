import glob
import os, sys
import subprocess

# subtract exons from already aged enhancers?
SHUF = 1


RE ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/syn/"

ExonP = "/dors/capra_lab/users/fongsl/data/ensembl/"
ExonF = "%sensGene_hg19_coding_exons.bed" % ExonP


if SHUF==0:

    EnhP = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/trimmed/"
    EnhFs = glob.glob("%sbreaks/trimmed-310-E*_age_breaks.bed" % EnhP)
    OutP = "%s" % EnhP

elif SHUF==1:
    EnhP = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/"
    EnhFs = glob.glob("%s/trimmed/shuffle/ages/shuf-trimmed-310_*_parallel_breaks.bed" % EnhP)
    OutP = "%strimmed/shuffle/breaks/non-genic/" % EnhP
    print(OutP)


print(len(EnhFs))


#%% function to subtract exons

def bed_subtract(inF, exonF, outP):

    if SHUF ==0:
        fid = (((inF.split("/")[-1]).split(".")[0]).split("_")[0])

        outF = "%sno-exon_%s_age_breaks.bed" % (outP, fid)


    elif SHUF ==1:
        fid = (((inF.split("/")[-1]).split(".")[0]).split("_")[1])
        outF = "%sno-exon_ROADMAP_%s_enh_and_shuf_age_arch_summary_matrix.bed" % (outP, fid)

    # use -v argument to subtract exons from shuffle file.
    cmd = "bedtools intersect -a %s -b %s -v > %s" % (inF, exonF, outF)

    subprocess.call(cmd, shell = True)
    no_exon =  len(open(outF).readlines(  ))

    outF_ex = "%sexonOverlap_%s.bed" % (outP, fid)
    cmd = "bedtools intersect -a %s -b %s > %s" %  (inF, exonF, outF_ex)
    print(cmd)
    subprocess.call(cmd, shell = True)
    exon =  len(open(outF_ex).readlines(  ))

    return no_exon, exon
#%% subtract enhancers
no_exon_list = {}
exon_list = {}
print(len(EnhFs))
#%%
for EnhF in EnhFs:
    fName = (EnhF.split("/")[-1]).split(".")[0]

    fid = fName.split("_")[1]

    #if fid in missing: # subtract exons from the missing files.
    print(fid)

    no_exon_count, exon_count = bed_subtract(EnhF, ExonF, OutP)
    no_exon_list[fName] = no_exon_count
    exon_list[fName] = exon_count


#%%


import numpy as np
exon_list.values()
#%%
print(np.mean(list(no_exon_list.values())))


print(np.mean(list(exon_list.values()))) # ~46% of roadmap enhancers overlap exons
#%%
"""
no_exon mean = 25933.577551020408
exon_mean = 16092.234693877552
"""
