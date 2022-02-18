import glob
import os, sys
import subprocess

# subtract exons from already aged enhancers?
AGE_BREAKS = 1


RE ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/syn/"

ExonP = "/dors/capra_lab/users/fongsl/data/ensembl/"
ExonF = "%sensGene_hg19_coding_exons.bed" % ExonP


EnhP = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/non-genic/breaks/to_delete/"
EnhFs = glob.glob("%sno-exon_ROADMAP_E*_enh_and_shuf_age_arch_summary_matrix.bed" % EnhP)


OutP = "%s" % EnhP
print(len(EnhFs))

#%% find out which datasets have already had exons subtracted

#%% function to subtract exons

def bed_subtract(inF, exonF, outP, AGED_BREAKS):

    fName = (((inF.split("/")[-1]).split(".")[0]).split("_")[2])

    outF = "%sno-exon_%s.bed" % (outP, fName)

    # use -v argument to subtract exons from shuffle file.
    cmd = "bedtools intersect -a %s -b %s -v > %s" % (inF, exonF, outF)

    subprocess.call(cmd, shell = True)
    no_exon =  len(open(outF).readlines(  ))

    outF_ex = "%sexonOverlap_%s.bed" % (outP, fName)
    cmd = "bedtools intersect -a %s -b %s > %s" %  (inF, exonF, outF_ex)
    print(cmd)
    subprocess.call(cmd, shell = True)
    exon =  len(open(outF_ex).readlines(  ))

    return no_exon, exon

#%% subtract enhancers
no_exon_list = {}
exon_list = {}
EnhFs
#%%
for EnhF in EnhFs:
    fName = (EnhF.split("/")[-1]).split(".")[0]

    if AGE_BREAKS ==1:

        fid = fName.split("_")[1]
    elif AGE_BREAKS ==0:
        fid = fName.split("_")[-1]


    #if fid in missing: # subtract exons from the missing files.
    print(fid)

    no_exon_count, exon_count = bed_subtract(EnhF, ExonF, OutP, AGE_BREAKS)
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
