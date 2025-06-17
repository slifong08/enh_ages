import glob
import os, sys
import subprocess

# subtract exons from already aged enhancers?
SHUF = 1


RE ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/syn/"

ExonP = "/dors/capra_lab/users/fongsl/data/ensembl/"
ExonF = "%sall_merged_exon.bed" % ExonP


if SHUF==0:

    EnhP = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/trimmed/"
    EnhFs = glob.glob("%sbreaks/trimmed-310-E*_summary_matrix.bed" % EnhP)
    OutP = "%s" % EnhP

elif SHUF==1:
    EnhP = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/"
    EnhFs = glob.glob("%strimmed/shuffle/old_breaks/shuf-trimmed-310_*_parallel_breaks.bed" % EnhP)
    OutP = "%strimmed/shuffle/breaks/non-genic/" % EnhP
    print(OutP)


print(len(EnhFs))
EnhFs[0]

#%% function to subtract exons

def bed_subtract(inF, exonF, outP, fid, SHUF):

    if SHUF ==0:


        outF = "%sno-exon_%s.bed" % (outP, fid)
        outF_ex = "%sexonOverlap_%s.bed" % (outP, fid)


    elif SHUF ==1:

        outF = "%sno-exon_shuf-%s.bed" % (outP, fid)
        outF_ex = "%sexonOverlap_shuf-%s.bed" % (outP, fid)

    # use -v argument to subtract exons from shuffle file.
    cmd = "bedtools intersect -a %s -b %s -v > %s" % (inF, exonF, outF)

    subprocess.call(cmd, shell = True)
    no_exon =  len(open(outF).readlines(  ))


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

    fid = fName.split("_")[0]

    #if fid in missing: # subtract exons from the missing files.
    print(fid)

    no_exon_count, exon_count = bed_subtract(EnhF, ExonF, OutP, fid, SHUF)
    no_exon_list[fName] = no_exon_count
    exon_list[fName] = exon_count


#%%


import numpy as np
exon_list.values()
#%%
print(np.mean(list(no_exon_list.values())))


print(np.mean(list(exon_list.values()))) # ~46% of roadmap enhancers overlap exons
#%% total trimmed = 30687
1612+29075

1612/30687 # % exon overlap = mean 4% for trimmed enhancers
#%% total trimmed shuffle = 2939438
110249 + 2829189
#%% 5% of shuffled enhancers overlap exons

110249/2939438
