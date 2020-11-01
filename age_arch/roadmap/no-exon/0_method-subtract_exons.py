import glob
import os, sys
import subprocess

RE ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/syn/"

ExonP = "/dors/capra_lab/users/fongsl/data/ensembl/"
ExonF = "%sensGene_hg19_coding_exons.bed" % ExonP


EnhP = "/dors/capra_lab/data/roadmap_epigenomics/release_v9/consolidated/histone/h3k27ac_plus_h3k4me3_minus_peaks/"
EnhFs = glob.glob("%sHsap_H3K27ac_plus_H3K4me3_minus_E*.bed" % EnhP)


OutP = "%snon-genic/" % EnhP
print(len(EnhFs))

already_done = glob.glob("%sno-exon_ROADMAP_E*_enh_and_shuf_age_arch_summary_matrix.bed" % OutP)
already_done_list = []
for f in already_done:
    sid = (f.split("/")[-1]).split("_")[2]
    already_done_list.append(sid)
len(already_done_list)
#%% define bedtools intersection command
missing = ["E003","E063", "E105", "E056", "E111", "E012", "E034", "E098", "E029", "E008", "E102", "E078", "E066", "E114", "E103", "E055"]
print(len(missing))
#%%
def bed_subtract(inF, exonF, outP):

    fName = (((inF.split("/")[-1]).split(".")[0]).split("_")[-1])

    outF = "%sno-exon_%s.bed" % (outP, fName)

    # use -v argument to subtract exons from shuffle file.
    cmd = "bedtools intersect -a %s -b %s -v > %s" % (inF, exonF, outF)

    subprocess.call(cmd, shell = True)
    print(cmd)
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

for EnhF in EnhFs:
    fName = (EnhF.split("/")[-1]).split(".")[0]
    fid = fName.split("_")[-1]
    if fid in missing:
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
