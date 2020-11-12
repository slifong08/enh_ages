import glob
import os, sys
import subprocess

# subtract exons from already aged enhancers?
REILLY = 0


RE ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/syn/"

ExonP = "/dors/capra_lab/users/fongsl/data/ensembl/"
ExonF = "%sensGene_hg19_coding_exons.bed" % ExonP


if REILLY ==0:

    EnhP = "/dors/capra_lab/projects/enhancer_ages/emera16/data/"
    EnhFs = glob.glob("%sbreaks/emera_2016_neocortical_dev_enhancers_hu_ms_parallel_breaks.bed" % EnhP)


elif REILLY==1:
    EnhP = "/dors/capra_lab/projects/enhancer_ages/reilly15/data/breaks/"
    EnhFs = glob.glob("%sHsap_brain_enhancers*.bed" % EnhP)


OutP = "%snon-genic/" % EnhP
print(len(EnhFs))
if os.path.exists(OutP) == False:
    cmd = "mkdir %s" % OutP
    subprocess.call(cmd, shell = True)


#%% define bedtools intersection command


#%% function to subtract exons

def bed_subtract(inF, exonF, outP):


    fName = (((inF.split("/")[-1]).split(".")[0]))

    outF = "%sno-exon_%s.bed" % (outP, fName)
    temp = "%sno-exon_%s_temp.bed" % (outP, fName)
    cmd = '''sed '1d' %s > %s & mv %s %s''' % (inF, temp, temp, inF)
    subprocess.call(cmd, shell = True)
    # use -v argument to subtract exons from shuffle file.
    cmd = "bedtools intersect -a %s -b %s -v > %s" % (inF, exonF, outF)
    print(cmd)
    subprocess.call(cmd, shell = True)
    no_exon =  len(open(outF).readlines(  ))

    outF_ex = "%sexonOverlap_%s.bed" % (outP, fName)
    cmd = "bedtools intersect -a %s -b %s > %s" %  (inF, exonF, outF_ex)

    subprocess.call(cmd, shell = True)

    exon =  len(open(outF_ex).readlines(  ))

    return no_exon, exon


#%% subtract enhancers


no_exon_list = {}
exon_list = {}
len(EnhFs)
#%%
for EnhF in EnhFs:
    fName = (EnhF.split("/")[-1]).split(".")[0]

    #if fid in missing: # subtract exons from the missing files.
    print(fName)

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
