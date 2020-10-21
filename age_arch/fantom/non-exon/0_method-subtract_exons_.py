import glob
import os, sys
import subprocess

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/syn/"

ExonP = "/dors/capra_lab/users/fongsl/data/ensembl/"
ExonF = "%sensGene_hg19_coding_exons.bed" % ExonP

ShufP = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"
ShufFs = glob.glob("%s*_age_breaks.bed" % ShufP)


EnhP = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"
EnhFs = glob.glob("%sFANTOM_enh_age_arch_*.tsv" % EnhP)
Shuf_collapsed = ["%sSHUFFLE_FANTOM_enh_age_arch_summary_matrix.tsv" % EnhP]

ShufOutP = "%snon-genic/" % EnhP
OutP = "%snon-genic/" % EnhP

#%% define bedtools intersection command


def bed_subtract(inF, exonF, outP):


        #fName = (inF.split("/")[-1]).split(".")[0]

    fName = (((inF.split("/")[-1]).split(".")[0]).split("noheader_")[0])
    outF = "%sno-exon_%s.bed" % (outP, fName)

    # use -v argument to subtract exons from shuffle file.
    cmd = "bedtools intersect -a %s -b %s -v > %s" % (inF, exonF, outF)
    subprocess.call(cmd, shell = True)
    no_exon =  len(open(outF).readlines(  ))
    print(cmd)

    outF_ex = "%sexonOverlap_%s.bed" % (outP, fName)
    cmd = "bedtools intersect -a %s -b %s > %s" %  (inF, exonF, outF_ex)
    subprocess.call(cmd, shell = True)
    exon = len(open(outF_ex).readlines(  ))
    print(no_exon, exon)
    return no_exon, exon

 #%% subtract shuffles
no_exon_list = []
exon_list = []
for ShufF in ShufFs:

    noex, ex = bed_subtract(ShufF, ExonF, ShufOutP)
    no_exon_list.append(noex)
    exon_list.append(ex)


#%%
import numpy as np
print(np.mean(no_exon_list)) # mean number of syntenic blocks not overlapping exons
#

print(np.mean(exon_list)) # mean number of syntenic blocks overlapping exons
# 4082 average # of overlaps

#%% subtract enhancers


for EnhF in EnhFs:
    fName = (EnhF.split("/")[-1]).split(".")[0]
    tempF = "%snoheader_%s.bed" % (OutP, fName)
    cmd = "sed '1d' %s > %s " % (EnhF, tempF)
    print(cmd)
    subprocess.check_call(cmd, shell = True)
    noex, ex = bed_subtract(tempF, ExonF, OutP)
    print(noex, ex)
    cmd = "rm %s" % tempF
    #subprocess.call(cmd, shell = True) # remove the temp file
#%%
for EnhF in Shuf_collapsed:
    fName = (EnhF.split("/")[-1]).split(".")[0]
    tempF = "%snoheader_%s.bed" % (OutP, fName)
    cmd = "sed '1d' %s > %s " % (EnhF, tempF)
    print(cmd)
    subprocess.check_call(cmd, shell = True)
    noex, ex = bed_subtract(tempF, ExonF, OutP)
    print(noex, ex)
    cmd = "rm %s" % tempF
    #subprocess.call(cmd, shell = True) # remove the temp file
