import glob
import os, sys
import subprocess

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/syn/"

ExonP = "/dors/capra_lab/users/fongsl/data/ensembl/"
ExonF = "%sensGene_hg19_coding_exons.bed" % ExonP

ShufP = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"
ShufFs = glob.glob("%s*_age_breaks.bed" % ShufP)

OutP = "%snon-genic/" % ShufP


#%% define bedtools intersection command


def bed_subtract(shufF, exonF, outP):

    fName = (shufF.split("/")[-1]).split(".")[0]
    outF = "%sno-exon_%s.bed" % (outP, fName)

    # use -v argument to subtract exons from shuffle file.
    cmd = "bedtools intersect -a %s -b %s -v > %s" % (shufF, exonF, outF)
    subprocess.call(cmd, shell = True)
    print(cmd)


#%%

for ShufF in ShufFs:

    bed_subtract(ShufF, ExonF, OutP)
