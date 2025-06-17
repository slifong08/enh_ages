import glob
import os, sys
import subprocess

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/syn/"

ExonP = "/dors/capra_lab/projects/enhancer_ages/te/data/"
ExonF = "%ssimple_repeats.bed" % ExonP

ShufP = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"
ShufFs = glob.glob("%s*_age_breaks.bed" % ShufP)


EnhP = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"
EnhFs = glob.glob("%sFANTOM_enh_age_arch_*.tsv" % EnhP)

ShufOutP = "%sno-repeats/" % EnhP
OutP = "%sno-repeats/" % EnhP

#%% define bedtools intersection command


def bed_subtract(shufF, exonF, outP):

    fName = (shufF.split("/")[-1]).split(".")[0]
    outF = "%sno-repeats_%s.bed" % (outP, fName)

    # use -v argument to subtract exons from shuffle file.
    cmd = "bedtools intersect -a %s -b %s -v > %s" % (shufF, exonF, outF)
    subprocess.call(cmd, shell = True)
    print(cmd)


#%% subtract shuffles


for ShufF in ShufFs:

    bed_subtract(ShufF, ExonF, ShufOutP)


#%% subtract enhancers


for EnhF in EnhFs:
    fName = (EnhF.split("/")[-1]).split(".")[0]
    tempF = "%stemp_%s.bed" % (OutP, fName)
    cmd = "sed '1d' %s > %s " % (EnhF, tempF)
    print(cmd)
    subprocess.check_call(cmd, shell = True)
    bed_subtract(tempF, ExonF, OutP)
