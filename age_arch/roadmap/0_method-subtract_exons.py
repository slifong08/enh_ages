import glob
import os, sys
import subprocess

RE ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/syn/"

ExonP = "/dors/capra_lab/users/fongsl/data/ensembl/"
ExonF = "%sensGene_hg19_coding_exons.bed" % ExonP


EnhP = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/breaks/"
EnhFs = glob.glob("%sROADMAP_*_enh_and_shuf_age_arch_summary_matrix.tsv" % EnhP)

OutP = "%snon-genic/" % EnhP

#%% define bedtools intersection command


def bed_subtract(inF, exonF, outP):

    fName = (((inF.split("/")[-1]).split(".")[0]).split("noheader_")[1])
    outF = "%sno-exon_%s.bed" % (outP, fName)

    # use -v argument to subtract exons from shuffle file.
    cmd = "bedtools intersect -a %s -b %s -v > %s" % (inF, exonF, outF)
    subprocess.call(cmd, shell = True)
    print(cmd)

#%% subtract enhancers


for EnhF in EnhFs:
    fName = (EnhF.split("/")[-1]).split(".")[0]
    tempF = "%snoheader_%s.bed" % (OutP, fName)
    cmd = "sed '1d' %s > %s " % (EnhF, tempF)
    print(cmd)
    subprocess.check_call(cmd, shell = True)
    bed_subtract(tempF, ExonF, OutP)
