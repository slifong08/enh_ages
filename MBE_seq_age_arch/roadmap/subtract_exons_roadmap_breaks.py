import glob
import os
import subprocess

def bed_subtract(inF, fid):

    inP = "/".join(inF.split("/")[:-1]) + "/"
    outP = "%snon-genic/" % inP
    outP = inP

    #mkdir(outP)

    ExonP = "/dors/capra_lab/users/fongsl/data/ensembl/"
    ExonF = "%sall_merged_exon.bed" % ExonP

    outF_noex = "%sno-exon_%s.bed" % (outP, fid)
    outF_ex = "%sexonOverlap_%s.bed" % (outP, fid)

    outP = "%s"
    # use -v argument to subtract exons from shuffle file.
    cmd = "bedtools intersect -a %s -b %s -v > %s" % (inF, ExonF, outF_noex)
    subprocess.call(cmd, shell = True)
    print(cmd)
    no_exon =  len(open(outF_noex).readlines(  ))


    cmd = "bedtools intersect -a %s -b %s > %s" %  (inF, ExonF, outF_ex)
    subprocess.call(cmd, shell = True)
    exon =  len(open(outF_ex).readlines(  ))

    return no_exon, exon
#%%

path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E*/cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_E*/breaks/"

fs = glob.glob("%s*summary_matrix.bed" % path)

for f in fs:
    fid = (f.split("/")[-4]).split("_")[-1]
    print(fid)
    bed_subtract(f, fid)

    
