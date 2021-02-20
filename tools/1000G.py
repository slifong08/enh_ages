### Summary ###


# intersect .bed with 1000G per chromosome.
# (1000G data is >1gb, so I split by chromosome to handle smaller data)
# ultimately, output file adds 11 columns to original .bed file
# 11 columns are
# [chr, pos-1, pos, qual, maf, 'EAS_AF_maf','AMR_AF_maf',\
# 'AFR_AF_maf', 'EUR_AF_maf', 'SAS_AF_maf','id']


#%%


# load modules
import argparse
import glob
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import subprocess

# Files and paths

arg_parser = argparse.ArgumentParser(description="Calculate linsight for .bed")

arg_parser.add_argument("bedfile", help='bed file w/ full path')

args = arg_parser.parse_args()

ENHANCERF = args.bedfile

ENHANCERPATH = "/".join(ENHANCERF.split("/")[:-1]) + "/"
ENHANCERFILE = ENHANCERF.split("/")[-1]
ENHANCER_ID = ".".join(ENHANCERFILE.split(".")[:-1])

THOUPATH = "/dors/capra_lab/projects/enhancer_ages/1000g/data/maf/"

OUTPATH = os.path.join(ENHANCERPATH, "1000g")


if not os.path.exists(OUTPATH):
    os.makedirs(OUTPATH)
    print("Directory '%s' created" %OUTPATH)

#%% FUNCTIONS


def chr_lst(): # make a list of autosomes

    chr_list = []

    for n in np.arange(1,23):

         chr_list.append("chr"+str(n))

    return chr_list


def bed_intersect(enhf, thouf, chr_num, outfile):

    bed_cmd = "bedtools intersect -a %s -b %s -wa -wb> %s" % (enhf, thouf, outfile)

    subprocess.call(bed_cmd, shell = True)

    print("finished", outfile)

def cleanup_dir(path, enh_id):

        # concat chromosome files
        cat_out = "%s/%s_x_1000g.bed" %(path, enh_id)

        cat_cmd = "cat %s/chr*.bed > %s" %(path, cat_out)

        subprocess.call(cat_cmd, shell = True)

        # clean up the chromosome files.
        cleanup_cmd = "rm %s/chr*.bed" %(path)

        subprocess.call(cleanup_cmd, shell = True)

        return cat_out


#%% run functions

chr_list = chr_lst()

for chr_num in chr_list:

    thoufile_chr = "trimmed.%s.phase3_shapeit2_mvncall_integrated_v5a.maf.bed.gz" % chr_num
    thouf_chr = os.path.join(THOUPATH, thoufile_chr)

    OUTFILE = "%s_%s_x_1000g.bed" % (chr_num, ENHANCER_ID)
    OUTF = os.path.join(OUTPATH, OUTFILE)

    bed_intersect(ENHANCERF, thouf_chr, chr_num, OUTF)

#%%

intersection_file = cleanup_dir(OUTPATH, ENHANCER_ID)
