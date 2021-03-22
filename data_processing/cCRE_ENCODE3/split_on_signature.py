# 20210311
# Sarah Fong

# take downloaded ENCODE3 Tier1 cCRE datasets and split on signatures to make separate bedfiles

#%%


import argparse
import os, sys
import subprocess

arg_parser = argparse.ArgumentParser(description="describe argument")

arg_parser.add_argument("bedfile", help = "bedfile with full path")

args = arg_parser.parse_args()

BEDFILE = args.bedfile



def split_signature(bedfile):

    CELL_LINE = bedfile.split("/")[-2]

    PATH = "/".join(bedfile.split("/")[:-1])

    os.chdir(PATH)

    cmd = '''awk '{print >$10"_%s.bed"}' %s''' % (CELL_LINE, bedfile)

    subprocess.call(cmd, shell = True)

#%%
split_signature(BEDFILE)
