import os, sys

import subprocess

CELL_LINES = ["K562", "GM12878", "HepG2", "MCF-7", "PANC1", "HCT116", "PC3", "H1"]
å
PATH = "/dors/capra_lab/data/encode/encode3_hg38/cCRE/"
OUTPATH = "/dors/capra_lab/projects/enhancer_ages/encode/data/"

def combine_files(cell_line, path, outpath):

    cell_path = os.path.join(path, cell_line)

    os.chdir(cell_path)

    outfile = f"{outpath}ELS_combined_{cell_line}.bed"

    cmd = f"cat dELS_{cell_line}.bed dELS,CTCF-bound_{cell_line}.bed pELS_{cell_line}.bed pELS,CTCF-bound_{cell_line}.bed > {outfile}"

    subprocess.call(cmd, shell = True)

    return outfile

def subexon(outfile, outpath, cell_line):

    no_exon_outfile = f"{outpath}no-exon_dELS_combined_{cell_line}.bed"

    cmd = f"exonsubhg38 {outfile} -v > {no_exon_outfile}"

    subprocess.call(cmd, shell = True)

#%%

for cell_line in CELL_LINES:

    outfile = combine_files(cell_line, PATH, OUTPATH)

    subexon(outfile, OUTPATH, cell_line)
