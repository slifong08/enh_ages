import os, sys

import subprocess

CELL_LINES = ["K562", "GM12878", "HepG2", "MCF-7", "PANC1", "HCT116", "PC3", "H1"]

PATH = "/dors/capra_lab/data/encode/encode3_hg38/cCRE/"
OUTPATH = "/dors/capra_lab/projects/enhancer_ages/encode/data/"
BUILD = 'hg38'

#%%
def combine_files(cell_line, path, outpath):

    cell_path = os.path.join(path, cell_line)

    os.chdir(cell_path)

    outfile = f"{outpath}ELS_combined_{cell_line}.bed"

    cmd = f"cat dELS_{cell_line}.bed dELS,CTCF-bound_{cell_line}.bed pELS_{cell_line}.bed pELS,CTCF-bound_{cell_line}.bed > {outfile}"

    subprocess.call(cmd, shell = True)

    return outfile

def subexon(outfile, outpath, cell_line):

    no_exon_outfile = f"{outpath}no-exon_ELS_combined_{cell_line}.bed"

    cmd = f"exonsubhg38 {outfile} -v > {no_exon_outfile}"

    #subprocess.call(cmd, shell = True)


def get_genomecov(file, build, outpath, sample_id):
    outf = f"{outpath}{sample_id}_{build}genomecov.bed"

    genome_file = f"/dors/capra_lab/data/dna/human/{build}/{build}_trim.chrom.sizes"
    cmd = f"bedtools genomecov -i {file} -g {genome_file} > {outf}"

    subprocess.call(cmd, shell = True)

#%%

for cell_line in CELL_LINES:

    outfile = combine_files(cell_line, PATH, OUTPATH)

    get_genomecov(outfile, BUILD, OUTPATH, cell_line)

    subexon(outfile, OUTPATH, cell_line)
