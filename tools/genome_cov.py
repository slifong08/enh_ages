import os, sys
import subprocess

#%%
def get_genomecov(file, build, outpath, sample_id):
    outf = f"{outpath}{sample_id}_{build}genomecov.bed"

    genome_file = f"/dors/capra_lab/data/dna/human/{build}/{build}_trim.chrom.sizes"
    cmd = f"bedtools genomecov -i {file} -g {genome_file} > {outf}"
    print(cmd)
    #subprocess.call(cmd, shell = True)

file = "/dors/capra_lab/projects/enhancer_ages/encode/data/ELS_combined_HepG2.bed"
build = "hg38"
outpath = "/".join(file.split("/")[:-1]) +"/"
sample_id = (file.split("/")[-1]).split(".bed")[0]


get_genomecov(file, build, outpath, sample_id)
