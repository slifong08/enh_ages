"""
### summary ###

# process GENCODE exon.bed coordinates downloaded from UCSC genome browser.
# handles hg19 and hg38 builds

# method
# remove autosomes, keep only bed coordinates from GENCODE_exon.bed
# strip white space
# sort and merge exon.bed coordinates
# concatenate with black_list regions
"""

import os
import pandas as pd
import subprocess

PATH = "/dors/capra_lab/projects/enhancer_ages/dna/"

GENPATH = os.path.join(PATH, "gencode")

# %% Functions


def get_genf(genpath, build):

    if build == "hg38":
        GENFILE = "GENCODE_v38_hg38"
    else:
        GENFILE = "GENCODE_v38_lift37_hg19"

    genf = os.path.join(genpath, GENFILE)

    return genf


# get the coding exon sequence boundaries
def parse_gencode(genf):

    # get the file name
    # genf = os.path.join(GENPATH, genf_name)

    # file to write coding exon boundaries
    newf = genf + "_cds.bed"
    newtxn = genf + "_tx.bed"

    # column names
    colnames = [
                "#chrom", "txStart", "txEnd",  "name", "score", "strand",
                "cdsStart", "cdsEnd",  "name2",
                "exonCount", "exonStarts", "exonEnds",
                ]

    # open dataframe
    df = pd.read_csv(genf, sep="\t", header=None)

    # rename data
    df.columns = colnames

    # make bedfile df, sum names
    ndf = df.groupby(["#chrom", "cdsStart", "cdsEnd"])[
                     "name"].sum().reset_index()
    txdf = df.groupby(["#chrom", "txStart", "txEnd"])[
                     "name"].sum().reset_index()

    # write .bed file
    ndf.to_csv(newf, sep="\t", index=False)
    txdf.to_csv(newtxn, sep="\t", index=False)

    return newf


# get autosomes
def get_unique_bed(f, genpath):

    # file for autosomes
    auto = (f.split("/")[-1]).split(".")[0] + "-autosome_only.bed"
    autof = os.path.join(genpath, auto)

    auto_sorted = (f.split("/")[-1]).split(".")[0] + "-autosome_sorted.bed"
    autos = os.path.join(genpath, auto_sorted)

    # file for uniq autosome coordinates.
    merged = (f.split("/")[-1]).split(".")[0] + "_autosome_merged.bed"
    mergedf = os.path.join(genpath, merged)

    # get the autosomes and print the bed file.
    cmd = '''awk ' !($1 ~ /chrX/ ) { print $1, "\t", $2, "\t", $3}' %s > %s'''\
        % (f, autof)
    subprocess.call(cmd, shell=True)

    # strip whitespace
    df = pd.read_csv(autof, header=None, sep='\t')
    df = df[[0, 1, 2]].astype(str).apply(lambda x: x.str.strip())
    df.to_csv(autof, sep='\t',  header=False, index=False)

    # sort bed
    cmd = "bedtools sort -i % s> %s" % (autof, autos)
    subprocess.call(cmd, shell=True)

    # merge bed coordinates
    cmd = "bedtools merge -i %s > %s" % (autos, mergedf)
    subprocess.call(cmd, shell=True)

    # clean up
    cmd = "rm %s %s" % (autof, autos)
    subprocess.call(cmd, shell=True)

    return mergedf


def combine_autosome_blacklist(f, path, build):

    black_list = "%s%s_blacklist_gap.bed" % (path, build)

    bl = pd.read_csv(black_list, sep='\t', header=None, usecols=[0, 1, 2])

    ex = pd.read_csv(f, sep='\t', header=None, usecols=[0, 1, 2])

    df = pd.concat([ex, bl])

    outf = "%s%s_blacklist_gap_gencode.bed" % (path, build)
    df.to_csv(outf, sep='\t', header=False, index=False)


# %%

BUILD = "hg19"

genf = get_genf(GENPATH, BUILD)  # get the gencode file

newf = parse_gencode(genf)  # slice the GENCODE FILES

mergedf = get_unique_bed(newf, GENPATH)

combine_autosome_blacklist(mergedf, PATH, BUILD)

# %%

BUILD = "hg38"

genf = get_genf(GENPATH, BUILD)  # get the gencode file

newf = parse_gencode(genf)  # slice the GENCODE FILES

mergedf = get_unique_bed(newf, GENPATH)

combine_autosome_blacklist(mergedf, PATH, BUILD)

genf
