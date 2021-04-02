### summary ###

# process GENCODE exon.bed coordinates downloaded from UCSC genome browser.
# handles hg19 and hg38 builds

# method
    # remove autosomes, keep only bed coordinates from GENCODE_exon.bed
    # strip white space
    # sort and merge exon.bed coordinates
    # concatenate with black_list regions



import os, sys
import pandas as pd
import subprocess

BUILD = "hg19"


PATH = "/dors/capra_lab/users/fongsl/data/"

ENSPATH =os.path.join(PATH, "ensembl")

if BUILD == "hg38":
    ENSFILE = "ensGene_v32_hg38_coding_exons.bed"
else:
    ENSFILE = "ensGene_v36lift37_hg19_coding_exons.bed"

ENSF = os.path.join(ENSPATH, ENSFILE)


#%%

def get_unique_bed(f, enspath):

    # file for autosomes
    auto = (f.split("/")[-1]).split(".")[0] + "-autosome_only.bed"
    autof = os.path.join(enspath, auto)

    auto_sorted = (f.split("/")[-1]).split(".")[0] + "-autosome_sorted.bed"
    autos = os.path.join(enspath, auto_sorted)

    # file for uniq autosome coordinates.
    merged = (f.split("/")[-1]).split(".")[0] + "_autosome_merged.bed"
    mergedf = os.path.join(enspath, merged)

    # get the autosomes and print the bed file.
    cmd = '''awk ' !($1 ~ /chrX/ ) { print $1, "\t", $2, "\t", $3}' %s > %s''' % (f, autof)
    subprocess.call(cmd, shell = True)

    # strip whitespace
    df = pd.read_csv(autof, header = None, sep ='\t')
    df = df[[0,1,2]].astype(str).apply(lambda x: x.str.strip())
    df.to_csv(autof, sep = '\t',  header = False, index = False)


    # sort bed
    cmd = "bedtools sort -i % s> %s" % (autof, autos)
    subprocess.call(cmd, shell = True)

    # merge bed coordinates
    cmd = "bedtools merge -i %s > %s" % (autos, mergedf)
    subprocess.call(cmd, shell = True)

    # clean up
    cmd = "rm %s %s" % (autof, autos)
    subprocess.call(cmd, shell = True)

    return mergedf


def combine_autosome_blacklist(f, path, build):

    black_list = "%s%s_blacklist_gap.bed" % (path, build)

    bl = pd.read_csv(black_list, sep ='\t', header = None, usecols = [0,1,2])

    ex = pd.read_csv(f, sep ='\t', header = None, usecols = [0,1,2])

    df = pd.concat([ex, bl])

    outf = "%s%s_blacklist_gap_ensemblexon.bed" % (path, build)
    df.to_csv(outf, sep ='\t', header = False, index = False)

#%%
mergedf = get_unique_bed(ENSF, ENSPATH)
mergedf
combine_autosome_blacklist(mergedf, PATH, BUILD)
