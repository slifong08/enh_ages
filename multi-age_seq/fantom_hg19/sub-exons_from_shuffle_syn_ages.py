import os
import pandas as pd
import subprocess

PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/first_round_breaks/"
SHUFFLENH_FILE = "breaks_shuf-all_fantom_enh_ages.bed"
SHUFFLEENHF = os.path.join(PATH, SHUFFLENH_FILE)


OUTFILE = "no-exon_breaks_shuf-all_fantom_enh_ages.bed"
OUTF = os.path.join(PATH, OUTFILE)

SHUFFLESYN_FILE = "syn_breaks_shuf-all_fantom_enh_ages.bed"
SHUFFLESYNF = os.path.join(PATH, SHUFFLESYN_FILE)

OUTSYN = "no-exon_syn_breaks_shuf-all_fantom_enh_ages.bed"
OUTSYNF = os.path.join(PATH, OUTSYN)
OUTSYNF

#%%
def subtract_exons(enh_file, outfile):
    exon_file = "/dors/capra_lab/users/fongsl/data/ensembl/ensGene_v36lift37_hg19_coding_exons.bed"

    cmd = "bedtools intersect -a %s -b %s -v > %s" % (enh_file, exon_file, outfile)
    subprocess.call(cmd, shell = True)

def get_enh_id(file):
    cols = ["chr", "start", "end", "enh_id"]
    df = pd.read_csv(file, sep = '\t', header = None, usecols =[0,1,2,3], names = cols)

    return df.enh_id.unique()

def keep_nonexons(syn, non_exon_ids, outsynf):
    cols = ["chr_syn", "start_syn", "end_syn",
        "enh_id",
        "chr", "start", "end",
        "seg_index", "core_remodeling", "core",
        "mrca",]
    df = pd.read_csv(syn, sep = '\t', header = None, names = cols)

    no_exon = df.loc[df.enh_id.isin(non_exon_ids)]
    no_exon.to_csv(outsynf, sep = '\t', header = None, index = None)
#%%

subtract_exons(SHUFFLEENHF, OUTF)
non_exon_ids = get_enh_id(SHUFFLEENHF)
keep_nonexons(SHUFFLESYNF, non_exon_ids, OUTSYNF)


#%%
