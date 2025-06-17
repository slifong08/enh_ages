import numpy as np
import os, sys
import pandas as pd
import subprocess

PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/first_round_breaks/"

EXONPATH = "/dors/capra_lab/users/fongsl/data/ensembl/"
EXONFILE = "ensGene_v36lift37_hg19_coding_exons_autosome_merged.bed"

EXONF = os.path.join(EXONPATH, EXONFILE)

#%%


def get_shuffile_names(iter, path):

    ENHFILE = "SHUFFLE_FANTOM_shuf-all_fantom_enh_112_tissues-%d_age_breaks_summary_matrix.bed" % iter
    ENHFILE_OUT = "shuf-all_fantom_enh_age_breaks_summary_matrix_noexon-%d.bed" % iter
    SYNFILE = "shuf-all_fantom_enh_112_tissues-%d_age_breaks.bed" % iter
    SYNFILE_OUT = "shuf-all_fantom_syn_agebreaks_noexon-%d.bed" % iter

    f_list = [ENHFILE, ENHFILE_OUT, SYNFILE, SYNFILE_OUT]
    f_out_list = []

    for f in f_list:
        f_formatted = os.path.join(path, f)
        f_out_list.append(f_formatted)

    return f_out_list


def bed_subtraction(bed, exonf, outf):

    cmd = "bedtools intersect -a %s -b %s -v > %s" % (bed, exonf, outf)
    subprocess.call(cmd, shell = True)


def crossref_noexon_enhids_in_syn(enhfile, enhfile_out, synfile, synfile_out, iter):

    # open noexon enhancer file.
    col = ["enh_id"]
    keep_ids = pd.read_csv(enhfile_out, sep = '\t',
    header = None, usecols = [3], names = col)
    enhdf_len = len(keep_ids)

    # cross reference syntenic blocks
    syn = pd.read_csv(synfile, sep = '\t',
    header = None)

    new_syn = syn.loc[syn[3].isin(keep_ids.enh_id)] # keep only non-exon enh_ids

    # check that there are the same number of enhancer ids
    #in the non-exon enhancer dataframe and the syn_df dataframe
    syndf_len = len(new_syn[3].unique())

    # write non-exon syn file if there are the same number of enhancers
    # and lengths of dfs are greater than zero
    if enhdf_len == syndf_len and syndf_len>0 and enhdf_len>0:

        new_syn.to_csv(synfile_out, sep = '\t', header = False, index = False)
        cmd = "rm %s %s" % (enhfile, synfile)
        subprocess.call(cmd, shell = True)
        print("removed exons from shuffle iteration", iter)

    else:
        print('something is wrong with these files', iter)


#%%

for iter in np.arange(30,100):

    fs = get_shuffile_names(iter, PATH )

    ENHFILE, ENHFILE_OUT, SYNFILE, SYNFILE_OUT = fs[0], fs[1], fs[2], fs[3]
    if os.path.exists(ENHFILE):
        bed_subtraction(ENHFILE, EXONF, ENHFILE_OUT)
        crossref_noexon_enhids_in_syn(ENHFILE, ENHFILE_OUT, SYNFILE, SYNFILE_OUT, iter)




    
