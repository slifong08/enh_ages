### summary ###

# from /dors/capra_lab/data/mpra/ernst16/*_SHARPR-MPRA_scores/
# BACKGROUND
# 3930 sequence tiles were designed from each H1, K562, hESC, and HepG2 enhancer datasets.
# These 16K tiles were tested for MPRA activity in 2 cell models (K562 and HepG2).
# MPRA activity was then estimated to the basepair level (./basepredictions_K562_ScaleUpDesign1and2_combinedP.bed)

# what I did:
# 1. merged basepredictions_K562_ScaleUpDesign1and2_combinedP.bed
#   to get genomic coordinates of each tile
#
# 2. split up tiles based on cell line enhancers that each "tile" was designed from.
#
# 3. age those enhancers, assemble architectures.
#
# 4. merge those architectures with the per basepair predictions.
#
# 5. evaluate distribution of activity scores for each enhancer architecture.

###

import os, sys
import pandas as pd
import subprocess


#%%


PATH = "/dors/capra_lab/projects/enhancer_ages/ernst16/new_data"

BASEPAIRFILE = "basepredictions_%s_ScaleUpDesign1and2_combinedP.bed" % CELL_MODEL

BASEPAIRF = os.path.join(PATH, BASEPAIRFILE)

#%%

def merge_bases(f, path, cell_model):

    # sort
    sort = os.path.join(path, "%s_sorted_bases.bed" % cell_model)

    cmd = "bedtools sort -i %s > %s" % (f, sort)
    subprocess.call(cmd, shell = True)
    print(cmd)

    # merge
    merged = os.path.join(path, "%s_tiles.bed" % cell_model)

    cmd = "bedtools merge -i %s > %s" % (sort, merged)
    subprocess.call(cmd, shell = True)
    print(cmd)

    #clean up
    cmd = "rm %s" % (sort)
    subprocess.call(cmd, shell = True)


    return merged


def age(merged, path):

    os.chdir(path)
    cmd = "python /dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/awk_breaks.py %s -i 0 -a 1 -b 1 -sh 0 -rt 1 -s hg19" % merged
    subprocess.call(cmd, shell = True)


def merge_age_with_bases(basepairf, path, cell_model):

    AGEPATH = os.path.join(path, "%s_tiles/ages" % cell_model)
    AGEFILE = "syn_breaks_%s_tiles_ages.bed" % cell_model

    AGEF = os.path.join(AGEPATH, AGEFILE)

    # combine base and tile information
    combined = os.path.join(path, "%s_combined_base_tiles_ages.bed" % cell_model)

    cmd = "bedtools intersect -a %s -b %s -wa -wb > %s" %(basepairf, AGEF, combined)
    subprocess.call(cmd, shell = True)

    return combined


#%%

CELL_MODEL = "HEPG2"
merged = merge_bases(BASEPAIRF, PATH, CELL_MODEL)

age(merged, PATH)

combined = merge_age_with_bases(BASEPAIRF, PATH, CELL_MODEL)


#%%

CELL_MODEL = "K562"
merged = merge_bases(BASEPAIRF, PATH, CELL_MODEL)

age(merged, PATH)

combined = merge_age_with_bases(BASEPAIRF, PATH, CELL_MODEL)
