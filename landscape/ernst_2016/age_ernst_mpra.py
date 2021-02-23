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
CELL_MODEL = "K562"

PATH = "/dors/capra_lab/data/mpra/ernst16/%s_SHARPR-MPRA_scores/" % CELL_MODEL
BASEPAIRFILE = "basepredictions_%s_ScaleUpDesign1and2_combinedP.bed" % CELL_MODEL

BASEPAIRF = os.path.join(PATH, BASEPAIRFILE)
