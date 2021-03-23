import pandas as pd
import os, sys

PATH = "/dors/capra_lab/users/fongsl/tyler/data/"

AGE = f"{PATH}GG-LL_species-specific_OCRs_rank/breaks/GG-LL_species-specific_OCRs_rank_enh_age_arch_summary_matrix.bed"

SPECIES_SPEC = f"{PATH}GG-LL_species-specific_OCRs_rank.bed"

OUTF = f"{PATH}GG-LL_species-specific_OCRs_rank/breaks/GG-LL_species-specific_OCRs_rank_enh_age_arch_summary_matrix_W_NAMES.bed"

agenames = ["chr", "start", "end", "id",
"sample", "seg_index", "core_remodeling",
"arch", "mrca", "taxon", "mrca_2", "taxon2"]

spenames = ["chr", "start", "end", "species_specific",]

age = pd.read_csv(AGE, sep = '\t', header = None)
age.head()
age.columns = agenames
spe = pd.read_csv(SPECIES_SPEC, sep = '\t', header = None)
spe.columns = spenames

#%%
age[["start", "end"]] = age[["start", "end"]].astype(int)
spe[["start", "end"]] = spe[["start", "end"]].astype(int)
age.head()
spe.head()
df = pd.merge(age, spe, how = "left")

df.to_csv(OUTF, sep = '\t', header = False, index = False)
