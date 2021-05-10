import os, sys
import subprocess


TE = "/dors/capra_lab/data/transposable_elements/repeatmasker/hg38.bed"
PATH = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/"
FS = {
"all": f"{PATH}all/GG-LL_all_OCRs.bed",
#"species-specific_10k": f"{PATH}GG-LL_species-specific_OCRs_rank/GG-LL_species-specific_OCRs_rank.bed",
"hu-specific": f"{PATH}hu_specific/GM12878_specific.bed",
"rhe-specific": f"{PATH}rhe_specific/LCL8664_specific.bed"
}


def subtract_te(key, f, te):

    f_path = "/".join(f.split("/")[:-1]) + "/" # get the file path
    te_less_path = os.path.join(f_path, "subtract_te/") # make a path for TE analysis

    if os.path.exists(te_less_path) == False: # make that directory if it does not exist
        os.mkdir(te_less_path)

    outf = f"{te_less_path}{key}_wo_te.bed"
    bedcmd = f"bedtools intersect -a {f} -b {te} -v > {outf}" # subtract the TE file.

    if os.path.exists(outf) == False:
        print(bedcmd)
        subprocess.call(bedcmd, shell = True)

    return outf

#%%


collection_dict = {}

for key, f in FS.items():
    outf = subtract_te(key, f, TE)
    collection_dict[key] = outf

#%%
collection_dict.values()
