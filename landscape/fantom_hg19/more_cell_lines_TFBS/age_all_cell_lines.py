import os
import subprocess

PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/ENCODE3/"
def get_cell_lines():
    sample_dict = {
    "GM12878": "GM12878_FANTOM5_tpm_hg19.bed",
    "HepG2": "HEPG2_FANTOM5_tpm_hg19.bed",
    "A549":"A549_FANTOM5_tpm_hg19.bed",
    "H1-ESC":"H1-esc_FANTOM5_tpm_hg19.bed",
    "PC-3": "PC-3_FANTOM5_tpm_hg19.bed",
    "K562":"CL_0000094_granulocyte_expressed_enhancers.bed",
    "MCF7":"CL_0002327_mammary_epithelial_cell_expressed_enhancers.bed",
    "liver":"UBERON_0002107_liver_expressed_enhancers.bed",
    }

    return sample_dict

def age_samples(file):

    cmd = "age %s 0 1 1 0 1 hg19" % file
    print(cmd)
    subprocess.call(cmd, shell = True)


#%%

 sample_dict = get_cell_lines()

 for file in sample_dict.values():
    os.chdir(PATH)
    age_samples(file)

#%%
print("done")
