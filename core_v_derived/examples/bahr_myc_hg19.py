import pandas as pd

path = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/enh_examples/"
age_path = "%sbahr_2018_hg19_myc/ages/" % path

enh_f = "%sbahr_2018_hg19_myc.bed" % path
age_f = "%ssyn_breaks_bahr_2018_hg19_myc_ages.bed" % age_path

enh_col = ["chr_enh", "start_enh", "end_enh", "enh"]
age_col = ["chr_syn", "start_syn", "end_syn", "enh_id", 
"chr_enh", "start_enh", "end_enh",
"seg_index", "core_remodeling", "core", "mrca"]

enh = pd.read_csv(enh_f, sep = '\t', header = None)
age = pd.read_csv(age_f, sep = '\t', header = None)

enh.columns = enh_col
age.columns = age_col

enh.head()
#%%

df = pd.merge(age, enh)
df["syn_len"] = df.end_syn - df.start_syn
df["enh_len"] = df.end_enh - df.start_enh

df
