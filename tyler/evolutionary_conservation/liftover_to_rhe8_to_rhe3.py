import matplotlib.pyplot as plt
import os, sys
import pandas as pd
import seaborn as sns
import subprocess

RE = "/dors/capra_lab/users/fongsl/tyler/results/liftover/"


def liftover(bedfile, path, from_build, to_build, key): # bedfile with full path

    # prepare annotations and files
    sid = (bedfile.split("/")[-1]).split(".")[0] # get the sample ID

    chainPath = "/dors/capra_lab/data/ucsc/liftOver/" # path to chain file
    chainf = f"{chainPath}{from_build}To{to_build}.over.chain.gz"

    #write to result files
    lifted = f"{path}{key}/{sid}.liftOver.to.{to_build}.bed" # name the liftover file
    notlifted = f"{path}{key}/{sid}.notlifted.to.{to_build}.bed" # name the notlifted file

    if os.path.exists(lifted) == False: # check that you haven't run this.

        ### sort the bedfile ###
        tempbed = f"{path}temp_{sid}.bed" # format the bed file into 5 columns

        # [[chr start end enh_id sample_id]] and sort by coordinates

        cmd = f"sort -k1,1 -k2,2 -k3,3 {bedfile} > {tempbed}"
        print("standardizing Bed format")
        subprocess.call(cmd, shell=True)

        ### liftover the formatted bedfile ###

        cmd = f"liftOver {tempbed} {chainf} {lifted} {notlifted}"
        print("liftingOver", sid)
        subprocess.call(cmd, shell=True)
        print("done lifting.")


        ### clean up temp ###
        if os.path.getsize(lifted) >0:
            os.remove(tempbed)
            print("cleaned up temp file")

    return lifted


def count_filelen(file):
    count = len(open(file).readlines(  ))
    return count


#%% ALL OCR
PATH = "/dors/capra_lab/users/fongsl/tyler/data/rhemac/"

file_dict = {
"all": f"{PATH}all/GG-LL_all_OCRs_rheMac10-from-hg38.bed",
"shared": f"{PATH}shared/GG-LL_shared_OCRs_rheMac10-from-hg38.bed",
"hu-specific": f"{PATH}hu-specific/GG-LL_GM12878-specific_OCRs_rheMac10-from-hg38.bed",
"rhe-specific":f"{PATH}rhe-specific/GG-LL_LCL8664-specific_OCRs_rheMac10-from-hg38.bed"}

results_dict = {}

for KEY, RM10 in file_dict.items():


    #RheMac10 to RheMac8
    FROM_BUILD = "rheMac10"
    TO_BUILD = "RheMac8"

    RM8 = liftover(RM10, PATH, FROM_BUILD, TO_BUILD, KEY)

    # RheMac8 to RheMac3
    FROM_BUILD = "rheMac8"
    TO_BUILD = "RheMac3"

    RM3 = liftover(RM8, PATH, FROM_BUILD, TO_BUILD, KEY)

    #  COUNTS
    RM10_count = count_filelen(RM10) # 149563
    RM8_count = count_filelen(RM8) #146902, lost 2661 from RheMac10
    RM3_count = count_filelen(RM3) # 136089, lost 13474 from RheMac10,


    RMdif_10_8 = RM10_count - RM8_count # lost 2661 from RheMac10
    RMdif_10_3 = RM10_count - RM3_count # lost 13474 from RheMac10,
    RMdif_8_3 = RM8_count - RM3_count # lost 10813 from Rhemac8
    key_count = KEY + f" (n= {str(RM10_count)})"
    results = pd.DataFrame({
    "dataset":[KEY, KEY, KEY],
    "comparison":["RheMac10 v. RheMac8", "RheMac10 v. RheMac3", "RheMac8 v. RheMac3"],
    "difference":[RMdif_10_8, RMdif_10_3, RMdif_8_3],
    "percent_overlap":[ (1-RMdif_10_8/RM10_count),  (1-RMdif_10_3/RM10_count), (1-RMdif_8_3/RM8_count)],
    "n_RheMac10":[RM10_count, RM10_count, RM10_count],
    "n_RheMac8":[RM8_count, RM8_count, RM8_count],
    "n_RheMac3":[RM3_count, RM3_count, RM3_count],
    })

    results_dict[KEY] = results
    print("Rhemac10", RM10_count, "Rhemac8", RM8_count, "Rhemac3", RM3_count)
    print("Rhemac10 - Rhemac3", RMdif_10_3, "Rhemac10 - Rhemac8", RMdif_10_8, "Rhemac8 - Rhemac3", RMdif_8_3)
    print("% retained Rhemac10 - Rhemac3 ", 1-RMdif_10_3/RM10_count, "% retained Rhemac10 - Rhemac8", 1-RMdif_10_8/RM10_count, "% retained Rhemac8 - Rhemac3", 1-RMdif_8_3/RM8_count)
#%%
df = pd.concat(results_dict.values())
df.head()
x = "comparison"
y = "percent_overlap"
data = df
hue = "dataset"
hue_order = ["all", "shared", "hu-specific", "rhe-specific"]
xlabs = ["10 v. 8", "10 v. 3", "8 v. 3"]
n = df[["dataset", "n_RheMac10"]].drop_duplicates()

n["legend"] = n.dataset + " n = " + n.n_RheMac10.map(str)

sns.set("poster")
fig, ax = plt.subplots(figsize = (9,9))

sns.barplot(
x =x,
y =y,
data = data,
hue = hue,
hue_order = hue_order
)

ax.set(ylim = (0.85, 1.0), xlabel = "RheMac build\n%s" % n["legend"].to_list())
ax.set_xticklabels(xlabs)

outf = f"{RE}liftOver_all_shared_specific.pdf"
plt.savefig(outf, bbox_inches = "tight")

outf = f"{RE}liftOver_all_shared_specific.tsv"
df.to_csv(outf, sep = '\t', index = False)
