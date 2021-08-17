import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy import stats
import seaborn as sns


branches = ["hg38", "rheMac8", "hg38-rheMac8"]
prunes = [True, False]
FDR = True
MSA_WAY = 30


# %% FUNCTIONS
def get_vars(branch, msa_way, FDR, prune):
    BASE = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/"
    REBASE = "/dors/capra_lab/projects/enhancer_ages/tyler/results/CON_ACC/"
    if prune is False:
        CONACC_F = f"{BASE}all/multiz30way_{branch}/all_con_acc.bed"
        PATH = "/".join(CONACC_F.split("/")[:-1]) + "/"  # the path
        if FDR is True:
            RE = f"{REBASE}{msa_way}way_{branch}_FDR/"
        else:
            RE = f"{REBASE}{msa_way}way_{branch}/"

        if os.path.exists(RE) is False:
            os.mkdir(RE)
    elif prune is True:
        CONACC_F = f"{BASE}prune/multiz30way_{branch}/all_con_acc.bed"
        PATH = "/".join(CONACC_F.split("/")[:-1]) + "/"  # the path

        if FDR is True:

            RE = f"{REBASE}prune/{msa_way}way_{branch}_FDR/"
        else:
            RE = f"{REBASE}prune/{msa_way}way_{branch}/"

        if os.path.exists(RE) is False:
            os.mkdir(RE)

    return CONACC_F, PATH, RE


def format_f(conacc_f):

    cols = ["#chr", "b", "bf", "start", "end", "conacc", "?", "??", "id"]
    test = pd.read_csv(conacc_f, sep='\t', nrows=5)

    if "start" not in list(test):  # do you need to rename the columns?
        df = pd.read_csv(conacc_f, sep='\t', header=None, names=cols)
        df = df.drop(["b", "bf", "?", "??"], axis=1)  # drop cols not needed
        df.to_csv(conacc_f, sep='\t', index=False)  # save this formatted file.

    else:  # or you've renamed them already and just need to read the df.
        df = pd.read_csv(conacc_f, sep='\t')

    return df


# %%


val = 0
for BRANCH in branches:
    for PRUNE in prunes:

        CONACC_F, PATH, RE = get_vars(BRANCH, MSA_WAY, FDR, PRUNE)

        df = format_f(CONACC_F)  # format the file

        # format the ID column
        df.id = df.id.apply(lambda x: x.split('"')[1])

        # rename the conacc col w/ info
        col_id = f"conacc_{BRANCH}_pruned-{PRUNE}"

        df = df.rename(columns={"conacc": col_id})

        # whittle down the df
        df_ = df[["id", col_id]].drop_duplicates()
        if val == 0:
            results = df_.copy()
        else:
            results = pd.merge(results, df_, how="left", on="id")
        val += 1
# %%

results.head()
# %%
sns.set("talk")
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

data = results

for ax in axes:
    ax.plot([-20, 20], [-20, 20], color="k",
            linestyle="--", transform=ax.transAxes)
    ax.axhline(0, color="grey", linestyle="--")
    ax.axvline(0, color="grey", linestyle="--")

ax = axes[0]
x, y = "conacc_hg38_pruned-False", "conacc_hg38_pruned-True"
results["hg38_delta"] = results[x] - results[y]

ax.plot([-20, 20], [-20, 20], color="k",
        linestyle="--", transform=ax.transAxes)
sns.regplot(x=x, y=y, data=data, ax=ax, scatter_kws={"s": 5})
ax.set(
    xlim=(-21, 21),
    ylim=(-21, 21)
)
print(stats.linregress(results[x], results[y]))

ax = axes[1]
x, y = "conacc_rheMac8_pruned-False", "conacc_rheMac8_pruned-True"
results["rheMac8_delta"] = results[x] - results[y]

ax.plot([-20, 20], [-20, 20], color="k",
        linestyle="--", transform=ax.transAxes)
ax.plot([-20, 20], [-20, 20], color="k",
        linestyle="--", transform=ax.transAxes)
sns.regplot(x=x, y=y, data=data, ax=ax, scatter_kws={"s": 5})
ax.set(
    xlim=(-21, 21),
    ylim=(-21, 21)
)

#  do rheMac8 linearegression w/o nans (n = 14276)
rna = results.dropna(how="any")  # drop nans
a, b = list(rna[x]), list(rna[y])
print(stats.linregress(b, a))

ax = axes[2]
x, y = "conacc_hg38-rheMac8_pruned-False", "conacc_hg38-rheMac8_pruned-True"
results["hg38-rheMac8_delta"] = results[x] - results[y]

ax.plot([-20, 20], [-20, 20], color="k",
        linestyle="--", transform=ax.transAxes)
#plt.plot([-20, 20], [-20, 20], ax=ax, color = "k", linestyle = "--")
sns.regplot(x=x, y=y, data=data, ax=ax, scatter_kws={"s": 5})
ax.set(
    xlim=(-21, 21),
    ylim=(-21, 21)
)
print(stats.linregress(results[x], results[y]))

outf = f"{RE}pruned_v_not_pruned.pdf"
plt.savefig(outf, bbox_inches="tight")
