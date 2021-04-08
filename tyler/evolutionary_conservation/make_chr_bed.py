import pandas as pd

PATH = "/dors/capra_lab/data/dna/human/hg38/"
F = f"{PATH}hg38_trim.chrom.sizes"

names = ["#chr", "end"]
df = pd.read_csv(F, sep ='\t', header = None, names = names)
df["start"] = 1
df = df[["#chr", "start", "end"]]
df.head()
#%%
for chr in df["#chr"].unique():
    outf = f"{PATH}chr_bed/{chr}.bed"
    test = df.loc[df["#chr"]== "chr"]
    test.to_csv(outf, sep = '\t', index = False)
    
