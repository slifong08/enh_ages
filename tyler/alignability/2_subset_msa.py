from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
import os, sys
import subprocess

CHR = "chr21"
PATH = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/all/"

NEW_PATH = "/dors/capra_lab/users/fongsl/tyler/data/alignment/"

#%%

def maf_x_bed(chr, chr_bed, path):

    os.chdir(path)

    MSA_PATH = "/dors/capra_lab/bin/"
    MAF_PATH = "/dors/capra_lab/data/ucsc/hg38/multiz30way/maf/"

    # args
    msa_func = f"{MSA_PATH}maf_parse"
    maf_arg = f"{MAF_PATH}{chr}.maf"
    feat_arg = f"--features {chr_bed}"
    order_arg = "--seqs hg38,rheMac8"
    out = f"{chr}_parse.maf"

    cmd = f" {msa_func} {maf_arg} {feat_arg} {order_arg}  > {out}"
    #subprocess.call(cmd, shell = True)
    print(cmd)


def make_chr_list():
    n = list(np.arange(1, 23))
    #n.append("X")

    chr_list = []
    for num in n:
        chrn = "chr" + str(num)
        chr_list.append(chrn)

    return chr_list
#%%

chrList = make_chr_list()

for CHR in chrList:
    F = f"{CHR}.bed"
    CHR_BED = os.path.join(PATH, F)
    maf_x_bed(CHR, CHR_BED, PATH)

#%%
cmd = f"mv {PATH}chr*_parse.maf {NEW_PATH}"
subprocess.call(cmd, shell = True)

#%%
