from Bio import Align
from Bio import AlignIO
from Bio import SeqIO
import glob
import numpy as np
import os
import subprocess


#%%

"""
1. make fasta files for human and rhesus sequences
2. per sequence, score percent identity.

"""
def bed_to_fasta(input_bed, seq_path, id):

    #script = "/dors/capra_lab/bin/bed2fasta.py"
    script = f"{DATA_PATH}bed2fasta.py"
    #INPUT_BED = sys.argv[1]
    #SEQ_PATH = sys.argv[2] # 2bit

    cmd = f"python2.7 {script} {input_bed} {seq_path} {id}"
    print(cmd)
    #subprocess.call(cmd, shell =True)

#%%

DATA_PATH = "/dors/capra_lab/users/fongsl/tyler/data/alignment/"

inputs_dict = {
"hg38":["/dors/capra_lab/users/fongsl/tyler/data/GG-LL_all_OCRs.bed",
"/dors/capra_lab/data/dna/human/hg38/hg38.2bit"],
"rheMac8": ["/dors/capra_lab/users/fongsl/tyler/data/rhemac/all/GG-LL_all_OCRs_rheMac10-from-hg38.liftOver.to.RheMac8.bed",
"/dors/capra_lab/data/dna/macaque/rheMac8/rheMac8.2bit"]
}

os.chdir(DATA_PATH)

#%%
# make the FASTA
for key, value in inputs_dict.items():
    INPUT_BED, SEQ_PATH = value[0], value[1]
    bed_to_fasta(INPUT_BED, SEQ_PATH, key)
    filepath = "/".join(INPUT_BED.split("/")[:-1])
