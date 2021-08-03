import glob
import os
import subprocess
import numpy as np
from Bio import Align
from Bio import AlignIO
from Bio import SeqIO
#%%

#%%
MSA = 30
ref = "hg38"
comp = "rheMac8"
msaway = str(MSA) + "way"

# the maf path
MAF_PATH = f"/dors/capra_lab/data/ucsc/hg38/multiz{msaway}/maf/"
OUTPATH = f'/dors/capra_lab/users/fongsl//'
#make a list of the MAF files separated by chromosome
MAF_PATH
maf_list = glob.glob(f"{MAF_PATH}chr*.maf")

maf_dict = {(i.split("/")[-1]).split(".")[0]: i for i in maf_list}

#%%
def get_percent_identity(refSeq, targetSeq):

    lenSeq = len(refSeq) # get the length of the sequence alignment.

    if targetSeq != None:

        aligner = Align.PairwiseAligner() # get the aligner
        aligner.mode = 'local'

        score = aligner.score(refSeq, targetSeq) # score how many positions match

        perc_ID = score/lenSeq # get the percent identity

    else: # Target sequence is not alignable

        perc_ID, score = -1, -1

    return perc_ID, lenSeq, score


#%%
for chr_num, MAF_FILE in maf_dict.items():

    outfile = f"{OUTPATH}_species_count_%s.bed" %(OUTPATH, chr_num)

    touch = "touch %s" %outfile

    os.system(touch)
    out_file = open(outfile, 'w')

    maf = AlignIO.parse(MAF_FILE, "maf")
    #print(maf)
    count = 0

    for block in maf:
        count += 1
        store_homolog_line = [] # write the identity of the homolog

        s_count = 0

        s_list = []
        #print(block)
        for row in block:
            #print(row.seq)
            ### this is where the parsing happens.
            species = row.id.split('.')[0]

            start = row.annotations['start']

            end =  row.annotations['start'] + row.annotations['size']

            # ref sequence, strand, chromosome
            if species == ref:

                refSeq = row.seq # sequence

                blockchr = row.id.split('.')[1] # chr

                if row.annotations['strand'] == 1: # strand
                    strand = "+"
                elif row.annotations['strand'] == -1:
                        strand = "-"
                else:
                    raise ValueError('strand parsing did not work')

            # target sequence
            elif species == target:
                targetSeq = row.seq

            # count species in block
            if species != ref:
                s_count +=1
                s_list.append(species)

        # get sequence identity

        perc_ID, lenSeq, score = get_percent_identity(refSeq, targetSeq)
        new_vals = [blockchr, start, end, strand, species, perc_ID, lenSeq, score]

        store_homolog_line.append(new_vals)

        store_homolog_line.append([str(s_count)])
        store_homolog_line = [value for v in store_homolog_line for value in v]
        store_homolog_line.append(s_list)

        out_file.write('\t'.join(map(str,store_homolog_line)) + '\n')

    out_file.close()

        break
    break
#%%
