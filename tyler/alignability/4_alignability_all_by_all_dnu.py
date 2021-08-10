from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import Align
from Bio import SeqIO
import numpy as np
import os, sys
import pandas as pd
import subprocess

PATH = "/dors/capra_lab/users/fongsl/tyler/data/alignment/"


#%%
def make_chr_list():
    n = list(np.arange(1, 23))
    #n.append("X")

    chr_list = []
    for num in n:
        chrn = "chr" + str(num)
        chr_list.append(chrn)

    return chr_list

def get_score(subjSeq, querySeq, lenSeq):

    subjSeq = subjSeq.upper()
    querySeq = querySeq.upper()

    aligner = Align.PairwiseAligner() # get the aligner
    alignment = aligner.align(subjSeq, querySeq)
    print(alignment[0])
    '''

    Pairwise aligner will apply
    Needleman-Wunsch, Smith-Waterman, Gotoh, and Waterman-Smith-Beyer global
    and local pairwise alignment algorithms
    to find the best-scoring between two sequences

    could also specify the aligner with:

    #aligner.mode = 'local'

    '''

    score = aligner.score(subjSeq.upper(), querySeq.upper()) # score how many positions match

    mismatch_score = (lenSeq-score)
    return mismatch_score

def get_seq_id_len(n, seq_list):

    seq = seq_list[n] # get the corresponding sequence at index
    key_lenSeq = len(seq)
    id = CHR + "-" + str(n) # make an id
    val = (id, seq) # make a val

    return key_lenSeq, val

def fill_dict(key, val, seq_dict):
    if key in seq_dict.keys():

        #print("similar length!", key)

        id_list = seq_dict[key] # append sequence id to list
        id_list.append(val)
        seq_dict[key] = id_list

    else:
        seq_dict[key] = [val]

    return seq_dict

#%% make a dictionary of hg38 sequences.

chr_list = make_chr_list()
hseqs = {}
rhseqs = {}
val = 0

for CHR in chr_list:


    outf = f"{PATH}{CHR}_seq_identity.tsv"

    df = pd.read_csv(outf, sep = '\t', usecols = ["hg38_seq", "rheMac8_seq"])


    h_ = list(df.hg38_seq) # make a list of the sequences
    r_ = list(df.rheMac8_seq) # make a list of the sequences


    for n in np.arange(len(df)):

        key_lenSeq, val = get_seq_id_len(n, h_)

        hseqs = fill_dict(key_lenSeq, val, hseqs)

        key_lenSeq, val = get_seq_id_len(n, r_) # make the rhesus dictionary

        rhseqs = fill_dict(key_lenSeq, val, rhseqs)

#%%
len(hseqs.keys())
len(rhseqs.keys())
#%%
for i in rhseqs.values():
    #print((i))
    print(len(i))

#%%

from random import choice

def String(length):

   DNA=""
   for count in np.arange(length):

      DNA+=choice("CGTA")
   return DNA

def Mutator(DNA,lenDNA):

    mutpos = choice(np.arange(lenDNA)) # pick a position to mutate
    nuc = DNA[mutpos] # get the nucleotide at that position
    nucs = ["C", "G", "T", "A"]
    nucs_ = set(nucs).difference(set(nuc)) # make  set of mutant nucleotides

    mutnuc = choice("".join(list(nucs_))) # choose a mutant nucleotide
    mutDNA = list(DNA) # make a mutant DNA list
    mutDNA[mutpos] = mutnuc # replace the position of the DNA
    mutDNA = "".join(mutDNA) # turn mutant DNA list into string

    return mutDNA

#%%
len = 160
pos = String(len)
pos1 = Mutator(pos, len)
pos2 = Mutator(pos1, len)
pos3 = Mutator(pos2, len)

#%%

#%%
rheseqs2 = {
160: [("a", pos), ("b", String(160))]
}
hseqs2 = {
160: [("a", pos), # perfect match, but identity is the same
("b", String(160)), # not a match at all
("c", pos), # perfect match, different identity
("d", pos1), # 1 nucleotide off, same length
("e", pos2), # 2 nucleotide off, same length
("f", pos3), # 2 nucleotide off, same length
]
}

#%%

potential_matches = {} # a dictionary of the potential sequence matches

for lenSeq, seq_list in rheseqs2.items(): # for each sequence in the rhesus sequence list

    hinfo = hseqs2[lenSeq] # get the human sequences with the matching sequence length.

    for tuple in seq_list: # for each rhesus id, sequence tuple

        rseq_id, rseq = tuple[0], tuple[1]

        for h_tuple in hinfo: # test each h_tuple (id, seq) with the same length.

            hseq_id, hseq = h_tuple[0], h_tuple[1]

            if rseq_id !=hseq_id: # if the rhesus seq id and the human seq id are not the same

                # if the first, second, or third nucleotides match
                #if rseq[0] == hseq[0] and rseq[-1] == hseq[-1]:
                    print(rseq_id, hseq_id,)
                    # calculate the mismatch_Score
                    mismatch_score = get_score(hseq, rseq, lenSeq)
                    print("mismatch score =", mismatch_score, "\n\n")

                    if mismatch_score <= 2: # if there are <= than two mismatches
                        print("close enough", rseq_id, hseq_id)

                        if rseq_id in potential_matches.keys(): # add new human sequence to the list
                            id_list = potential_matches[rseq_id]

                            id_list.append(hseq_id)
                            potential_matches[rseq_id] = id_list

                        else:
                            # or make a new entry and add the human sequence to the list
                            # add new human sequence to the list
                            potential_matches[rseq_id] = [hseq_id]
#%%


#%%
for key, value in potential_matches.items():
    print(key, value)
#%%


np.unique(h, return_counts = True)
