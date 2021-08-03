from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
import os, sys
import subprocess

CHR = "chr21"
PATH = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/all/"


maf_f = "/dors/capra_lab/users/fongsl/tyler/data/alignment/chr21_parse.maf"
val = 0
maf = AlignIO.parse(maf_f, "maf")



for block in maf:
    align = []
    align2 = []
    last_end = 0
    align

    for row in block:
        #print(row.seq, row.annotations['size'])

        species = row.id.split('.')[0]

        start = row.annotations['start']

        end =  row.annotations['start'] + row.annotations['size']

        seq = row.seq
        print(species, start, end )
        #print(len(align))
        if species == "hg38":
            if len(align) == 0:
                #print("before", align)
                align.append(seq)
                #print("after", align)
            else:
                print("hg38", align)
                align = align + seq

        else:
            if len(align2) == 0:
                align2.append(seq)
            else:
                align2 = align2 + seq


        if last_end ==0 and species == "hg38":
            last_end = end
            print("new last end", last_end)

        elif start == last_end and species == "hg38": # when there is a continuous block
            last_end = end
            print("updated last end", last_end)

        elif start != last_end and species == "hg38": # when there is a continuous block
            print("STARTAGAIN", last_end)
            last_end = 0
            break
