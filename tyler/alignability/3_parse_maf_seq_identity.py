from Bio import AlignIO
from Bio import Align
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

def get_percent_identity(subjSeq, querySeq):

    lenSeq = len(subjSeq) # get the length of the sequence alignment.
    subjSeq = subjSeq.upper()
    querySeq = querySeq.upper()
    if querySeq != None:

        aligner = Align.PairwiseAligner() # get the aligner

        '''

        Pairwise aligner will apply
        Needleman-Wunsch, Smith-Waterman, Gotoh, and Waterman-Smith-Beyer global
        and local pairwise alignment algorithms
        to find the best-scoring between two sequences

        could also specify the aligner with:

        #aligner.mode = 'local'

        '''

        score = aligner.score(subjSeq.upper(), querySeq.upper()) # score how many positions match

        '''
        score - assign 1 to nucleotide matches.

        percent identity hg38 = score / hg38 sequence length
        percent identity rheMac8 = score / rheMac8 sequence length
        '''


    else: # Target sequence is not alignable

        perc_ID, score = -1, -1

    return score

def parse_rows(block):


    coor =[]

    for row in block:

        species = row.id.split('.')[0]

        seq = row.seq
        #print(species, start, end, row.annotations['size'])

        if species == "hg38":

            chr = row.id.split('.')[1]

            start = row.annotations['start']

            end =  row.annotations['start'] + row.annotations['size']

            subjSeq = seq

            hg38_size = row.annotations['size']

            coor.extend([chr, start, end, hg38_size])

        elif species == "rheMac8":
            rh_size = row.annotations['size']

            querySeq = seq

            coor.append(rh_size)

        else:
            print("problems")

    if len(coor)<5:
        print("No rhe alignment")
        querySeq = "".join(list(["-"]*hg38_size))
        rh_size = 0
        coor.append(rh_size)

    return coor, subjSeq, querySeq

def make_region_df(info, h, r, husize, rhsize, score):

    chr, locus_start, locus_end = info[0], info[1], info[2]

    df = pd.DataFrame({
    "#chr" :[chr],
    "start":[locus_start],
    "end":[locus_end],
    "hg38_Seqlen": [husize],
    "rheMac8_Seqlen": [rhsize],
    "hg38_seq":["".join(list(h))],
    "rheMac8_seq":["".join(list(r))],
    "score":[score],
    "per_id_hg38":[(score/husize)],
    "per_id_rheMac8":[(score/rhsize)]
    })

    return df

def make_block_df(results_dict, chr, path):

    # concat the dictionary
    re = pd.concat(results_dict.values()).drop_duplicates()
    re = re[[
            '#chr', 'start', 'end',
            'hg38_Seqlen', 'rheMac8_Seqlen',
            'n_matched_bases', 'identity',
            'hg38_seq', 'rheMac8_seq'
            ]]
    # write an outfile
    outf = f"{chr}_seq_identity.tsv"
    out = os.path.join(path, outf)

    # write the file
    re.to_csv(out, sep = '\t', index = False)

    return re

#%%

chrlist = make_chr_list()

for CHR in chrlist:
    outf = f"{PATH}{CHR}_seq_identity.tsv"
    if os.path.exists(outf) == False or os.path.getsize(outf) ==0:
        maf_f = f"{PATH}{CHR}_parse.maf"

        maf = AlignIO.parse(maf_f, "maf") # open the maf file

        results_dict ={} # dictionary to collect results for the chromosome

        locus_start, locus_end = 0, 0 # record the start, end of regulatory region

        last_end, val = 0,0  # value for recording the number of regulatory regions,
        husize, rhsize = 0, 0 # empty values to add size of contiguous msa blocks
        h, r = "", "" # empty string to collect sequence records

        for n, block in enumerate(maf):

            info, subjSeq, querySeq = parse_rows(block)

            chr, start, end = info[0], info[1], info[2]

            if last_end == 0: # for the first block in the msa file

                last_end = end # create a variable to store and compare last end with next start
                # if the last_end == next start, we're piecing together contiguous msa blocks

                locus_start = start # set the start of the region

                h += subjSeq # add the human sequence
                r += querySeq # add the rhesus sequence

                husize += info[3] # add the size of the human sequence to the contiguous sequence
                rhsize += info[4] # add the size of the rhe sequence to the contiguous sequence

                print("START", info)

            elif last_end == start: # for every other block in the msa file
                last_end = end # update the last end value

                h += subjSeq # add the human sequence
                r += querySeq # add the rhe sequence
                husize += info[3] # add the size of the human sequence to the contiguous sequence
                rhsize += info[4] # add the size of the rhe sequence to the contiguous sequence

                print("CONTINUE", info)

            elif last_end != start: # if the last_end != next start, we need to
            # (1) summarize the last sequence region.
            #(2) reset last end, locus start, h, r, husize, rhsize for next region

                locus_end = last_end
                score = get_percent_identity(h, r)
                region_info = [chr, locus_start, locus_end]
                df = make_region_df(info, h, r, husize, rhsize, score)

                results_dict[n] = df

                print("new locus!")
                print("NEW LOCUS START", info)
                # after
                last_end = end # update the last end value
                locus_start = start

                del h, r, husize, rhsize
                h = subjSeq # add the human sequence
                r = querySeq # add the rhe sequence
                husize = info[3] # add the size of the human sequence to the contiguous sequence
                rhsize = info[4] # add the size of the rhe sequence to the contiguous sequence

                val +=1
                print(val)

            #if val ==5:
            #    break
        re = make_block_df(results_dict, CHR, PATH)
    else:
        print("made this already", CHR)
