#180601
#20180614 update

# sarahfong
# Tony Capra drafted this script. Sarah edited.
# biopython is loaded in the sf_test env with 'conda install -c bioconda biopython'

# the species count does not count hg19. Therefore, any human-specific region = 0

import sys, os
import glob

#SARAH- ADD THE PATH WITH BIOPYTHON
bio_package = '/home/fongsl/.conda/envs/sf_test/lib/python3.6/site-packages'
if bio_package not in sys.path:

    sys.path.append('/home/fongsl/.conda/envs/sf_test/lib/python3.6/site-packages')
sys.path

from Bio import AlignIO
from Bio import SeqIO

# path to data

DATAPATH = '/dors/capra_lab/data/ucsc/hg19/maf/'

OUTPATH = '/dors/capra_lab/data/ucsc/hg19/maf/synteny_age_hg19/'

#make a list of the MAF files separated by chromosome

maf_list = glob.glob("%schr*.maf"%DATAPATH)
maf_dict = {(i.split("/")[-1]).split(".")[0]: i for i in maf_list}

for chr_num, MAF_FILE in maf_dict.items():

# this grafted from Abin's script. Credits to him.

    # make a file for each chromosome
    outfile = "%shg19_species_count_%s.bed" %(OUTPATH, chr_num)

    touch = "touch %s" %outfile

    os.system(touch)
    out_file = open(outfile, 'w')

    maf = AlignIO.parse(MAF_FILE, "maf")
    count = 0

    for block in maf:
        count += 1
        store_homolog_line = []

        s_count = 0

        s_list = []
        for row in block:
            ### this is where the parsing happens.
            species = row.id.split('.')[0]

            blockchr = row.id.split('.')[1]

            start = row.annotations['start']

            end =  row.annotations['start'] + row.annotations['size']

            if row.annotations['strand'] == 1:
                strand = "+"
            elif row.annotations['strand'] == -1:
                    strand = "-"
            else:
                raise ValueError('strand parsing did not work')

            if species != "hg19":
                s_count +=1
                s_list.append(species)

            else:
                store_homolog_line.append([blockchr, start, end, strand, species]) # only write the row with the h19 information. Not all the species information.

        store_homolog_line.append([str(s_count)])
        store_homolog_line = [value for v in store_homolog_line for value in v]
        store_homolog_line.append(s_list)

        out_file.write('\t'.join(map(str,store_homolog_line)) + '\n')

    out_file.close()
