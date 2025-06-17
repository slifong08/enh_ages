import argparse
import configparser
from datetime import datetime
import glob
import itertools
import numpy as np
import os
import pandas as pd
import subprocess
import sys

sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
import config_readwrite as crw
import chr_functions

###
# arguments
###

arg_parser = argparse.ArgumentParser(description= "compute sequence identity between two species")

arg_parser.add_argument("-b","--bedfile", help='.bed file')
arg_parser.add_argument("-g","--genome_build", help='repeatmasker genome build (e.g. hg38)')
arg_parser.add_argument("-o","--outdirectory", help='out directory to save results')


args = arg_parser.parse_args()

BEDF = args.bedfile
BUILD = args.genome_build
OUTDIR = args.outdirectory
CHR_PATH = os.path.join(OUTDIR, "chr")


"""
INPUT = config["INPUTS"]["BEDFILE"]
DISTANCES = config["INPUTS"]["species_46way_hg19_distances"]

DATA_PATH = config["SEQ_ID"]["PATH"]
DATA_PATHS = config["MKPATHS"]["LIST"]

CHR_PATH = config["CHR"]["PATH"]

MAF_PATH = config["DORS_PATH"]["MSA_46WAY_hg19"]
MSA_SPLIT_BIN = config["SEQ_ID"]["MSA_SPLIT_BIN"]

"""

def load_constants(build):
    constants_dict = {"MSA_SPLIT_BIN":"/dors/capra_lab/bin/./msa_split",
                      "MAF_PATH":f"/dors/capra_lab/data/ucsc/{build}/maf", 

    }
    return constants_dict


def make_paths(path):
 
        if os.path.exists(path) is False:

            os.mkdir(path)
            print("made", i, path)

def split_chr(tiling, chr_path):
    chr_functions.split_into_chr_bed(tiling, chr_path)

def msa_split_x_bed(chr_, path, msa_split_bin, maf_path, chr_raw_path):
    
    """
    1. write chr info to a sequence identity file. 
    2. prepare input chrN.bed with full path
    3. change dir to write dir
    4. prepare input chrN.maf
    5. check if chrN.maf is zipped. (5A) if zipped, unzip before proceeding. 
    6. write feature arg string. This will split the maf file based on .bed coordintes!
    7. arg for writing to dir
    8. write MSA split command string
    9. check whether the chrN.bed file (n_lines) has already been split into n_lines*.fa files
    10.run MSA split command if not already split into n_lines*.fa files
    11. return list of split .fa files for chrN.bed
    
    """
    print(chr_)
    
    #1
    outf = os.path.join(path, f"{chr_}_seq_identity.tsv")  # chr-seq identity file to write to
    
    #2
    chrF = f"{chr_}.bed"  # regions to do msasplit on. 
    chr_bed = os.path.join(chr_raw_path, chrF)  # with full path. 
    
    #3
    os.chdir(path)
    
    #4
    maf_arg = os.path.join(maf_path, f"{chr_}.maf")
    
    #5
    zipped = maf_arg + ".gz"
    if os.path.exists(zipped) is True:
        cmd = f"gunzip {zipped}"
        print("\ngunzip\n", cmd, "\n")
        subprocess.call(cmd, shell = True)
    
    #6
    feat_arg = f"--features {chr_bed} --for-features"  # --for-features will split on each row of bed file. 
    
    #7
    out_root_arg = f"--out-root {chr_}"
    
    #8
    cmd = f"{msa_split_bin} {maf_arg} --in-format MAF {feat_arg} {out_root_arg}"
    
    #9
    already_split = len(glob.glob(f"{chr_}*.fa"))
    n_lines = sum(1 for line in open(chr_bed))
    
    #10
    if already_split !=n_lines:
        print(cmd)
        subprocess.call(cmd, shell = True)
        
    else:
        print("split .bed to .fa already")

    #11
    #already_split = glob.glob(f"{chr_}*.fa")
    
    #return already_split

def get_percent_identity(subjSeq, querySeq):

    lenSeq = len(subjSeq) # get the length of the sequence alignment.

    count_identical = 0
    count_gap = 0
    count_non_identical = 0
    
    # parse through sequence and ask if alignments match. 
    for a,b in zip(subjSeq,querySeq):

        if a==b:
            count_identical+=1  # count identical bases

        elif a != b:
            count_non_identical +=1  # count non-identical bases
            
        if a == "-" or b == "-":
            count_gap +=1  # count gap bases
            
    percent = count_identical/lenSeq  # return percent identity

    return count_identical, count_gap, percent

def id_age_dict(input_bed):
    """
    make dict[region_id]:mrca_age
    1. load file w specific columns, rename col names
    2. adjust start to add 1. idk why msa_split does this. 
    3. create id for region to match .fa file id
    4. make dictionary of id:age
    """
    #1
    cols = ["#chr", 'start', "end", "mrca"]
    srcAge = pd.read_csv(input_bed, sep = '\t', header = None, usecols =[0,1,2,10], names = cols)
    #2
    srcAge["start"] = srcAge["start"]+1
    #3
    srcAge["id"] = srcAge["#chr"]+":"+srcAge["start"].map(str) + "-" + srcAge["end"].map(str)
    #4
    id_dict = dict(zip(srcAge["id"], srcAge["mrca"].round(3)))

    return id_dict

def load_species_distances(distance):
    
    # make dictionary of species names and distances. 
    
    distances = pd.read_csv(distance)
    distances["taxon2"]= distances["taxon2"].apply(lambda x: x.replace("'", ""))
    dist = dict(zip(distances["taxon2"], distances["taxon1_mrca_dist"].round(3)))

    return dist

def extract_pairwise_fa_data(fa_handle, dist, id_dict, chr_):
    region_id = chr_+":"+(fa_handle.split("chr")[1].split(".")[1])
    age = id_dict[region_id]
    print(region_id)
    
    with open(fa_handle, "r") as fa_reader:
        """
        (1) set empty values for collecting species' sequence and sequence size
        (2) if species is human, set species variable to human
        (3) else, set species variable to rhesus
        (4) if neither hg38 or rheMac10 annotation, use species variable to recore sequence, size
        """
        
        #(1)

        hg19_seq, current_seq,  = "", "",  # keep track of hg19 and current sequences

        husize, current_size,  = 0, 0

        current_species = None  # keep track of current species

        current_mrca = 0  # keep track of the oldest MRCA age

        seqid=0

        n = 0

        ref_species = [
            "panTro2",
            "gorGor1",
            "ponAbe2",
            "rheMac2",
            "calJac1",
            "tarSyr1",
            "otoGar1",# mouse lemur on same MRCA. 
            "tupBel1",
            "mm9",
            "canFam2",
            "loxAfr3",
            "monDom5",
            "ornAna1","galGal3",
            "xenTro2",
            "danRer6",
            "petMar1"
            ]

        for i, line in enumerate(fa_reader):

            if ">" in line:

                # check whether last species has alignment
                if current_mrca == age and current_species in ref_species:
                    if list(set(current_seq)) != ['*']:

                        # get sequence identity
                        count_identical, count_gap, percent = get_percent_identity(hg19_seq, current_seq)

                        seqid=percent
                    else:
                        seqid=-1


                # begin new species
                current_species = ((line.split(" ")[1]).split("\n")[0]) # get the species
                current_seq = ""  # reset current sequence
                current_mrca = dist[current_species]  # get current age

                n+=1

            # get hg19 sequence
            elif current_species == "hg19":
                hg19_seq += line.strip("\n")  # add to hg19seq str
                husize += len(line.strip("\n"))

            # get other species sequence

            elif current_species in ref_species:
                current_seq += line.strip("\n")
                current_size += len(line.strip("\n"))

        # make dataframe of ages and sequence identity scores
        newdf = pd.DataFrame({
                            "region_id":[region_id],
                            "MRCA":[age],
                            "percent_id":[seqid]
                        })

    return newdf, region_id

def quantify_fa_seqids(chr_, data_path, distance, id_dict):
    
    out = os.path.join(data_path, f"{chr_}_seqid.tsv")
    fa_handles = glob.glob(os.path.join(data_path, f"{chr_}.*.fa"))
    print("len_fa_handles", len(set(fa_handles)))
    #if os.path.exists(out) is False and len(fa_handles)>0:

    #fa_handles = glob.glob(os.path.join(data_path, f"{chr_}*.fa"))

    results = {}

    dist = load_species_distances(distance)

    for f in fa_handles:
        seqid_df, region_id= extract_pairwise_fa_data(f, dist, id_dict, chr_)
        results[region_id] = seqid_df

    # concat sequence ids
    print(len(results.keys()))
    outdf = pd.concat(results.values())

    # save file
    outdf.to_csv(out, sep = '\t', index = False)

    for f in fa_handles:
        os.remove(f)
    """
    #elif os.path.exists(out) is True and len(fa_handles)>0:
        for f in fa_handles:
            os.remove(f)
    
    elif os.path.exists(out) is False and len(fa_handles)==0:
        print("need to do", chr_, "\n\n")
    else:
        print("already calculated sequence id for", chr_, "\n\n")
    """
    return out, outdf



def main(argv):
    dors = constants_dict(BUILD)

    MSA_SPLIT_BIN, MAF_PATH = dors["MSA_SPLIT_BIN"], dors["MAF_PATH"]
    
    chrList = chr_functions.make_chr_list()  # get chromosomes

    """
    (0) Make output paths
    
    (1) split file by chromosome number

    (2) perform msa_splits in parallel per chromosome 
    """

    #(0) make paths for results
    make_paths(CHR_PATH)

    #(1) split input file into chromosomes
    split_chr(BEDF, CHR_PATH)
    
    for chr_ in chrList:
        out = os.path.join(OUTDIR, f"{chr_}_seqid.tsv")
        if os.path.exists(out) is True:

            chrList.remove(chr_)
            
    print(chrList)

    
    #(2) perform msa_splits
    
    for chr_ in chrList:

        msa_split_x_bed(chr_, OUTDIR, MSA_SPLIT_BIN, MAF_PATH, CHR_PATH)
        id_dict = id_age_dict(BEDF)
   
        out, outdf = quantify_fa_seqids(chr_, OUTDIR, DISTANCES, id_dict)
        print(f"\n\n finished splitting {chr_}\n\n")

        
if __name__ == "__main__":
    main(sys.argv[1:])