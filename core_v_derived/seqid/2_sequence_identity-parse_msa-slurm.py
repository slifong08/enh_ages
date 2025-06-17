#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import configparser
from datetime import datetime
import glob

import itertools
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import subprocess
import sys

sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
import config_readwrite as crw
import chr_functions


# In[ ]:


arg_parser = argparse.ArgumentParser(description="Calculate sequence identity")

arg_parser.add_argument("input", help='input bed')

args = arg_parser.parse_args()

# CONSTANTS

INPUT = args.input


# In[ ]:


#CHR = "chr22"
#INPUT = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/shuffle_syn_breaks/no-exon_shuf-all_fantom_enh_112_tissues-58_age_breaks.bed"

# handle naming dir strategy for inputs
if "shuf" in INPUT:
    print("shuf")
    INPUT_NAME = (INPUT.split("/")[-1]).split("_")[5]

else:
    print("target")
    INPUT_NAME = (INPUT.split("/")[-1]).split(".bed")[0]


# In[ ]:


# load the sequence identity config
NAME = os.path.join(os.getcwd(),"config_seqarch")
config, configfile_name = crw.read_config(NAME)  # config is an object, configfile_name is a string


# In[ ]:


#INPUT = config["INPUTS"]["BEDFILE"]
DISTANCES = config["INPUTS"]["species_46way_hg19_distances"]

DATA_PATH = config["SEQ_ID"]["PATH"]+"/"+ INPUT_NAME # make a path for each input
DATA_PATHS = config["MKPATHS"]["LIST"]

CHR_PATH = config["CHR"]["PATH"] + INPUT_NAME # make a path for each input

DATA_PATHS = DATA_PATHS + "," + CHR_PATH + "," + DATA_PATH # add new input path

MAF_PATH = config["DORS_PATH"]["MSA_46WAY_hg19"]
MSA_SPLIT_BIN = config["SEQ_ID"]["MSA_SPLIT_BIN"]


# In[ ]:


def make_paths(data_paths):
    for i, path in enumerate(data_paths.split(",")):
 
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


def load_species_distances(distance):
    
    # make dictionary of species names and distances. 
    
    distances = pd.read_csv(distance)
    distances["taxon2"]= distances["taxon2"].apply(lambda x: x.replace("'", ""))
    dist = dict(zip(distances["taxon2"], distances["taxon1_mrca_dist"]))
    
    return dist


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


def extract_pairwise_fa_data(fa_handle, dist, id_dict, chr_):
    
    region_id = chr_+":"+(fa_handle.split("chr")[1].split(".")[1])
    age = id_dict[region_id] # get age to query. 
    
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

        seqid = 0

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

                        # assign seq id variable
                        seqid = percent
                    else:
                        seqid = -1


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


    #print(region_id, age, seqid,"\n\n")
    # make dataframe of ages and sequence identity scores
    
    newdf = pd.DataFrame({
                        "region_id":[region_id],
                        "MRCA":[age],
                        "percent_id":[seqid]
                         })
    
    return newdf, region_id


def quantify_fa_seqids(chr_, data_path, distance, IdAgeDict):
    
    out = os.path.join(data_path, f"{chr_}_seqid.tsv")
    fa_handles = glob.glob(os.path.join(data_path, f"{chr_}.*.fa"))
    
    if os.path.exists(out) is False and len(fa_handles)>0:
        
        results = {}
        
        dist = load_species_distances(distance)
        
        for f in fa_handles:
            seqid_df, region_id= extract_pairwise_fa_data(f, dist, IdAgeDict, chr_)
            results[region_id] = seqid_df
        
        # concat sequence ids
        outdf = pd.concat(results.values())
        
        # save file
        outdf.to_csv(out, sep = '\t', index = False)
        """
        for f in fa_handles:
            os.remove(f)
            
    elif os.path.exists(out) is True and len(fa_handles)>0:
        for f in fa_handles:
            os.remove(f)
    """
    elif os.path.exists(out) is False and len(fa_handles)==0:
        print("\n\nneed to do", chr_, "\n\n")
        
    else:
        print("\n\nalready calculated sequence id for", chr_, "\n\n")

    return out


# # main 

# In[ ]:


def main(argv):

    chrList = chr_functions.make_chr_list()  # get chromosomes

    """
    (0) Make output paths
    
    (1) split file by chromosome number

    (2) perform msa_splits in parallel per chromosome 
    
    (3) compute seq_id
    """

    #(0) make paths for results
    make_paths(DATA_PATHS)

    IdAgeDict = id_age_dict(INPUT) # dictionary of sequence ids and ages. 
    
    for CHR in chrList:
        test_chr = os.path.join(CHR_PATH, f"{CHR}.bed")
        if os.path.exists(test_chr) is False:
            
            #(1) split input file into chromosomes
            split_chr(INPUT, CHR_PATH)

        out = os.path.join(DATA_PATH, f"{CHR}_seqid.tsv")
        fa_handles = glob.glob(os.path.join(DATA_PATH, f"{CHR}.*.fa"))
        if os.path.exists(out) is False:     
            
            if len(fa_handles) ==0:
                #(2) perform msa_splits
              
                msa_split_x_bed(CHR, DATA_PATH, MSA_SPLIT_BIN, MAF_PATH, CHR_PATH)

            #(3) quantify sequence identity for .fa files
            out = quantify_fa_seqids(CHR, DATA_PATH, DISTANCES, IdAgeDict)

            print(f"\n\n finished splitting {CHR}\n\n")
        else:
            print(f"already split {INPUT_NAME} {CHR}\n\n")

        
if __name__ == "__main__":
    main(sys.argv[1:])





