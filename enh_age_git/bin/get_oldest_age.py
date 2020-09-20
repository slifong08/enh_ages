import os
import sys, traceback
import argparse
import csv
from datetime import datetime
from functools import reduce
import glob
from itertools import groupby
import numpy as np
import pandas as pd
from functools import partial
from multiprocessing import Pool
import subprocess

from joblib import Parallel, delayed
import multiprocessing

# TO RUN

# python script input_file sample_id

# python /dors/capra_lab/fongsl/enh_age/bin/age_enhancers.py UBERON_0002372_tonsil_expressed_enhancers.bed UBERON0002372


###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate enhancer age.")

arg_parser.add_argument("region_file_1", help='bed file 1 (enhancers to age) w/ full path')


arg_parser.add_argument("-s", "--species", type=str, default='hg19', choices=['hg19', 'hg38', 'mm10'],
                        help='species and assembly; default=hg19')

arg_parser.add_argument("-n", "--num_threads", type=int, default=100,
                        help='number of threads; default=100')

arg_parser.add_argument("-a", "--age", type=int, default=1,
                        help="age sequence w/ syntenic blocks")

arg_parser.add_argument("-o", "--oldest", type=int, default=1,
                        help="assemble breaks from aged sequence w/ syntenic blocks")

arg_parser.add_argument("-rt", "--run_testbed", type=int, default=1,
                        help="shuffle and calculate ages/breaks for input file")

arg_parser.add_argument("-sh", "--shuffles", type=int, default=0,
                        help="shuffle and calculate ages/breaks for input file")

args = arg_parser.parse_args()

# CONSTANTS
TEST_ENH = args.region_file_1
SPECIES = args.species
RUN_TEST_ENH = args.run_testbed
ANALYZE_AGE = args.age
ANALYZE_OLDEST = args.oldest
SHUFFLES = args.shuffles

if args.num_threads:
    NUM_THREADS = args.num_threads
else:
    NUM_THREADS = 100
print(ANALYZE_AGE, ANALYZE_AGE, RUN_TEST_ENH, "\nnum_threads", NUM_THREADS)

# EXTRACT OTHER CONSTANTS
TEST_PATH = "/".join(TEST_ENH.split("/")[:-1])+"/"
print("TEST PATH", TEST_PATH)
SAMPLE_ID = (TEST_ENH.split("/")[-1]).split(".")[0]

AGE_OUTFILE ="%s%s_enh_ages.bed" %(TEST_PATH, SAMPLE_ID)
SHUFFLE_PATH = "%sshuffle/" % TEST_PATH # make a shuffle path file
print("SHUFFLE PATH", SHUFFLE_PATH)
###
#   functions
###

def loadConstants(species):  # note chrom.sizes not used in current implementation | 2018.10.29
    return {'hg19': ("/dors/capra_lab/users/bentonml/data/dna/hg19/hg19_blacklist_gap.bed", "/dors/capra_lab/data/dna/human/hg19/hg19_trim.chrom.sizes"),
            'hg38': ("/dorvs/capra_lab/users/bentonml/data/dna/hg38/hg38_blacklist_gap.bed", "/dors/capra_lab/data/dna/human/hg38/hg38_trim.chrom.sizes"),
            'mm10': ("/dors/capra_lab/users/bentonml/data/dna/mm10/mm10_blacklist_gap.bed", "/dors/capra_lab/data/dna/mouse/mm10/mm10_trim.chrom.sizes")
            }[species]

def mkdir(path):
    if os.path.isdir(path) == False:
        cmd = "mkdir %s" % path
        subprocess.call(cmd, shell = True)
        
def os_remove(files):
    cmd = "rm %s" % files
    subprocess.call(cmd, shell = True)

def stripTabs(infile, tempfile):
    cmd = '''tr -s " \t" < %s > %s && mv %s %s''' % (infile, tempfile, tempfile, infile)
    subprocess.call(cmd, shell = True)
    
    rm_cmd = "rm %s" % tempfile
    subprocess.call(cmd, shell = True)
    
def enh_age(test_enh, sample_id, test_path, species):
    
    print("aging", sample_id)
    outpath = "%sages/" % test_path # mkdir ./ages/
    mkdir(outpath)
    
    os.chdir(outpath)
    
    chr_cmd = '''awk '{print >$1"_%s_temp.bed"}' %s''' % (sample_id, test_enh) # split test into chrN.bed files
    subprocess.call(chr_cmd, shell = True)
    
    
    enh_chr_list = glob.glob("%schr*_%s_temp.bed" % (outpath, sample_id)) # glob chromosomes
    sex_chr = ["chrX", "chrY", "chrM", "chr"]
    
    for f in enh_chr_list:
        chr_num = (f.split("/")[-1]).split("_")[0]
        if chr_num not in sex_chr: # filter out sex chromosomes. 
            print(chr_num)
            
            syn_path = "/dors/capra_lab/data/ucsc/%s/synteny_age_bkgd_%s/" %(species, species)
            syn_file= ("%s%s_syn_age.bed" % (syn_path, chr_num))
            
            outfile = "%s%s_%s_ages.bed" %(outpath, chr_num, sample_id)
            
            BEDcmd = "bedtools intersect -a %s -b %s -wo > %s"  % (f, syn_file, outfile)
            subprocess.call(BEDcmd, shell = True)
            
            # write any lines where enhancer does not overlap syntenic break, does not have mrca. 
            no_syn = "%s%s-%s_no_syn_alignment_tiles.txt" % (outpath, chr_num, sample_id) # write non-overlapping regions

            NOSYNcmd = "bedtools intersect -a %s -b %s -v -wo > %s"  % (f, syn_file, no_syn)
            subprocess.call(NOSYNcmd, shell = True)  
            
        else:
            os_remove(f)
            
    concat_file = "%s%s_enh_ages.bed" % (outpath, sample_id)        
    concat_cmd = "cat %schr*_%s_ages.bed > %s" % (outpath, sample_id, concat_file)

    print("concatenating aged enhancer chromosome files")
    subprocess.call(concat_cmd, shell = True) 
  
    syn_concat_file = "%s%s_no_syn_alignment_tiles.txt" % (outpath, sample_id)        
    syn_concat_cmd = "cat %schr*-%s_no_syn_alignment_tiles.txt > %s" % (outpath, sample_id, syn_concat_file)#concatenate chromosomes

    subprocess.call(syn_concat_cmd, shell = True) 
    os_remove("%schr*-%s_no_syn_alignment_tiles.txt"%(outpath, sample_id))
    os_remove("%schr*_%s_*.bed" %(outpath, sample_id))       
    return concat_file


def assemble_breaks(mrca_list, test_id):
    
    ages, arch = np.unique(mrca_list, return_inverse = True) # get the unique ages, architecture order
    arch = [k for k, g in groupby(arch)] # get the architecture w/ itertools groupby

    # create core_remodeling binary
    if len(ages) ==1: 
        core_remodeling = 0
    elif len(ages)>1:
        core_remodeling = 1

    newdf = pd.DataFrame({"test_id": [test_id],
                      "max_seg": [len(arch)], 
                      "core_remodeling": [core_remodeling],
                      "arch":[arch],
                      "max_age": [max(ages).round(3)]})
    return newdf


def format_df(df, agepath, sid):    

    #drop axis
    df = df.drop(["strand","ref","num_species", "patr"], axis = 1)

    # drop syn_blocks > 5 bp in len
    df = df[df.len_syn_overlap > 5]

    # add a test_id
    df["test_id"] = df["chr_enh"] + ":" + df["start_enh"].map(str) + "-" + df["end_enh"].map(str) 

    # figure out difference between enh and syn start
    # positive = syn block inside enh
    # negative = syn block starts before enh - needs to be trimmed!
    df["startdif"] = df.start_enh - df.start_syn

    # trim the start positions (when syn start is before enh start)
    df.loc[df.startdif>0, "start_syn"] = df.start_enh.copy() 

    # figure out difference between enh and syn start
    # positive = syn inside enh
    # negative = syn ends after enh - needs to be trimmed!
    df["enddif"] = df.end_enh - df.end_syn

    # trim the end positions (when syn end is after enh end)
    df.loc[df.enddif<0, "end_syn"] = df.end_enh.copy() 

    collection_dict = {} # collect unique enhancer architectures, max ages, core_remodeling

    # iterate through test_ids
    for test_id in df.test_id.unique():

        if test_id not in collection_dict.keys():

            mrca_list = df.loc[df.test_id==test_id, "mrca"].to_list() # get a list of ages.
            newdf = assemble_breaks(mrca_list, test_id) # assemble the breaks!!!
            collection_dict[test_id] = newdf
    
    archdf = pd.concat(collection_dict.values())
    originaldf = df[["chr_enh", "start_enh", "end_enh", "id", "test_id"]].drop_duplicates()
    results = pd.merge(originaldf, archdf, how = "left")
    chrnum =results.chr_enh.unique().item()
    print('results! %s' % chrnum)
    results.to_csv("%s%s_%s_parallel.txt" % (agepath, chrnum, sid), sep ='\t', header = True, index =False)
    return results

### SPLIT BY CHROMOSOME, TEST PARALLEL ###
def run_parallel_breaks(F):
    chr_dict = {}
    df = pd.read_csv(F, sep ='\t', header = None)

    df.columns = ["chr_enh", "start_enh", "end_enh", "id",  
                            "chr_syn", "start_syn","end_syn","strand",
                            "ref","num_species", "len_syn","mrca",
                            "patr","len_syn_overlap"]

    for chrnum in df.chr_enh.unique():
        chrdf = df.loc[df.chr_enh ==chrnum]
        chr_dict[chrnum] = chrdf
    print(chr_dict.keys())

    ### RUN PARALLEL ###

    num_cores = multiprocessing.cpu_count()
    num_cores = len(chr_dict.values())

    # print start time
    start = datetime.now()
    print("start parallel", start)

    # run parallel jobs
    results = Parallel(n_jobs=num_cores, verbose=100, prefer="threads")(delayed(format_df)(i, agepath, sid) for i in chr_dict.values())

    # how long did the run take?
    end = datetime.now()
    print("end parallel", datetime.now())
    print("total time", end - start)

    catf = "%sshuf-trimmed-310_%s_parallel_breaks.bed" %(agepath, (sid.split("-")[1]).split(".")[0])
    cmd = "cat %schr*_%s_parallel.txt > %s" % (agepath, sid, catf)
    os.system(cmd)

    rm_cmd = "rm %schr*_%s_parallel.txt" % (agepath, sid) # remove all the parallel files
    os.system(rm_cmd)

    rmF_cmd = "rm %s" % (F) # remove the input file
    os.system(rmF_cmd)

def get_oldest(concat_file, sample_id, test_path, num_threads):
    
    inpath = "%sages/" % test_path # mkdir ./ages/
    outpath = "%sbreaks/" % test_path # mkdir ./break/
    mkdir(outpath) 

    oldest_out = "%s%s_oldest_ages.bed" % (outpath, sample_id)
    print(oldest_out)
    
    age_temp = "%stemp_%s.bed" % (inpath, sample_id)
    
    stripTabs(concat_file, age_temp) # strip any double tabs

    run_parallel_breaks(concat_file) # get the breaks
    
def formatBedfile(test_enh, out_test_enh, sample_id):
    
    cmd = '''awk 'OFS=" " {print $1"\t", $2"\t", $3"\t", $4}' %s | tr -d " "| sort -k1,1 -k2,2 -k3,3 > %s''' % (test_enh, out_test_enh) 
    print("standardizing Bed format") 
    subprocess.call(cmd, shell=True)
    print(out_test_enh)
    
    df = pd.read_csv(out_test_enh, sep ='\t', header = None)
    df[3] = sample_id # label the interation
    df.to_csv(out_test_enh, sep = '\t', header = False, index = False)

def shuffleBedfile(f, species, sample_id, path, iteration):
    
    bl, g = loadConstants(species)  # note CHROM_SZ not used
    outf = "%sshuf-%s-%s.bed" % (path, sample_id, iteration) # outfile

    #bedtools cmd hg19, exclude blacklist, same chromosome, noOverlapping bits
    cmd = "bedtools shuffle -i %s -g %s -excl %s -chrom -noOverlapping > %s" % (f, g, bl, outf) 

    os.system(cmd)
    
    # label the shuffle file
    df = pd.read_csv(outf, sep ='\t', header = None)
    df[3] = df[3] + "-shuf-"+ str(iteration) # label the interation
    #df = df.sample(frac =0.1) # sample 10% of the dataframe x 100 times
    df.to_csv(outf, sep = '\t', header = False, index = False)

    return outf

###
#   main
###

def main(argv):
    print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.now())[:20]))
    
    ### standardize the enhancer file ###
    TEST_ENH_CUT = "%scut-%s.bed" % (TEST_PATH, SAMPLE_ID) 

    
    if RUN_TEST_ENH ==1:

        if ANALYZE_AGE == 1:
            formatBedfile(TEST_ENH, TEST_ENH_CUT, SAMPLE_ID) # format the enhancer bed file and sort
            age_file = enh_age(TEST_ENH_CUT, SAMPLE_ID, TEST_PATH, SPECIES)
            
        if ANALYZE_OLDEST == 1:
            if ANALYZE_AGE==1:
                get_oldest(age_file, SAMPLE_ID, TEST_PATH, NUM_THREADS)
                os_remove(age_file)
                
            else:
                AGE_F = "%s%s.bed" % (TEST_PATH, SAMPLE_ID)    
                get_oldest(AGE_F, SAMPLE_ID, TEST_PATH, NUM_THREADS)
                os_remove(AGE_F)

    if SHUFFLES >0:
        mkdir(SHUFFLE_PATH)

        formatBedfile(TEST_ENH, TEST_ENH_CUT, SAMPLE_ID) # format the enhancer bed file and sort
        
        shuffleList = []
        
        for i in range(SHUFFLES):
            
            new_shuffle = shuffleBedfile(TEST_ENH_CUT, SPECIES, SAMPLE_ID, SHUFFLE_PATH, i)
            shuffleList.append(new_shuffle) # make a list of the new shuffle files
        
        concat_file = "%sshuf-%s_concat.bed" % (SHUFFLE_PATH, SAMPLE_ID) # concat all the shuffles together.        
        concat_cmd = "cat %sshuf-%s-*.bed > %s" % (SHUFFLE_PATH, SAMPLE_ID, concat_file)
        os.system(concat_cmd)
        
        for shuffle in shuffleList: # remove the extra shuffles
            os_remove(shuffle)
        
        if ANALYZE_AGE == 1:

            age_file = enh_age(concat_file, SAMPLE_ID, SHUFFLE_PATH, SPECIES)
        
        if ANALYZE_OLDEST == 1:
            if ANALYZE_AGE==1:
                get_oldest(age_file, SAMPLE_ID, SHUFFLE_PATH, NUM_THREADS)
                os_remove(age_file)

            else:
                AGE_F = "%s%s.bed" % (TEST_PATH, SAMPLE_ID)    
                get_oldest(AGE_F, SAMPLE_ID, TEST_PATH, NUM_THREADS)
                os_remove(AGE_F)
        
    os_remove(TEST_ENH_CUT)
#    if "shuf" in SAMPLE_ID:
 #       os_remove(TEST_ENH)
    
if __name__ == "__main__":
    main(sys.argv[1:])    