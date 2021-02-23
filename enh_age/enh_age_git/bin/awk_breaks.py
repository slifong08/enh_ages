import os
import sys, traceback
import argparse
import csv
from datetime import datetime
import glob
import numpy as np
import pandas as pd
from functools import partial

from multiprocessing.pool import Pool, TimeoutError


import subprocess


# TO RUN

# python script input_file sample_id

# python /dors/capra_lab/fongsl/enh_age/bin/age_enhancers.py UBERON_0002372_tonsil_expressed_enhancers.bed UBERON0002372


###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate enhancer age.")

arg_parser.add_argument("bedfile_to_age", help='bed file (enhancers to age) w/ full path')

arg_parser.add_argument("-i", "--iters", type=int, default=100,
                        help='number of simulation iterations; default=100')

arg_parser.add_argument("-s", "--species", type=str, default='hg19', choices=['hg19', 'hg38', 'mm10'],
                        help='species and assembly; default=hg19')

arg_parser.add_argument("-a", "--age", type=int, default=1,
                        help="age sequence w/ syntenic blocks")

arg_parser.add_argument("-b", "--breaks", type=int, default=1,
                        help="assemble breaks from aged sequence w/ syntenic blocks")

arg_parser.add_argument("-t", "--tfbs_den", type=int, default=0,
                        help="calculate tfbs density using 161 ChIP-Seq datasets in aged sequence w/ syntenic blocks")

arg_parser.add_argument("-sh", "--shuffle", type=int, default=0,
                        help="shuffle and calculate ages/breaks for input file")

arg_parser.add_argument("-rt", "--run_testbed", type=int, default=1,
                        help="shuffle and calculate ages/breaks for input file")



args = arg_parser.parse_args()

# CONSTANTS

TEST_ENH = args.bedfile_to_age
ITERATIONS = args.iters
SPECIES = args.species
RUN_TEST_ENH = args.run_testbed
ANALYZE_AGE = args.age
ANALYZE_BREAKS = args.breaks
ANALYZE_TFBS_DEN = args.tfbs_den
ANALYZE_SHUFFLE = args.shuffle

NUM_THREADS = 100


print("\n", ANALYZE_AGE, ANALYZE_BREAKS, ANALYZE_TFBS_DEN, ANALYZE_SHUFFLE, RUN_TEST_ENH)

# EXTRACT OTHER CONSTANTS
TEST_PATH = "/".join(TEST_ENH.split("/")[:-1])
SAMPLE_ID = (TEST_ENH.split("/")[-1]).split(".")[0]
print("\nSAMPLE_ID:",SAMPLE_ID)
print("\nFILE TO RUN:",TEST_ENH)


###
#   functions
###


def loadConstants(species):  # note chrom.sizes not used in current implementation | 2018.10.29
    return {'hg19': ("/dors/capra_lab/users/fongsl/data/hg19_blacklist_gap_ensemblexon.bed", "/dors/capra_lab/data/dna/human/hg19/hg19_trim.chrom.sizes"),
            'hg38': ("/dors/capra_lab/users/fongsl/data/hg38_blacklist_gap_ensemblexon.bed", "/dors/capra_lab/data/dna/human/hg38/hg38_trim.chrom.sizes"),
            'mm10': ("/dors/capra_lab/users/bentonml/data/dna/mm10/mm10_blacklist_gap.bed", "/dors/capra_lab/data/dna/mouse/mm10/mm10_trim.chrom.sizes")
            }[species]

# make directory
def mkdir(path):
    if os.path.isdir(path) == False:

        cmd = "mkdir %s" % path
        print("MKDIR", cmd)
        os.system(cmd)

# remove files
def os_remove(files):

    cmd = "rm %s" % files
    print("REMOVE FILE", cmd)
    os.system(cmd)

# remove extra tabs
def stripTabs(infile, tempfile):
    print("STRIPTABS")
    cmd = '''tr -s " \t" < %s > %s''' % (infile, tempfile)
    os.system(cmd)

    rename_cmd = "mv %s %s" % (tempfile, infile) # the infile is now stripped, no need to return file.
    os.system(rename_cmd)

# sort bedfile
def sort_bed(infile):
    sid = (infile.split("/")[-1]).split(".")[0]
    path = "/".join(infile.split("/")[:-1])
    temp = "%s%s_sorted.bed" % (path, sid)
    cmd = "sort -k1,1 -k2,2 -k3,3 %s > %s && mv %s %s" % (infile, temp, temp, infile)
    print("SORTING FILE")
    os.system(cmd)


# bedtools intersection w/ chr-specific syntenic block.
def bedintersect_syn(f, species, chr_num, sample_id, outpath):

    print("SYNTENY INTERSECTION", chr_num, f, species, chr_num, sample_id, outpath)
    syn_path = "/dors/capra_lab/data/ucsc/%s/synteny_age_%s/" %(species, species)
    syn_file= ("%s%s_syn_age.bed.gz" % (syn_path, chr_num))

    outfile = "%s/%s_%s_ages.bed" %(outpath, chr_num, sample_id)

    BEDcmd = "bedtools intersect -a %s -b %s -wao > %s"  % (f, syn_file, outfile)
    os.system(BEDcmd) # do the intersection


    return outfile

# get the enhancer ages
def age_enh(test_enh, sample_id, test_path, species):

    print("AGING", sample_id, test_enh, test_path, species)
    outpath = "%s/ages" % test_path # mkdir ./ages/
    mkdir(outpath)

    os.chdir(outpath)

    sex_chr = ["chrX", "chrM", "chrY"]


    # for files w/ regions from multiple chromosomes.
    if "chr" not in sample_id:


        chr_cmd = '''awk '{print >$1"_%s_temp.bed"}' %s''' % (sample_id, test_enh) # split test into chrN.bed files
        os.system(chr_cmd)


        enh_chr_list = glob.glob("%s/chr*_%s_temp.bed" % (outpath, sample_id)) # glob chromosomes

        for f in enh_chr_list:

            chr_num = (f.split("/")[-1]).split("_")[0]

            if chr_num not in sex_chr: # filter out sex chromosomes.
                sort_bed(f) # sort the file
                bedintersect_syn(f, species, chr_num, sample_id, outpath) # syntenic blcok intersection
                os_remove(f) # remove the chromosome file temp

            else:
                os_remove(f) # delete the sex chromosomes

        # concatenate all the chromosomes together again
        concat_file = "%s/%s_ages.bed" % (outpath, sample_id)
        concat_cmd = "cat %s/chr*_%s_ages.bed > %s" % (outpath, sample_id, concat_file)

        print("concatenating aged enhancer chromosome files")
        subprocess.call(concat_cmd, shell = True)

        os_remove("%s/chr*_%s_ages.bed" %(outpath, sample_id))

        return concat_file

    else: # for files w/ regions from one chromosome
        chr_num = ((test_enh.split("/")[-1]).split("_")[0]).split("-")[1]
        outfile = bedintersect_syn(test_enh, species, chr_num, sample_id, outpath) # syntenic blcok intersection

        print("nothing to delete - single chromosome analysis")

        return outfile

def break_scripts(age_file, sample_id, test_path):

    outpath = "%s/breaks/" % test_path

    if os.path.exists(outpath) == False:
        mkdir(outpath)


    # sort bed file
    print("sorting age file")
    sort_bed(age_file)

    # remove lines that do not overlap syntenic blocks, overlap X chromosome
    cleanup_file = "%sclean_ages_%s.bed" % (outpath, sample_id)
    cleanup = '''awk '!($5 ~ /[.]/ ) && !($1 ~ /chrX/ ) && ($14 > 5 ) { print $1, "\t", $2,  "\t",$3, "\t", $4, "\t", $12, "\t", $14 }' %s > %s''' % (age_file, cleanup_file)
    print("cleanup\n")
    subprocess.call(cleanup, shell = True)

    # write all the lines that do not overlap syntentic blocks.
    no_overlap_file = "%sno_syn_alignment_%s.txt" % (outpath, sample_id)
    no_overlap_cmd = '''awk '($5 ~ /[.]/ )' %s > %s''' % (age_file, no_overlap_file)
    print("no_overlap")
    subprocess.call(no_overlap_cmd, shell = True)


    # add enhancer id column
    temp = "%stemp_%s.bed" % (outpath, sample_id)
    add_enh_id = '''awk '{$(NF+1)=$1":"$2"-"$3 ; print}' %s > %s && mv %s %s''' % (cleanup_file, temp, temp, cleanup_file)
    print("add enh_id")
    subprocess.call(add_enh_id, shell = True)


    # add tabs to cleanup file
    tab_cmd = '''awk '{$1=$1}1' OFS="\t" %s > %s && mv %s %s''' % (cleanup_file, temp, temp, cleanup_file)
    print("tabs")
    subprocess.call(tab_cmd, shell = True)


    # get max mrca per enhancer from cleanup file
    mrca_file= "%smax_mrca_%s.bed" % (outpath, sample_id)
    mrca_cmd = '''awk '$5>max[$7]{max[$7]=$5; row[$7]=$0} END{for (i in row) print row[i]}' %s > %s '''% (cleanup_file, mrca_file)
    print("get max mrca")
    subprocess.call(mrca_cmd, shell = True)


    # index syntenic blocks
    syn_index_cmd = "python /dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/syntenic_assembly.py %s" % age_file
    print("syn arch")
    subprocess.call(syn_index_cmd, shell = True)

    # get number of segments per enhancer from outpath file
    seg_count_file = "%sage_seg_count_%s.bed" %(outpath, sample_id)
    seg_count_cmd = '''cut -f 5,7 %s | uniq | cut -f2 | uniq -c > %s''' % (cleanup_file, seg_count_file)
    print("count age segments")
    subprocess.call(seg_count_cmd, shell = True)


    # add tabs to seg_index_count
    tab_cmd = '''awk '{$1=$1}1' OFS="\t" %s > %s && mv %s %s''' % (seg_count_file, temp, temp, seg_count_file)
    print("more tabs")
    subprocess.call(tab_cmd, shell = True)


    ### reference MRCA file ###

    syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t',
    usecols = ["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]) # read the file

    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages


    # open up and merge a bunch of dataframes

    print("opening files")

    mrca_df = pd.read_csv(mrca_file, sep='\t', header = None)
    mrca_df.columns = ["chr_enh", "start_enh", "end_enh", "sample_id", "mrca", "syn_len", "enh_id"]
    mrca_df.sort_values(by = "enh_id")

    seg_df = pd.read_csv(seg_count_file, sep='\t', header = None)
    seg_df.columns = ["seg_index", "enh_id"]
    seg_df.sort_values(by = "enh_id")

    breaksdf = pd.merge(mrca_df, seg_df, how = "left", on = "enh_id")

    # round the mrca value
    breaksdf.mrca = breaksdf.mrca.round(3)

    # add binary for simple/complex architecture based on median break value
    breaksdf["core_remodeling"] = 0
    breaksdf.loc[breaksdf[ "seg_index"] >= breaksdf[ "seg_index"].median(), "core_remodeling"] = 1

    # add annotation based on simple/complex architecture
    breaksdf["arch"] = "simple"
    breaksdf.loc[breaksdf["core_remodeling"] ==1, "arch"] = "complexenh"

    # reorder columns
    breaksdf = breaksdf[["chr_enh", "start_enh", "end_enh", "enh_id", "sample_id",
                "seg_index", "core_remodeling", "arch", "mrca"]]

    # merge with other age, taxon annotations
    breaksdf = pd.merge(breaksdf,
    syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]],
    how = "left", on = "mrca")


    # save all this information.

    headerf = "%ssummary_matrix_header.txt" % (outpath)

    out_summarized_bed = "%s%s_enh_age_arch_summary_matrix.bed" % (outpath, sample_id)

    breaksdf.to_csv(out_summarized_bed, sep = '\t', index = False, header = False)

    if os.path.exists(headerf) ==False:

        cols = list(breaksdf.columns) # make the columns into a list
        np.savetxt(headerf, cols, delimiter="\t", fmt = '%s') # save the column headers

    count = len(open(out_summarized_bed).readlines(  )) # makesure the file wrote before deleting ages.

    if count >0:
        os_remove(cleanup_file)


# generate shuffles
def calculateExpected(test_enh, sample_id, shuffle_path, species, analyze_age, analyze_breaks, analyze_tfbs_den, iters):


    print("shuffling", shuffle_path, iters)

    BLACKLIST, CHROM_SZ = loadConstants(species)  # note CHROM_SZ not used

    exp_sum = 0

    if os.path.exists(shuffle_path) == False:

        mkdir(shuffle_path)

    shuffle_id = "shuf-%s-%s" % (sample_id, iters)

    unformatted_rand_out = '%s/%s-%s_unformatted.bed'% (shuffle_path, shuffle_id, iters) # make shuffle file

    BEDshuf = "bedtools shuffle -i %s -g %s -excl %s -chrom -noOverlapping -maxTries 5000 > %s" % (test_enh, CHROM_SZ, BLACKLIST, unformatted_rand_out)

    os.system(BEDshuf)

    #sid = '%s-%s'% (shuffle_id, iters) # make shuffle file

    rand_out = preformatBedfile(unformatted_rand_out, shuffle_id, shuffle_path) #format the shuffle file again

    rm = "rm %s" % unformatted_rand_out
    os.system(rm)

    return rand_out



def preformatBedfile(test_enh, sample_id, test_path):
# before analysis, format bed file. Allow only 4 cols (chr, start, end, sample_id)

    if "shuf" in sample_id:
        test_enh_cut = "%s/%s.bed" % (test_path, sample_id) # prepare to format test_enh
    else:
        test_enh_cut = "%s/cut-%s.bed" % (test_path, sample_id) # prepare to format test_enh

    cmd = '''awk '{print $1"\t", $2"\t", $3"\t", $4}' %s | tr -d " "| sort -k1,1 -k2,2 -k3,3 > %s''' % (test_enh, test_enh_cut)

    print("standardizing Bed format")
    subprocess.call(cmd, shell=True)

    temp = "%ssid_%s.bed" %(test_path, sample_id)
    awk_cmd = '''awk '{$4="%s"}1' FS="\t" OFS="\t"'' %s > %s && mv %s %s''' %(sample_id, test_enh_cut, temp, temp,  test_enh_cut)
    subprocess.call(awk_cmd, shell=True)


    return test_enh_cut


# concatenate shuffles before aging concatenated shuffles w/ bedintersect_syn.
def concat_shuffles(shuffle_path, shuffle_id):

    print("CONCATENATING SHUFFLES")

    catf = "%s/cat_%s_enh_ages.bed" % (shuffle_path, shuffle_id)
    cattemp = "%s/%s_cat_temp.bed" % (shuffle_path, shuffle_id)
    cmd = "cat %s/%s*.bed > %s" %(shuffle_path, shuffle_id, catf)

    subprocess.call(cmd, shell = True)

    stripTabs(catf, cattemp)

    # add enhancer id column
    add_enh_id = '''awk '{$(NF+1)=$1":"$2"-"$3 ; print}' %s > %s && mv %s %s''' % (catf, cattemp, cattemp, catf)
    print("add enh_id")
    subprocess.call(add_enh_id, shell = True)

    # make sure file is tab separated
    tab_cmd = '''awk '{$1=$1}1' OFS="\t" %s > %s && mv %s %s''' % (catf, cattemp, cattemp, catf)
    print("add tabs")
    subprocess.call(tab_cmd, shell = True)

    rm = 'rm %s/%s*-*.bed' % (shuffle_path, shuffle_id)
    subprocess.call(rm, shell = True)

    return catf


# put the pipeline together
def runscripts(TEST_ENH, SAMPLE_ID, TEST_PATH, SPECIES, ANALYZE_AGE, ANALYZE_BREAKS):

    outpath = "%s/%s" % (TEST_PATH, SAMPLE_ID)

    if os.path.exists(outpath) == False and "shuffle" not in outpath:
        mkdir(outpath)

        TEST_PATH = outpath

    else:

        TEST_PATH = TEST_PATH

    TEST_ENH_CUT = "%s/cut-%s.bed" % (TEST_PATH, SAMPLE_ID) # prepare to format test_enh

    if ANALYZE_AGE ==1:
        print("AGING")

        test_enh_formatted = preformatBedfile(TEST_ENH, SAMPLE_ID, TEST_PATH) # format the enhancer bed file and sort

        age_file = age_enh(test_enh_formatted, SAMPLE_ID, TEST_PATH, SPECIES) # age the enhancer file

        temp = "%s/ages/%s-temp.bed" % (TEST_PATH, SAMPLE_ID)
        stripTabs(age_file, temp)

        if ANALYZE_BREAKS ==1: # assemble age architecture
            print("BREAKS", age_file)
            break_file = break_scripts(age_file, SAMPLE_ID, TEST_PATH)

        os_remove(test_enh_formatted)

    elif ANALYZE_BREAKS ==1: # you've already aged the enhancers, just assemble architecture


        AGE_F = TEST_ENH

        print("NO AGING, JUST", AGE_F)

        break_file = break_scripts(AGE_F, SAMPLE_ID, TEST_PATH)



###
#   main
###


def main(argv):


    print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.now())[:20]))


    if RUN_TEST_ENH ==1:

        runscripts(TEST_ENH, SAMPLE_ID, TEST_PATH,\
        SPECIES, ANALYZE_AGE, ANALYZE_BREAKS)


    if ANALYZE_SHUFFLE != 0:

        # create pool and run simulations in parallel

        shuffle_id =  "shuf-"+(SAMPLE_ID).split("_enhancers")[0]
        shuffle_path = "%s/%s/shuffle" % (TEST_PATH, SAMPLE_ID)

        test_enh_formatted = preformatBedfile(TEST_ENH, SAMPLE_ID, TEST_PATH) # format the enhancer bed file and sort
        """
        if ITERATIONS !=0:

            exp_sum_list= []
            val = 0
            for i in range(ITERATIONS):

                shuffled_f = calculateExpected(test_enh_formatted, SAMPLE_ID,\
                TEST_PATH, SPECIES, ANALYZE_AGE, \
                ANALYZE_BREAKS, ANALYZE_TFBS_DEN, val)
                exp_sum_list.append(shuffled_f)
                                val +=1
         """
        pool = Pool(NUM_THREADS)
        partial_calcExp = partial(calculateExpected,\
                                      test_enh_formatted, SAMPLE_ID,\
                                      shuffle_path, SPECIES, ANALYZE_AGE, \
                                      ANALYZE_BREAKS, ANALYZE_TFBS_DEN)


        exp_sum_list = pool.map(partial_calcExp, [i for i in range(ITERATIONS)])
        pool.close()
        pool.join()

        """
            for i in exp_sum_list:

                shuffle_path = "/".join(i.split("/")[:-1])
                shuffle_iter = (i.split("/")[-1]).split(".")[0]
                shuffle_id = 'shuf-' + SAMPLE_ID +"-" + shuffle_iter
        """

        if os.path.exists(shuffle_path) == False:
            mkdir(shuffle_path)

        if "enh_ages.bed" in TEST_ENH and "cat" not in TEST_ENH:
            print("SHUFFLE_ID", shuffle_id)
            runscripts(TEST_ENH, shuffle_id, TEST_PATH,\
            SPECIES, ANALYZE_AGE, ANALYZE_BREAKS)

        elif "cat" not in TEST_ENH:
            catf = concat_shuffles(shuffle_path, shuffle_id)

            runscripts(catf, shuffle_id, shuffle_path,\
            SPECIES, ANALYZE_AGE, ANALYZE_BREAKS)
        else:
            print("sarah, address these problems with shuffle not running")



        rm_cmd = "rm %s/cut-*%s*.bed" %(TEST_PATH, SAMPLE_ID)
        os.system(rm_cmd)

if __name__ == "__main__":
    main(sys.argv[1:])
