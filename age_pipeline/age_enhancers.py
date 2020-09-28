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
#from pybedtools import BedTool
import subprocess

# TO RUN

# python script input_file sample_id

# python /dors/capra_lab/fongsl/enh_age/bin/age_enhancers.py UBERON_0002372_tonsil_expressed_enhancers.bed UBERON0002372


###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate enhancer age.")

arg_parser.add_argument("region_file_1", help='bed file 1 (enhancers to age) w/ full path')

arg_parser.add_argument("-i", "--iters", type=int, default=100,
                        help='number of simulation iterations; default=100')

arg_parser.add_argument("-s", "--species", type=str, default='hg19', choices=['hg19', 'hg38', 'mm10'],
                        help='species and assembly; default=hg19')

arg_parser.add_argument("-n", "--num_threads", type=int, default=100,
                        help='number of threads; default=100')

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
TEST_ENH = args.region_file_1
ITERATIONS = args.iters
SPECIES = args.species
RUN_TEST_ENH = args.run_testbed
ANALYZE_AGE = args.age
ANALYZE_BREAKS = args.breaks
ANALYZE_TFBS_DEN = args.tfbs_den
ANALYZE_SHUFFLE = args.shuffle
if args.num_threads:
    NUM_THREADS = args.num_threads
else:
    NUM_THREADS = 100
print(ANALYZE_AGE, ANALYZE_BREAKS, ANALYZE_TFBS_DEN, ANALYZE_SHUFFLE, RUN_TEST_ENH, "\nnum_threads", NUM_THREADS)

# EXTRACT OTHER CONSTANTS
TEST_PATH = "/".join(TEST_ENH.split("/")[:-1])
SAMPLE_ID = (TEST_ENH.split("/")[-1]).split(".")[0]

AGE_OUTFILE ="%s/%s_enh_ages.bed" %(TEST_PATH, SAMPLE_ID)
SHUFFLE_PATH = "%s/shuffle" % TEST_PATH # make a shuffle path file

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
        os.system(cmd)

def os_remove(files):
    cmd = "rm %s" % files
    os.system(cmd)

def stripTabs(infile, tempfile):
    cmd = '''tr -s " \t" < %s > %s && mv %s %s''' % (infile, tempfile, tempfile, infile)
    os.system(cmd)

    rm_cmd = "rm %s" % tempfile
    os.system(cmd)

def enh_age(test_enh, sample_id, test_path, species):

    print("aging", sample_id)
    outpath = "%s/ages/" % test_path # mkdir ./ages/
    mkdir(outpath)

    os.chdir(outpath)

    chr_cmd = '''awk '{print >$1"_%s_temp.bed"}' %s''' % (sample_id, test_enh) # split test into chrN.bed files
    os.system(chr_cmd)

    enh_chr_list = glob.glob("%schr*_%s_temp.bed" % (outpath, sample_id)) # glob chromosomes
    sex_chr = ["chrX", "chrY", "chrM", "chr"]

    for f in enh_chr_list:
        chr_num = (f.split("/")[-1]).split("_")[0]
        if chr_num not in sex_chr: # filter out sex chromosomes.
            print(chr_num)

            syn_path = "/dors/capra_lab/data/ucsc/%s/synteny_age_bkgd_%s/" %(species, species)
            syn_file= ("%s%s_syn_age.bed.gz" % (syn_path, chr_num))

            outfile = "%s%s_%s_ages.bed" %(outpath, chr_num, sample_id)

            BEDcmd = "bedtools intersect -a %s -b %s -wo > %s"  % (f, syn_file, outfile)
            os.system(BEDcmd)

            # write any lines where enhancer does not overlap syntenic break, does not have mrca.
            no_syn = "%s%s-%s_no_syn_alignment_tiles.txt" % (outpath, chr_num, sample_id) # write non-overlapping regions
            NOSYNcmd = "bedtools intersect -a %s -b %s -v -wo > %s"  % (f, syn_file, no_syn)
            os.system(NOSYNcmd)

            os_remove(f)
        else:
            os_remove(f) # delete the sex chromosomes

    concat_file = "%s%s_enh_ages.bed" % (outpath, sample_id)
    concat_cmd = "cat %schr*_%s_ages.bed > %s" % (outpath, sample_id, concat_file)

    print("concatenating aged enhancer chromosome files")
    os.system(concat_cmd)

    syn_concat_file = "%s%s_no_syn_alignment_tiles.txt" % (outpath, sample_id)
    syn_concat_cmd = "cat %schr*-%s_no_syn_alignment_tiles.txt > %s" % (outpath, sample_id, syn_concat_file)#concatenate chromosomes

    os.system(syn_concat_cmd)
    os_remove("%schr*-%s_no_syn_alignment_tiles.txt"%(outpath, sample_id))
    os_remove("%schr*_%s_ages.bed" %(outpath, sample_id))
    return concat_file

def assemble_break(df, id_dict, iteration):

    i = id_dict[iteration] # get the enh_id from the id_dict

    test = df[["chr_enh",  # filter df to only enh_id
                   "start_enh", "end_enh","id",
               "chr_syn","start_syn","end_syn",
                   "mrca","len_syn_overlap","test_id"]]\
        .loc[df["test_id"]==i]\
        .drop_duplicates()\
        .sort_values(["chr_syn","start_syn","end_syn"]).reset_index()

    test = test.loc[test.len_syn_overlap>5] # remove segments less than 5bp in length

    new_data = test[["chr_enh","start_enh", "end_enh","id","test_id", "chr_syn"]]\
    .drop_duplicates()# new dataframe for reassembled breaks

    collect_df = pd.DataFrame()

    age_seg = [(round(age,3), sum(1 for i in rows)) for age,rows in groupby(test["mrca"])] # round the age and sum the number of rows per age.

    start, end = test["start_enh"].astype(int), test["end_enh"].astype(int) # change start and end

    df_index = 0 # placeholder - for syntenic segment index
    seg_count = 0 # placeholder - for syntenic segment count
    core_age = max(i for i, k in age_seg) # GET OLDEST/MAX AGE

    for tup in age_seg: # unpack the list of tuples
        age, idx = tup

        if len(age_seg)== 1: # ONE AGE SEGMENT - SIMPLE ARCHITECTURE

            new_data["start_syn"] = start
            new_data["end_syn"] = end
            new_data["mrca_seg"] = age
            new_data["seg_index"] = 0 # index syntenic segments
            new_data["core"] = 1 # binary core measure
            new_data["core_remodeling"] = 0 # binary - core remodeling measure
            collect_df= collect_df.append(new_data)

        else: # COMPLEX

            new_data["seg_index"]= seg_count # index syntenic segments,
            new_data["core_remodeling"] = 1 # binary - core remodeling measure

            if age == core_age: # OLDEST AGE
                new_data["core"] = 1

            else:
                new_data["core"] = 0

            if seg_count == 0: # trim first syntenic block to start
                new_data["start_syn"] = start
                new_data["end_syn"] = test.loc[idx-1, "end_syn"]

                new_data["mrca_seg"] = round(age, 3)
                collect_df= collect_df.append(new_data)

            elif seg_count == len(age_seg)-1: # trim last syntenic block
                new_data["mrca_seg"] = age
                new_data["start_syn"] = test.loc[df_index, "start_syn"]
                new_data["end_syn"] = end
                collect_df= collect_df.append(new_data)

            else: # deal with all the blocks in between first and last syntenic block
                new_data["mrca_seg"] = age
                new_data["start_syn"]= test.loc[df_index, "start_syn"]
                new_data["end_syn"]= test.loc[df_index + idx -1, "end_syn"]
                collect_df= collect_df.append(new_data)

        df_index +=idx # go to next index
        seg_count +=1 # count age segments

    return collect_df

def break_tfbs(concat_file, sample_id, test_path, num_threads):

    inpath = "%s/ages/" % test_path # mkdir ./ages/
    outpath = "%s/breaks/" % test_path # mkdir ./break/
    mkdir(outpath)

    age_breaks_out = "%s%s_age_breaks.bed" % (outpath, sample_id)

    age_temp = "%stemp_%s.bed" % (inpath, sample_id)

    stripTabs(concat_file, age_temp) # strip any double tabs

    indf = pd.read_csv(concat_file, sep = '\t', header = None , low_memory=False) # open up the bed file and assemble breaks.

    indf.columns = ["chr_enh", "start_enh", "end_enh", "id",
                    "chr_syn", "start_syn","end_syn","strand",
                    "ref","num_species", "len_syn","mrca",
                    "patr","len_syn_overlap"]# rename the columns

    indf["test_id"] = indf["chr_enh"] + ":" + indf["start_enh"].map(str) + "-" + indf["end_enh"].map(str) # add a test_id

    df = indf.loc[indf["len_syn_overlap"]>5] # remove all the records that do not overlap syntenic blocks

    ENH_COUNT = len(df["test_id"].unique())
    print("unique enhancers =", ENH_COUNT) # count the number of enhancers

    print("# rows to reduce and assemble breaks = ", len(df.drop_duplicates())) # drop duplicates
    df = df.drop_duplicates()

    #####
    # prepare to join breaks
    #####

    df["mrca"] = df["mrca"].astype(float).round(3) # round MRCA distance

    id_list = df["test_id"].unique() # get enh_ids
    id_dict = dict(enumerate(df["test_id"].unique())) # dictionary of enhancer ids

    #####
    # start the break assembly
    #####

    start = datetime.now()
    print("Start assembling breaks", start)

    val = 0
    end_val = len(id_list) # batches of

    while val <  end_val:

        print(val, end_val)

        if end_val - val > num_threads:
            new_range = np.arange(num_threads) + val
        else:
            num_threads = abs(end_val - val)

            new_range = np.arange(num_threads) + val
            print(end_val,"minus", val, "equals", num_threads)

        pool = Pool(num_threads)
        partial_calcExp = partial(assemble_break, df, id_dict)
        results = pool.map(partial_calcExp, [i for i in new_range])
        pool.close()
        pool.join()

        temp = pd.concat(results, sort = True)
        temp = temp[['chr_syn','start_syn','end_syn',"id",'test_id',
                         'chr_enh','start_enh','end_enh','seg_index',
                         'core_remodeling','core','mrca_seg']]
        with open(age_breaks_out, 'a') as f:
                temp.to_csv(f, sep = '\t', header = False, index = False) # write last enhancers
        val +=num_threads

    clean_up = "rm %s" % concat_file # cleanup the input file
    os.system(clean_up)

    print("Finished assembling breaks", datetime.now())

    return age_breaks_out

def tfbs_density(core_breaks_file, sample_id, test_path):

    inpath = "%s/breaks/" % test_path

    outpath = "%s/tfbs/"% test_path

    #####
    # prepare workspace
    #####

    outfile_wao = "%s%s_x_raw_tfbs_midpeak.bed"% (outpath, sample_id)
    outfile_waoTemp = "%stemp.bed"% (outpath)

    enh_out = "%s%s_enh_tfbs_density.bed" % (outpath, sample_id)
    syn_out = "%s%s_syn_tfbs_density.bed" % (outpath, sample_id)

    tfbs_file = "/dors/capra_lab/data/encode/midpeaks_wgEncodeRegTfbsClusteredV3.bed" # TFBS data, 30bp trimmed ChIP-Peaks

    mkdir(outpath) # make directory

    ########
    # PART 1 - sort bed file
    ########
    #c = BedTool(core_breaks_file).sort()

    ########
    # PART 2 - Bed intersect file with TFBS
    ########
    #t = BedTool(tfbs_file).sort()
    #c.intersect(t, wao =True).saveas(outfile_wao)

    BEDcmd = "bedtools intersect -a %s -b %s -wao > %s" % ( core_breaks_file, tfbs_file, outfile_wao)
    os.system(BEDcmd)
    stripTabs(outfile_wao, outfile_waoTemp)

    ########
    # open the sample_id file
    ########

    enh = pd.read_csv(outfile_wao, header= None, sep= '\t',low_memory = False)
    print(len(enh), "raw df len")

    ########
    # reformat file to get columns/data I want
    ########

    print("the number syntenic blocks that did not overlap TFBS" ,\
          len(enh.iloc[:,-1].loc[enh[20]==0]))

    #select the columns you want for analyzing break density
    enh = enh[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 20]].drop_duplicates()

    # rename the columns
    enh.columns = ["chr_syn", "start_syn", "end_syn", "enh_id",
                   "chr_enh", "start_enh", "end_enh","seg_index",
                   "core_remodeling", "core", "mrca","chr_tf",
                   "start_tf", "end_tf", "tf", "len_tfbs_overlap"]

    # assign unique tf_id

    enh["tf_id"] = enh["chr_tf"] + ":" + enh["start_tf"].map(str) + "-" + enh["end_tf"].map(str) +"_"+ enh["tf"]

    print("the number of unique enhancers", len(enh["enh_id"].unique()))

    enh_tf_list = enh["enh_id"].unique() # make a list of unique enhancer ids

    enh_tf_list = enh_tf_list
    val = 0

    # int transform
    enh[["start_enh", "end_enh", "start_tf", "start_syn","end_syn","end_tf", "seg_index", "core_remodeling",       "core","len_tfbs_overlap"]] = enh[["start_enh",  "end_enh", "start_tf","start_syn", "end_syn",
         "end_tf", "seg_index", "core_remodeling", "core","len_tfbs_overlap"]].astype(int)

    enh["mrca"] = enh["mrca"].astype(float).round(3)

    enh["enh_len"]= enh["end_enh"]- enh["start_enh"] # enhancer lengths
    enh["syn_len"]= enh["end_syn"]- enh["start_syn"] # syntenic lengths
    enh["syn_id"] = enh["chr_syn"]+ ":" + enh["start_syn"].map(str) + "-" + enh["end_syn"].map(str)

    ########
    # calculate TFBS density
    #
    # count tfs associated with enhancers
    # count uniq tfs associated with enhancers
    # count syn tfs associated with syntenic blocks
    # count uniq syn tfs associated with syntenic blocks
    ########

    enh["tf_bin"] = 0
    enh["tf_bin"].loc[enh.len_tfbs_overlap >5] = 1

    # syn_tf_den

    syntf_count = enh.groupby(['enh_id','syn_id', 'mrca', 'seg_index',
                               'core', 'core_remodeling', 'syn_len'])["tf_bin"].sum().reset_index()

    syntf_count["syn_tf_den"] = syntf_count.tf_bin.divide(syntf_count.syn_len) # calculate syn_density

    syntf_count.columns = ['enh_id','syn_id', 'mrca', 'seg_index',
                               'core', 'core_remodeling', 'syn_len',
                       'syn_tf_count', 'syn_tf_den'] # rename columns

    # syn_tf_u_den

    syntfu_count = enh[['enh_id','syn_id', 'seg_index', 'enh_len', 'syn_len', 'tf', 'tf_bin']]\
    .drop_duplicates().reset_index() # get unique tf density

    syntfu_count = syntfu_count.groupby(['enh_id','syn_id', 'seg_index', 'enh_len', 'syn_len'])["tf_bin"]\
    .sum().reset_index()

    syntfu_count["syn_tf_u_den"] = syntfu_count.tf_bin.divide(syntfu_count.syn_len) # calculate syn_density

    syntfu_count.columns = ['enh_id','syn_id', 'seg_index', 'enh_len', 'syn_len', 'syn_tf_u_count', 'syn_tf_u_den'] # rename columns

    # enh_tf_den

    enhtf_count = enh.groupby(['enh_id','core_remodeling', 'enh_len'])["tf_bin"].sum().reset_index()

    enhtf_count["enh_tf_den"] = enhtf_count.tf_bin.divide(enhtf_count.enh_len) # calculate syn_density

    enhtf_count.columns = ['enh_id','core_remodeling', 'enh_len', 'enh_tf_count', 'enh_tf_den'] # rename columns

    enh_mrca = enh.groupby(['enh_id'])['mrca'].max().reset_index()

    enhtf_count = pd.merge(enhtf_count, enh_mrca, how = "left", on = "enh_id") # get max age

    # enh_tf_u_den

    enhtfu_count = enh[['enh_id','enh_len', 'tf', "tf_bin"]].drop_duplicates().reset_index() # get unique tf density

    enhtfu_count = enhtfu_count.groupby(['enh_id','enh_len'])["tf_bin"].sum().reset_index()

    enhtfu_count["enh_tf_u_den"] = enhtfu_count.tf_bin.divide(enhtfu_count.enh_len) # calculate syn_density

    enhtfu_count.columns = ['enh_id', 'enh_len', 'enh_tf_u_count', 'enh_tf_u_den'] # rename columns

    enh = pd.merge(enhtf_count, enhtfu_count, how = "left", on = (['enh_id','enh_len']))

    syn = pd.merge(syntf_count, syntfu_count, how = "left", on =(['enh_id','syn_id', 'seg_index', 'syn_len']))

    with open(syn_out, 'a') as f:
        syn.to_csv(f, sep = '\t', header = True, index = False) # write last enhancers

    with open(enh_out, 'a') as f:
        enh.to_csv(f, sep = '\t', header = True, index = False) # write last enhancers

    print("Finished calculating SYN TFBS density", datetime.now())

def calculateExpected(test_enh, sample_id, test_path, species, analyze_age, analyze_breaks, analyze_tfbs_den, iters):
    print("shuffling", iters)

    BLACKLIST, CHROM_SZ = loadConstants(species)  # note CHROM_SZ not used

    exp_sum = 0

    shuffle_path = "%s/shuffle/" % test_path # make a shuffle path file
    mkdir(shuffle_path)

    shuffle_id = "shuf-%s" % (sample_id)

    unformatted_rand_out = '%s%s-%s_unformatted.bed'% (shuffle_path, shuffle_id, iters) # make shuffle file

    BEDshuf = "bedtools shuffle -i %s -g %s -excl %s -chrom -noOverlapping -maxTries 50000 > %s" % (test_enh, CHROM_SZ, BLACKLIST, unformatted_rand_out)
    os.system(BEDshuf)

    #rand_file = b(test_enh).shuffle(genome='hg19', excl=BLACKLIST, chrom=True, noOverlapping=True, maxTries = 50000) # shuffle bed
    #rand_file.saveas(unformatted_rand_out)# save shuffle

    rand_out = '%s%s-%s.bed'% (shuffle_path, shuffle_id, iters) # make shuffle file

    formatBedfile(unformatted_rand_out, rand_out) #format the shuffle file again

    rm = "rm %s" % unformatted_rand_out
    os.system(rm)

    return rand_out

def formatBedfile(test_enh, out_test_enh):

    cmd = '''awk 'OFS=" " {print $1"\t", $2"\t", $3"\t", $4}' %s | tr -d " "| sort -k1,1 -k2,2 -k3,3 > %s''' % (test_enh, out_test_enh)
    print("standardizing Bed format")
    subprocess.call(cmd, shell=True)

###
#   main
###

def main(argv):
    print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.now())[:20]))

    ### standardize the enhancer file ###
    TEST_ENH_CUT = "%s/cut-%s.bed" % (TEST_PATH, SAMPLE_ID)

    formatBedfile(TEST_ENH, TEST_ENH_CUT) # format the enhancer bed file and sort

    if RUN_TEST_ENH ==1:

        if ANALYZE_AGE != 0:
            age_file = enh_age(TEST_ENH_CUT, SAMPLE_ID, TEST_PATH, SPECIES)

        if ANALYZE_BREAKS != 0:
            if ANALYZE_AGE!=0:
                break_file = break_tfbs(age_file, SAMPLE_ID, TEST_PATH, NUM_THREADS)

            else:
                AGE_F = "%s/ages/%s_enh_ages.bed" % (TEST_PATH, SAMPLE_ID)
                break_file = break_tfbs(AGE_F, SAMPLE_ID, TEST_PATH, NUM_THREADS)

        if ANALYZE_TFBS_DEN != 0:
            if ANALYZE_BREAKS != 0: # if break tfbs was run.
                tfbs_density(break_file, SAMPLE_ID, TEST_PATH)

            else:  # if break tfbs was not run.
                BREAK_F = "%s/breaks/%s_age_breaks.bed" % (TEST_PATH, SAMPLE_ID)
                tfbs_density(TEST_ENH, SAMPLE_ID, TEST_PATH)

    # create pool and run simulations in parallel

    if ANALYZE_SHUFFLE != 0:

        pool = Pool(NUM_THREADS)

        partial_calcExp = partial(calculateExpected,\
                                  #BedTool(TEST_ENH_CUT),\
                                  TEST_ENH_CUT, SAMPLE_ID,\
                                  TEST_PATH, SPECIES, ANALYZE_AGE, \
                                  ANALYZE_BREAKS, ANALYZE_TFBS_DEN)

        exp_sum_list = pool.map(partial_calcExp, [i for i in range(ITERATIONS)])

        for i in exp_sum_list:
            print(i)

            shuffle_path = "/".join(i.split("/")[:-1])
            shuffle_id = (i.split("/")[-1]).split(".")[0]
            print(shuffle_id)

            if ANALYZE_AGE !=0:
                shuf_age_file = enh_age(i, shuffle_id, shuffle_path, SPECIES) # age the shuffles


            if ANALYZE_BREAKS !=0:
                if ANALYZE_AGE!=0: # if continous with ANALYZE_AGE
                    shuf_break_file = break_tfbs(shuf_age_file, shuffle_id, shuffle_path, NUM_THREADS)

                else: # not continous with ANALYZE_AGE
                    AGE_F = "%s/ages/%s_enh_ages.bed" % (shuffle_path, shuffle_id)
                    break_file = break_tfbs(AGE_F, shuffle_id, TEST_PATH, NUM_THREADS)


            if ANALYZE_TFBS_DEN !=0:
                if ANALYZE_BREAKS != 0: # if break tfbs was run.
                    tfbs_density(shuf_break_file, shuffle_id, shuffle_path)

                else:  # if break tfbs was not run.
                    BREAK_F = "%s/breaks/%s_age_breaks.bed" % (TEST_PATH, SAMPLE_ID)
                    tfbs_density(BREAK_F, shuffle_id, TEST_PATH)

    os_remove(TEST_ENH_CUT)

if __name__ == "__main__":
    main(sys.argv[1:])
