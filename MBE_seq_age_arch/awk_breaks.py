
import argparse
from datetime import datetime
from functools import partial
import glob
from multiprocessing.pool import Pool
import numpy as np
import os
import pandas as pd
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
ANALYZE_SHUFFLE = args.shuffle

NUM_THREADS = 100


print("\n", ANALYZE_AGE, ANALYZE_BREAKS, ANALYZE_SHUFFLE, RUN_TEST_ENH)

# EXTRACT OTHER CONSTANTS
TEST_PATH = "/".join(TEST_ENH.split("/")[:-1])
SAMPLE_ID = (TEST_ENH.split("/")[-1]).split(".")[0]
print("\nSAMPLE_ID:",SAMPLE_ID)
print("\nFILE TO RUN:",TEST_ENH)


###
#   functions
###


def loadConstants(species):  # note chrom.sizes not used in current implementation | 2018.10.29
    return {'hg19': ("/dors/capra_lab/projects/enhancer_ages/dna/hg19_blacklist_gap_gencode.bed", "/dors/capra_lab/data/dna/human/hg19/hg19_trim.chrom.sizes"),
            'hg38': ("/dors/capra_lab/projects/enhancer_ages/dna/hg38_blacklist_gap_gencode.bed", "/dors/capra_lab/data/dna/human/hg38/hg38_trim.chrom.sizes"),
            'mm10': ("/dors/capra_lab/users/bentonml/data/dna/mm10/mm10_blacklist_gap.bed", "/dors/capra_lab/data/dna/mouse/mm10/mm10_trim.chrom.sizes")
            }[species]



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
    temp = os.path.join(path, f"{sid}_sorted.bed")
    cmd = f"sort -k1,1 -k2,2 -k3,3 {infile} > {temp} && mv {temp} {infile}"
    print("SORTING FILE")
    os.system(cmd)


# create a list of autosome chromosomes
def make_chr_list():
    n = list(np.arange(1, 23))
    #n.append("X")

    chr_list = []
    for num in n:
        chrn = "chr" + str(num)
        chr_list.append(chrn)

    return chr_list


# bedtools intersection w/ chr-specific syntenic block.
def bedintersect_syn(f, species, chr_num, sample_id, outpath):

    print("SYNTENY INTERSECTION", chr_num, f, species, chr_num, sample_id, outpath)
    syn_path = f"/dors/capra_lab/data/ucsc/{species}/synteny_age_{species}/"
    syn_file= (f"{syn_path}{chr_num}_syn_age.bed.gz")

    outfile = "%s/%s_%s_ages.bed" %(outpath, chr_num, sample_id)

    BEDcmd = "bedtools intersect -a %s -b %s -wao > %s"  % (f, syn_file, outfile)
    os.system(BEDcmd) # do the intersection


    return outfile

# get the enhancer ages
def age_enh(test_enh, sample_id, test_path, species):

    print("AGING", sample_id, test_enh, test_path, species)
    outpath = os.path.join(test_path, "ages")

    if os.path.exists(outpath) == False:
        os.mkdir(outpath)

    os.chdir(outpath)

    autosomes = make_chr_list()


    # for files w/ regions from multiple chromosomes.
    if "chr" not in sample_id:

        chr_cmd = '''awk '{print >$1"_%s_temp.bed"}' %s''' % (sample_id, test_enh) # split test into chrN.bed files
        print(chr_cmd)
        os.system(chr_cmd)


        enh_chr_list = glob.glob(f"{outpath}/chr*_{sample_id}_temp.bed") # glob chromosomes

        for f in enh_chr_list:

            chr_num = (f.split("/")[-1]).split("_")[0]

            if chr_num in autosomes: # only include autosomes
                sort_bed(f) # sort the file
                bedintersect_syn(f, species, chr_num, sample_id, outpath) # syntenic blcok intersection
                os.remove(f) # remove the chromosome file temp

            else:
                os.remove(f) # delete the sex chromosomes

        # concatenate all the chromosomes together again
        concat_file = os.path.join(outpath, f"{sample_id}_ages.bed")
        concat_cmd = f"cat {outpath}/chr*_{sample_id}_ages.bed > {concat_file}"

        print("concatenating aged enhancer chromosome files")
        subprocess.call(concat_cmd, shell = True)

        rm_cmd = f"rm {outpath}/chr*_{sample_id}_ages.bed"
        subprocess.call(rm_cmd, shell = True)

        return concat_file

    else: # for files w/ regions from one chromosome
        chr_num = ((test_enh.split("/")[-1]).split("_")[0]).split("-")[1]
        outfile = bedintersect_syn(test_enh, species, chr_num, sample_id, outpath) # syntenic blcok intersection

        print("nothing to delete - single chromosome analysis")

        return outfile

def reEval_PrimComplex(df):

    # get all the complex enhancers w/ primate core ages
    prComEnhID = df.loc[(df.core ==1) &
    (df.core_remodeling ==1) &
    (df.taxon2.str.contains("Primate"))]["enh_id"].unique()

    # get all the complex enhancer ids where there is a real human derived sequence
    pr_complex = df.loc[(df.enh_id.isin(prComEnhID)) &
    (df.core_remodeling == 1) &
    (df.core ==0) &
    (df.mrca ==0),
    ]["enh_id"]


    # i'm going to reassign any primate complex enhancer
    # where derived regions are from other primates
    # get the set of primate complex enhancers w/ primate derived sequences
    # and rename them as simple enhancers
    pr_simple = set(prComEnhID) - set(pr_complex)

    # reassign core and core remodeling columns
    df.loc[df.enh_id.isin(pr_simple), "core"] = 1
    df.loc[df.enh_id.isin(pr_simple), "core_remodeling"] = 0
    return df

def assemble_break_data(mrca_file, seg_count_file, syn_gen_bkgd, outpath, sample_id):

    mrca_df = pd.read_csv(mrca_file, sep='\t')
    mrca_df.sort_values(by = "enh_id")

    seg_df = pd.read_csv(seg_count_file, sep='\t', header = None)
    seg_df.columns = ["seg_index", "enh_id"]
    seg_df.sort_values(by = "enh_id")

    breaksdf = pd.merge(mrca_df, seg_df, how = "left", on = "enh_id")

    # round the mrca value
    breaksdf.mrca = breaksdf.mrca.round(3)
    breaksdf["len"] = breaksdf.end_enh - breaksdf.start_enh
    breaksdf = breaksdf.loc[breaksdf.len >5] # exclude any enhancers less than 5bp

    # add binary for simple/complex architecture based on median break value
    breaksdf["core_remodeling"] = 0

    relative_simple =  breaksdf["seg_index"].median() # get median number of age segments.

    if relative_simple == 1: # if the median is one, then complex enhancers >= 2 age segments
        relative_simple+=1

    breaksdf.loc[breaksdf[ "seg_index"] >= relative_simple, "core_remodeling"] = 1

    # add annotation based on simple/complex architecture
    breaksdf["arch"] = "simple"
    breaksdf.loc[breaksdf["core_remodeling"] ==1, "arch"] = "complexenh"

    # reorder columns
    breaksdf = breaksdf[["#chr_enh", "start_enh", "end_enh", "enh_id", "sample_id",
                "seg_index", "core_remodeling", "arch", "mrca"]]

    # merge with other age, taxon annotations
    breaksdf = pd.merge(breaksdf,
    syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]],
    how = "left", on = "mrca").drop_duplicates()


    # save all this information.

    out_summarized_bed = os.path.join(outpath, f"{sample_id}_enh_age_arch_summary_matrix.bed")

    breaksdf.to_csv(out_summarized_bed, sep = '\t', index = False)


    return out_summarized_bed


def break_scripts(age_file, sample_id, test_path, species):

    outpath = os.path.join(test_path, "breaks")

    if os.path.exists(outpath) == False:
        os.mkdir(outpath)


    # sort bed file
    print("sorting age file")
    sort_bed(age_file)

    # remove lines that do not overlap syntenic blocks, overlap X chromosome
    cleanup_file = os.path.join(outpath, f"clean_ages_{sample_id}.bed")

    # remove any lines where
    # !($5 ~ /[.]/ )  syntenic block doesn't overlap the coordinates
    # !($1 ~ /chrX/ ) no sex chromosomes
    # ($14 > 5 )  length of syntenic overlap is greater than 5bp long.

    # print only enhancer coordinates ($1, $2, $3), sample_id ($4), mrca ($12), syntenic lengths ($14)
    cleanup = '''awk '!($5 ~ /[.]/ ) && !($1 ~ /chrX/ ) && ($13 > 5 ) { print $1, "\t", $2,  "\t",$3, "\t", $4, "\t", $11, "\t", $13 }' %s > %s''' % (age_file, cleanup_file)

    print("cleanup\n")
    subprocess.call(cleanup, shell = True)

    # write all the lines that do not overlap syntentic blocks.
    no_overlap_file = os.path.join(outpath, "no_syn_alignment_%s.txt")
    no_overlap_cmd = '''awk '($5 ~ /[.]/ )' %s > %s''' % (age_file, no_overlap_file)
    print("no_overlap")
    subprocess.call(no_overlap_cmd, shell = True)


    # add enhancer id column
    temp = os.path.join(outpath, f"temp_{sample_id}.bed")
    add_enh_id = '''awk '{$(NF+1)=$1":"$2"-"$3 ; print}' OFS="\t" %s > %s && mv %s %s''' % (cleanup_file, temp, temp, cleanup_file)
    print("add enh_id")
    subprocess.call(add_enh_id, shell = True)


    # get max mrca per enhancer from cleanup file
    mrca_file=  os.path.join(outpath, f"max_mrca_{sample_id}.bed")
    #mrca_cmd = '''awk '$5>max[$7]{max[$7]=$5; row[$7]=$0} END{for (i in row) print row[i]}' %s > %s '''% (cleanup_file, mrca_file)
    # first format the file
    mrca_cmd = f"cut -f 1,2,3,4,5,7 {cleanup_file} | sort | uniq > {mrca_file}"
    subprocess.call(mrca_cmd, shell = True)

    # then use python script to get max val.
    mrca_cmd = f"python /dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/get_max_mrca.py {mrca_file}"
    print("get max mrca")
    subprocess.call(mrca_cmd, shell = True)


    # get number of segments per enhancer from outpath file
    # produces # of age segments per enhancer id.
    seg_count_file =  os.path.join(outpath,f"age_seg_count_{sample_id}.bed")
    seg_count_cmd = '''cut -f 5,7 %s | uniq | cut -f2 | uniq -c > %s''' % (cleanup_file, seg_count_file)
    print("count age segments")
    subprocess.call(seg_count_cmd, shell = True)


    # add tabs to seg_index_count
    tab_cmd = '''awk '{$1=$1}1' OFS="\t" %s > %s && mv %s %s''' % (seg_count_file, temp, temp, seg_count_file)
    print("more tabs")
    subprocess.call(tab_cmd, shell = True)


    ### reference MRCA file ###

    syn_gen_bkgd_file = f"/dors/capra_lab/projects/enhancer_ages/{species}_syn_taxon.bed"
    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t',
    usecols = ["mrca", "taxon", "mrca_2", "taxon2"]) # read the file

    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages


    # open up and merge a bunch of dataframes

    print("opening files")
    breakdf = assemble_break_data(mrca_file, seg_count_file, syn_gen_bkgd, outpath, sample_id)
    count = len(open(breakdf).readlines(  )) # makesure the file wrote before deleting ages.

    if count >0:
        os.remove(cleanup_file)
        os.remove(seg_count_file)
        os.remove(mrca_file)

    # index syntenic blocks last.
    syn_index_cmd = f"python /dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/syntenic_assembly.py {age_file}"
    print("syn arch")
    subprocess.call(syn_index_cmd, shell = True)

# generate shuffles
def calculateExpected(test_enh, sample_id, shuffle_path, species, iters):


    print("shuffling", shuffle_path, iters)

    BLACKLIST, CHROM_SZ = loadConstants(species)  # note CHROM_SZ not used

    exp_sum = 0

    shuffle_id = f"shuf-{sample_id}-{iters}"
    unformatted_rand_out = os.path.join(shuffle_path,f'{shuffle_id}-{iters}_unformatted.bed')

    # make shuffle file
    BEDshuf = f"bedtools shuffle -i {test_enh} -g {CHROM_SZ} -excl {BLACKLIST} -chrom -noOverlapping -maxTries 5000 > {unformatted_rand_out}"

    subprocess.call(BEDshuf, shell = True)

    #format the shuffle file again
    rand_out = preformatBedfile(unformatted_rand_out, shuffle_id, shuffle_path)

    #rm
    os.remove(unformatted_rand_out)

    return rand_out


# keep only the first 4 columns of the bed file.
def preformatBedfile(test_enh, sample_id, test_path):

# before analysis, format bed file. Allow only 4 cols (chr, start, end, sample_id)

    if "shuf" in sample_id:
        test_enh_cut = os.path.join(test_path, f"{sample_id}.bed") # prepare to format test_enh
    else:
        test_enh_cut = os.path.join(test_path, f"cut-{sample_id}.bed") # prepare to format test_enh

    # cut and sort the first 4 columns.
    cmd = '''awk '{print $1"\t", $2"\t", $3"\t", $4}' %s | tr -d " "| sort -k1,1 -k2,2 -k3,3 > %s''' % (test_enh, test_enh_cut)
    print("standardizing Bed format")
    subprocess.call(cmd, shell=True)

    # add sample_id as 4th column
    temp = os.path.join(test_path, f"sid_{sample_id}.bed")
    awk_cmd = '''awk '{$4="%s"}1' FS="\t" OFS="\t"'' %s > %s && mv %s %s''' %(sample_id, test_enh_cut, temp, temp,  test_enh_cut)
    #subprocess.call(awk_cmd, shell=True)

    return test_enh_cut


# put the pipeline together
def runscripts(TEST_ENH, SAMPLE_ID, TEST_PATH, SPECIES, ANALYZE_AGE, ANALYZE_BREAKS):

    outpath = os.path.join(TEST_PATH, SAMPLE_ID)

    if os.path.exists(outpath) == False and "shuffle" not in outpath:
        os.mkdir(outpath)

        TEST_PATH = outpath

    else:

        TEST_PATH = TEST_PATH

    TEST_ENH_CUT = os.path.join(TEST_PATH, f"/cut-{SAMPLE_ID}.bed") # prepare to format test_enh

    if ANALYZE_AGE ==1:
        print("AGING")

        if "shuf" in SAMPLE_ID:
            test_enh_formatted = TEST_ENH
        else:
            test_enh_formatted = preformatBedfile(TEST_ENH, SAMPLE_ID, TEST_PATH) # format the enhancer bed file and sort

        age_file = age_enh(test_enh_formatted, SAMPLE_ID, TEST_PATH, SPECIES) # age the enhancer file

        temp = os.path.join(TEST_PATH, f"/ages/{SAMPLE_ID}-temp.bed")
        stripTabs(age_file, temp)

        if ANALYZE_BREAKS ==1: # assemble age architecture
            print("BREAKS", age_file)

            break_file = break_scripts(age_file, SAMPLE_ID, TEST_PATH, SPECIES)

        os.remove(test_enh_formatted)

        print("syn arch")
    #subprocess.call(syn_index_cmd, shell = True)

    elif ANALYZE_BREAKS ==1: # you've already aged the enhancers, just assemble architecture


        AGE_F = TEST_ENH

        print("NO AGING, JUST", AGE_F)
        test_path = "/".join(TEST_PATH.split("/")[:-1]) +"/"


        break_file = break_scripts(AGE_F, SAMPLE_ID, test_path, SPECIES)



###
#   main
###


def main(argv):


    print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.now())[:20]))


    if RUN_TEST_ENH ==1:

        runscripts(TEST_ENH, SAMPLE_ID, TEST_PATH,\
        SPECIES, ANALYZE_AGE, ANALYZE_BREAKS)


    if ANALYZE_SHUFFLE == 1:

        # create pool and run simulations in parallel

        shuffle_id =  "shuf-"+(SAMPLE_ID).split("_enhancers")[0]
        shuffle_path = os.path.join(TEST_PATH, SAMPLE_ID, "shuffle")

        if os.path.exists(shuffle_path) == False:
            os.mkdir(shuffle_path)

        test_enh_formatted = preformatBedfile(TEST_ENH, SAMPLE_ID, TEST_PATH) # format the enhancer bed file and sort

        pool = Pool(NUM_THREADS)
        partial_calcExp = partial(calculateExpected,\
                                      test_enh_formatted, SAMPLE_ID,\
                                      shuffle_path, SPECIES)


        exp_sum_list = pool.map(partial_calcExp, [i for i in range(ITERATIONS)])
        pool.close()
        pool.join()

        if os.path.exists(shuffle_path) == False:
            os.mkdir(shuffle_path)

        if "enh_ages.bed" in TEST_ENH:

            print("SHUFFLE_ID", shuffle_id)
            runscripts(TEST_ENH, shuffle_id, TEST_PATH,\
            SPECIES, ANALYZE_AGE, ANALYZE_BREAKS)

        elif "enh_ages.bed" not in TEST_ENH:

            shuf_fs = glob.glob(f"{shuffle_path}/{shuffle_id}*.bed") # get all the shuffle files

            val = 0

            for shuf_f in shuf_fs: # age each shuffle file.

                iter_id = shuffle_id + "-" + str(val)
                print("iter_id", iter_id)
                runscripts(shuf_f, iter_id, shuffle_path,\
                SPECIES, ANALYZE_AGE, ANALYZE_BREAKS)

                val +=1
        else:
            print("sarah, address these problems with shuffle not running")



        rm_cmd = f"rm {TEST_PATH}/cut-*{SAMPLE_ID}*.bed"
        os.system(rm_cmd)

if __name__ == "__main__":
    main(sys.argv[1:])
