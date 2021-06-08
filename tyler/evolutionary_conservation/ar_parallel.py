import argparse
import glob
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing.pool import Pool, TimeoutError
import numpy as np
import os, sys
import subprocess

#%% ARGPARSE arguments.

arg_parser = argparse.ArgumentParser(description=" describe argparse")

arg_parser.add_argument("bedfile", help ='bed file w/ full path')
arg_parser.add_argument("-br","--branches", help ='hg38, rheMac3')
arg_parser.add_argument("-m", "--multiz", help ='20-, 30-, 100-way multiple sequence alignments in hg38')


# PARSE THE ARGUMENTS
args = arg_parser.parse_args()

F = args.bedfile # the bedfile
BRANCH = args.branches # the branches to test.
PATH = "/".join(F.split("/")[:-1]) + "/" # the path
MSA_WAY = args.multiz # multiple sequence alignment.

random_seed = 42


#%% FUNCTIONS

# run all the chromosomes
def make_chr_list():
    n = list(np.arange(1, 23))
    #n.append("X")

    chr_list = []
    for num in n:
        chrn = "chr" + str(num)
        chr_list.append(chrn)

    return chr_list


# in case you need to split file on chromosome number before getting started
def split_by_chr(f):

    if "chr" not in f:
        cmd = '''awk '{print>$1".bed"}' %s ''' %f
        subprocess.call(cmd, shell = True)


# run phylop

def run_phylop(msa, ocr, chrnum, path, random_seed, branch):

    PHAST_PATH = "/dors/capra_lab/bin/"

    msaway = str(msa) + "way"

    # for multiple branches, need to remove whitespace
    # before generating outpath

    if " " in branch:
        branch_path = "_".join(branch.split(" "))
    else:
        branch_path = branch

    # the neutral tree
    mod = f"/dors/capra_lab/data/ucsc/hg38/multiz{msaway}/hg38.phastCons{msaway}.mod"

    # the multiple sequence alignment file
    maf_zipped = f"/dors/capra_lab/data/ucsc/hg38/multiz{msaway}/maf/{chrnum}.maf.gz"
    maf_unzipped = maf_zipped.split(".gz")[0]

    # if maf needs to be unzipped
    if os.path.exists(maf_unzipped) == False:
        cmd = f"gunzip {maf_zipped}"
        subprocess.call(cmd, shell = True)

    # make the outpath, outfile
    outpath = f"{path}multiz{msaway}_{branch_path}/"

    print("\n\nOUTPATH:", branch_path, "\n\n")

    if os.path.exists(outpath) == False:
        os.mkdir(outpath)

    outf = f"{outpath}{chrnum}_con_acc.bed"

    # check to see that you haven't done the phyloP analysis on this file.
    if os.path.exists(outf) == False:

        # run phyloP!
        cmd = f"{PHAST_PATH}./phyloP --features {ocr} --msa-format MAF --method LRT --branch {branch} --mode CONACC -d {random_seed} -g {mod} {maf_unzipped}> {outf}"

        print(f"starting {msaway}")

        subprocess.call(cmd, shell = True)

        print(f"done with {msaway}")

    # if this has already been run -
    else:
        print("already processed", outf)

    return outf


# reduce the amount of information in input file

def cut_file(path, chrnum):
    ocr = f"{path}{chrnum}.bed"
    temp = f"{path}temp_{chrnum}.bed"
    cmd = ''' awk '{print $1, $2, $3, $4}' FS="\t" OFS="\t" %s > %s''' %(ocr, temp)

    subprocess.call(cmd, shell = True)

    return temp


#%% MAIN


def main(argv):
    chr_list = make_chr_list() # generate a list of chr numbers

    split_by_chr(F) # split the file into chromosome files.

    chr_dict = {}
    for CHRNUM in chr_list:
        ocr = cut_file(PATH, CHRNUM) # format the file

        chr_dict[CHRNUM] = ocr

    # prepare to run parallel jobs
    num_cores = multiprocessing.cpu_count()
    print("number of cores", num_cores)

    # run parallel jobs
    print(MSA_WAY, PATH, random_seed, BRANCH)
    results = Parallel(n_jobs=num_cores, verbose=100, prefer="threads")(delayed(run_phylop)(MSA_WAY, ocr, chrnum, PATH, random_seed, BRANCH) for chrnum, ocr in chr_dict.items())

    for r in results:

        # check that there are results first.
        if os.path.getsize(r) > 0:

            chrnum = "chr" + (r.split("chr")[1]).split("_con_acc")[0] # extract the chromosome number

            temp = f"{PATH}temp_{chrnum}.bed"
            os.remove(temp)
            print("removed", temp)


if __name__ == "__main__":
    main(sys.argv[1:])
#%%
"""
example command:

 /dors/capra_lab/bin/./phyloP --features /dors/capra_lab/users/fongsl/tyler/data/CON_ACC/chr21.bed --msa-format MAF --method LRT --mode CONACC -d 42 -g /dors/capra_lab/data/ucsc/hg38/multiz20way/hg38.phastCons20way.mod /dors/capra_lab/data/ucsc/hg38/multiz20way/maf/chr21.maf > /dors/capra_lab/users/fongsl/tyler/data/CON_ACC/multiz20way/chr21_cons_acc.bed

# run phyloP_results
f"{phast_path}./phyloP --features {phastcons} --msa-format MAF --method LRT --mode ACC --subtree {species_of_interest} -d {random_seed} -g {auto_neutral_model} {species_maf}"

From Kathleen's pipeline:s
#${params.phast_path}./phyloP --features ${phastcons} --msa-format MAF --method LRT --mode ACC --subtree ${params.species_of_interest} -d ${params.random_seed} -g ${params.nonauto_neutral_model}${chrom}.neutral.mod ${species_maf}

From PhyloP documentation
 PHYLOP --help
    --features, -f <file>
        Read features from <file> (GFF or BED format) and output a
        table of p-values and related statistics with one row per
        feature.  The features are assumed to use the coordinate frame
        of the first sequence in the alignment.  Not for use with
        --null or --posterior.  See also --gff-scores.

    --msa-format, -i FASTA|PHYLIP|MPM|MAF|SS
        Alignment format (default is to guess format from file contents).

   --subtree, -s <node-name>
        (Not available in GERP mode) Partition the tree into the subtree
        beneath the node whose name is given and the complementary
        supertree, and consider conservation/acceleration in the subtree
        given the supertree.  The branch above the specified node is
        included with the subtree.  Thus, given the tree
        "((human,chimp)primate,(mouse,rat)rodent)", the option "--subtree
        primate" will create one partition consisting of human, chimp, and
        the branch leading to them, and another partition consisting of the
        rest of the tree; "--subtree human" will create one partition
        consisting only of human and the branch leading to it and another
        partition consisting of the rest of the tree.  In 'SPH' mode, a

   --gff-scores, -g
        (For use with features)  Instead of a table, output a GFF and
        assign each feature a score equal to its -log p-value.

   --seed, -d <seed>
        Provide a random number seed, should be an integer >=1.  Random
        numbers are used in some cases to generate starting values for
        optimization.  If not specified will use a seed based on the
	current time.



auto_neutral_model = results of phastcons run
"""
