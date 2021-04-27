import os
import subprocess

#%%
msa_ways = [20, 30, 100]
CHRNUM = "chr21"
outfs =[]

random_seed = 42

PATHS = [
"/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/species_specific_10k/",
"/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/all/",
"/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/hars/",
"/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/hu_specific/",
"/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/rhe_specific/",
]

PATHS[3:]

#%% FUNCTIONS

# in case you need to split file on chromosome number before getting started

def split_by_chr(f):
    cmd = '''awk '{print>$1".bed}' %s ''' %f
    subprocess.call(cmd, shell = True)


# run phylop
def run_phylop(msa, chrnum, path, random_seed):

    PHAST_PATH = "/dors/capra_lab/bin/"

    msaway = str(msa) + "way"

    # get the files
    ocr = f"{path}{chrnum}.bed"

    mod = f"/dors/capra_lab/data/ucsc/hg38/multiz{msaway}/hg38.phastCons{msaway}.mod"

    maf = f"/dors/capra_lab/data/ucsc/hg38/multiz{msaway}/maf/{chrnum}.maf"

    # make the outpath, outfile
    outpath = f"{path}multiz{msaway}/"

    if os.path.exists(outpath) == False:
        os.mkdir(outpath)

    outf = f"{outpath}{chrnum}_con_acc.bed"

    # run phyloP
    cmd = f"{PHAST_PATH}./phyloP --features {ocr} --msa-format MAF --method LRT --branch hg38 --mode CONACC -d {random_seed} -g {mod} {maf}> {outf}"
    print(f"starting {msaway}")
    subprocess.call(cmd, shell = True)
    print(f"done with {msaway}")

    return outf


# reduce the amount of information in file for PhyloP run.
def cut_file(path, chrnum):
    ocr = f"{path}{chrnum}.bed"
    temp = f"{path}temp_{chrnum}.bed"
    cmd = ''' awk '{print $1, $2, $3, $4}' FS="\t" OFS="\t" %s > %s && mv %s %s''' %(ocr, temp, temp, ocr)
    subprocess.call(cmd, shell = True)

#%%

for msa in msa_ways:
    for PATH in PATHS[3:]:

        cut_file(PATH, CHRNUM)
        outf = run_phylop(msa, CHRNUM, PATH, random_seed)
        outfs.append(outf)
#%%
"""
example command:

 /dors/capra_lab/bin/./phyloP --features /dors/capra_lab/users/fongsl/tyler/data/CON_ACC/chr21.bed --msa-format MAF --method LRT --mode CONACC -d 42 -g /dors/capra_lab/data/ucsc/hg38/multiz20way/hg38.phastCons20way.mod /dors/capra_lab/data/ucsc/hg38/multiz20way/maf/chr21.maf > /dors/capra_lab/users/fongsl/tyler/data/CON_ACC/multiz20way/chr21_cons_acc.bed

# run phyloP_results
f"{phast_path}./phyloP --features {phastcons} --msa-format MAF --method LRT --mode ACC --subtree {species_of_interest} -d {random_seed} -g {auto_neutral_model} {species_maf}"

From Kathleen's pipeline:
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
