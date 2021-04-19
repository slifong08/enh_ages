import os
import subprocess

#%%
msa_ways = [20, 30, 100]
CHRNUM = "chr21"
outfs =[]

random_seed = 42

PATH = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/species_specific_10k/"
PATH = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/all/"

#%%

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
    cmd = f"{PHAST_PATH}./phyloP --features {ocr} --msa-format MAF --method LRT --mode CONACC -d {random_seed} -g {mod} {maf}> {outf}"
    print(f"starting {msaway}")
    subprocess.call(cmd, shell = True)
    print(f"done with {msaway}")

    return outf

def cut_file(path, chrnum):
    ocr = f"{path}{chrnum}.bed"
    temp = f"{path}temp_{chrnum}.bed"
    cmd = ''' awk '{print $1, $2, $3, $4}' FS="\t" OFS="\t" %s > %s && mv %s %s''' %(ocr, temp, temp, ocr)
    subprocess.call(cmd, shell = True)

#%%

for msa in msa_ways:
    cut_file(PATH, CHRNUM)
    outf = run_phylop(msa, CHRNUM, PATH, random_seed)
    outfs.append(outf)
#%%
# prune tree
#f"{phast_path}./tree_doctor -P ${species_list} ${tree} > pruned_tree.nh"


# parse species of interest and mask this out from MAF. I think we want to separate hu chain_nets and rhe chain_nets?
#${params.phast_path}./maf_parse --features ${chrom_bed_path}/${chrom}.bed --mask-features ${params.species_of_interest} ${species_maf} > ${fname}_${params.species_of_interest}_masked.maf
 /dors/capra_lab/bin/./phyloP --features /dors/capra_lab/users/fongsl/tyler/data/CON_ACC/chr21.bed --msa-format MAF --method LRT --mode CONACC -d 42 -g /dors/capra_lab/data/ucsc/hg38/multiz20way/hg38.phastCons20way.mod /dors/capra_lab/data/ucsc/hg38/multiz20way/maf/chr21.maf > /dors/capra_lab/users/fongsl/tyler/data/CON_ACC/multiz20way/chr21_cons_acc.bed

# run phyloP_results
f"{phast_path}./phyloP --features {phastcons} --msa-format MAF --method LRT --mode ACC --subtree {species_of_interest} -d {random_seed} -g {auto_neutral_model} {species_maf}"

#${params.phast_path}./phyloP --features ${phastcons} --msa-format MAF --method LRT --mode ACC --subtree ${params.species_of_interest} -d ${params.random_seed} -g ${params.nonauto_neutral_model}${chrom}.neutral.mod ${species_maf}

""" PHYLOP --help
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
"""

"""
auto_neutral_model = results of phastcons run
