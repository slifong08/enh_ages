species_of_interest = "Homo_sapiens"
phastcons = "GG_LL_OCRs.bed" # I think these need to be the OCRs.
species_maf = "hg38_RheMac.maf" # just a maf w/ hu and rhesus net chains.
auto_neutral_model = "hg38_tree.mod" # mod file

random_seed = 42
#%%


# prune tree
f"{phast_path}./tree_doctor -P ${species_list} ${tree} > pruned_tree.nh"


# parse species of interest and mask this out from MAF. I think we want to separate hu chain_nets and rhe chain_nets?
${params.phast_path}./maf_parse --features ${chrom_bed_path}/${chrom}.bed --mask-features ${params.species_of_interest} ${species_maf} > ${fname}_${params.species_of_interest}_masked.maf


# run phyloP_results

#${params.phast_path}./phyloP --features ${phastcons} --msa-format MAF --method LRT --mode ACC --subtree ${params.species_of_interest} -d ${params.random_seed} -g ${params.nonauto_neutral_model}${chrom}.neutral.mod ${species_maf}

f"{phast_path}./phyloP --features {phastcons} --msa-format MAF --method LRT --mode CONACC # later --subtree {species_of_interest} -d {random_seed} -g ${auto_neutral_model} ${species_maf}"



""" PHYLOP --help
    --features, -f <file>
        Read features from <file> (GFF or BED format) and output a
        table of p-values and related statistics with one row per
        feature.  The features are assumed to use the coordinate frame
        of the first sequence in the alignment.  Not for use with
        --null or --posterior.  See also --gff-scores.

    --msa-format, -i FASTA|PHYLIP|MPM|MAF|SS
        Alignment format (default is to guess format from file contents).

    --method, -m SPH|LRT|SCORE|GERP
        Method used to compute p-values or conservation/acceleration scores
        (Default SPH).  The likelihood ratio test (LRT) and score test
        (SCORE) compare an alternative model having a free scale parameter
        with the given neutral model, or, if --subtree is used, an
        alternative model having free scale parameters for the supertree
        and subtree with a null model having a single free scale parameter.
        P-values are computed by comparing test statistics with asymptotic
        chi-square null distributions.  The GERP-like method (GERP)
        estimates the number of "rejected substitutions" per base by
        comparing the (per-site) maximum likelihood expected number of
        substitutions with the expected number under the neutral model.
        Currently LRT, SCORE, and GERP can be used only with
        --base-by-base, --wig-scores, or --features.

   --mode, -o CON|ACC|NNEUT|CONACC
        (For use with --wig-scores, --base-by-base, or --features) Whether
        to compute one-sided p-values so that small p (large -log p)
        indicates unexpected conservation (CON; the default) or
        acceleration (ACC); or two-sided p-values such that small p
        indicates an unexpected departure from neutrality (NNEUT).  The
        fourth option (CONACC) uses positive values (p-values or scores) to
        indicate conservation and negative values to indicate acceleration.
        In GERP mode, CON and CONACC both report the number of rejected
        substitutions R (which may be negative), while ACC reports -R, and
        NNEUT reports abs(R).

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
