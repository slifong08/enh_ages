#20180809
#sarahfong

# the purpose of this graph is to take bed files of synteny blocks from 46-way multiZ HG19 multiple sequence alignments
# and assign ages to each of those synteny blocks based on the maximum patristic phylogenetic distance between nodes,
# the edge length (not weighted) to the Most Recent Common Ancestor (MRCA)
# I want to use this analysis as a control to compare the whole genome species overlap versus

# I'm thinking about the informaction calculating the average patristic and MRCA would give us. This would describe how strong the max signal is.
# If there is a large difference between the average age and the max age, this means the max distance may be undersampled in the average.
# I guess the average would give me an idea of how much the max age is an outlier, but it might also be a reflection of how undersampled the genomes are on the branch containing the max age alignment

import os, sys
import pandas
import timeit
import glob
import numpy as np

############################################################
# select the synteny block samples to age.
############################################################
dist_path = "/dors/capra_lab/users/fongsl/synteny_blocks/data/synteny_species_hg19/"

file_list = glob.glob("%shg19_species_count_chr*.bed" % dist_path)

############################################################
# get the phylogenetic distances from hg38 to each 100-way multiz species using the neutral tree.
############################################################

phylo_file = "/dors/capra_lab/users/fongsl/enh_age/data/hg19_46way_dist.csv"

phdf = pandas.read_csv(phylo_file, sep = ',')
phdf.columns = ["taxon1", "taxon2", "patr_dist", "edge_dist", "taxon1_mrca_dist", "taxon1_mrca_edge_dist", "taxon2_mrca_dist", "taxon2_mrca_edge_dist"]

phdf=phdf.sort_values(by = 'patr_dist')

species_list = list(phdf["taxon2"].unique())

###### Function to manipulate age####
def enh_age(file):

    sample =((file.split("/")[-1]).split("_")[-1]).split(".")[0]

    df = pandas.read_csv(file, sep='\t', header = -1)

    sanity_check = df[4].unique()

    df.columns = ["chr", "start", "end", "strand", "taxon1", "species_count", "taxon2"]
    df["len"] =df["end"]- df["start"]

    # parse through the species-overlap list, created a matrix of species overlap and phylogenetic values in order to find the oldest enhancer age.

    for i in species_list:

        p = "%spatr"%i

        m= "%smrca"%i

        patr= float(phdf["patr_dist"].loc[phdf["taxon2"]==i])
        mrca= float(phdf["taxon1_mrca_dist"].loc[phdf["taxon2"]==i]) # this is the distance from the hg38 node to the mrca for taxon2


        df[p] = np.where(df["taxon2"].str.contains(i), patr, 0)
        df[m] = np.where(df["taxon2"].str.contains(i), mrca, 0)


    # find the Max phylogenetic distance for each enhancer fragement
    df["max_patr"] =df.filter(like='patr', axis=1).max(axis=1)
    df["max_mrca"] =df.filter(like='mrca', axis=1).max(axis=1)


    # find the Max phylogenetic distance for each enhancer fragement
#    df["mean_patr"] =df.filter(like='patr', axis=1).mean(axis=1)
#    df["mean_mrca"] =df.filter(like='mrca', axis=1).mean(axis=1)


    df= df[["chr", "start", "end", "strand", "taxon1", "species_count", "len", "max_mrca","max_patr"]] # chr, start, end, overlap number, sample_id, multiz start, multiz end, number of species, length, furtherst phylogenetic distance.

    # save the dataframe
    # print("/dors/capra_lab/users/fongsl/enh_age/data/synteny_age_bkgd/%s_syn_age.bed" % sample)
    df.to_csv("/dors/capra_lab/users/fongsl/synteny_blocks/data/synteny_age_bkgd_hg19/%s_syn_age.bed" % sample, index = False, sep = '\t', header = False)

for file in file_list:
    print(file)
    start = timeit.default_timer()
    print("working on %s" % file)
    enh_age(file)
    stop = timeit.default_timer()
    print(stop-start)
