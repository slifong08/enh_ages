# created

# 2019-07-03
# sarah fong

# inputs - region_file_1, region_file_2, interations , species

# outputs - overlap count, fold change, p-values, bootstrapped CI from 10000 bootstraps of shuffled fold-change overlaps. 

import os
import sys, traceback
import argparse
import csv
import datetime
import numpy as np
from functools import partial
import matplotlib.pyplot as plt
from multiprocessing import Pool
import pandas
from pybedtools import BedTool


###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate enrichment between bed files.")

arg_parser.add_argument("region_file_1", help='bed file 1 (shuffled)')

arg_parser.add_argument("region_file_2", help='bed file 2 (not shuffled)')

arg_parser.add_argument("-i", "--iters", type=int, default=100,
                        help='number of simulation iterations; default=100')

arg_parser.add_argument("-s", "--species", type=str, default='hg19', choices=['hg19', 'hg38', 'mm10'],
                        help='species and assembly; default=hg19')

arg_parser.add_argument("-n", "--num_threads", type=int,
                        help='number of threads; default=SLURM_CPUS_PER_TASK or 1')

arg_parser.add_argument("-o", "--outfile", type=str, default=None,
                        help="print expected counts to file")

args = arg_parser.parse_args()

# save parameters
ANNOTATION_FILENAME = args.region_file_1
TEST_FILENAME = args.region_file_2
COUNT_FILENAME = args.outfile
ITERATIONS = args.iters
SPECIES = args.species

# calculate the number of threads
if args.num_threads:
    num_threads = args.num_threads
else:
    num_threads = int(os.getenv('SLURM_CPUS_PER_TASK', 1))

###
#   functions
###

def loadConstants(species):  # note chrom.sizes not used in current implementation | 2018.10.29
    return {'hg19': ("/dors/capra_lab/users/bentonml/data/dna/hg19/hg19_blacklist_gap.bed", "/dors/capra_lab/data/dna/human/hg19/hg19_trim.chrom.sizes"),
            'hg38': ("/dors/capra_lab/users/bentonml/data/dna/hg38/hg38_blacklist_gap.bed", "/dors/capra_lab/data/dna/human/hg38/hg38_trim.chrom.sizes"),
            'mm10': ("/dors/capra_lab/users/bentonml/data/dna/mm10/mm10_blacklist_gap.bed", "/dors/capra_lab/data/dna/mouse/mm10/mm10_trim.chrom.sizes")
            }[species]

def calculateObserved(annotation, test):
    obs_sum = 0

    obs_intersect = annotation.intersect(test, wo=True) # intersect obs regions w/ test file (e.g. enh x gwas snps)

    for line in obs_intersect:
        obs_sum += int(line[-1]) # sum snps in obs regions

    return obs_sum

def calculateExpected(annotation, test, species, iters):
    
    BLACKLIST, CHROM_SZ = loadConstants(species)  # note CHROM_SZ not used
    exp_sum = 0

    rand_file = annotation.shuffle(genome='hg19', excl=BLACKLIST, chrom=True, noOverlapping=True) # shuffle obs regions

    exp_intersect = rand_file.intersect(test, wo=True) # intersect shuffled obs regions w/ test file

    for line in exp_intersect:
        exp_sum += int(line[-1]) # sum snps in shuffled exp regions

    return exp_sum

def calculateEmpiricalP(obs, exp_sum_list):
    
    mu = np.mean(exp_sum_list) # 20190523 mean is now median
    sigma = np.std(exp_sum_list)
    dist_from_mu = [exp - mu for exp in exp_sum_list]
    p_sum = sum(1 for exp_dist in dist_from_mu if abs(exp_dist) >= abs(obs - mu))

    fold_change = (obs + 1.0) / (mu + 1.0) # fold chage of the median
    fold_changes = list((obs + 1.0) / (m + 1.0) for m in exp_sum_list)
    fold_change_95 = (obs + 1.0) / (np.percentile(exp_sum_list, 95) + 1.0) # 20190527 95th percentile value from exp_sum_list
    fold_change_5 = (obs + 1.0) / (np.percentile(exp_sum_list, 5) + 1.0) # # 20190527 5th percentile value from exp_sum_list
    p_val = (p_sum + 1.0) / (len(exp_sum_list) + 1.0)

    info = [obs, mu, sigma, 
    fold_change, p_val, 
    fold_change_95, fold_change_5, str(datetime.datetime.now())]
    return info, fold_changes

def bootstrapCI(fold_changes_list): # bootstrap C.I.s
    
    n = len(fold_changes_list) # size of distribution to bootstrap

    xbar = np.mean(fold_changes_list) # get the real median fold-change from 500 shuffles

    nboot = 10000 # resample 10000 times
    val = 0
    bs_means = []

    while val < nboot:

        bs_dist = np.random.choice(fold_changes_list, replace = True, size = n)
        bsmean = np.mean(bs_dist)
        bs_means.append(bsmean)
        val +=1
    
    bs = pandas.DataFrame(data = bs_means, index = np.arange(nboot), columns = ["bs_means"]) # make dataframe of bootstraps

    bs["deltas"] = bs.bs_means - xbar

    bs = bs.sort_values(by = "deltas", ascending= False)

    low = bs.deltas.quantile(0.025) # get 95th CI
    high = bs.deltas.quantile(0.975)
    print("CI", low, high)
    
    ci = xbar - [high, low]
    return ci

def main(argv):
    
    print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())))

    sid = (ANNOTATION_FILENAME.split("/")[-1]).split("_")[0] # name
    arch = ((ANNOTATION_FILENAME.split("/")[-1]).split("_")[-1]).split(".")[0] # name
    target = (TEST_FILENAME.split("/")[-1]).split(".")[0]
    print("outfile:", COUNT_FILENAME)

    # run initial intersection and save overlap counts
    obs_sum = calculateObserved(BedTool(ANNOTATION_FILENAME), BedTool(TEST_FILENAME))
    print("obs_sum", obs_sum)

    # create pool and run simulations in parallel
    pool = Pool(num_threads)
    partial_calcExp = partial(calculateExpected, BedTool(ANNOTATION_FILENAME), BedTool(TEST_FILENAME), SPECIES)
    exp_sum_list = pool.map(partial_calcExp, [i for i in range(ITERATIONS)])

    # wait for results to finish before calculating p-value
    pool.close()
    pool.join()

    # calculate empirical p value
    obs_p, fold_changes= calculateEmpiricalP(obs_sum, exp_sum_list)

    # bootstrap CI
    ci = bootstrapCI(fold_changes)

    # gather write data
    obs_p.extend(ci)
    info =[arch, sid, target]
    obs_p.extend(info)

    header = ['Observed', 'Expected', 'StdDev',
              'FoldChange_mean', 'p-value', 'fold_change_95_mean',
              'fold_change_5_mean', 'date_time',"ci_975", "ci_025", "arch", "sid", "target"]
    print(header)
    print(obs_p)

    csvfile = open(COUNT_FILENAME, "a", newline='')
    csv_writer = csv.writer(csvfile, delimiter = '\t')

    csv_writer.writerow(obs_p)
    csvfile.close()
    
if __name__ == "__main__":
    main(sys.argv[1:])

