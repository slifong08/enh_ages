import os
import sys, traceback
import argparse
import datetime
import glob
import numpy as np
import pandas as pd

from functools import partial
import matplotlib.pyplot as plt
from multiprocessing import Pool

from statsmodels.stats import multitest

import datetime
import subprocess

last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())


### parse args ###
arg_parser = argparse.ArgumentParser(description="Calculate GWAS enrichment.")

arg_parser.add_argument("region_file_1", help='bed file 1 (enhancers to calculate enrichment) w/ full path')

arg_parser.add_argument("-g", "--gwas", type = str, default = "2018", choices=["2018", "2018LDex", "2019", "2019LDex"], help='choose GWAS bed file from 20180723 or 20190924')

args = arg_parser.parse_args()

FILE = args.region_file_1

GWAS_TEST = args.gwas
print("GWAS version:", GWAS_TEST)

PATH = "/".join(FILE.split("/")[:-1]) +"/"
SAMPLE_ID = (FILE.split("/")[-1]).split(".")[0]

print(FILE, "/n", PATH, "/n", SAMPLE_ID)


###
#   functions
###

def os_remove(files):
    cmd = "rm %s" % files
    os.system(cmd)

def mkdir(path):
    if os.path.isdir(path) == False:
        cmd = "mkdir %s" % path
        os.system(cmd)
        
def loadConstants(species):  # note chrom.sizes not used in current implementation | 2018.10.29
    return {'hg19': ("/dors/capra_lab/users/bentonml/data/dna/hg19/hg19_blacklist_gap.bed", "/dors/capra_lab/data/dna/human/hg19/hg19_trim.chrom.sizes"),
            'hg38': ("/dors/capra_lab/users/bentonml/data/dna/hg38/hg38_blacklist_gap.bed", "/dors/capra_lab/data/dna/human/hg38/hg38_trim.chrom.sizes"),
            'mm10': ("/dors/capra_lab/users/bentonml/data/dna/mm10/mm10_blacklist_gap.bed", "/dors/capra_lab/data/dna/mouse/mm10/mm10_trim.chrom.sizes")
            }[species]

def loadGWAS(gwas_test):  # note chrom.sizes not used in current implementation | 2018.10.29
    return {'2018': ("/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/gwasCatalog_2018-07-23_hg19_unique_cleaned_p5e-8.bed"),
            '2018LDex': ("/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/gwasCatalog_2018-07-23_hg19_unique_cleaned_LDEx_p5e-8.bed"),
            '2019': ("/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/gwasCatalog_2019-09-24_hg19_unique_cleaned_p5e-8.bed"),
            '2019LDex': ("/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/gwasCatalog_2019-09-24_hg19_unique_cleaned_LDex_p5e-8.bed")
            }[gwas_test]

def calculateObserved(f, test, sample_id, path):
    obs_out =  "%s%s_x_GWAS.bed" % (path, sample_id) # outfile
    
    obs_intersect = "bedtools intersect -a %s -b %s -wo > %s" % (f, test, obs_out)

    subprocess.call(obs_intersect, shell = True) 
    
    
    count = len(open(obs_out).readlines(  )) 
    os_remove(obs_out)
    
    return count
    
def calculateExpected(f, test, species, sample_id, path, iteration):
    bl, g = loadConstants(species)  # note CHROM_SZ not used
        
    outf = "%sshuf-%s-%s.bed" % (path, sample_id, iteration) # outfile

    #bedtools cmd hg19, exclude blacklist, same chromosome, noOverlapping bits
    cmd = "bedtools shuffle -i %s -g %s -excl %s -chrom -noOverlapping > %s" % (f, g, bl, outf) 
        
    subprocess.call(cmd, shell = True)
    
    exp_out =  "%sshuf-%s-%s_x_GWAS.bed" % (path, sample_id, iteration) # outfile
    
    exp_intersect = "bedtools intersect -a %s -b %s -wo > %s" % (outf, test, exp_out)

    subprocess.call(exp_intersect, shell = True) 
    
    count = len(open(exp_out).readlines(  )) 
    
    os_remove(exp_out)
    os_remove(outf)
    
    return count


def calculateEmpiricalP(obs, exp_sum_list):
    
    mu = np.median(exp_sum_list) # 20190523 median
    sigma = np.std(exp_sum_list)
    dist_from_mu = [exp - mu for exp in exp_sum_list]
    p_sum = sum(1 for exp_dist in dist_from_mu if abs(exp_dist) >= abs(obs - mu))

    fold_change = (obs + 1.0) / (mu + 1.0) # fold chage of the median
    fold_changes = list((obs + 1.0) / (m + 1.0) for m in exp_sum_list)

    p_val = (p_sum + 1.0) / (len(exp_sum_list) + 1.0)

    info = [obs, mu, sigma, 
    fold_change, p_val, 
     str(datetime.datetime.now())]
    return info, fold_changes

def bootstrapCI(fold_changes_list): # bootstrap C.I.s
    
    n = len(fold_changes_list)

    xbar = np.mean(fold_changes_list) # get the real median fold-change from 500 shuffles

    nboot = 10000 # resample 10000 times
    val = 0
    bs_means = []

    while val < nboot:

        bs_dist = np.random.choice(fold_changes_list, replace = True, size = n)
        bsmean = np.mean(bs_dist)
        bs_means.append(bsmean)
        val +=1
    
    bs = pd.DataFrame(data = bs_means, index = np.arange(nboot), columns = ["bs_means"])

    bs["deltas"] = bs.bs_means - xbar

    bs = bs.sort_values(by = "deltas", ascending= False)

    low = bs.deltas.quantile(0.025)
    high = bs.deltas.quantile(0.975)
    print(low, high)
    
    ci = xbar - [high, low]
    return ci

def split_simple_complex(f, sample_id, path):

    os.chdir(path)
    cmd = "awk '{print >$6\"_%s_ROADMAP_complexenh.bed\"}' %s " % ( sample_id, f, )
    print(cmd)
    subprocess.call(cmd, shell=True)
    
    complexenhf, simplef = "%s1_%s_ROADMAP_complexenh.bed" % (path, sample_id), "%s0_%s_ROADMAP_complexenh.bed"% (path, sample_id,)
    f_list = [complexenhf, simplef]
    
    return f_list

def get_results(f_list, test, sample_id, path, species, iterations, gwas):
    f_dict = {}
    for f in f_list:

        count = len(open(f).readlines(  ))
        print(count)

        print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())))

        
        arch = ((f.split("/")[-1]).split("_")[0])

        target = (test.split("/")[-1]).split(".")[0]
        key = "_".join([sample_id, arch])
        print(sample_id, arch)

        # run initial intersection and save overlap counts
        obs_sum = calculateObserved(f, test, sample_id, path)
        print("obs_sum", obs_sum)

        # create pool and run simulations in parallel
        pool = Pool(num_threads)
        partial_calcExp = partial(calculateExpected, f, test, SPECIES, sample_id, path)
        exp_sum_list = pool.map(partial_calcExp, [i for i in range(iterations)])

        # wait for results to finish before calculating p-value
        pool.close()
        pool.join()

        # remove iterations that throw bedtools exceptions
        final_exp_sum_list = [x for x in exp_sum_list if x >= 0]
        exceptions = exp_sum_list.count(-999)

        # calculate empirical p value
        if exceptions / iterations <= .1:
            obs_p, fold_changes = calculateEmpiricalP(obs_sum, exp_sum_list)

            ci = bootstrapCI(fold_changes)
            obs_p.extend(ci)

            info =[arch, sample_id, target]
            obs_p.extend(info)

            header = ['Observed', 'Expected', 'StdDev',
                      'FoldChange_Med', 'p-value',  'date_time',"ci_975", "ci_025", "arch", "sid", "target"]

            enh = pd.DataFrame(data = [obs_p], columns = header)
            enh["arch"] = arch

            f_dict[key] = enh
        else:
#                print(f'iterations not completed: {exceptions}\nresulted in nonzero exit status', file=sys.stderr)
            cleanup()
            sys.exit(1)
        os_remove(f)
    df = pd.concat(f_dict.values())
    outf = "/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/ROADMAP_%s_dataset_GWAS_%s_enrichment.tsv" % (sample_id, gwas)
    df.to_csv(outf, sep ='\t', header = True, index = False)
    return outf

###
#   main
###
def main(argv):
    print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())[:20]))
    
SPECIES = "hg19"
ITERATIONS = 100

num_threads = int(os.getenv('SLURM_CPUS_PER_TASK', 1))


TEST_FILENAME = loadGWAS(GWAS_TEST)

# split roadmap files into simple-only and complex-only bed files

f_list = split_simple_complex(FILE, SAMPLE_ID, PATH)

outf = get_results(f_list, TEST_FILENAME, SAMPLE_ID, PATH, SPECIES, ITERATIONS, GWAS_TEST)

if __name__ == "__main__":
    main(sys.argv[1:])    