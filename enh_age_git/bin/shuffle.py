import os
import sys, traceback
import argparse
from datetime import datetime
from functools import reduce
from itertools import groupby
from functools import partial
from multiprocessing import Pool



###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate enhancer age.")

arg_parser.add_argument("region_file_1", help='bed file 1 (enhancers to age) w/ full path')

arg_parser.add_argument("sample_id", help='str label for files')

arg_parser.add_argument("-i", "--iters", type=int, default=100,
                        help='number of simulation iterations; default=100')

arg_parser.add_argument("-s", "--species", type=str, default='hg19', choices=['hg19', 'hg38', 'mm10'],
                        help='species and assembly; default=hg19')

arg_parser.add_argument("-n", "--num_threads", type=int,
                        help='number of threads; default=SLURM_CPUS_PER_TASK or 1')

args = arg_parser.parse_args()

TEST_ENH = args.region_file_1
SAMPLE_ID = args.sample_id
ITERATIONS = args.iters
SPECIES = args.species
TEST_PATH = "/".join(TEST_ENH.split("/")[:-1])

SHUFFLE_ID = "shuf-%s" % SAMPLE_ID

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

def calculateExpected(test_enh, shuffle_id, test_path, species, iters):
    
    a = datetime.now()
    time = "%s%s%s" % (a.hour,  a.minute, a.microsecond)
    BLACKLIST, CHROM_SZ = loadConstants(species)  # note CHROM_SZ not used
    
    exp_sum = 0

    shuffle_path = "%s/shuffle" % test_path # make a shuffle path file
    if os.path.exists(shuffle_path) == False:
        os.system("mkdir %s" % shuffle_path)

    rand_out = '%s/rand_file_%s_%s.bed'% (shuffle_path, shuffle_id, time) # make shuffle file
    
    #bedtools cmd hg19, exclude blacklist, same chromosome, noOverlapping bits
    cmd = "bedtools shuffle -i %s -g %s -excl %s -chrom -noOverlapping > %s" % (test_enh, CHROM_SZ, BLACKLIST, rand_out) 

    os.system(cmd)
    
    
###
# main
###

def main(argv):
    print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.now())[:20]))

    a = datetime.now()
    TIME = "%s%s%s" % (a.hour,  a.minute, a.second)

    # create pool and run simulations in parallel
    pool = Pool(num_threads)
    partial_calcExp = partial(calculateExpected, TEST_ENH, SHUFFLE_ID, TEST_PATH, SPECIES)
    exp_sum_list = pool.map(partial_calcExp, [i for i in range(ITERATIONS)])

if __name__ == "__main__":
    main(sys.argv[1:])


