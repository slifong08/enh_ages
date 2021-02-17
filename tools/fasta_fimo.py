import os, sys, traceback
import argparse
import glob
import pandas as pd
from collections import Counter

# TO RUN

# python script [input_file] [sample_id] [-i]

# python /dors/capra_lab/fongsl/enh_age/bin/age_enhancers.py UBERON_0002372_tonsil_expressed_enhancers.bed UBERON0002372


###
#   arguments
###
arg_parser = argparse.ArgumentParser(description="Calculate enhancer age.")

arg_parser.add_argument("bedfile", help='.bed file w/ full path')

arg_parser.add_argument("sample_id", help='str label for files')

arg_parser.add_argument("-i", "--iters", type=int, default=100,
                        help='number of simulation iterations; default=100')

TEST_ENH = args.region_file_1
SAMPLE_ID = args.sample_id
ITERATIONS = args.iters
TEST_PATH = "/".join(TEST_ENH.split("/")[:-1])
RESULTS_PATH = "%s/tfbs_motif/%s/" % (TEST_PATH, SAMPLE_ID)
mkdir = "mkdir %s" % RESULTS_PATH
os.system(mkdir)

#%%
def getfasta(test_enh, sample_id, test_path, results_path):
    os.chdir(test_path)

    chr_cmd = "awk '{print >$1\"_%s_temp.bed\"}' %s" % (sample_id, test_enh) # split test into chrN.bed files
    #os.system(chr_cmd)

    enh_chr_list = glob.glob("%s/chr*_%s_temp.bed" % (test_path, sample_id)) # glob chromosomes

    for enh_chr in enh_chr_list: # intersect test chrN.bed w/ syntenic block

        chr_num = (enh_chr.split("/")[-1]).split("_")[0]

        dna_path = "/dors/capra_lab/data/dna/human/hg19/"

        fasta_in = ("%s%s.fa.gz" % (dna_path, chr_num)) #gunzip -fi file
        gunzip = "gunzip %s" % fasta_in
        os.system(gunzip)

        fasta_in = ("%s%s.fa" % (dna_path, chr_num))
        fasta_out = "%s%s_%s.fa" % (results_path, sample_id, chr_num)

        cmd = "bedtools getfasta -fi %s -bed %s -fo %s" %(fasta_in, enh_chr, fasta_out)
        #os.system(cmd)

        cleanup_cmd = "rm %s" % enh_chr # clean up chr temp file
        os.system(cleanup_cmd)


    cat_fasta = "%s%s.fa" % (results_path, sample_id)
    cat_cmd = "cat %s%s_chr*.fa > %s" % (results_path, sample_id, cat_fasta)
    os.system(cat_cmd)

    cleanup_cmd = "rm %s%s_chr*.fa" % (results_path, sample_id)
    os.system(cleanup_cmd)

    return cat_fasta

def fimo(fasta, sample_id, results_path):
    motif_file = "/dors/capra_lab/data/tf_motif/meme/meme_motif_databases.12.17/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
    out_dir = "%sfimo_%s" % (results_path, sample_id)
    fimo_cmd = "fimo -o %s %s %s" % (out_dir, motif_file, fasta)
    os.system(fimo_cmd)
    print(fimo_cmd)

    return out_dir

def countgc(sequence):
    gc = []
    total = []
    letters = ["G", "C"]
    counts = Counter(sequence)

    for letter in letters:
        gc.append(int(counts[letter]))

    gc_sum = sum(gc)
    return gc_sum

def fimo2bed(fimo_outdir, sample_id, results_path):
    df = pd.read_csv("%s/fimo.txt" % fimo_outdir, sep = '\t')
    df["tf"] = df["#pattern name"].apply(lambda x: x.split("_")[0])

    df["chr_enh"] = df["sequence name"].apply(lambda x: x.split(":")[0])
    df["start_enh"] = df["sequence name"].apply(lambda x: (x.split(":")[1]).split("-")[0])
    df["end_enh"] = df["sequence name"].apply(lambda x: (x.split(":")[1]).split("-")[1])

    df["start_motif"] = df["start_enh"].astype(int) + df["start"].astype(int)
    df["end_motif"] = df["end_enh"].astype(int) + df["stop"].astype(int)

    df["GC"] = df["matched sequence"].apply(lambda x: countgc(x))
    df["motif_len"] = df["matched sequence"].apply(lambda x: len(x))

    df = df[['chr_enh',"start_motif", "end_motif", 'start_enh', 'end_enh', 'tf', '#pattern name',
         'sequence name', 'start', 'stop', 'strand', 'score',
         'p-value', 'q-value', 'GC', 'motif_len', 'matched sequence']]

    outfile = "%s/%s_motifs.bed" % (results_path, sample_id)
    df.to_csv(outfile, sep = "\t", header = False, index = False)


###
#   main
###

def main(argv):

    FASTA = getfasta(TEST_ENH, SAMPLE_ID, TEST_PATH, RESULTS_PATH) # make fasta
    FIMO_OUTDIR = fimo(FASTA, SAMPLE_ID, RESULTS_PATH) # make fimo
    fimo2bed(FIMO_OUTDIR, SAMPLE_ID, RESULTS_PATH) # fimo text to bed

if __name__ == "__main__":
    main(sys.argv[1:])
