# 20210318
# sarahfong
# this needs to be optimized for anyscript!
import argparse
import os, sys
import pandas as pd
"""
arg_parser = argparse.ArgumentParser(description=" describe argparse")

arg_parser.add_argument("fasta", help='fasta file w/ full path')
arg_parser.add_argument("-b", "--build", type=str, default='hg19', choices=['hg19', 'hg38'],
                        help='human genome assembly; default=hg19')


args = arg_parser.parse_args()

F = args.fasta
BUILD = args.build
PATH = "/".join(F.split("/")[:-1]) + "/"
"""
#%%


BUILD = "hg38"
PATH = "/dors/capra_lab/data/transposable_elements/repeatmasker/"
F = f"{PATH}{BUILD}.fa.out"


#%%


def replace_w_tabs(f, path, build):

    # replace white spaces with tabs
    outf = f"{path}{build}.txt"

    cmd = '''awk -v OFS="\t" '$1=$1' %s > %s''' % (f, outf)
    os.system(cmd)

    return outf


def make_bed(outf, path, build):
    # make a bedfile
    os.chdir(path)

    bedf = f"{path}{build}.bed"

    cmd_bed = f"cut -f 5,6,7,10,11 {outf} > {bedf}"

    os.system(cmd_bed)

    return bedf


def fasta_2_bed(f, path, build):

    outf = replace_w_tabs(f, path, build) # replace w tabs
    bedf = make_bed(outf, path, build) # make a bedfile

    return bedf


#%%

bedf = fasta_2_bed(F, PATH, BUILD)

df = pd.read_csv(bedf, sep = '\t', skiprows = 2, header = None)

df.head()
#%%
simple = df.loc[df[4] == "Simple_repeat"]
print(df.shape, simple.shape)
#%%
outsimplef = "/dors/capra_lab/projects/enhancer_ages/te/data/simple_repeats.bed"
simple.to_csv(outsimplef, sep = '\t', header = False, index = False)
