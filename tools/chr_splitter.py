import argparse
import subprocess


arg_parser = argparse.ArgumentParser(description="split .bed coordinates by chr")

arg_parser.add_argument("bedfile", help='bed file w/ full path')

args = arg_parser.parse_args()

F = args.bedfile

def chr_splitter(F):
    FNAME = (F.split("/")[-1]).split(".")[0]
    PATH = "/".join(F.split("/")[:-1])

    os.chdir(PATH)
    cmd = '''awk '{print >$1"_%s.bed"}' %s''' % FNAME
    subprocess.call(cmd, shell = True)
