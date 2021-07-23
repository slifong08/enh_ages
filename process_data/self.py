# 20120712
# sarahfong
# massage self bed to remove duplicate rows, count the number of copies.
# save as a clean file
import os, sys
import pandas as pd
import subprocess

PATH = "/dors/capra_lab/data/ucsc/hg38/self/"
SELF_F = f"{PATH}hg38_repeats_self.bed"
CLEAN_F = f"{PATH}hg38_repeats_self_coor.bed"

cmd = f"cut -f 1,2,3 {SELF_F} | sort | uniq > {CLEAN_F}"
subprocess.call(cmd, shell = True)
