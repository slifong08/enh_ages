# 20201007
# sarahfong
# this needs to be optimized for anyscript!

import os, sys
import pandas as pd

path = "/dors/capra_lab/data/transposable_elements/repeatmasker/"
f = "%shg19.fa.out" % path
outf = "%shg19.txt" % path


# replace white spaces with tabs
cmd = '''awk -v OFS="\t" '$1=$1' %s > %s''' % (f, outf)
os.system(cmd)

# make a bedfile
os.chdir(path)
bedf = "%shg19.bed" %path
cmd_bed = "cut -f 5,6,7,10,11 %s > %s" % (outf, bedf)
os.system(cmd_bed)

#%%
df = pd.read_csv(bedf, sep = '\t', skiprows = 2, header = None)
#%%
df.head()
simple = df.loc[df[4] == "Simple_repeat"]
print(df.shape, simple.shape)
#%%
outsimplef = "/dors/capra_lab/projects/enhancer_ages/te/data/simple_repeats.bed"
simple.to_csv(outsimplef, sep = '\t', header = False, index = False)
