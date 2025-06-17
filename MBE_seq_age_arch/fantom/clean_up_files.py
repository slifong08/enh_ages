import os
import glob
import pandas as pd
import subprocess

path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/"
fs = glob.glob("%s*.bed" % path)
len(fs)
print(fs[0])

#%%

def mkdir(outDir):
    if os.path.exists(outDir) == False:
        subprocess.call("mkdir %s" % outDir, shell = True)


#%% # remove very long enhancers, chrx
def clean_up(f, fid, sid):

    fname = "/cleaned_" + fid + ".bed"

    outpath = "/".join(f.split("/")[:-1]) +"/" + fid

    mkdir(outpath)

    outf = outpath + fname

    # now process the file
    df = pd.read_csv(f, sep = '\t', header = None)

    df = df.loc[df[0]!= "chrX"] # remove chrX
    df = df.loc[df[0]!= "chrY"] # remove chrY

    df[3] = df[2] - df[1]

    df = df.loc[df[3]<10000] # remove enhancers longer than 10kb

    df.to_csv(outf, sep ='\t', header = False, index = False)

    return outf


#%% subtract exons

def bed_subtract(inF, fid):

    inP = "/".join(inF.split("/")[:-1]) + "/"
    outP = "%snon-genic/" % inP

    mkdir(outP)

    ExonP = "/dors/capra_lab/users/fongsl/data/ensembl/"
    ExonF = "%sall_merged_exon.bed" % ExonP

    outF_noex = "%sno-exon_%s.bed" % (outP, fid)
    outF_ex = "%sexonOverlap_%s.bed" % (outP, fid)

    outP = "%s"
    # use -v argument to subtract exons from shuffle file.
    cmd = "bedtools intersect -a %s -b %s -v > %s" % (inF, ExonF, outF_noex)
    subprocess.call(cmd, shell = True)
    no_exon =  len(open(outF_noex).readlines(  ))


    cmd = "bedtools intersect -a %s -b %s > %s" %  (inF, ExonF, outF_ex)
    subprocess.call(cmd, shell = True)
    exon =  len(open(outF_ex).readlines(  ))

    return no_exon, exon

#%%

for f in fs:
    fid = (f.split("/")[-1]).split(".")[0]
    sid = "_".join(((f.split("/")[-1]).split(".")[0]).split("_")[:2])

    cleanf = clean_up(f, fid, sid)

    bed_subtract(cleanf, sid)

#%%
for f in fs:
    fid = (f.split("/")[-1]).split(".")[0]

    fP = "%s%s/" %(path, fid)

    oldP = "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/breaks/"
    oldF = "%s%s_enh_ages_enh_age_arch_summary_matrix.bed" % (oldP, fid)
    cmd = "mv %s %s" % (oldF, fP)
    subprocess.call(cmd, shell = True)
