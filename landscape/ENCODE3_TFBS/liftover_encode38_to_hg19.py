### SUMMARY ###

# makes directory.
# sorts encode3 TFBS clustered with cell lines file
# liftOvers the entire file between hg builds
# write liftOver and notlifted files.
import os

TO = 19
FROM = 38

ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/"
ENCODEFILE = "trimmed_encRegTfbsClusteredWithCells.hg38.bed"

ENCODE = os.path.join(ENCODEPATH, ENCODEFILE)

CHAINPATH = "/dors/capra_lab/data/ucsc/liftOver/"
CHAINFILE = "hg%sToHg%s.over.chain.gz" % (FROM, TO)

CHAIN = os.path.join(CHAINPATH, CHAINFILE)

OUTDIR = "liftOver_hg%s" % TO

OUTPATH = os.path.join(ENCODEPATH, OUTDIR)

if os.path.exists(OUTPATH) == False:
    os.mkdir(OUTPATH)

#%%


def liftover(bedfile, bedpath, outpath, chain, to):

    sid = (bedfile.split("/")[-1]).split(".")[0]
    temp = "%s/temp_%s.bed" % (outpath,sid)

    cmd = '''awk 'OFS=" " {print $1"\t", $2"\t", $3"\t",\
    $4"\t", $5"\t", $7"\t", $9}' %s \
    | tr -d " "| sort -k1,1 -k2,2 -k3,3 > %s''' % (bedfile, temp)
    print(cmd)

    subprocess.call(cmd, shell = True)

    lifted = "%s/%s.liftOver.to.hg%s.bed" % (outpath, sid, to)
    notlifted = "%s/%s.notlifted.to.hg%s.bed" % (outpath, sid, to)
    cmd = "liftOver %s %s %s %s" % (temp, chain, lifted, notlifted)
    print(cmd)
    subprocess.call(cmd, shell = True)

    cmd = "rm %s" % temp
    subprocess.call(cmd, shell = True)
    return lifted

#%%

lifted = liftover(ENCODE, ENCODEPATH, OUTPATH, CHAIN, TO)

print(lifted)
