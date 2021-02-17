import os

TO = 19
FROM = 38

ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/"
ENCODEFILE = "HepG2.bed"
ENCODEFILE = "trimmed_encRegTfbsClusteredWithCells.hg38.bed"

ENCODE = os.path.join(ENCODEPATH, ENCODEFILE)

CHAINPATH = "/dors/capra_lab/data/ucsc/liftOver/"
CHAINFILE = "hg%sToHg%s.over.chain.gz" % (FROM, TO)

CHAIN = os.path.join(CHAINPATH, CHAINFILE)


#%%


def liftover(bedfile, bedpath, chain, to):

    sid = (bedfile.split("/")[-1]).split(".")[0]
    temp = "%stemp_%s.bed" % (bedpath,sid)

    cmd = '''awk 'OFS=" " {print $1"\t", $2"\t", $3"\t",\
    $4"\t", $5"\t", $6"\t", $7}' %s \
    | tr -d " "| sort -k1,1 -k2,2 -k3,3 > %s''' % (bedfile, temp)
    print(cmd)

    subprocess.call(cmd, shell = True)

    lifted = "%s%s.liftOver.to.hg%s.bed" % (bedpath, sid, to)
    notlifted = "%s%s.notlifted.to.hg%s.bed" % (bedpath, sid, to)
    cmd = "liftOver %s %s %s %s" % (temp, chain, lifted, notlifted)
    print(cmd)
    subprocess.call(cmd, shell = True)

    cmd = "rm %s" % temp
    subprocess.call(cmd, shell = True)
    return lifted

#%%

lifted = liftover(ENCODE, ENCODEPATH, CHAIN, TO)

print(lifted)
