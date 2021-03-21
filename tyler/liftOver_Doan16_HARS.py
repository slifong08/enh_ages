import os
import subprocess

HARPATH = "/dors/capra_lab/data/human_accelerated_regions/"
HARFILE = "Doan2016_HARs.bed"
HARF = os.path.join(HARPATH, HARFILE)

#%%

def liftover(bedfile, path): # bedfile with full path

    sid = "Doan2016_HARs" # get the sample ID
    tempbed = "%stemp_%s.bed" %(path, sid) # format the bed file into 5 columns

    cmd = '''awk -F'\t' -v OFS="\t" '{ print $1, $2, $3, $4, $5, $6, $7}'\
    %s | sort -k1,1 -k2,2 -k3,3 > %s''' % (bedfile, tempbed)
    subprocess.call(cmd, shell=True)

    chainPath = "/dors/capra_lab/data/ucsc/liftOver/" # path to chain file
    chainf = "%shg19ToHg38.over.chain.gz" % chainPath # Hg19 to Hg38 chain
    lifted = "%s%s.liftOver.to.hg38.bed" % (path, sid) # name the liftover file    
    notlifted = "%s%s.notlifted.to.hg38.bed" % (path, sid) # name the notlifted file

    cmd = "liftOver %s %s %s %s" % (tempbed, chainf, lifted, notlifted)

    subprocess.call(cmd, shell=True)
    print("liftedOver", sid)
#%%
liftover(HARF, HARPATH)
