


def liftover(bedfile, path): # bedfile with full path

    sid = (bedfile.split("/")[-1]).split(".")[0] # get the sample ID
    
    ### format the bedfile ###
    
    tempbed = "%stemp_%s.bed" %(path, sid) # format the bed file into 5 columns
    
    # [[chr start end enh_id sample_id]] and sort by coordinates
    
    cmd = '''awk 'OFS=" " {print $1"\t", $2"\t", $3"\t", $10"\t", $11}'\
    %s | tr -d " "| sort -k1,1 -k2,2 -k3,3 > %s''' % (bedfile, tempbed)
    
    print("standardizing Bed format")
    
    subprocess.call(cmd, shell=True)
    
    ### liftover the formatted bedfile ###
    
    chainPath = "/dors/capra_lab/data/ucsc/liftOver/" # path to chain file
    
    chainf = "%shg19ToHg38.over.chain.gz" % chainPath # Hg19 to Hg38 chain



    lifted = "%s%s.liftOver.to.hg38.bed" % (path, sid) # name the liftover file
    
    notlifted = "%s%s.notlifted.to.hg38.bed" % (path, sid) # name the notlifted file
    
    cmd = "liftOver %s %s %s %s" % (tempbed, chainf, lifted, notlifted)
    
    subprocess.call(cmd, shell=True)
    
    print("liftedOver", sid)



    ### clean up temp ###

    

    cmd = "rm %s" % tempbed



    subprocess.call(cmd, shell=True)
    
    print("cleaned up temp file")



    return lifted
