#%% function to subtract exons

def bed_subtract(inF, fid):

    inP = "/".join(inF.split("/")[:-1]) + "/"
    outP = "%snon-genic/" % inP

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
