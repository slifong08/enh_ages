#20201021

# these files are originally from Ling's directory: /dors/capra_lab/projects/EnhancerCodeConservation/data/reilly15/
#%%
import glob
import os
import subprocess


path = "/dors/capra_lab/projects/enhancer_ages/reilly15/data/"
ms_bed = "%sMmus_brain_enhancers_reilly15_gapexcluded.bed4" % path
rh_bed = "%sMmul_brain_enhancers_reilly15_gapexcluded.bed4" % path

#%% I tested a bunch of other mouse files to see how much overlap I should expect from
# these mouse coordinates. I was shocked by how little the mouse coordinates overlapped.

# one ChIP-seq experiment at one time point, very raw data, no consensus.
other_ms = "%sGSM1489701_Mm_e11_H3K27ac_rep1_regions.bed" %path

# consensus peaks w/ gaps included.
other2_ms = "%sMmus_brain_enhancers_reilly15_andgap.bed" %path

# all peaks at all time points, regardless of consensus between replicates at timepoints
all_Mm_H3K27ac = '/dors/capra_lab/data/enhancers/reilly15/all_Mm_H3K27ac.bed'


#%% get the chain files. 


chainPath = "/dors/capra_lab/data/ucsc/liftOver/" # path to chain file

Mmus_chainf = "%smm9ToHg19.over.chain.gz" % chainPath # rhesus chain file
Mmul_chainf = "%srheMac2ToHg19.over.chain.gz" % chainPath # rhesus chain file


#%%


def liftover(bedfile, path, chainf):
    sid = ((bedfile.split("/")[-1]).split(".")[0])
    tempbed = "%stemp_%s.bed" % (path, sid)

    # [[chr start end enh_id sample_id]] and sort by coordinates

    cmd = '''awk 'OFS=" " {print $1"\t", $2"\t", $3"\t", $10"\t", $11}'\
    %s | tr -d " "| sort -k1,1 -k2,2 -k3,3 > %s''' % (bedfile, tempbed)

    print("standardizing Bed format")

    subprocess.call(cmd, shell=True)

    ### liftover the formatted bedfile ##

    build = (chainf.split("To")[1]).split(".")[0] # get the target build

    lifted = "%s%s.liftOver.to.%s.bed" % (path, sid, build) # name the liftover file

    notlifted = "%s%s.notlifted.to.%s.bed" % (path, sid, build) # name the notlifted file

    cmd = "liftOver %s %s %s %s" % (tempbed, chainf, lifted, notlifted)


    subprocess.call(cmd, shell=True)

    print("liftedOver", sid)

    ### clean up temp ###

    cmd = "rm %s" % tempbed

    subprocess.call(cmd, shell=True)

    print("cleaned up temp file")

    return lifted
#%%
liftover(ms_bed, path, Mmus_chainf)
#%%
liftover(other_ms, path, Mmus_chainf)

#%%
liftover(other2_ms, path, Mmus_chainf)

#%%
liftover(rh_bed, path, Mmul_chainf)
#%%
liftover(all_Mm_H3K27ac, path, Mmus_chainf)
