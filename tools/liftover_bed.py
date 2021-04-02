import os, sys
import subprocess



def liftover(bedfile, path, from_build, to_build): # bedfile with full path

    # prepare annotations
    sid = (bedfile.split("/")[-1]).split(".")[0] # get the sample ID

    ### sort the bedfile ###
    tempbed = f"{path}temp_{sid}.bed" # format the bed file into 5 columns

    # [[chr start end enh_id sample_id]] and sort by coordinates

    cmd = f"sort -k1,1 -k2,2 -k3,3 {bedfile} > {tempbed}"

    print("standardizing Bed format")
    subprocess.call(cmd, shell=True)

    ### liftover the formatted bedfile ###

    chainPath = "/dors/capra_lab/data/ucsc/liftOver/" # path to chain file
    chainf = f"{chainPath}{from_build}To{to_build}.over.chain.gz"

    #write to result files
    lifted = f"{path}{sid}. liftOver.to.{to_build}.bed" # name the liftover file
    notlifted = f"{path}{sid}.notlifted.to.{to_build}.bed" # name the notlifted file

    cmd = f"liftOver {tempbed} {chainf} {lifted} {notlifted}"
    print("liftingOver", sid)
    subprocess.call(cmd, shell=True)
    print("done lifting")


    ### clean up temp ###
    if os.path.getsize(lifted) >0:
        os.remove(tempbed)
        print("cleaned up temp file")

    return lifted
