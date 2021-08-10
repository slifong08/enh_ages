import os, sys
import subprocess

# add the path with the file containing all the paths to other files.
PATH = '/dors/capra_lab/users/fongsl/enh_ages/tyler/evolutionary_conservation'
sys.path.append(PATH)

import config # import the config.py file with all the files for tyler's project

all_ = config.all

#%%
"""

make a test file

"""
testDir = "/dors/capra_lab/users/fongsl/enh_ages/tyler/alignability/test/"
nShuf = 10

if os.path.exists(testDir) == False:
    cmd = f"mkdir {testDir}"

    subprocess.call(cmd, shell = True)

testFile = f"{testDir}test_{nShuf}.bed"

if os.path.exists(testFile) == False:
    cmd = f"shuf -n {nShuf} {all_} > {testFile}"
    subprocess.call(cmd, shell = True)
#%%

"""
liftOver and allow multiple mappings

"""
RE = "/dors/capra_lab/users/fongsl/tyler/results/liftover/"


def sort(bedfile, path, sid):
    """
    cut the first 5 fields of the bedfile

    sort the bedfile

    format the bed file into 5 columns [[chr start end enh_id sample_id]]
    and sort by coordinates
    """

    os.chdir(path)
    tempbed = os.path.join(path, f"temp_{sid}.bed")

    cut = f"cut -f 1-5 {bedfile} >  {tempbed}"
    subprocess.call(cut, shell = True)



    cmd = f"sort -k1,1 -k2,2 -k3,3 {tempbed} > t && mv t {tempbed}"
    print("standardizing Bed format")
    subprocess.call(cmd, shell=True)
    print(cut, "\n\n", cmd)
    return tempbed


def liftover(bedfile, path, from_build, to_build): # bedfile with full path

    # prepare annotations and files
    sid = (bedfile.split("/")[-1]).split(".")[0] # get the sample ID

    chainPath = "/dors/capra_lab/data/ucsc/liftOver/" # path to chain file
    chainf = f"{chainPath}{from_build}To{to_build}.over.chain.gz"

    #write to result files

    lifted = os.path.join(path, f"{sid}.liftOver.to.{to_build}.bed")  # name the liftover file
    notlifted = os.path.join(path, f"{sid}.notlifted.to.{to_build}.bed")

    if os.path.exists(lifted) == False: # check that you haven't run this.

        tempbed = sort(bedfile, path, sid)

        """
        liftover the formatted bedfile
        """

        cmd = f"liftOver {tempbed} {chainf} {lifted} {notlifted} -multiple"
        print("liftingOver", sid, "\n\n", cmd)
        subprocess.call(cmd, shell=True)
        print("done lifting.")

        """
        clean up temp
        """

        if os.path.getsize(lifted) >0:
            os.remove(tempbed)
            print("cleaned up temp file")

    return lifted


def get_multi_maps(lifted, path, from_build, to_build):

    sid = f"{from_build}_to_{to_build}"

    multi_counts = f"{sid}_multimap.txt" # make a new file name

    os.chdir(path) # go to the dir

    cmd = f"cut -f 4 {lifted} | sort | uniq -c > {multi_counts}" # find the number of overlaps

    print(cmd)
    subprocess.call(cmd, shell=True) # run in cmdline

    return multi_counts


#%%

PATH = "/dors/capra_lab/users/fongsl/tyler/data/liftover/"

file_dict = {
"hg38-RheMac8": os.path.join(PATH, "all_ocr_sorted.bed"),
"rheMac8-Hg38": os.path.join(PATH, "all_ocr_sorted.liftOver.to.RheMac8.bed")
}

results_dict = {}

#%%
for KEY, BED in file_dict.items():


    #RheMac10 to RheMac8
    FROM_BUILD = KEY.split("-")[0]
    TO_BUILD = KEY.split("-")[1]

    lifted = liftover(BED, PATH, FROM_BUILD, TO_BUILD)
    maps = get_multi_maps(lifted, PATH, FROM_BUILD, TO_BUILD)


#%%
lifted = '/dors/capra_lab/users/fongsl/tyler/data/liftover/all_ocr_sorted.liftOver.to.RheMac8.bed'

maps