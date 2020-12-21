import os, sys
import subprocess
import glob

# arguments
iters, age, breaks, tfbs, shuf, run_test = 0, 0, 1, 0, 0, 1
run = 0

# function for running aging/age segment pipeline

def run_python(infile, iters, age, breaks, tfbs, shuf, run_test, run):

    cmd = "python /dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/awk_breaks.py\
    %s -i %d -a %d -b %d -t %d -sh %d -rt %d" % (infile, iters, age, breaks, tfbs,
    shuf, run_test)
    print(cmd)

    if run ==1:
        subprocess.call(cmd, shell = True)


#%%
# gather the files to run


path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/trimmed/shuffle/ages/"
fs = glob.glob("%s*.bed" % path)
print(len(fs))

missing = [ 'E050', 'E108', 'E074', 'E071', 'E019','E041', 'E042', 'E062',
 'E047', 'E029', 'E045', 'E127', 'E069', 'E115', 'E123', 'E032', 'E058', 'E055',
 'E039', 'E044', 'E072', 'E126', 'E048', 'E040', 'E124', 'E122', 'E046', 'E118' ]#, 'E116']


#%%
#run the pipeline

for f in fs:

    run_python(f, iters, age, breaks, tfbs, shuf, run_test, run)

#%%
# post-pipeline clean up

fs = glob.glob("/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E*/shuffle/cat*.bed")
for f in fs:
    path = "/".join(f.split("/")[:-1])

    temp = "%st.bed" % path

    add_enh_id = '''awk '{$(NF+1)=$1":"$2"-"$3 ; print $4"\t"$5}' %s > %s && mv %s %s''' % (f, temp, temp, f)
    print("add enh_id")
    subprocess.call(add_enh_id, shell = True)

#%%
