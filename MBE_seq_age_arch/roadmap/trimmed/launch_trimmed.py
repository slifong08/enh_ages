import glob
import os, sys
import subprocess
import time

PATH = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/"

FS = glob.glob(f"{PATH}/Hsap_H3K27ac_plus_H3K4me3_minus_E*/non-genic/trimmed/trim_310_E*.bed")
done = ['/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E092/non-genic/trimmed/trim_310_E092.bed',
 '/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E016/non-genic/trimmed/trim_310_E016.bed']
#%%
val = 10
for F in FS[10:]:
    if F not in done:
        slurm_script = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/awk_breaks.slurm"
        cmd = f"sbatch {slurm_script} {F} 10 1 1 1 1 hg19"
        subprocess.call(cmd, shell = True)
        print(F)
        val +=1
    else:
        done.append(F)
    if val%10 == 0 and val !=0:
        time.sleep(1200)

#%%
FS = glob.glob(f"{PATH}/Hsap_H3K27ac_plus_H3K4me3_minus_E*/non-genic/no-exon_E*.bed")

new_done = []
val = 0
for F in FS[20:]:
    if F not in done and "summary" not in F:
        slurm_script = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/awk_breaks.slurm"
        cmd = f"sbatch {slurm_script} {F} 10 1 1 1 1 hg19"
        subprocess.call(cmd, shell = True)
        print(F)
        val +=1
    else:
        new_done.append(F)
    #if val%20 == 0 and val !=0:
    #    time.sleep(1200)
