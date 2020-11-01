#20201028
# evaluate roadmap enhancer features without exons
import glob
import os
import subprocess

path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/breaks/non-genic/"
fs = glob.glob("%sno-exon_ROADMAP_E*_enh_and_shuf_age_arch_summary_matrix.bed" % path)
fs[0]
#%%
for F in fs:
    script = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/break_stats_summary_matrix.py"
    cmd = "python %s %s" % (script, F)
    print(cmd)
    os.system(cmd)


#%%
cmd
