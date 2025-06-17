#20201028
# evaluate roadmap enhancer features without exons
import glob
import os
import subprocess

path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/non-genic/breaks/"
fs = glob.glob("%sno-exon_E*_parallel_breaks_enh_age_arch_summary_matrix.bed" % path)
fs[0]
#%%

for F in fs:
    sid = (F.split("/")[-1]).split("_")[1]

    print(sid)
    script = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/break_stats_summary_matrix.py"
    cmd = "python %s %s" % (script, F)
    print(cmd)
    os.system(cmd)



#%%
shufpath = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/non-genic/shuffle/breaks/"
shuffs = glob.glob("%sno-exon_*.bed" % shufpath)
len(shuffs)
#%%
redo_list = ['E098', 'E102', 'E103', 'E111', 'E055', 'E008', 'E034', 'E063',
       'E012', 'E105', 'E078', 'E029', 'E114', 'E003', 'E066', 'E056']
#%%
for F in shuffs:
    sid = (F.split("/")[-1]).split("_")[1]



    print(sid)
    script = "/dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/break_stats_summary_matrix.py"
    cmd = "python %s %s" % (script, F)
    print(cmd)
    os.system(cmd)
#%%
print(sid)
