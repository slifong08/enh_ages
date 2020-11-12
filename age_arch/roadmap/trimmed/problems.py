import os
import pandas as pd
import glob
import subprocess

path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/trimmed/shuffle/ages/"

sids = ['E032', 'E074', 'E118', 'E127',
 'E124', 'E122', 'E072', 'E069', 'E042', 'E047','E019', 'E048', 'E041', 'E126',
 'E045', 'E115', 'E046', 'E034','E071', 'E108', 'E029', 'E040', 'E062','E039','E123']
for sid in sids:
    fs = glob.glob("%sshuf-trimmed-310_%s_*_parallel_breaks.bed.gz" % (path, sid))

    for f in fs:
        unzipped = f.split(".gz")[0]
        cmd = "gunzip %s" % f
        subprocess.call(cmd, shell = True)

    out = "%sshuf-trimmed-310_%s_parallel_breaks.bed" %(path, sid)
    cat_cmd = "cat %sshuf-trimmed-310_%s_*_parallel_breaks.bed.gz > %s" % (sid, path, out)
    subprocess.call(cat_cmd, shell = True)
    print(sid)
    break
