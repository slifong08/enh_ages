import glob
import os
import pandas as pd

path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/"
samples = glob.glob("%sshuffle/breaks/shuf-*_age_breaks-*_parallel_breaks_enh_age_arch_summary_matrix.bed"%path)

sid_done = []

#%% Remove the 4th column. provides peak reads that I don't need and will mess up analysis
for sample in samples:

    df = pd.read_csv(sample, sep = '\t', header = None,
    usecols=[0,1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17])
    df.to_csv(sample, sep = '\t', header = False, index = False)


#%%
for sample in samples:
    sid = ((sample.split("/")[-1]).split("_")[1]).split("-")[0]
    print(sid)

    if sid not in sid_done:
        sid_done.append(sid)

        basef = "%sbreaks/no-exon_ROADMAP_%s_enh_and_shuf_age_arch_summary_matrix.bed" % (path, sid)
        shuffs = "%snon-genic/shuffle/breaks/shuf-no-exon_%s-*_parallel_breaks_enh_age_arch_summary_matrix.bed"%(path, sid)
        CONCAT_TEST = "%snon-genic/breaks/CONCAT.TXT" % path
        concat_cmd = "cat %s %s > %s && mv %s %s"  %(basef, shuffs, CONCAT_TEST, CONCAT_TEST, basef)
        os.system(concat_cmd)
#%%
sid_done

#%%
sid
