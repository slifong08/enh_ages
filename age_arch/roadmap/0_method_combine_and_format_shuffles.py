import glob
import os

path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/"
samples = glob.glob("%sshuffle/breaks/shuf-E*_age_breaks-*_parallel_breaks_enh_age_arch_summary_matrix.bed"%path)

sid_done = []

for sample in samples:
    sid = ((sample.split("/")[-1]).split("-")[1]).split("_")[0]
    print(sid)
    if sid not in sid_done:
        sid_done.append(sid)

        basef = "%sbreaks/ROADMAP_%s_enh_and_shuf_age_arch_summary_matrix.tsv" % (path, sid)
        shuffs = "%sshuffle/breaks/shuf-%s_age_breaks-*_parallel_breaks_enh_age_arch_summary_matrix.bed"%(path, sid)
        CONCAT_TEST = "%sbreaks/CONCAT.TXT" % path
        concat_cmd = "cat %s %s > %s && mv %s %s"  %(basef, shuffs, CONCAT_TEST, CONCAT_TEST, basef)
        os.system(concat_cmd)
