import glob
import os, sys

path = "/dors/capra_lab/projects/enhancer_ages/reilly15/data/breaks/"

fs = glob.glob("%snon-genic/no-exon_Hsap_brain_enhancers_reilly15_gapexcluded_parallel_breaks_enh_age_arch_summary_matrix.bed" % path)
fraction_overlap = 0.5
fs

all_beds = glob.glob("%sMmu*_brain_enhancers_reilly15_gapexcluded_parallel_breaks_enh_age_arch_summary_matrix.bed" % path)
fs = fs + all_beds
print(len(fs))
outpath = "/dors/capra_lab/projects/enhancer_ages/reilly15/data/multiintersect/"
#%%

# Entire Enhancer #
# -a

infile = "%snon-genic/no-exon_Hsap_brain_enhancers_reilly15_gapexcluded_parallel_breaks_enh_age_arch_summary_matrix.bed" % (path)

# ALL OTHER ROADMAP TISSUES #
    # -b

all_beds_minus_one = []
mask_file = "%sHsap_brain_enhancers_reilly15_gapexcluded_3000_1set_negs_parallel_breaks_enh_age_arch_summary_matrix.bed" % (path)
for all_bed_f in fs: # remove CORE file from intersection
    if all_bed_f != mask_file:
        all_beds_minus_one.append(all_bed_f)

all_beds_str = ' '.join(str(bed) for bed in all_beds_minus_one)

# Enhancer + ROADMAP Bedtools intersection
multi_cmd = "bedtools intersect -a %s -b %s -f %s -c > %strim-%s_multiintersect_hu_count.bed" %\
(infile, all_beds_str, fraction_overlap,  outpath, fraction_overlap)
print(multi_cmd)
os.system(multi_cmd)
