# 2019-05-27
# sarahfong

# using Mary Lauren's script to calculate GTEx eQTL enrichment using the script calculate_enrichment.py
# FANTOM enhancer architectures were used in this analysis. 
# I copied this script from /dors/capra_lab/users/bentonml/resources/bin/calculate_enrichment.py

import glob
import os, sys

tissue = sys.argv[1] # full path
target = sys.argv[2] # full path
outpath = sys.argv[3] # full path

tissue_path = '/'.join((tissue.split("/")[:-1]))

tissue_id = (tissue.split("/")[-1]).split("*")[0]
# sample tissue: "/dors/capra_lab/projects/enhancer_ages/fantom/data/fantom_enh_age/architecture_coordinates/fibroblast-of-lymphatic-vessel_FANTOM_*.bed"

target_id = (target.split("/")[-1]).split(".")[0]
# sample target id: "/dors/capra_lab/projects/enhancer_ages/eqtl/gtex_v7/Whole_Blood_v7_signif_variant_gene_pairs.bed"
# gwas = "/dors/capra_lab/projects/enhancer_ages/gwas_catalog/gwasCatalog_2018-07-23_hg19_unique_cleaned.bed" # 49405 unique snps

# Mary Lauren's script (copy)
ml_script = "/dors/capra_lab/users/fongsl/enh_age/bin/calculate_enrichment.py"

# Glob FANTOM enhancer architecture files
arch_list = glob.glob("%s" % (tissue))

# Make an output file
outfile = "%s%s_x_%s.txt" % (outpath, tissue_id, target_id )
touch = "touch %s" % outfile
os.system(touch) # make a new output file

for file in arch_list:
    file_name = (file.split("/")[-1]).split(".")[0]
    cmd = "python %s %s %s --print_counts_to %s -i 500" %(ml_script, file, target, outfile)
    print(file_name, target_id)
    os.system(cmd)