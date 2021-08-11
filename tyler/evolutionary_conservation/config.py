# all the files I have used for this human acceleration analysis:


#"/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/species_specific_10k/"

"""
# all elements
"""

all = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/all/GG-LL_all_OCRs.bed"

hu_specific = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/hu_specific/GM12878_specific.bed"
rhe_specific = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/rhe_specific/LCL8664_specific.bed"

"""
# subtracting TEs

"""

all_noTE = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/all/subtract_te/all_wo_te.bed"
hu_specific_noTE = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/hu_specific/subtract_te/hu-specific_wo_te.bed"
rhe_specific_noTE = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/rhe_specific/subtract_te/rhe-specific_wo_te.bed"

"""
# positive control
"""

hars = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/hars/hars_merged_PMIDs_hg38.bed"

"""
# negative control
"""

phastCons = "/dors/capra_lab/users/fongsl/tyler/data/CON_ACC/phastCons/phastConsElements100way_hg38.bed.gz"


"""
models and maf files (chr sep)
"""

mod_maf = {
"20way":["/dors/capra_lab/data/ucsc/hg38/multiz20way/hg38.phastCons20way.mod",
"/dors/capra_lab/data/ucsc/hg38/multiz20way/maf/*.maf.gz"],
"30way":["/dors/capra_lab/data/ucsc/hg38/multiz30way/hg38.phastCons30way.mod",
"/dors/capra_lab/data/ucsc/hg38/multiz30way/maf/*.maf.gz"],
"100way":["/dors/capra_lab/data/ucsc/hg38/multiz100way/hg38.phastCons100way.mod",
"/dors/capra_lab/data/ucsc/hg38/multiz20way/maf/*.maf.gz"],
}

"""
multiple mappings using liftOver
"""

multiple_maps = {
"Hg38-rheMac8": "/dors/capra_lab/users/fongsl/tyler/data/liftover/hg38_to_RheMac8_multimap.txt",
"RheMac8-hg38": "/dors/capra_lab/users/fongsl/tyler/data/liftover/rheMac8_to_Hg38_multimap.txt"
}


"""
software
"""
PHAST_PATH = "/dors/capra_lab/bin/"
phylop = f"{PHAST_PATH}./phyloP"
treedoctor = f"{PHAST_PATH}./tree_doctor"
