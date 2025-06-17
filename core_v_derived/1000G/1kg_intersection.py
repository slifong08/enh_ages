import argparse
import glob
import os, sys
import subprocess


arg_parser = argparse.ArgumentParser(description= "intersect bed with 1KG global population variants")

arg_parser.add_argument("-b", "--bed", help='.bed file w/ full path')

arg_parser.add_argument("-o", "--outdirectory", help='directory to dump gc information')

args = arg_parser.parse_args()

BEDF = args.bed
OUTPATH = args.outdirectory

###
# functions
###

def get_1kg_path():
    
    thou_path = "/dors/capra_lab/projects/enhancer_ages/1000g/data/maf/"
    common_var_beds = os.path.join(thou_path, f"trimmed.chr*.phase3_shapeit2_mvncall_integrated_v5a.maf.bed.gz")
    var_chrs = glob.glob(common_var_beds) ## get chr files

    return thou_path, var_chrs

# Functions

def intersect_1kg(file, outpath):
    
    """
    return .bed file intersected w/ 1kg alleles
    1. get name of file
    2. get thousand genomes bed files
    3. make a list of files to return
    4. do 1kg x bed intersection if the outfile does not already exist
    5. concatenate files
    6. 
    
    """
    #1
    sid = (file.split("/")[-1]).split(".")[0]
    
    #2
    thou_path, common_var_beds = get_1kg_path()
    
    #3 
    outfs = []
    
    #4
    for chr_var in common_var_beds:
        
        chr_num = (chr_var.split("/")[-1]).split(".")[1]
        
        outfile = os.path.join(outpath, f"{sid}_x_1kg_trimmed-{chr_num}.bed")
        
        bed_cmd = f"bedtools intersect -a {file} -b {chr_var}  -wa -wb> {outfile}" 
    
        if os.path.exists(outfile) is False:

            subprocess.call(bed_cmd, shell = True)

            print("intersecting", bed_cmd)

        else:
            print("already intersected", outfile)
        outfs.append(outfile)

    #5        
    to_cat =  os.path.join(outpath, f"{sid}_x_1kg_trimmed-chr*.bed")
    out_cat = os.path.join(outpath, f"{sid}_x_1kg_trimmed.bed")
    cmd = f"cat {to_cat} > {out_cat}"
    print("cat!", cmd)
    subprocess.call(cmd, shell = True)
    
    #6 
    #subprocess.call(f"rm {to_cat}", shell = True)
    
    
    return outfs


def main(argv):
                          
    intersect_1kg(BEDF, OUTPATH)
    
if __name__ == "__main__":
    main(sys.argv[1:])
