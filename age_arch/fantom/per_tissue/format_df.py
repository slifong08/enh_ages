import os
import glob
import pandas as pd

#%%

syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t',
usecols = ["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]) # read the file

syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

#%%
def postformat_breaksdf(df, datatype, syn_gen_bkgd):

    print("POST BREAK ASSEMBLY FORMATTING")

    df = df.loc[df.chr_enh !="chrX"] # autosomes only

    df[[ "start_enh", "end_enh", "seg_index", "core_remodeling"]]\
     = df[[ "start_enh", "end_enh", "seg_index", "core_remodeling"]].astype(float).astype(int)

    df["mrca"] = df["mrca"].astype(float).round(3) # round mrca value

    df["enh_len"]= df.end_enh - df.start_enh
    df = df.loc[df.enh_len <10000] # remove all enhancer longer than 10kb

    if "chr_syn" in list(df):

        df["syn_id"] = df["chr_syn"] + ":" + df["start_syn"].map(str) + "-" + df["end_syn"].map(str)
        df["syn_len"]= df.end_syn - df.start_syn
        df["seg_rep"] = df["syn_len"].divide(df["enh_len"]) # percent of the enhancer occupied by each syntenic block

        #####FILTER out short syntenic blocks#####

        df = df.loc[df["syn_len"]>=6]

        #####ARCHITECTURE#####

        df["code"] = "complex_core"
        df.loc[(df["core"] ==0),"code"] = "derived"
        df.loc[(df["core_remodeling"] ==0), "code"] = "simple"

        df["seg_index"] =  (df["seg_index"] +1) # pseudocount for measuring break density

    df["arch"] = 0
    df.loc[(df["core_remodeling"] ==0),"arch"] = "simple"
    df.loc[(df["core_remodeling"] ==1), "arch"] = "complexenh"

    df = pd.merge(df, syn_gen_bkgd,\
     how = "left", on = "mrca" ).sort_values(by="mrca_2")# add species annotation
    df["mrca_2"] = df["mrca_2"].round(3) # round mrca value

    df = df.drop_duplicates()

    # get the max breaks and oldest mrca

    breakDF = df.groupby([ "chr_enh","start_enh", "end_enh","enh_id", "fourth_col", "id", \
         "core_remodeling", "arch"])["seg_index", "mrca", "enh_len"].max().reset_index()

    breakDF["mrca"] = breakDF["mrca"].round(3) # round mrca value

    breakDF = pd.merge(breakDF, syn_gen_bkgd,\
    how = "left", on = "mrca").sort_values(by="mrca_2")# add species annotation

    breakDF["seg_den"] = (breakDF["seg_index"]).divide(breakDF["enh_len"])

    breakDF["datatype"] = datatype

    breakDF = breakDF.drop_duplicates()

    return df, breakDF

#%%

path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/shuffle/breaks/"

fs = glob.glob("%s*_expressed_parallel_breaks.bed"% path)
fs
#%%
done = []
for f in fs:
    if f not in done:
        fid = (f.split("/")[-1]).split(".")[0]
        print(fid, f)

        breaksdf = pd.read_csv(f, sep='\t')

        num_cols = len(list(breaksdf))

        breaksdf.columns = ["chr_enh", "start_enh", "end_enh", "fourth_col", "enh_id",
                "seg_index", "core_remodeling", "arch", "mrca"]

        breaksdf = breaksdf.loc[breaksdf.mrca != "max_age"] # these are header columns that need to get removed.



        breaksdf["id"] = fid # add sample_id to column

        formatted_df, formatted_df_breaks = postformat_breaksdf(breaksdf, fid, syn_gen_bkgd)

        # save all this information.

        #out_full = "%s%s_enh_age_arch_full_matrix.tsv" % (path, sample_id)

        headerf = "%s/summary_matrix_header.txt" % (path)

        out_summarized_bed = "%s/%s_enh_age_arch_summary_matrix.bed" % (path, fid)
        formatted_df_breaks.to_csv(out_summarized_bed, sep = '\t', index = False, header = False)

        done.append(f)

#%%
formatted_df_breaks.head()
