import glob
import pandas as pd
import os


PATH = "/dors/capra_lab/data/fantom/fantom5/fantom5_phase1-2_enhancers"
FILE = "human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz"
NAMEFILE = "Human.sample_name2library_id.txt"

FULL = os.path.join(PATH, FILE)
NAMEF = os.path.join(PATH, NAMEFILE)

OUTPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/ENCODE3"

ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/liftOver_hg19"
ENCODEFILES = glob.glob("%s/*.bed" % ENCODEPATH)

ENCODEFANTOM_MAPPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data"
ENCODEFANTOM_MAPFILE = "encode_cl_fantom_cl.tsv"
ENCODEFANTOM_MAPF = os.path.join(ENCODEFANTOM_MAPPATH, ENCODEFANTOM_MAPFILE)

#%%


def get_encode_cell_lines(encodefiles):
    cell_lines = []

    for f in encodefiles:
        if "trimmed" not in f:
            cell_line = (f.split("/")[-1]).split(".bed")[0]
            cell_lines.append(cell_line)
    return cell_lines


def get_fantom_df(namef):

    cols = ["info", "sample_id"]
    df = pd.read_csv(namef, sep ='\t')
    df.columns = cols

    return df

def get_mapfile(mapfile):
    cols = ["ENCODE_cl", "FANTOM_cl", "FANTOM_facet"]
    df = pd.read_csv(mapfile, sep = '\t', header = None, usecols=[0,1,2])
    df.columns = cols

    return df


def get_Ids_sets(ids, threshold, df):
    set_list = []

    for i in ids:
        if i != "Id":
            test_set = set(df.loc[df[i] >=threshold, "Id"]) # get set of enh_id with reads over threshold
            set_list.append(test_set)

    return set_list


#%%

encode_cl = get_encode_cell_lines(ENCODEFILES)
encode_df = pd.DataFrame({"cell_line":encode_cl}) # encode cl dataframe
map_df = get_mapfile(ENCODEFANTOM_MAPF) # map file
fandf = get_fantom_df(NAMEF) # fantom df


#%%

#%%

sample_dict = {}
fandf["FANTOM_cl"] = "" # make a column to merge mapfile with FANTOM

#%% H1 ESC is tricky, so I manually coded this.

fandf.loc[fandf["info"].str.contains('H1 embryonic stem cells'), "cl"] = "H1-esc"
IDs = fandf.loc[fandf["info"].str.contains('H1 embryonic stem cells'), "sample_id"].to_list()
IDs.append("Id")
sample_dict["H1-esc"] = IDs # ugh. stupid str searches


#%%


for cl in map_df.ENCODE_cl.unique():
    fantom_cl = map_df.loc[map_df["ENCODE_cl"] == cl, "FANTOM_cl"].iloc[0]

    if fantom_cl != "FANTOM_cl" and fantom_cl != "X" and str(fantom_cl) != "nan":
        print(fantom_cl)
        fandf.loc[fandf["info"].str.contains(fantom_cl), "FANTOM_cl"] = cl
        IDs = fandf.loc[fandf["FANTOM_cl"] == cl, "sample_id"].to_list()

        if len(IDs) >0:
            IDs.append("Id")

            sample_dict[cl] = IDs

sample_dict.keys()
#%%

THRESHOLD = 0.2 # any erna with >=5 reads will be included.

for cl, ids in sample_dict.items():
    print(cl)


    #load dataframe

    df = pd.read_csv(FULL, sep = '\t', usecols = ids)

    df = df.loc[~df["Id"].str.contains("chrX")] # remove chrX


    # make the bed file

    df["chr"] = df["Id"].apply(lambda x: x.split(":")[0])
    df["start"] = df["Id"].apply(lambda x: (x.split(":")[1]).split("-")[0])
    df["end"] = df["Id"].apply(lambda x: (x.split(":")[1]).split("-")[1])


    # get df with more than 5 reads in each replicate, n = 728

    set_list = get_Ids_sets(ids, THRESHOLD, df)
    u = set.intersection(*set_list) # intersect all the set_list ids
    morethan5_df = df.loc[df["Id"].isin(u)]


    # rearrange the columns in bed format
    morethan5_df = morethan5_df[[
    "chr", "start", "end",
    ]]

    morethan5_df.info()


    #write the file
    OUTFILE = "%s_FANTOM5_tpm_hg19.bed" % cl
    FULL_OUT = os.path.join(OUTPATH, OUTFILE)
    morethan5_df.to_csv(FULL_OUT, sep = '\t', header = False, index = False)
