import os, sys
import pandas as pd
from scipy import stats
import statsmodels
import statsmodels.api as sm
import subprocess

#%%

PATH = "/dors/capra_lab/projects/enhancer_ages/landscape/results/cCRE_x_tfbs_encode3/K562/data/"

CL = "ELS_combined_K562"
BUILD = "hg38"
MIN_GENES = 5 # min genes in group

# possible pairwise comparisons to run GO annotation enrichment on.
ANALYSIS_list = ["der_v_core", "der_v_bkgd", "simple_v_core", "core_v_bkgd", "simple_v_bkgd"]

ANALYSIS = ANALYSIS_list[0]

#%% Functions

# load GO dataframe
def load_GO():


    GOPATH = "/dors/capra_lab/data/gene_ontology/"
    GOHU = f"{GOPATH}goa_human.tsv"
    GOSLIM = f"{GOPATH}PANTHERGOslim_gene_annot.tsv"

    go_ = pd.read_csv(GOHU, sep = '\t')

    #Reduce the go annotation dataframe to key elements
    go = go_[["DB_Object_Symbol", "name", "namespace",
    "Taxon", "GO_ID", "Aspect"]].drop_duplicates()

    return go_, go

# load the cell line TF enrichment file_name
def get_f_dict(cl):

    cl_name = cl.split("ELS_combined_")[1]

    F_DICT = {
    "core_v_bkgd" : f"{cl_name}_core_v_bkgdOR_per_MRCA.tsv",
    "der_v_bkgd":f"{cl_name}_der_v_bkgdOR_per_MRCA.tsv",
    "der_v_core":f"{cl_name}_der_v_coreOR_per_MRCA.tsv",
    "simple_v_bkgd": f"{cl_name}_simple_v_bkgdOR_per_MRCA.tsv",
    "simple_v_core":f"{cl_name}_simple_v_coreOR_per_MRCA.tsv",
    "simple_v_der":f"{cl_name}_simple_v_derOR_per_MRCA.tsv"
    }

    return F_DICT

# get the cell line TF enrichment dataframe
def get_df(path, fdict_key):
    F = os.path.join(path, fdict_key) # join the path
    df = pd.read_csv(F, sep = '\t')
    return df

# perform an FDR correction
def fdr_correction(collection_dict, alpha):

    df = pd.concat(collection_dict.values()) # concat the dictionary into a df

    pvals = df["P"] # get a vector of pvals from df

    # run the multitest
    df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=alpha)

    return df

# get the GO ids linked to TF genes
# filter GO ids to:
#(1) TF, (2) specific process, and (3) GO id gene group >=  min_genes

def get_GO_IDs(tf_ids, df, go, process, min_genes_in_group):


    # (1) get the go ids linked to enriched TFs and (2) process
    tf_go = go.loc[(go.DB_Object_Symbol.isin(tf_ids)) &(go.Aspect == process)].drop_duplicates()

    # count how many genes are in GO annotation related to TF
    names = tf_go.groupby(["Aspect", "name", "GO_ID"])["Taxon"].count().reset_index()

    # (3) filter to GO groups with more than 3 TF genes in the group
    names = names.loc[names.Taxon > min_genes_in_group]

    return names, tf_go

# get the TF genes that are enriched in comp1
def get_tf_ids(comp1, df):

    # get only sig values
    sig = df.loc[df.reject_null == True].drop_duplicates()
    # make a table of log2 values per age
    summed_log2 = pd.pivot(sig, index = "tf", columns = "mrca_2", values = "log2").sum(axis = 1)

    # drop the na, and reset the index
    summed_log2 = summed_log2.dropna().reset_index()

    # log2 > 0 = enriched in comp1
    # GT = "Greater Than"
    GTzero = summed_log2.loc[summed_log2[0]>0, "tf"].to_list()

    # log2 < 0 = enriched in comp2
    # LT = "Lesser Than"
    LTzero = summed_log2.loc[summed_log2[0]<0, "tf"].to_list()

    print(comp1, GTzero, len(GTzero))

    return GTzero, LTzero

def get_OR_go(comp1, comp2, df, process, min_genes_in_group):

    # make a name for this analysis
    comparison_name = f"{comp1}_v_{comp2}_{process}"

    # get the list of tfs for the comparison
    comp1_tf, comp2_tf = get_tf_ids(comp1, df)

    # get TF genes, GO annotations associated with comp1
    comp1_go_ids, comp1df, = get_GO_IDs( comp1_tf, df, go, process, min_genes_in_group)

    # get TF genes, GO annotations associated with comp2
    comp2_go_ids, comp2df, = get_GO_IDs(comp2_tf, df, go, process, min_genes_in_group)

    # combine the genes enriched in comp1 and comp2
    test_ids = pd.concat([comp1_go_ids, comp2_go_ids]) # concat arch1 and arch2 ids

    # only run tests if there are annotations to run in the first comparison
    if len(comp1_go_ids)>0:

        collection_dict = {} # collect the OR results per GO ID

        for goid in test_ids.GO_ID.unique():

            test = comp1df.loc[comp1df.GO_ID == goid].drop_duplicates() # get df overlapping go term
            test_genes = list(test.DB_Object_Symbol.unique()) # get the genes that overlap the term in test.
            test_obs = test.shape[0] # count the number of genes overlapping go term in test
            test_all = len(comp1df.DB_Object_Symbol.unique()) - test_obs # count total number of genes # w/o overlap?


            bkgd_test = comp2df.loc[~comp2df.DB_Object_Symbol.isin(test_genes)] # subtract genes that are in test set, don't want to test them twice.
            bkgd_test = bkgd_test[bkgd_test.GO_ID == goid].drop_duplicates() # count number of genes that overlap go term

            bkgd_test_genes = list(bkgd_test.DB_Object_Symbol.unique()) # get the genes that overlap the term in test.
            bkgd_overlap = bkgd_test.shape[0]
            bkgd_all = len(comp2df.DB_Object_Symbol.unique()) - bkgd_overlap# count total number of genes # exlude test set?

            obs = [[test_obs,test_all], [bkgd_overlap,bkgd_all]]
            print(obs)
            OR, P = stats.fisher_exact(obs)
            #table = sm.stats.Table2x2(obs) # get confidence interval
            #odds_ci = table.oddsratio_confint()


            newdf = pd.DataFrame({"comparison_name":comparison_name,
                                  "a":obs[0][0], "b":obs[0][1],
                                  "c":obs[1][0], "d":obs[1][1],
                                  "OR":[OR], "P":[P],
                                  #"ci_lower" :[odds_ci[0]],
                                  #"ci_upper" :[odds_ci[1]],
                                  "GO_ID": goid
                                })

            collection_dict[goid] = newdf

        if len(collection_dict.keys())>0:
            alpha = 0.1
            resultsdf = fdr_correction(collection_dict, alpha)
            return resultsdf
        else:
            print("no GO annotation enrichment for", comparison_name, process)
    else:
        print("no GO annotations for", comparison_name, process)

def get_annot(resultsdf, go):

    sig_ids = resultsdf.loc[resultsdf.reject_null == True, "GO_ID"] # list of significant ids

    annots_der = go.loc[go.GO_ID.isin(sig_ids), ["GO_ID","name"]].drop_duplicates()

    annot = pd.merge(resultsdf, annots_der)

    return annot

def get_comps(analysis):

    comp1 = analysis.split("_")[0] #"der_v_bkgd"
    comp2 = analysis.split("_")[2]
    print(comp1, comp2)

    return comp1, comp2

def run_GO_enrichment(comp1, comp2, df, go, min_genes):
    result_dict = {}
    for process in go.Aspect.unique():
        print("testing process", process)
        result_df = get_OR_go(comp1, comp2, df, process, MIN_GENES)

        if result_df is not None:
            annot_results = get_annot(result_df, go)
            result_dict[process] = annot_results
            results = pd.concat(result_dict.values())
            results.sort_values(by = "OR")

            return results
        else:
            return None


    # concatenate results


#%% set some constants

F_DICT =  get_f_dict(CL)

#%%
## Load GO annotation file, TF enrichment in CL enhancers

go_, go = load_GO()
df = get_df(PATH, F_DICT[ANALYSIS])

#%% test w/ min number of genes in group

comp1, comp2 = get_comps(ANALYSIS)
MIN_GENES = 10
results = run_GO_enrichment(comp1, comp2, df, go, MIN_GENES)

#%%
results
for n in results["name"].unique():
    print(n)

#%%
# RUN AGE-SPECIFIC ANALYSIS W/ FOR LOOP HERE
mrca_results = {}
MIN_GENES = 5
for mrca_2 in df.mrca_2.unique():
    test = df.loc[df.mrca_2 == mrca_2]
    results = run_GO_enrichment(comp1, comp2, test, go, MIN_GENES)
    if results is not None:
        results["mrca_2"] = mrca_2
        mrca_results[mrca_2] = results
mrca_resultsdf = pd.concat(mrca_results.values())

for n in mrca_resultsdf.name.unique():
    print(n)
#%%
mrca_resultsdf
syn
