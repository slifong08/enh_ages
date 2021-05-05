import os, sys
import pandas as pd
from scipy import stats
import statsmodels
import statsmodels.api as sm
import subprocess

#%%

PATH = "/dors/capra_lab/projects/enhancer_ages/landscape/results/cCRE_x_tfbs_encode3/HepG2/data/"


#%% Functions


def get_f_dict():
    F_DICT = {
    "core_v_bkgd" : "HepG2_core_v_bkgdOR_per_MRCA.tsv",
    "der_v_bkgd":"HepG2_der_v_bkgdOR_per_MRCA.tsv",
    "der_v_core":"HepG2_der_v_coreOR_per_MRCA.tsv",
    "simple_v_bkgd": "HepG2_simple_v_bkgdOR_per_MRCA.tsv",
    "simple_v_core":"HepG2_simple_v_coreOR_per_MRCA.tsv",
    "simple_v_der":"HepG2_simple_v_derOR_per_MRCA.tsv"
    }
    return F_DICT

def fdr_correction(collection_dict, alpha):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=alpha)

    return df

def get_GO_IDs(arch, comp1, df, go, process, min_genes_in_group):

    # get the list of tfs
    tf_ids = get_tf_ids(arch, df, comp1)

    # get the go ids linked to tfs
    tf_go = go.loc[(go.DB_Object_Symbol.isin(tf_ids)) &(go.Aspect == process)].drop_duplicates()

    # reduce down the go dataframe to the relevant to the TFs
    names = tf_go.groupby(["Aspect", "name", "GO_ID"])["Taxon"].count().reset_index()
    names = names.loc[names.Taxon > min_genes_in_group]

    return names, tf_go

def get_tf_ids(arch, df, comp1):
    if comp1 == "der":
        comp_dict = {
        "core":list(df.loc[(df.reject_null == True) & (df.OR < 1), "tf"].unique()),
        "der":list(df.loc[(df.reject_null == True) & (df.OR > 1),  "tf"].unique()),
        "simple":list(df.loc[(df.reject_null == True) & (df.OR < 1),  "tf"].unique()),
        "bkgd":list(df["tf"].unique())
        }
    elif comp1 == "core":
        comp_dict = {
        "core":list(df.loc[(df.reject_null == True) & (df.OR >1), "tf"].unique()),
        "der":list(df.loc[(df.reject_null == True) & (df.OR < 1),  "tf"].unique()),
        "simple":list(df.loc[(df.reject_null == True) & (df.OR < 1),  "tf"].unique()),
        "bkgd":list(df["tf"].unique())
        }
    elif comp1 == "simple":
        comp_dict = {
        "core":list(df.loc[(df.reject_null == True) & (df.OR <1), "tf"].unique()),
        "der":list(df.loc[(df.reject_null == True) & (df.OR < 1),  "tf"].unique()),
        "simple":list(df.loc[(df.reject_null == True) & (df.OR > 1),  "tf"].unique()),
        "bkgd":list(df["tf"].unique())
        }
    tf_ids = comp_dict[arch]
    print(tf_ids, arch)

    return tf_ids

def get_OR_go(comp1, comp2, df, process, min_genes_in_group):

    comparison_name = f"{comp1}_v_{comp2}_{process}"

    comp1_go_ids, comp1df, = get_GO_IDs(comp1, comp1, df, go, process, min_genes_in_group)

    comp2_go_ids, comp2df, = get_GO_IDs(comp2, comp1, df, go, process, min_genes_in_group)


    test_ids = pd.concat([comp1_go_ids, comp2_go_ids]) # concat arch1 and arch2 ids

    # only run tests if there are annotations to run

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

            OR, P = stats.fisher_exact(obs)
            table = sm.stats.Table2x2(obs) # get confidence interval
            odds_ci = table.oddsratio_confint()


            newdf = pd.DataFrame({"comparison_name":comparison_name,
                                  "a":obs[0][0], "b":obs[0][1],
                                  "c":obs[1][0], "d":obs[1][1],
                                  "OR":[OR], "P":[P],
                                  "ci_lower" :[odds_ci[0]],
                                  "ci_upper" :[odds_ci[1]],
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


#%% set some constants

F_DICT =  get_f_dict()
BUILD = "hg38"
CL = "ELS_combined_HepG2"

#%%
## Open the GO annotation file

GOPATH = "/dors/capra_lab/data/gene_ontology/"
GOHU = f"{GOPATH}goa_human.tsv"
GOSLIM = f"{GOPATH}PANTHERGOslim_gene_annot.tsv"

go_ = pd.read_csv(GOHU, sep = '\t')


go_.head()


#%% Reduce the go annotation dataframe to key elements


go = go_[["DB_Object_Symbol", "name", "namespace", "Taxon", "GO_ID", "Aspect"]].drop_duplicates()

#%% initialize results dictionary

result_dict = {}

#%% get the dataframe


ANALYSIS = "der_v_core"
"""
ANALYSIS = "der_v_bkgd"
ANALYSIS = "simple_v_core"
ANALYSIS = "core_v_bkgd"
ANALYSIS = "simple_v_bkgd"
"""
F = os.path.join(PATH, F_DICT[ANALYSIS])
df = pd.read_csv(F, sep = '\t')
df.head()

#%% test w/ min number of genes in group

comp1 = ANALYSIS.split("_")[0] #"der_v_bkgd"
comp2 = ANALYSIS.split("_")[2]
min_genes_in_group = 3


for process in go.Aspect.unique():
    result_df = get_OR_go(comp1, comp2, df, process, min_genes_in_group)
    annot_results = get_annot(result_df, go)
    result_dict[process] = annot_results
#%%

results = pd.concat(result_dict.values())
results.sort_values(by = "OR")
