import os, sys
import subprocess
import pandas as pd

"""
Example of the file format:


[Term]
id: GO:0000132
name: establishment of mitotic spindle orientation
namespace: biological_process
alt_id: GO:0030607
alt_id: GO:0030609
def: "A cell cycle process that sets the alignment of mitotic spindle relative to other cellular structures." [GOC:ems]
synonym: "establishment of spindle orientation during mitosis" RELATED [GOC:dph, GOC:tb]
synonym: "establishment of spindle orientation involved in mitotic cell cycle" EXACT [GOC:dph, GOC:tb]
synonym: "mitotic spindle orientation" EXACT []
synonym: "orienting of mitotic spindle" EXACT []
is_a: GO:0040001 ! establishment of mitotic spindle localization
is_a: GO:0051294 ! establishment of spindle orientation
intersection_of: GO:0051294 ! establishment of spindle orientation
intersection_of: part_of GO:0000278 ! mitotic cell cycle


"""

PATH = "/dors/capra_lab/data/gene_ontology/"
F = f"{PATH}PANTHERGOslim.obo"


# make a bunch of lists to collect the results
go_ids, names, namespaces = [], [], []

# open the file
with open(F, "r") as file:
    lines = file.readlines() # read the lines
    for i in range(0, len(lines)): # parse through lines
        line = lines[i] # figure out the line index

        if line == "[Term]\n": # if the line is a new term

            # get the go_id fields (1 below), name (2 lines below), and namespace (3 lines below)
            go_id = lines[i + 1].split("id: ")[1].rstrip() # do some str processing to get data you want.
            name = lines[i + 2].split("name: ")[1].rstrip()
            namespace = lines[i + 3].split("namespace: ")[1].rstrip()

            # append fields to list
            go_ids.append(go_id)
            names.append(name)
            namespaces.append(namespace)
            #print(go_id, name, namespace,'\n')


#%% make a dataframe the only way you know how

df = pd.DataFrame({
"GO_ID": go_ids,
"name": names,
"namespace": namespace
})

outf = f"{PATH}PANTHERGOslim.tsv"

df.to_csv(outf, sep = '\t', index = False)
#%%


GOPATH = "/dors/capra_lab/data/gene_ontology/"
GOHU = f"{GOPATH}goa_human.tsv"


go_ = pd.read_csv(GOHU, sep = '\t')

go_.head()


#%% Reduce the go annotation dataframe to key elements


go = go_[["DB_Object_Symbol", "name", "namespace", "Taxon", "GO_ID", "Aspect"]].drop_duplicates()


#%% reduce GO dataframe down to terms that overlap SLIM annotations. Get genes.


merged=pd.merge(go, df, how = "left", on = ["name", "namespace"])
slim = merged.loc[~merged["GO_ID_y"].isna()]
outf = f"{PATH}PANTHERGOslim_gene_annot.tsv"

slim = slim[["DB_Object_Symbol", "name", "namespace", "Aspect", "Taxon", "GO_ID_y"]].drop_duplicates()
slim = slim.rename(columns={'GO_ID_y': 'GO_ID'})
slim.shape
slim.to_csv(outf, sep = '\t', index = False)

slim.groupby("name")["namespace"].count()
gg = go.groupby("name")["namespace"].count().reset_index()
gg.loc[gg.name == "digestion"]
