import os, sys
import pandas as pd
import numpy as np
from scipy import stats

FID = "ENCFF325VXG"
PATH = "/dors/capra_lab/projects/enhancer_ages/Tippens/data/"
F = f"{FID}/ages/syn_breaks_{FID}_ages.bed"


AGE_FILE = os.path.join(PATH, F)
AGE_FILE
ACT = "ENCFF591EIJ_tippens_activity.csv"
ACT_FILE = os.path.join(PATH, ACT)

ID_ = f"{FID}.bed"
ID_FILE = os.path.join(PATH, ID_)


#%% open the activity file
def formatACT_FILE(act_file):

    df = pd.read_csv(act_file, sep = ",")

    df.loc[df["ID"].str.contains("EK"), "extended"] = 1
    df.loc[(df["ID"].str.contains("A")) & (df["ID"].str.contains("TRE")), "whole"] = 1
    df.loc[(df["ID"].str.contains("B")) & (df["ID"].str.contains("TRE")), "whole"] = 5
    df.loc[(df["ID"].str.contains("c")) & (df["ID"].str.contains("TRE")), "whole"] = 3

    df["ID2"] = df.loc[df["ID"].str.contains("EK")]["ID"].apply(lambda x: "K" + x.split("EK")[1])

    df.loc[df["ID2"].isna(), "ID2"] = df["ID"]
    return df

def formatAGE_FILE(age_file):

    cols = ["chr_s", "start_s", "end_s", "enh_id",
     "chr_e", "start_e", "end_e", "seg_index", "core_remodeling", "core", "mrca"]
    #open the enhancer age file
    df = pd.read_csv(age_file, sep = "\t", header = None)
    df.columns = cols

    # calculate syntenic length
    df["syn_len"] = df.end_s - df.start_s

    # syn_id
    df["syn_id"] = df.chr_s + ":" + df.start_s.map(str) + "-" + df.end_s.map(str)

    return df

def formatID_FILE(id_file):

    cols = ["chr_e", "start_e", "end_e", "ID"]
    df = pd.read_csv(id_file, sep = "\t", header = None)
    df.columns = cols

    return df

#%%

act = formatACT_FILE(ACT_FILE)
age = formatAGE_FILE(AGE_FILE)
id = formatID_FILE(ID_FILE)
age.head()
age_id = pd.merge(id, age, how = "left")
complex_age = age_id.loc[age_id.core_remodeling == 1]
complex_id = complex_age["ID"].unique() # get the ids for complex enhancers
print(len(complex_id)) # 94 regions are complex


complex_act = act.loc[act.ID.isin(complex_id)]

complex = pd.merge(complex_age, complex_act)
#%%

todrop = ['chr_e', 'start_e', 'end_e', 'forward_q', 'reverse_q', "forward_logFC", "reverse_logFC",
'chr_s', 'start_s', 'end_s', 'extended', "whole"]
complex = complex.drop(todrop, axis = 1).drop_duplicates()
#%%
complex.loc[complex.ID.str.contains("B")]


#%% CATCH ALL THE ENHANCERS where extending the boundaries

not_same_change_act = []
not_same = []
same = []
for ID_ in complex.ID2.unique():
    re = complex.loc[complex.ID.str.contains(ID_)]
    val = 0
    call = None
    for uid in re.ID.unique():
        test_seg = re.loc[re.ID == uid, "seg_index"].max()
        test_call = re.loc[re.ID == uid, "call"].iloc[0]
        #print(uid, test)
        if val ==0:
            val = test_seg
            call = test_call
        elif test_seg == val:
            print("same number of segments tested", ID_)
            same.append(ID_)
        elif test_seg != val:
            print("NOT number of segments tested", ID_)
            if test_call != call:
                not_same_change_act.append(ID_)
            else:
                not_same.append(ID_)
#%%
len(same) # 60 elements that have the same length
same_change_act =[] # cases where there is no addition of new sequence age, but extension of same sequence ages that change activity.
val = 0
for ID_ in same:
    re = complex.loc[complex.ID.str.contains(ID_)]
    if len(re.ID.unique())>1:
        val_ = 0
        call = None
        for uid in re.ID.unique():

            test = re.loc[re.ID == uid, "core"].to_list()
            testAct = re.loc[re.ID == uid, ["ID", "logFC", "call", 'Size', "syn_id"]].drop_duplicates()
            test_call = re.loc[re.ID == uid, "call"].iloc[0]

            if val_ == 0:
                call = test_call
            elif call != test_call:
                print(ID_, test, re)
                same_change_act.append(ID_)
            elif call == test_call:
                continue

            val_ +=1

    val +=1
    print("\n\n")
#%%
len(not_same) # 53 enhancers
val = 0
for ID_ in not_same_change_act:
    re = complex.loc[complex.ID.str.contains(ID_)]
    for uid in re.ID.unique():
        test = re.loc[re.ID == uid, "core"].to_list()
        testAct = re.loc[re.ID == uid, ["ID", "logFC", "call", 'Size', "syn_id"]].drop_duplicates()
        print(test, testAct)
    val +=1
    print("\n\n")


  #KUUAE0005, KUUAE0024, KUUAT0049 enhancer core to enhancer
  #EKNGAT0016, KUUAE0010, KUUAE0051 inactive core  to enhancer
#%%
inact_2_act = ["KNGAT0016", "KUUAE0010", "KUUAE0051"]
act_2_act = ["KUUAE0005", "KUUAE0024", "KUUAT0049"]

outf = f"{PATH}_inact_2_act.tsv"
i2a = complex.loc[complex.ID2.isin(not_same_change_act)]
i2a.to_csv(outf, sep = '\t', index = False)

outf = f"{PATH}_act_2_act.tsv"
a2a = complex.loc[complex.ID2.isin(act_2_act)]
a2a.to_csv(outf, sep = '\t', index = False)

outf = f"{PATH}_act_change_act.tsv"
a2a = complex.loc[complex.ID2.isin(same_change_act)]
a2a.to_csv(outf, sep = '\t', index = False)
#%%
for ID_ in inact_2_act:
    re = complex.loc[complex.ID.str.contains(ID_)]
    print(re, "\n\n")

#%%
for ID_ in act_2_act:
    re = complex.loc[complex.ID.str.contains(ID_)]
    print(re, "\n\n")

#%%
for ID_ in complex.ID2.unique():
    if "TRE" not in ID_:
        test = complex
        print(act.loc[act.ID.str.contains(ID_), ["ID", "logFC", "call", "Size"]])




#%%KUUAT0042
"""
Both active Full>> core
EKSUUT0025
EKUUAE0147
EKUUAE0104
EKNGAT0067
EKUUAE0034
EKSUUT0003
EKUUAE0047

Full active, truncated not active
EKUUAE0046
EKUUAT0030
EKUUUT0012
EKUUAE0010
EKUUAE0051
EKSUUT0007* added derived?
EKSSAT0026
EKNGAE0056
EKUUAT0005
EKSUAT0013
KUUAE0066
EKSUAT0156
EKNGAE0132
EKUUAE0063
KSUUT0116
KSUUT0055
EKNGAT0016
EKUUAE0005

Core is active, full is in
KSUAE0003 (core is active, full and truncated are not )

# no pairs?
KUUAT0015
KUUAT0042
"""
#%%
simp_age = age_id.loc[age_id.core_remodeling != 1]
simp_id = simp_age["ID"].unique() # get the ids for simp enhancers
print(len(simp_id)) # 76 regions are simp


simp_act = act.loc[act.ID.isin(simp_id)]

simp = pd.merge(simp_age, simp_act)
#%%

todrop = ['chr_e', 'start_e', 'end_e', 'forward_q', 'reverse_q', "enh_id", "forward_logFC", "reverse_logFC",
"core_remodeling", 'chr_s', 'start_s', 'end_s', 'extended', "whole"]
simp = simp.drop(todrop, axis = 1).drop_duplicates()
#%%

for ID_ in simp.ID2.unique():
    if "TRE" not in ID_:

        print(act.loc[act.ID.str.contains(ID_), ["ID", "logFC", "call", "Size"]])
