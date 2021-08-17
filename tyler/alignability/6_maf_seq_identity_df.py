import glob
import os, sys
import pandas as pd
import subprocess

PATH = "/dors/capra_lab/users/fongsl/tyler/data/alignment/"

FS = glob.glob(f"{PATH}chr*_seq_identity.tsv")

data_dict = {}
for n, F in enumerate(FS):
    use_cols = [
    "#chr",
    "start",
    "end",
    "hg38_Seqlen",
    "rheMac8_Seqlen",
    "score",
    "per_id_hg38",
    "per_id_rheMac8",
    #"hg38_seq",
    #"rheMac8_seq"
    ]
    df = pd.read_csv(F, sep = '\t', usecols = use_cols)

    data_dict[n] = df

df = pd.concat(data_dict.values())

df.head()


#%% describe the dataframe


df[["hg38_Seqlen",
"rheMac8_Seqlen",
"score",
"per_id_hg38",
"per_id_rheMac8"
]].describe()
"""
	hg38_Seqlen	rheMac8_Seqlen	score	per_id_hg38	per_id_rheMac8
count	145202.000000	145202.000000	145202.000000	145202.000000	143773.000000
mean	466.575743	459.266346	431.618724	0.921215	0.940185
std	500.260567	493.822968	464.647624	0.101226	0.022513
min	39.000000	0.000000	0.000000	0.000000	0.621818
25%	185.000000	182.000000	171.000000	0.916505	0.927273
50%	291.000000	287.000000	270.000000	0.934132	0.941476
75%	548.000000	542.000000	509.000000	0.949749	0.954545
max	17078.000000	16459.000000	15487.000000	1.000000	1.030303
"""

#%% evaluate the regions where per_id_rheMac8 > 1 (higher sequence identity score than length?)
# one part of this is due to half-open bedtools coordinates and the way maf size


df.loc[df["per_id_rheMac8"] >1]

"""
	#chr	start	end	hg38_Seqlen	rheMac8_Seqlen	score	per_id_hg38	per_id_rheMac8
14	chr5	262051	262107	4458	865	870.0	0.194587	1.005780
3582	chr17	50802069	50802233	510	98	99.0	0.193738	1.010204
5250	chr7	101798313	101798457	184	33	34.0	0.183784	1.030303
3025	chr16	66874410	66874442	212	58	59.0	0.276995	1.017241
3694	chr10	72259488	72259717	551	48	49.0	0.088768	1.020833
2640	chr22	46862674	46862730	296	51	52.0	0.175084	1.019608
12208	chr2	241863108	241863171	496	83	85.0	0.170341	1.024096
"""
#%%

df.per_id_hg38.value_counts(bins =10, normalize = True)

"""
(0.9, 1.0]       0.884141
(0.8, 0.9]       0.099406
(-0.002, 0.1]    0.010069
(0.7, 0.8]       0.002486
(0.6, 0.7]       0.001302
(0.5, 0.6]       0.000778
(0.3, 0.4]       0.000613
(0.4, 0.5]       0.000523
(0.1, 0.2]       0.000365
(0.2, 0.3]       0.000317
Name: per_id_hg38, dtype: float64
"""

df.per_id_rheMac8.value_counts(bins =10, normalize = True)

"""
(0.908, 0.949]    0.560467
(0.949, 0.989]    0.346242
(0.867, 0.908]    0.069035
(0.989, 1.03]     0.009869
(0.826, 0.867]    0.003884
(0.785, 0.826]    0.000468
(0.744, 0.785]    0.000096
(0.704, 0.744]    0.000062
(0.663, 0.704]    0.000021
(0.62, 0.663]     0.000014
Name: per_id_rheMac8, dtype: float64
"""
