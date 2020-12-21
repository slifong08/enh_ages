import pandas as pd

from scipy import stats
import statsmodels
import statsmodels.api as sm

collection_dict = {}

def get_2x2(val, df, rand):


    obs = [[a,b], [c,d]]
    OR, P = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()
    newdf = pd.DataFrame({"seg_index":[seg_index], "a":[a], "b":[b], "c":[c], "d":[d],
                         "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                        "ci_upper" :[odds_ci[1]]})

    print(seg_index, obs, OR, P)

    return newdf

def fdr_correction(collection_dict):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05)

    return df
