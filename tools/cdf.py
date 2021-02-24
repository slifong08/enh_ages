import numpy as np
import pandas as pd

# input dataframe and column w/ values for calculating cdf
def get_cumsum(df, col):

    cdf= np.cumsum(df.col)/1 # calculate cdf
    newdf = pd.DataFrame({"cdf": cdf}) # create df of cdf

    #merge cdf values w/ original dataframe
    testdf = pd.merge(df, newdf, left_index = True, right_index = True)

    return testdf
