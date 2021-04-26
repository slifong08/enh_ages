
def assign_architecture(df):
    df["code"] = ""
    df.code.loc[(df.core_remodeling == 0)& (df.core == 1)] = "simple"
    df.code.loc[(df.core_remodeling == 1)& (df.core == 1)] = "complex_core"
    df.code.loc[(df.core_remodeling == 1)& (df.core == 0)] = "derived"

    df["arch"] = ""
    df.arch.loc[(df.core_remodeling == 0)] = "simple"
    df.arch.loc[(df.core_remodeling == 1)] = "complexenh"

    return df 
