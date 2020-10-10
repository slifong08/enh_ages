def custom_round(x, base=5):
    return int(base * round(float(x)/base))

#df = pd.Series([11,16,21]).apply(lambda x: custom_round(x, base=5))
