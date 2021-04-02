import numpy as np

def make_chr_list():
    n = list(np.arange(1, 23))
    #n.append("X")

    chr_list = []
    for num in n:
        chrn = "chr" + str(num)
        chr_list.append(chrn)

    return chr_list

#%%
chr_list = []
n = list(np.arange(1, 23))

(chr_list.append(("chr" + str(num) for num in n)))

a = (list(chr_list))

for num in n:
    if num in a:
        print(num)
a
num
