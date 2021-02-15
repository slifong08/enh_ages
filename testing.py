lol = {"hello": "gdbye", "yes":"no", "p":"u"}
print(lol)
print(len(lol))
#%%
i = list(lol.keys())

for n, key in enumerate(i):
    val =  lol.pop(key)

    print(n, key, val, lol)
#%%
import numpy as np

val = 0
max = 10
for i in np.arange(max):
    val +=1
    if val % 4 == 0:
        print(val)
    elif val == max:
        print("done")
