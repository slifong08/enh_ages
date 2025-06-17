import matplotlib

import matplotlib.colors as mcolors
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap

import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import seaborn as sns

#colors = [ "amber", "dusty purple", "windows blue"]
#PAL = sns.xkcd_palette(colors)
#sns.palplot(PAL)

#%%
yg = '#f7fcb9'
stemg = '#addd8e'
kellyg = '#31a354'

y = mcolors.to_rgba(yg)
lg = mcolors.to_rgba(stemg)
g = mcolors.to_rgba(kellyg)


colors = [g, lg, y]  # kelly (core) -> light (derived) - > yellow (simple)
cmap_name = 'yg_'
CM = LinearSegmentedColormap.from_list(cmap_name, colors, N=3)
print("COLORMAP - CM", CM)


#%%
colors = [ "ecru"]  # simple
y_ = sns.xkcd_palette(colors)
print("y_")
sns.palplot(y_)

colors = [ "light sage"]  # derived
stemg = sns.xkcd_palette(colors)
print("stemg")
sns.palplot(stemg)

colors = [ "shamrock green"]  # core
kellyg = sns.xkcd_palette(colors)
print("kellyg")
sns.palplot(kellyg)

colors = ["ecru","shamrock green", "light sage"]
yg = sns.xkcd_palette(colors)
print("yg")
sns.palplot(yg)

colors = ["shamrock green", "light sage"]
cd = sns.xkcd_palette(colors)
print("cd")
sns.palplot(cd)

colors = ["ecru",  "putty","shamrock green", "cool grey","light sage", "lichen"]
yg_shuf = sns.xkcd_palette(colors)
print("yg_shuf")
sns.palplot(yg_shuf)


colors = ["shamrock green", "cool grey","light sage", "lichen"]
gg_shuf = sns.xkcd_palette(colors)
print("gg_shuf")
sns.palplot(gg_shuf)

### FONTS ###

matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.size'] = 18