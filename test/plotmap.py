'''
Figure 8. in original draft

'''

from utils import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl  
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, FormatStrFormatter
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import gridspec
import matplotlib.ticker as ticker

import matplotlib.ticker as mticker
f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
import os

mpl.rc('font',family='Times New Roman')
fontsz = 17

font = {'family': 'Times New Roman',
    'color':  'k',
    'weight': 'normal',
    'size': fontsz,
    }
legend_tick_size = fontsz
font2 = {'family': 'Times New Roman',
    'weight': 'normal',
    'size': fontsz,
    }


cmap1 = 'seismic'
datapath = "./data/"

lens_params = read_lens_system_triple(datapath+"lens_system_triple.dat")
xsCenter, ysCenter = lens_params["xsCenter"], lens_params["ysCenter"]
Xs0, Ys0 , mags0 = read_cppmap(datapath+"magmap0.05.dat")

# mapraysht = np.load(datapath+"rayshtmap_5e-2_128_raynum2e4.dat")
mapraysht = np.load(datapath+"rayshtmap_5e-2_128_raynum5e4.dat")

filename = "rayshtmap_data"
dirpath = "/Users/anything/THU/astro/softwares/gravlens/triplelens/data/response_figure/data/"
if not os.path.exists(dirpath+filename):
    print("saving raysht mag map {}".format(dirpath+filename))
    np.savez(dirpath+filename, Xs0=Xs0, Ys0=Ys0, mapraysht = mapraysht)
else:
    print("{} already exists, do you sure to overwrite?".format(dirpath+filename))


# mapraysht = np.load("/Users/anything/THU/astro/softwares/aeroastro/gravlen/critical_and_caustics/pys/totahiti191210/muRayshoot_raynum1.0e+04.npy")
# mapraysht = mapraysht.T

Imgsize = 128
mags0 = mags0.reshape(Imgsize, Imgsize)

resmap = (mags0 - mapraysht)/mapraysht

mags0 = np.log10(mags0)
mapraysht = np.log10(mapraysht)



# fig, ax = plt.subplots(figsize=(8,8))

fig = plt.subplots(figsize=(14,6), dpi=100)
gs = gridspec.GridSpec(1,2)
plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.9, left = 0.1, hspace = 0, wspace = 0.2)
# plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.9, left = 0.1, hspace = 0, wspace = 0.35)


main = plt.subplot(gs[0])
# ax.yaxis.set_minor_locator(AutoMinorLocator(4))
# ax.xaxis.set_minor_locator(AutoMinorLocator(4))
main.tick_params(axis='both', labelsize = legend_tick_size, direction="in")
main.set_xlabel(r"$x/ \theta_E$", fontdict = font)
main.set_ylabel(r"$y/ \theta_E$", fontdict = font)

# mapfig = main.imshow(mags0.T, extent=(-0.2,1,-0.6,0.6), origin='lower',cmap=cmap1)
mapfig = main.imshow(mags0.T, extent=(-0.2,1,-0.6,0.6), origin='lower',cmap="Reds")

# mapfig = plt.imshow(mapraysht.T, extent=(-0.2,1,-0.6,0.6), origin='lower',cmap=cmap1)

cbar = add_colorbar(mapfig)
cbar.set_label(r"log($\mu$)",fontdict=font)
cbar.ax.tick_params(labelsize=legend_tick_size)

secnd = plt.subplot(gs[1])
# ax.yaxis.set_minor_locator(AutoMinorLocator(4))
# ax.xaxis.set_minor_locator(AutoMinorLocator(4))
secnd.tick_params(axis='both', labelsize = legend_tick_size, direction="in")
secnd.set_xlabel(r"$x/ \theta_E$", fontdict = font)
# ax.set_ylabel(r"$y(\theta_E)$", fontdict = font)
# mapfig = secnd.imshow(mapraysht.T, extent=(-0.2,1,-0.6,0.6), origin='lower',cmap=cmap1)

mapfig = secnd.imshow(resmap.T, extent=(-0.2,1,-0.6,0.6), origin='lower',cmap=cmap1, vmin = -7e-4, vmax = 7e-4)

cbar = add_colorbar(mapfig)#format=ticker.FuncFormatter(fmt),format=OOMFormatter(-2, mathText=True)
cbar.set_label(r"Residual",fontdict=font)
cbar.ax.tick_params(labelsize=legend_tick_size)
cbar.ax.yaxis.get_offset_text().set_fontsize(legend_tick_size)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()





# if input("save or not? [y/n]") == "y":

# inset_axes = inset_axes(ax,
#                     width="32%", # width = 30% of parent_bbox
#                     height="32%", # height : 1 inch
#                     loc=1)
# inset_axes.plot(cx, cy, '-', color='red', markersize=1)
# inset_axes.plot(sx, sy, 'k')
# inset_axes.tick_params(axis='both', labelsize = 12, direction="in")
# inset_axes.set_ylim(-0.1,0.1)
# inset_axes.set_xlim(-0.1,0.1)


# plt.savefig("./data/magmap.png".format(xsCenter,ysCenter), dpi=300)
plt.show()
