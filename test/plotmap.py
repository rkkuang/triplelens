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

Imgsize = 128
mags0 = np.log10(mags0.reshape(Imgsize, Imgsize))

cx, cy = readFile(datapath+"caustics.dat", 0, 1, expected_elem_each_row=2)

fig, ax = plt.subplots(figsize=(8,8))
# fig = plt.subplots(figsize=(14,6), dpi=100)
gs = gridspec.GridSpec(1,1)
main = plt.subplot(gs[0])
main.plot(cx, cy, '-', color='yellow', markersize=1)
# ax.yaxis.set_minor_locator(AutoMinorLocator(4))
# ax.xaxis.set_minor_locator(AutoMinorLocator(4))
main.tick_params(axis='both', labelsize = legend_tick_size, direction="in")
main.set_xlabel(r"$x/ \theta_E$", fontdict = font)
main.set_ylabel(r"$y/ \theta_E$", fontdict = font)

mapfig = main.imshow(mags0.T, extent=(-0.2,1,-0.6,0.6), origin='lower',cmap="Reds")
cbar = add_colorbar(mapfig)
cbar.set_label(r"log($\mu$)",fontdict=font)
cbar.ax.tick_params(labelsize=legend_tick_size)
plt.show()
