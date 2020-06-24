from utils import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl  
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, FormatStrFormatter

from mpl_toolkits.axes_grid.inset_locator import inset_axes


mpl.rc('font',family='Times New Roman')

font = {'family': 'Times New Roman',
    'color':  'k',
    'weight': 'normal',
    'size': 17,
    }
legend_tick_size = 17
font2 = {'family': 'Times New Roman',
    'weight': 'normal',
    'size': 17,
    }

datapath = "./data/"

lens_params = read_lens_system_triple(datapath+"lens_system_triple.dat")
xsCenter, ysCenter = lens_params["xsCenter"], lens_params["ysCenter"]

fig, ax = plt.subplots(figsize=(8,8))


# ax.yaxis.set_minor_locator(AutoMinorLocator(4))
# ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.tick_params(axis='both', labelsize = legend_tick_size, direction="in")

cx, cy = readFile(datapath+"caustics.dat", 0, 1, expected_elem_each_row=2)
ax.plot(cx, cy, '-', color='red', markersize=1)

# inset_axes.plot(x, y, '-', color='red', markersize=1)


x, y = readFile(datapath+"critical_curves.dat", 0, 1, expected_elem_each_row=2)
ax.plot(x, y, '--', color='r', markersize=1)


lensx, lensy = readFile(datapath+"lens_system.dat", 1, 2, expected_elem_each_row=3)
lensm, _ = readFile(datapath+"lens_system.dat", 0, 2, expected_elem_each_row=3)

for i in range(len(lensm)-1):
    plt.plot(lensx[1+i], lensy[1+i], '+', color='k')#markersize=5*lensm[i+1]
    # plt.text(lensx[1+i], lensy[1+i],"lens{}".format(i+1))

print("plot pureImgPoints.dat")

xpureImg, ypureImg = readFile(datapath+"pureImgPoints.dat", 0, 1, expected_elem_each_row=2)
ax.plot(xpureImg, ypureImg, '.', color='cyan', markersize=1)
print("len imgs: ", len(x))

xfalseimg, yfalseimg = readFile(datapath+"pureImgPoints_falseimg.dat", 0, 1, expected_elem_each_row=2)
ax.plot(xfalseimg, yfalseimg, '.', color='b', markersize=1)

f = open(datapath+"lens_system.dat", "r")
#read in the source information
full_line = f.readline()
f.close()
line=full_line.split()

xs=np.float(line[0])
ys=np.float(line[1])
rs=np.float(line[2])
print("xs, ys, rs in py", xs,ys,rs)

nphi=150
# phi=arange(0.0, 2*pi, 2.0*pi/(nphi-1))
phi=np.linspace(0.0, 2*np.pi, nphi)
sx=xs+rs*np.cos(phi)
sy=ys+rs*np.sin(phi)
ax.plot(sx, sy, 'k')
# plt.text(xs,ys,"src")
# ax.set_ylim(-1.25,1.75)
# ax.set_xlim(-1.1,2.1)
plt.axis('equal')
# 
# 
ax.annotate('({}, {})'.format(xsCenter,ysCenter), xy=(0.2, 0.9), xycoords='axes fraction', fontsize=17,
                horizontalalignment='right', verticalalignment='bottom')
# if input("save or not? [y/n]") == "y":

# inset_axes = inset_axes(ax,
#                     width="32%", # width = 30% of parent_bbox
#                     height="32%", # height : 1 inch
#                     loc=1) #1 top right, 
# inset_axes.plot(cx, cy, '-', color='red', markersize=1)
# inset_axes.plot(sx, sy, 'k')
# inset_axes.plot(xpureImg, ypureImg, '.', color='cyan', markersize=1)
# inset_axes.plot(xfalseimg, yfalseimg, '.', color='b', markersize=1)
# inset_axes.plot(cx, cy, '-', color='red', markersize=1)
# inset_axes.tick_params(axis='both', labelsize = 12, direction="in")
# inset_axes.set_xlim(-0.06,0.04)
# inset_axes.set_ylim(-0.05,0.05)


# plt.savefig("./data/topo_{}_{}_rs{}.png".format(xsCenter,ysCenter,lens_params["rs"]), dpi=300)
plt.show()
