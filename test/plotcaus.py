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
ax.scatter(cx, cy,marker = '.', color='red', s=1) # caustics

x, y = readFile(datapath+"critical_curves.dat", 0, 1, expected_elem_each_row=2)
ax.scatter(x, y, marker = '.', color='g', s=1) # critical curves


lensx, lensy = readFile(datapath+"lens_system.dat", 1, 2, expected_elem_each_row=3)
lensm, _ = readFile(datapath+"lens_system.dat", 0, 2, expected_elem_each_row=3)

for i in range(len(lensm)-1):
    plt.plot(lensx[1+i], lensy[1+i], '+', color='k')
    # plt.text(lensx[1+i], lensy[1+i],"lens{}".format(i+1))

print("plot pureImgPoints.dat")

xpureImg, ypureImg = readFile(datapath+"pureImgPoints.dat", 0, 1, expected_elem_each_row=2)
# ax.plot(xpureImg, ypureImg, '.', color='cyan', markersize=1)
ax.plot(xpureImg, ypureImg, '.', color='b', markersize=1)
print("len imgs: ", len(x))

xfalseimg, yfalseimg = readFile(datapath+"pureImgPoints_falseimg.dat", 0, 1, expected_elem_each_row=2)
ax.plot(xfalseimg, yfalseimg, '.', color='gray', markersize=1)

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
phi=np.linspace(0.0, 2*np.pi, nphi)
sx=xs+rs*np.cos(phi)
sy=ys+rs*np.sin(phi)
ax.plot(sx, sy, 'k')
plt.axis('equal') 

ax.annotate('(${:.1f}$, ${:.1f}$)'.format(xsCenter,ysCenter), xy=(0.3, 0.9), xycoords='axes fraction', fontsize=17,
                horizontalalignment='right', verticalalignment='bottom')
plt.show()
