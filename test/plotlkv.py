from utils import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl  
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, FormatStrFormatter

from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import gridspec

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


times0, Xs0, Ys0 , mags0 = read_cpplkv(datapath+"ftrilkv_rs_0.dat")
times1, Xs1, Ys1 , mags1 = read_cpplkv(datapath+"ftrilkv_rs_0.15.dat")
times2, Xs2, Ys2 , mags2 = read_cpplkv(datapath+"ftrilkv_rs_5e-2.dat")
times3, Xs3, Ys3 , mags3 = read_cpplkv(datapath+"ftrilkv_rs_1e-2.dat")
# times4, Xs4, Ys4 , mags4 = read_cpplkv(datapath+"ftrilkv_rs_2.2e-4.dat")
mags0 = np.log10( mags0 )
mags1 = np.log10( mags1 )
mags2 = np.log10( mags2 )
mags3 = np.log10( mags3 )
# mags4 = np.log10( mags4 )

fig = plt.subplots(figsize=(13,7), dpi=100)
# gs = gridspec.GridSpec(2,1,height_ratios=[5,1])
# plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.9, left = 0.1, hspace = 0, wspace = 0)
# main = plt.subplot(gs[0])
main = plt.subplot()

main.set_ylabel(r"log($\mu$)", fontdict=font)
main.set_xlabel('HJD - 2450000', fontdict=font)
main.tick_params(axis='both', labelsize = legend_tick_size, direction = "in")#labelcolor=color
# main.yaxis.set_minor_locator(AutoMinorLocator(5))

colors = ["k","b", "g", "r", "cyan"]

main.plot(times0, mags0, color=colors[0],linewidth = 2)
main.plot(times0, mags1, color=colors[1],linewidth = 2)
main.plot(times0, mags2, color=colors[2],linewidth = 2)
main.plot(times0, mags3, color=colors[3],linewidth = 2)
# main.plot(times0, mags4, color=colors[4],linewidth = 2)

# "t0,u0,tE,s2,q2,alpha,s3,q3,psi,rs,xsCenter,ysCenter".split(",")
msg = r"""
    $t_0$ = {}
    $u_0$ = {}
    $t_E$ = {} d
    $s_2$ = {}
    $q_2$ = {}
    $s_3$ = {}
    $q_3$ = {}
    $\alpha$ = {}
    $\psi$ = {}
    """.format(lens_params["t0"], lens_params["u0"], lens_params["tE"] ,lens_params["s2"],lens_params["q2"],lens_params["s3"], lens_params["q3"],lens_params["alpha"],lens_params["psi"])
# print("msg", msg) t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs
main.text(0.0, 0.5, msg, transform=main.transAxes, fontdict = font)


plt.savefig("./data/lkv.png".format(xsCenter,ysCenter), dpi=300)

plt.show()

fig, ax = plt.subplots(figsize=(8,8))
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

line = ax.plot(Xs1, Ys1, c="g")[0] # source tragectory
add_arrow(line)

# circle1 = plt.Circle((Xs1[0], Ys1[0]), 2.2e-1, color=colors[1])
# circle2 = plt.Circle((Xs1[0], Ys1[0]), 2.2e-2, color=colors[2])
# circle3 = plt.Circle((Xs1[0], Ys1[0]), 2.2e-3, color=colors[3])
# circle4 = plt.Circle((Xs1[0], Ys1[0]), 2.2e-4, color=colors[4])
# ax.add_artist(circle1)
# ax.add_artist(circle2)
# ax.add_artist(circle3)
# ax.add_artist(circle4)

draw_circle(ax, Xs1[0], Ys1[0], 0.15, color =colors[1])
draw_circle(ax, Xs1[0], Ys1[0], 0.05, color =colors[2])
draw_circle(ax, Xs1[0], Ys1[0], 0.01, color =colors[3])
# draw_circle(ax, Xs1[0], Ys1[0], 2.2e-4, color =colors[4])
ax.scatter(Xs1[0], Ys1[0], color=colors[0], s = 0.5)
ax.set_xlabel(r"$x/ \theta_E $", fontdict = font)
ax.set_ylabel(r"$y/ \theta_E $", fontdict = font)



f = open(datapath+"lens_system.dat", "r")
#read in the source information
full_line = f.readline()
f.close()
line=full_line.split()

xs=np.float(line[0])
ys=np.float(line[1])
rs=np.float(line[2])
# print("xs, ys, rs in py", xs,ys,rs)

# nphi=150
# phi=np.linspace(0.0, 2*np.pi, nphi)
# sx=xs+rs*np.cos(phi)
# sy=ys+rs*np.sin(phi)
# ax.plot(sx, sy, 'k')

# ax.set_ylim(-1.25,1.75)
# ax.set_xlim(-1.1,2.1)
plt.axis('equal')

# ax.annotate('({}, {})'.format(xsCenter,ysCenter), xy=(0.2, 0.9), xycoords='axes fraction', fontsize=17,
                # horizontalalignment='right', verticalalignment='bottom')
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


plt.savefig("./data/lkv_geo.png", dpi=300)
plt.show()
