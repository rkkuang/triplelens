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


iffinite_tris = []
iffinite_vbbs = []
magsts = []
magsbs = []
times = []
Xss = []
Yss = []
Residuals = []
timeratios = []

def add_1group_data(trifile, binfile, timfile):
    times0, Xs0, Ys0 , mags0t, iffinite_tri = read_cpplkv(datapath+trifile, raws = 5)
    times0, Xs0, Ys0 , mags0b, iffinite_vbb = read_cpplkv(datapath+binfile, raws = 5)
    Residual0 = (mags0t - mags0b)/mags0b

    dttri, dtvbb, dtratio = read_timratio(datapath+timfile)

    magsts.append(mags0t)
    magsbs.append(mags0b)

    iffinite_tris.append(iffinite_tri)
    iffinite_vbbs.append(iffinite_vbb)

    times.append(times0)
    Xss.append(Xs0)
    Yss.append(Ys0)
    Residuals.append(Residual0) 
    timeratios.append(dtratio)   

datapath = "./data/"

lens_params = read_lens_system_triple(datapath+"lens_system_triple.dat")
# xsCenter, ysCenter = lens_params["xsCenter"], lens_params["ysCenter"]





# u0s = []
# add_1group_data("flkvcmp_tri_rs1e-2_u08e-1.dat", "flkvcmp_bin_rs1e-2_u08e-1.dat", "flkvcmp_tim_rs1e-2_u08e-1.dat")
# u0s.append(0.8)
# add_1group_data("flkvcmp_tri_rs1e-2_u05e-1.dat", "flkvcmp_bin_rs1e-2_u05e-1.dat", "flkvcmp_tim_rs1e-2_u05e-1.dat")
# u0s.append(0.48)
# add_1group_data("flkvcmp_tri_rs1e-2_u0-3e-1.dat", "flkvcmp_bin_rs1e-2_u0-3e-1.dat", "flkvcmp_tim_rs1e-2_u0-3e-1.dat")
# u0s.append(-0.24)

def getnames(strrs, stru0, limb = False):
    if limb:
        fnametri = "limblkvcmp_tri_rs";
        fnameb = "limblkvcmp_bin_rs";
        fnametim = "limblkvcmp_tim_rs";  
    else:      
        fnametri = "flkvcmp_tri_rs";
        fnameb = "flkvcmp_bin_rs";
        fnametim = "flkvcmp_tim_rs";
    final_name1 = fnametri+strrs+"u0"+stru0+".dat";
    final_name2 = fnameb+strrs+"u0"+stru0+".dat";
    final_name3 = fnametim+strrs+"u0"+stru0+".dat";
    return final_name1, final_name2, final_name3

limb = True
limba1 = 0.51
strrs = "1e-3";

u0s = []
stru0 = "-0.24";
final_name1, final_name2, final_name3 = getnames(strrs, stru0, limb = limb)
print(final_name1)
print(final_name2)
print(final_name3)
add_1group_data(final_name1, final_name2, final_name3)
u0s.append(-0.24)

stru0 = "0.8";
final_name1, final_name2, final_name3 = getnames(strrs, stru0, limb = limb)
add_1group_data(final_name1, final_name2, final_name3)
u0s.append(0.8)

stru0 = "0.48";
final_name1, final_name2, final_name3 = getnames(strrs, stru0, limb = limb)
add_1group_data(final_name1, final_name2, final_name3)
u0s.append(0.48)



# times0, Xs0, Ys0 , mags0t = read_cpplkv(datapath+"flkvcmp_tri_rs1e-2_u05e-1.dat")
# times0, Xs0, Ys0 , mags0b = read_cpplkv(datapath+"flkvcmp_bin_rs1e-2_u05e-1.dat")
# Residual0 = (mags0t - mags0b)/mags0b

# magsts.append(mags0t)
# magsbs.append(mags0b)
# times.append(times0)
# Xss.append(Xs0)
# Yss.append(Ys0)
# Residuals.append(Residual0)



# times0, Xs0, Ys0 , mags0t = read_cpplkv(datapath+"flkvcmp_tri_rs1e-2_u0-3e-1.dat")
# times0, Xs0, Ys0 , mags0b = read_cpplkv(datapath+"flkvcmp_bin_rs1e-2_u0-3e-1.dat")
# Residual0 = (mags0t - mags0b)/mags0b

# magsts.append(mags0t)
# magsbs.append(mags0b)
# times.append(times0)
# Xss.append(Xs0)
# Yss.append(Ys0)
# Residuals.append(Residual0)




Nexp = len(magsts)


# times4, Xs4, Ys4 , mags4 = read_cpplkv(datapath+"ftrilkv_rs_2.2e-4.dat")
# mags0t = np.log10( mags0t )
# mags0b = np.log10( mags0b )
# mags1t = np.log10( mags1t )
# mags1b = np.log10( mags1b )


for i in range(Nexp):
    magsts[i] = np.log10( magsts[i] )
    magsbs[i] = np.log10( magsbs[i] )


fig = plt.subplots(figsize=(13,7), dpi=100)
gs = gridspec.GridSpec(2,1,height_ratios=[5,1])
plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.9, left = 0.1, hspace = 0, wspace = 0)
main = plt.subplot(gs[0])

main.set_ylabel(r"log($\mu$)", fontdict=font)
main.set_xlabel('HJD - 2450000', fontdict=font)
main.tick_params(axis='both', labelsize = legend_tick_size, direction = "in")#labelcolor=color
# main.yaxis.set_minor_locator(AutoMinorLocator(5))

colors = ["k","b", "g", "cyan"]

# main.plot(times0, mags0t, color=colors[0],linewidth = 2)
# main.plot(times0, mags0b, color=colors[0],linewidth = 2)
# main.plot(times1, mags1t, color=colors[1],linewidth = 2)
# main.plot(times1, mags1b, color=colors[1],linewidth = 2)

for i in range(Nexp):
    main.scatter(times[i], magsts[i], color=colors[i], marker = ".")
    main.scatter(times[i], magsbs[i], color=colors[i], marker = "*")

    main.plot(times[i], magsts[i], color=colors[i],linestyle = "-.",linewidth = 2)
    main.plot(times[i], magsbs[i], color=colors[i],linewidth = 2, label = r"$u_0 = {}$".format(u0s[i]))

# plt.legend(loc = 2, fontsize = 17, frameon=False)
leg = main.legend(prop = font2, frameon=False, loc="upper left")
# for text,color in zip(leg.get_texts(), colors):
#     text.set_color(color)



# main.plot(times0, mags4, color=colors[4],linewidth = 2)

# "t0,u0,tE,s2,q2,alpha,s3,q3,psi,rs,xsCenter,ysCenter".split(",")
msg = r"""
    $t_0 = {}$
    $t_E = {}$ d
    $s = {}$
    $q = {}$
    $\alpha = {}$
    $\rho = {}$
    {}
    """.format(lens_params["t0"], lens_params["tE"] ,lens_params["s2"],lens_params["q2"],lens_params["alpha"], lens_params["rs"],
        "limb a = {}".format(limba1) if limb else ""
        )
# print("msg", msg) t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs
main.text(0.25, 0.6, msg, transform=main.transAxes, fontdict = font)



# ax_traj = plt.subplot(gs[1])
ax_traj = plt.subplot(gs[1],sharex=main)
ax_traj.axhline(0,color='k')
ax_traj.set_xlabel('HJD - 2450000', fontdict=font, color='k')
ax_traj.set_ylabel('Residual',fontdict=font, color='k')
ax_traj.xaxis.set_tick_params(labelsize=legend_tick_size)
ax_traj.yaxis.set_tick_params(labelsize=legend_tick_size)
ax_traj.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# yi = 10**yi + 1 # interpolation
# combinemags = 10**combinemags + 1

for i in range(Nexp):
    ax_traj.plot( times[i], Residuals[i], c=colors[i], linewidth = 2 )
    ax_traj.scatter( times[i], Residuals[i],marker = ".", s=10, facecolors='none', edgecolors=colors[i])

# ax_traj.plot( times1, Residual1, c=colors[1], linewidth = 2 )
# ax_traj.scatter( times1, Residual1,marker = ".", s=10, facecolors='none', edgecolors=colors[1])


temprelerr = abs(Residuals[0])
maxidx = np.argmax(temprelerr)
print("idx = {}, maxerr: {}, maxerr time: {:.15f}".format(maxidx, temprelerr[maxidx], times[0][maxidx]))



# fig, ax = plt.subplots(figsize=(8,8))
ax2 = fig[0].add_axes([.54, .54, .4, .4], aspect=1)
# ax2.plot(caustic_curves[:, 0], caustic_curves[:, 1], 'k-')
# ax2.plot(y1, y2, 'k--')

# ax2.grid(True)
# ax2.set_xlim(-1,1)
# ax2.set_ylim(-1,1)

ax2.tick_params(axis='both', labelsize = legend_tick_size, direction="in")
cx, cy = readFile(datapath+"caustics.dat", 0, 1, expected_elem_each_row=2)
ax2.plot(cx, cy, '-', color='r', markersize=1)
x, y = readFile(datapath+"critical_curves.dat", 0, 1, expected_elem_each_row=2)
ax2.plot(x, y, '--', color='r', markersize=1)

lensx, lensy = readFile(datapath+"lens_system.dat", 1, 2, expected_elem_each_row=3)
lensm, _ = readFile(datapath+"lens_system.dat", 0, 2, expected_elem_each_row=3)

for i in range(len(lensm)-1):
    ax2.plot(lensx[1+i], lensy[1+i], '+', color='k', markersize=10)#
    # plt.plot(lensx[1+i], lensy[1+i], '+', color='k', markersize=15)#
# ax2.axis("equal")

# print(np.shape(iffinite_vbbs[0]))

for i in range(Nexp):
    line0 = ax2.plot(Xss[i], Yss[i], c=colors[i])[0] # source tragectory
    add_arrow(line0)
    ax2.text(Xss[i][0]-0.2, Yss[i][0], r"${:.1f}$".format(timeratios[i]))#transform=main.transAxes,
    if not limb:
        ax2.text(Xss[i][-1]-0.2, Yss[i][-1], r"${:.1f}$".format( np.sum(iffinite_vbbs[i])/np.shape(iffinite_vbbs[i])[0] ))
        ax2.text(Xss[i][-1]-0.2, Yss[i][-1]-0.2, r"${:.1f}$".format(np.sum(iffinite_tris[i])/np.shape(iffinite_tris[i])[0]))

    # line1 = ax2.plot(Xs1, Ys1, c=colors[1])[0] # source tragectory
# add_arrow(line1)


plt.savefig("./data/lkvcomp_rs{}_{}.png".format(strrs, "limb" if limb else ""), dpi=300)

plt.show()







'''

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
    plt.plot(lensx[1+i], lensy[1+i], '+', color='k', markersize=15)#markersize=5*lensm[i+1]

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
'''
