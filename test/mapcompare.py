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


'''####
# compare raysht result using differnet total ray number
# cpp code is npM_lens_inverse_ray_shoot.cpp in xxx/aeroastro/gravlen/critical_and_caustics/cpps/
maprayshtcppimgplane = read_cppmap1c("/Users/anything/THU/astro/softwares/aeroastro/gravlen/critical_and_caustics/cpps/datafromdpcg/magmap_imgplane_1280_1e5.dat")
maprayshtcppsrcplane = read_cppmap1c("/Users/anything/THU/astro/softwares/aeroastro/gravlen/critical_and_caustics/cpps/datafromdpcg/magmap_1280_1e5.dat")
print("maprayshtcppimgplane: ", maprayshtcppimgplane[:10])
print("maprayshtcppimgplane: ", maprayshtcppimgplane[-10:])

maprayshtcppimgplane = maprayshtcppimgplane.reshape(1280, 1280)
maprayshtcppsrcplane = maprayshtcppsrcplane.reshape(1280, 1280)

maprayshtcpp = np.log10(maprayshtcppsrcplane/maprayshtcppimgplane)
# maprayshtcpp1e5 = maprayshtcpp.reshape(1280, 1280)

maprayshtcppimgplane2e5 = read_cppmap1c("/Users/anything/THU/astro/softwares/aeroastro/gravlen/critical_and_caustics/cpps/datafromdpcg/magmap_imgplane_1280_2e5.dat")
maprayshtcppsrcplane2e5 = read_cppmap1c("/Users/anything/THU/astro/softwares/aeroastro/gravlen/critical_and_caustics/cpps/datafromdpcg/magmap_1280_2e5.dat")
print("maprayshtcppimgplane: ", maprayshtcppimgplane2e5[:10])
print("maprayshtcppimgplane: ", maprayshtcppimgplane2e5[-10:])

maprayshtcppimgplane2e5 = maprayshtcppimgplane2e5.reshape(1280, 1280)
maprayshtcppsrcplane2e5 = maprayshtcppsrcplane2e5.reshape(1280, 1280)

maprayshtcpp = np.log10(maprayshtcppsrcplane2e5/maprayshtcppimgplane2e5)
# maprayshtcpp2e5 = maprayshtcpp.reshape(1280, 1280)


maprayshtcpp128_1e5 = cutsubimg4finiteSource([1280, 1280], 0.05, (-2,2),(-2,2), 128, maprayshtcppsrcplane, maprayshtcppimgplane,-0.1,1,-0.6,0.6 )
maprayshtcpp128_2e5 = cutsubimg4finiteSource([1280, 1280], 0.05, (-2,2),(-2,2), 128, maprayshtcppsrcplane2e5, maprayshtcppimgplane2e5, -0.1,1,-0.6,0.6 )


mapraysht_relerr = (maprayshtcpp128_2e5 - maprayshtcpp128_1e5)/maprayshtcpp128_1e5
fig = plt.subplots(figsize=(16,8), dpi=100)
gs = gridspec.GridSpec(1,2)
plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.9, left = 0.1, hspace = 0, wspace = 0.3)
main = plt.subplot(gs[1])
main.tick_params(axis='both', labelsize = legend_tick_size, direction="in")
main.set_xlabel(r"$x/ \theta_E$", fontdict = font)
main.set_ylabel(r"$y/ \theta_E$", fontdict = font)
mapfig = main.imshow(mapraysht_relerr, extent=(-0.2,1,-0.6,0.6), origin='lower',cmap=cmap1)
# mapfig = main.imshow(np.log10(maprayshtcpp128_2e5), extent=(-0.2,1,-0.6,0.6), origin='lower',cmap=cmap1)
# mapfig = main.imshow(mapraysht_relerr, extent=(-2,2,-2,2), origin='lower',cmap=cmap1)
cbar = add_colorbar(mapfig)
cbar.set_label(r"Residual",fontdict=font)
cbar.ax.tick_params(labelsize=legend_tick_size)
cbar.ax.yaxis.get_offset_text().set_fontsize(legend_tick_size)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()
plt.title("Relative error between total ray number ${}$ and ${}$".format( f._formatSciNotation('%1.10e' %(1e10)), f._formatSciNotation('%1.10e' %(4e10)) ) , fontsize = 14)

main = plt.subplot(gs[0])
main.tick_params(axis='both', labelsize = legend_tick_size, direction="in")
main.set_xlabel(r"$x/ \theta_E$", fontdict = font)
main.set_ylabel(r"$y/ \theta_E$", fontdict = font)
# mapfig = main.imshow(mapraysht_relerr.T, extent=(-0.2,1,-0.6,0.6), origin='lower',cmap=cmap1)
mapfig = main.imshow(np.log10(maprayshtcpp128_2e5), extent=(-0.2,1,-0.6,0.6), origin='lower',cmap=cmap1)
# mapfig = main.imshow(mapraysht_relerr, extent=(-2,2,-2,2), origin='lower',cmap=cmap1)
cbar = add_colorbar(mapfig)
cbar.set_label(r"Magnification [Log]",fontdict=font)
cbar.ax.tick_params(labelsize=legend_tick_size)
cbar.ax.yaxis.get_offset_text().set_fontsize(legend_tick_size)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()

plt.show()
input(">>>>>>")


####'''

#speed comparison:
# magmapTimeRatio4 10 3
# magmapTimeRatio6_2 30 5
# magmapTimeRatio6_3 30 5

TIMRATIO = 0
TIMRATIOMAP = 1

if TIMRATIO:
    rhos, timvbbls, timtrils = readTimeRarioData(datapath+"magmapTimeRatio6_10_2.dat")
    timVBBL = [np.mean(i) for i in timvbbls]
    # print("timeVBBL: ", timVBBL, "\n\n")
    timTRIL = [np.mean(i) for i in timtrils]
    print(timvbbls[-1])
    print(timtrils[-1])
    stdVBBL = [np.std(i) for i in timvbbls]
    stdTRIL = [np.std(i) for i in timtrils]
    # print("stdTRIL", stdTRIL)



    timratio = [x/y for (x,y) in zip(timTRIL, timVBBL)]
    print("time ratio: ", timratio)


    stdRatio = [ timratio[i]*( (stdVBBL[i]/timVBBL[i])**2 + (stdTRIL[i]/timTRIL[i])**2 )**0.5  for i in range(len(rhos))]

    # print("rhos:", rhos, len(rhos))
    fig, ax= plt.subplots(figsize = (16,9))
    plt.subplots_adjust(top = 0.95, bottom = 0.15, right = 0.9, left = 0.1, hspace = 0, wspace = 0.35)
    # ax.plot(rhos, timratio)
    # ax.plot(rhos, timratio, "or")
    l1 = ax.errorbar(rhos, timratio, yerr=stdRatio, fmt="bo:", label = "Time ratio")
    ax.set_xscale('log')
    ax.tick_params(axis='both', labelsize = legend_tick_size, direction = "in", which='both')#direction="in", , width=1.5, length=4.5
    ax.set_xlabel(r"$\rho$", fontdict = font)
    ax.set_ylabel("Time ratio", fontdict = font, color = "b")
    # ax.legend(loc = 2)


    # ax_t = ax.secondary_xaxis('top')
    # ax_t.tick_params(axis='x', direction='in')



    ax2 = ax.twinx()
    ax2.set_ylabel('Time/s', fontdict = font)  # we already handled the x-label with ax1
    l2 = ax2.errorbar(rhos, timVBBL, yerr=stdVBBL, fmt="go:", label="Time (VBBL)", ms = 5)
    l3 = ax2.errorbar(rhos, timTRIL, yerr=stdTRIL, fmt="ko:", label="Time (Triple code)", ms = 5)
    ax2.tick_params(axis='y', labelsize = legend_tick_size, direction = "in")#, width=1.5, length=4.5
    # ax2.legend(loc = 9)

    # plt.legend([l1, l2, l3])

    # added these three lines
    # lns = l1[0]+l2[0]+l3[0]
    # labs = [l.get_label() for l in lns]
    plt.legend([l2, l3, l1], ["Time (VBBL)", "Time (Triple code)", "Time ratio"], loc=0, fancybox=0, framealpha=0.1, fontsize = 17)

    plt.savefig("./data/newimg_havexy/timeratio.png".format(xsCenter,ysCenter), dpi=300)
    # plt.show()
    # input(">>>")


if TIMRATIOMAP:
    # Xs0, Ys0 , mags0 = read_cppmap(datapath+"magmapTRIL_s0.8q0.1rs1e-2_new5.dat")
    # Xs0, Ys0 , magsVBB = read_cppmap(datapath+"magmapvbbl_s0.8q0.1rs1e-2_new5.dat")

    Xs0, Ys0 , mags0 = read_cppmap(datapath+"magmapTRIL_s0.8q0.1rs1e-1_new6.dat")
    Xs0, Ys0 , magsVBB = read_cppmap(datapath+"magmapvbbl_s0.8q0.1rs1e-1_new6.dat")
    
    filename = "vbblmap_data"
    dirpath = "/Users/anything/THU/astro/softwares/gravlens/triplelens/data/response_figure/data/"
    if not os.path.exists(dirpath+filename):
        print("saving vbbl mag map {}".format(dirpath+filename))
        np.savez(dirpath+filename, Xs0=Xs0, Ys0=Ys0, magsVBB = magsVBB)
    else:
        print("{} already exists, do you sure to overwrite?".format(dirpath+filename))

    # mapraysht0 = np.load(datapath+"rayshtmap_5e-2_128_raynum2e4.dat")
    # mapraysht = np.load(datapath+"rayshtmap_5e-2_128_raynum5e4.dat")
    # mapraysht = np.log10(mapraysht)


    Imgsize = 128
    mags0 = mags0.reshape(Imgsize, Imgsize)
    magsVBB = magsVBB.reshape(Imgsize, Imgsize)

    resmap = (mags0 - magsVBB)/magsVBB

    magsVBB = np.log10(magsVBB)
    mags0 = np.log10(mags0)


    # fig, ax = plt.subplots(figsize=(8,8))

    fig = plt.subplots(figsize=(14,6), dpi=100)
    gs = gridspec.GridSpec(1,2)
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.9, left = 0.1, hspace = 0, wspace = 0.35)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.9, left = 0.1, hspace = 0, wspace = 0.2)


    main = plt.subplot(gs[0])
    # ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    # ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    main.tick_params(axis='both', labelsize = legend_tick_size, direction="in")
    main.set_xlabel(r"$x/ \theta_E$", fontdict = font)
    main.set_ylabel(r"$y/ \theta_E$", fontdict = font)

    # mapfig = main.imshow(mags0.T, extent=(-0.5,0.5,-0.5,0.5), origin='lower',cmap=cmap1)
    mapfig = main.imshow(mags0.T, extent=(-0.5,0.5,-0.5,0.5), origin='lower',cmap="Reds")


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
    mapfig = secnd.imshow(resmap.T, extent=(-0.5,0.5,-0.5,0.5), origin='lower',cmap=cmap1, vmin=-7.5e-4,vmax=7.5e-4)
    # mapfig = secnd.imshow(resmap.T, extent=(-0.5,0.5,-0.5,0.5), origin='lower',cmap=cmap1)

    cbar = add_colorbar(mapfig)#format=ticker.FuncFormatter(fmt),format=OOMFormatter(-2, mathText=True)
    cbar.set_label(r"Residual",fontdict=font)
    cbar.ax.tick_params(labelsize=legend_tick_size)
    cbar.ax.yaxis.get_offset_text().set_fontsize(legend_tick_size)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.update_ticks()

    maxres = np.argmax( np.abs(resmap) )
    maxj = maxres % Imgsize
    maxi = maxres // Imgsize
    # print("max err place: xsCenter = {}; ysCenter = {}; ".format( Xs0[maxres], Ys0[maxres]  ) )
    print("max err place: xsCenter = {}; ysCenter = {}; err = {:.3e}".format( (Xs0.reshape(Imgsize, Imgsize))[maxi][maxj], (Ys0.reshape(Imgsize, Imgsize))[maxi][maxj], resmap[maxi][maxj] ))



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


    # plt.savefig("./data/newimg_havexy/timeratiomap_rs1e-2.png".format(xsCenter,ysCenter), dpi=300)


plt.show()
