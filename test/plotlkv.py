from utils import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl  
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, FormatStrFormatter

from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import gridspec

mpl.rc('font',family='Times New Roman')

legend_tick_size = 25
font = {'family': 'Times New Roman',
    'color':  'k',
    'weight': 'normal',
    'size': legend_tick_size,
    }

font2 = {'family': 'Times New Roman',
    'weight': 'normal',
    'size': legend_tick_size,
    }

datapath = "./data/"

lens_params = read_lens_system_triple(datapath+"lens_system_triple.dat")
xsCenter, ysCenter = lens_params["xsCenter"], lens_params["ysCenter"]

times0, Xs0, Ys0 , mags0 = read_cpplkv(datapath+"fulltrilkv_rs5e-3.dat")
mags0 = np.log10( mags0 )
rho = 5e-3

legend_tick_size = 17
font = {'family': 'Times New Roman',
    'color':  'k',
    'weight': 'normal',
    'size': legend_tick_size,
    }

font2 = {'family': 'Times New Roman',
    'weight': 'normal',
    'size': legend_tick_size,
    }

fig = plt.subplots(figsize=(13,7), dpi=100)
main = plt.subplot()
main.set_ylabel(r"log($\mu$)", fontdict=font)
main.set_xlabel('HJD - 2450000', fontdict=font)
main.tick_params(axis='both', labelsize = legend_tick_size, direction = "in")#labelcolor=color
main.plot(times0, mags0, color="r",linewidth = 2, label = r"$\rho={}$".format(rho))
leg = main.legend(prop = font2, frameon=False, loc="upper right")
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
main.text(0.0, 0.5, msg, transform=main.transAxes, fontdict = font)
plt.show()

if 1:
    times0, Xs0, Ys0 , mags0 = read_cpplkv(datapath+"unilkv.dat")
    mags0 = np.log10( mags0 )

    times1, Xs1, Ys1 , mags1 = read_cpplkv(datapath+"limblkv.dat")
    mags1 = np.array(mags1)
    mags1 = np.log10( mags1)


    fig = plt.subplots(figsize=(13,7), dpi=100)
    main = plt.subplot()
    main.set_ylabel(r"log($\mu$)", fontdict=font)
    main.set_xlabel('HJD - 2450000', fontdict=font)
    main.tick_params(axis='both', labelsize = legend_tick_size, direction = "in")#labelcolor=color

    main.plot(times0, mags0, color="r",linewidth = 2, label = "Without limb-darkening")
    main.plot(times0, mags1, color="b",linewidth = 2, label = "With limb-darkening")
    leg = main.legend(prop = font2, frameon=False, loc="upper right")
    plt.show()











