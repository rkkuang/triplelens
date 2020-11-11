from utils import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl  
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, FormatStrFormatter

from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import gridspec
from scipy import interpolate
import matplotlib.ticker as mticker
from matplotlib.ticker import FuncFormatter   ### 今天的主角
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, FormatStrFormatter

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


f = mticker.ScalarFormatter(useOffset=False, useMathText=True)


def plot_1D_adap_samp():
    lens_params = read_lens_system_triple("./data/lens_system_triple.dat")

    # times, _, _, mags = read_cpplkv( "./data/ftrilkv_rs2.2e-4_adap.dat")
    # adaptive_times, adaptive_mags = read_cpplkv_adap( "./data/ftrilkv_adaptive_rs2.2e-4.dat")

    times, _, _, mags = read_cpplkv( "./data/ftrilkv_rs5e-3_adap_v2.dat")
    adaptive_times, adaptive_mags = read_cpplkv_adap( "./data/ftrilkv_adaptive_rs5e-3_v2.dat")
    lens_params["rs"] = 5e-3

    # times, _, _, mags = read_cpplkv( "./data/ftrilkv_rs1e-3_adap_v2.dat")
    # adaptive_times, adaptive_mags = read_cpplkv_adap( "./data/ftrilkv_adaptive_rs1e-3_v2.dat")

    print("full lightcurve: t0 = {}, tend = {}".format(times[0], times[-1]))

    mags = np.array(mags)
    mags = np.log10( mags - 1 )
    

    combinetimes = np.concatenate((adaptive_times, times))
    combinemags = np.concatenate((adaptive_mags, mags))
    stidx = np.argsort(combinetimes)
    combinetimes = combinetimes[stidx]
    combinemags = combinemags[stidx]

    k = 1 # The degree of the spline fit
    tck, fp, ier, msg = interpolate.splrep(adaptive_times, (adaptive_mags), s=0,k=k, full_output = 1)# Given degree of the spline (k=0) is not supported. (1<=k<=5)
    yi = interpolate.splev(combinetimes, tck, der=0)

    print("len times:", len(times))
    fig = plt.subplots(figsize=(14,7), dpi=100)
    gs = gridspec.GridSpec(2,1,height_ratios=[5,1])
    plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0, wspace = 0)
    main = plt.subplot(gs[0])

    main.set_ylabel(r"log($\mu$ - 1)", fontdict=font)
    main.tick_params(axis='y', labelsize = legend_tick_size)#labelcolor=color
    main.yaxis.set_minor_locator(AutoMinorLocator(5))

    colors = ['tab:green', 'r',"b", "k"]#orange
    main.plot(combinetimes, (combinemags), color=colors[0],linewidth = 2, label = "Full ${}$ points".format(f._formatSciNotation('%1.10e' %len(mags))))
    main.plot(adaptive_times, (adaptive_mags), '.', color=colors[1], label = "Adaptive {} points".format(len(adaptive_mags)))
    main.plot(combinetimes, (yi), color=colors[2], linewidth = 1, label = "Interpolation, order = {}".format(k).format(len(mags)))
    leg = main.legend(prop = font2, frameon=False, loc="upper right") #, ncol=2
    for text,color in zip(leg.get_texts(), colors):
        text.set_color(color)
    # plt.show()
    # add text, t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs
    msg = r"""
    $t_0$ = {}
    $u_0$ = {}
    $t_E$ = {} d
    $\rho$ = ${}$, $2\rho t_E$ = ${}$ m
    full light curve: cadence = {:.2f} m (${}$ points)
    adaptive sampling: ${}$ points
    """.format(lens_params["t0"], lens_params["u0"], lens_params["tE"] ,f._formatSciNotation('%1.10e' % lens_params["rs"]) ,f._formatSciNotation('%1.1e' %  (2*lens_params["tE"]*lens_params["rs"]*24*60) ),(times[-1]-times[0])/len(times)*24*60, f._formatSciNotation('%1.10e' %len(times)), len(adaptive_times))
    print("msg", msg)
    main.text(0.0, 0.6, msg, transform=main.transAxes, fontdict = font)

    # ax_traj = plt.subplot(gs[1])
    ax_traj = plt.subplot(gs[1],sharex=main)
    ax_traj.axhline(0,color='k')
    ax_traj.set_xlabel('HJD - 2450000', fontdict=font, color='k')
    ax_traj.set_ylabel('Residual',fontdict=font, color='k')
    ax_traj.xaxis.set_tick_params(labelsize=legend_tick_size)
    ax_traj.yaxis.set_tick_params(labelsize=legend_tick_size)
    ax_traj.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    yi = 10**yi + 1 # interpolation
    combinemags = 10**combinemags + 1
    ax_traj.plot( combinetimes, (yi-combinemags)/combinemags, c="k", linewidth = 2 )
    ax_traj.scatter( combinetimes, (yi-combinemags)/combinemags,marker = ".", s=10, facecolors='none', edgecolors='k')

    temprelerr = abs((yi-combinemags)/combinemags)
    maxidx = np.argmax(temprelerr)
    print("idx = {}, maxerr: {}, maxerr time: {:.15f}".format(maxidx, temprelerr[maxidx], combinetimes[maxidx]))

    # ax_traj.set_ylim(-5e-3, 5e-3)
    ax_traj.set_ylim(-6e-4, 6e-4)
    # https://blog.csdn.net/monotonomo/article/details/83826621
    from matplotlib.ticker import FuncFormatter   ### 今天的主角
    def formatnum(x, pos):
        if x != 0:
            # return '$%.0f$x$10^{-3}$' % (x/1e-3)
            return '$%.0f$x$10^{-4}$' % (x/1e-4)
        else:
            return 0
    formatter = FuncFormatter(formatnum)
    ax_traj.yaxis.set_major_formatter(formatter)


plot_1D_adap_samp()

# plt.savefig("./data/newimg_havexy/lkvadav1e-3.png", dpi=300)
plt.show()
