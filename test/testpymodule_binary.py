"""
Test triple code using binary lens system

To compute binary lens magnification, please change the "NLENS" parameter in TripleLensingLibrary.h
and reinstall the TripleLensing python module by executing:

sudo python3 setup install

Compared with MulensModel package
"""

import TripleLensing
from MulensModel.binarylens import BinaryLens

def getbi_light_curve(bilens, xs, ys, rs, u_limb_darkening = None, accuracy = 1e-3):
    mus = []
    for x,y in zip(xs, ys):
        tempmag = bilens.vbbl_magnification(x, y, rs, u_limb_darkening=u_limb_darkening, accuracy=accuracy)
        mus.append(tempmag)
    return np.array(mus)

# from mulens model, https://github.com/rpoleski/MulensModel/tree/master/source
"""
Linear limb-darkening parameters. Both *gamma* and *u* conventions
can be used.  The *u* convention is more frequently used in
studies other than microlensing.  It has fixed flux at the center.
`An et al. 2002 (ApJ 572, 521)
<https://ui.adsabs.harvard.edu/abs/2002ApJ...572..521A/abstract>`_
introduced the *gamma* convention:
gamma = (2 * u) / (3 - u)
u = 3 * gamma / (2 + gamma)
Note that the gamma convention has fixed total flux.
"""

# The two functions below implement convention introduced by:
# An et al. 2002 (ApJ 572, 521)
# https://ui.adsabs.harvard.edu/abs/2002ApJ...572..521A/abstract
def gamma_to_u(gamma):
    """
    Transform gamma limb darkening coefficient to u in convention
    introduced by `An et al. 2008 (ApJ 681, 1593)
    <https://ui.adsabs.harvard.edu/abs/2002ApJ...572..521A/abstract>`_.
    Parameters :
        gamma: *float*
            Limb darkening coefficient in gamma convention.
    Returns :
        u: *float*
            Limb darkening coefficient in u convention.
    """
    return 3. * gamma / (2. + gamma)

def u_to_gamma(u):
    """
    Transform u limb darkening coefficient to gamma in convention
    introduced by `An et al. 2008 (ApJ 681, 1593)
    <https://ui.adsabs.harvard.edu/abs/2002ApJ...572..521A/abstract>`_.
    Parameters :
        u: *float*
            Limb darkening coefficient in u convention.
    Returns :
        gamma: *float*
            Limb darkening coefficient in gamma convention.
    """
    return (2. * u) / (3. - u)

s = 0.8
q = 0.1
mlens = [1.0/(1.0+q), q/(1.0 + q)]
zlens = [-q * s / (1.0 + q), 0, s / (1.0 + q), 0]
rs = 0.01





bilens = BinaryLens(mass_1=1/(1.0+q), mass_2=q/(1.0 + q), separation=s)
TRIL = TripleLensing.TripleLensing()

# mag = bilens.vbbl_magnification(xsCenter, ysCenter, rs, gamma=0.51, accuracy=0.001)
# print(mag)



# #plot lens system topology, critical curves and caustics
from utils import *
# xsCenter=0.145626;
# ysCenter = 0.007150
# plot_critcaus_srcimgs(mlens, zlens, xsCenter, ysCenter, rs)
# plt.show()


# test light curve generation:
# triple lens event parameters
t0 = 7494.153;
u0 = 0.021;
tE = 74.62;
# s2 = 1.396;
# q2 = 0.029;
alpha = 2.948; #//rad
# s3 = 1.168;
# q3 = 3.27e-3;
# psi = 5.332; #//rad
# rs = 0.01#0.22e-3 #0.005#
salpha = np.sin(alpha)
calpha = np.cos(alpha)
# params = [t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs]
# source position
ts = np.linspace(7470, 7510, 200) # 1000个点时，在index = 334 的点会卡住
tn = (ts - t0) / tE;
y1s = u0 * salpha + tn * calpha;
y2s = u0 * calpha - tn * salpha;
start, end = 0, -1
# start, end = 166, 168
# start, end = 334, 336







########################################
########################################
########################################
# compare with vbbl using binary system, without limb darkening

# print("generating light curve triple-code ...")
# mus = TRIL.TriLightCurve(mlens, zlens, y1s[start:end], y2s[start:end], rs, 45, 2, 1e-5, 1e-3)

# print("generating light curve vbbl ...")
# bimags = getbi_light_curve(bilens, y1s[start:end], y2s[start:end], rs, u_limb_darkening = 0)
# main, gs = pltlkv(ts[start:end], mus, label = "triple")
# main.scatter( ts[start:end], np.log10(mus) ,marker = ".", s=30, facecolors='none', edgecolors='r')

# main.plot(ts[start:end], np.log10( bimags ), color="b",linewidth = 1, label = "binary")
# plt.legend()
# msg = r"""
#     without limb darkening
#     $t_0$ = {}
#     $u_0$ = {}
#     $t_E$ = {} d
#     $s$ = {}
#     $q$ = {}
#     $\alpha$ = {}
#     $\rho$ = {:.1e}
#     """.format(t0, u0, tE, s, q, alpha, rs)
# main.text(0.0, 0.4, msg, transform=main.transAxes, fontdict = font)
# ax_traj = plt.subplot(gs[1],sharex=main)
# ax_traj.axhline(0,color='k')
# ax_traj.set_xlabel('HJD - 2450000', fontdict=font, color='k')
# ax_traj.set_ylabel('Residual',fontdict=font, color='k')
# ax_traj.xaxis.set_tick_params(labelsize=17)
# ax_traj.yaxis.set_tick_params(labelsize=17)
# ax_traj.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax_traj.plot( ts[start:end], (mus-bimags)/bimags, c="k", linewidth = 2 )
# ax_traj.scatter( ts[start:end], (mus-bimags)/bimags ,marker = ".", s=10, facecolors='none', edgecolors='k')

# temprelerr = abs((mus-bimags)/bimags)
# maxidx = np.argmax(temprelerr)
# print("idx = {}, maxerr: {}, maxerr time: {:.15f}".format(maxidx, temprelerr[maxidx], ts[start:end][maxidx]))
# plt.show()
# print("muvbbl = {:.15f}, mutri = {:.15f}, err = {:.5e}".format( bimags[maxidx], mus[maxidx], (mus[maxidx]-bimags[maxidx])/bimags[maxidx] ))
# input("without limb darkening finished >>>")





# mus = TRIL.TriLightCurve(mlens, zlens, y1s[start:end], y2s[start:end], 0.0000809248218210688790913, 45, 2, 1e-5, 1e-3)
# mus = TRIL.TriLightCurve(mlens, zlens, [0.145275], [0.007081], 0.0000809248218210688790913, 45, 2, 1e-5, 1e-3)

# print(rs*0.36784009918667670558, mus)
# plot_critcaus_srcimgs(mlens, zlens, y1s[335], y2s[335], 0.0000809248218210688790913)
# plot_critcaus_srcimgs(mlens, zlens, 0.145275, 0.007081, 0.000081)
# plt.show()

# print(
# """
# rs = {:.25f}
# xs = {:.25f}
# ys = {:.25f}
# mlens = {}
# zlens = {}
# """.format(rs, y1s[335], y2s[335], mlens, zlens)
#     )

# input("Stop")

##################################################
##################################################
##################################################
'''
special case:
rs = 0.0002200000000000000078323
xs = 0.1452752184254198497548316
ys = 0.0070811503507445944238796
mlens = [0.968738798957637, 0.028093425169771, 0.003167775872591]
zlens = [-0.039343051506317, 0, 1.356656948493683, 0, 0.6389361960108, -0.950873946634155]


'''


########################################
########################################
########################################
# compare with vbbl using binary system, with limb darkening
RelTolLimb = 1e-4
AbsTolLimb = 1e-5
Gamma = 0.51#0.51 #0.51
limbu = gamma_to_u(Gamma)

print("generating light curve with limb darkening ...")
# mus = TRIL.TriLightCurve(mlens, zlens, y1s[start:end], y2s[start:end], rs, 45, 2, 1e-5, 1e-3)
mus = TRIL.TriLightCurveLimb(mlens, zlens, y1s[start:end], y2s[start:end], rs, 45, 3, 1e-5, 1e-3, RelTolLimb, AbsTolLimb, limbu)
print("generating vbbl light curve ...")
bimags = getbi_light_curve(bilens, y1s[start:end], y2s[start:end], rs, u_limb_darkening = limbu, accuracy = 1e-5)
bimagsnolimb = getbi_light_curve(bilens, y1s[start:end], y2s[start:end], rs, u_limb_darkening = 0)


main, gs = pltlkv(ts[start:end], mus, label = "triple")
main.scatter( ts[start:end], np.log10(mus) ,marker = ".", s=30, facecolors='none', edgecolors='r')

main.plot(ts[start:end], np.log10( bimags ), color="b",linewidth = 1, label = "binary")
main.plot(ts[start:end], np.log10( bimagsnolimb ), color="g",linewidth = 1, label = "without limb")


plt.legend()
msg = r"""
    with limb darkening
    $t_0$ = {}
    $u_0$ = {}
    $t_E$ = {} d
    $s$ = {}
    $q$ = {}
    $\alpha$ = {}
    $\rho$ = {:.1e}
    Gamma = {}
    """.format(t0, u0, tE, s, q, alpha, rs, Gamma)
main.text(0.0, 0.4, msg, transform=main.transAxes, fontdict = font)
ax_traj = plt.subplot(gs[1],sharex=main)
ax_traj.axhline(0,color='k')
ax_traj.set_xlabel('HJD - 2450000', fontdict=font, color='k')
ax_traj.set_ylabel('Residual',fontdict=font, color='k')
ax_traj.xaxis.set_tick_params(labelsize=17)
ax_traj.yaxis.set_tick_params(labelsize=17)
ax_traj.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_traj.plot( ts[start:end], (mus-bimags)/bimags, c="k", linewidth = 2 )
ax_traj.scatter( ts[start:end], (mus-bimags)/bimags ,marker = ".", s=10, facecolors='none', edgecolors='k')

temprelerr = abs((mus-bimags)/bimags)
maxidx = np.argmax(temprelerr)
print("idx = {}, maxerr: {}, maxerr time: {:.15f}".format(maxidx, temprelerr[maxidx], ts[start:end][maxidx]))
print("muvbbl = {:.15f}, mutri = {:.15f}, mu_no_limb = {:.15f}, err = {:.5e}".format( bimags[maxidx], mus[maxidx],bimagsnolimb[maxidx] , (mus[maxidx]-bimags[maxidx])/bimags[maxidx] ))
plt.show()
input("with limb darkening finished >>>")


plt.show()






