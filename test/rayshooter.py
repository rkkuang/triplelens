# Testing by ray shooting

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import math, os
from matplotlib import gridspec
from numba import jit

from numpy import logical_and, where, sqrt, load, save, array
from numpy import abs as npabs
from numpy import sum as npsum
from matplotlib.ticker import MaxNLocator, AutoMinorLocator, FormatStrFormatter
from mpl_toolkits.axes_grid.inset_locator import inset_axes as iax

# https://blog.csdn.net/monotonomo/article/details/83826621
from matplotlib.ticker import FuncFormatter   ### 今天的主角
import time

M_PI = math.pi
M_PI2 = 2*M_PI

import TripleLensing

font = {'family': 'Times New Roman',
    'color':  'k',
    'weight': 'normal',
    'size': 17,
    }
legend_tick_size = 17
font2 = {'family': 'Times New Roman',
    'weight': 'normal',
    'size': 16,
    }


datapath = "./data/rayshtdata/"

def getbi_light_curve(bilens, xs, ys, rs, u_limb_darkening = None, accuracy = 1e-3):
    mus = []
    for x,y in zip(xs, ys):
        tempmag = bilens.vbbl_magnification(x, y, rs, u_limb_darkening=u_limb_darkening, accuracy=accuracy)
        mus.append(tempmag)
    return np.array(mus)


@jit(nopython=True)
def genxy(x, y):
    X,Y = np.meshgrid(x,y)
    return X.reshape(1,-1)[0], Y.reshape(1,-1)[0]

@jit(nopython=True)
def inverse_ray_shooting_use_numba(thetax, thetay, mlens, xs, ys):
    betax, betay = thetax.copy(), thetay.copy()
    for i in range(len(xs)):
        Ri2s = (thetax-xs[i])**2 + (thetay-ys[i])**2
        betax -= mlens[i]*(thetax - xs[i])/Ri2s
        betay -= mlens[i]*(thetay - ys[i])/Ri2s
    return betax, betay

@jit(nopython=True)
def findidxs(betax, betay , x, y, rho, rho2):
    # https://stackoverflow.com/questions/13869173/numpy-find-index-of-the-elements-within-range
    idxs = where(logical_and( npabs(betax - x) <= rho, npabs(betay - y) <= rho))[0]
    bxs = betax[idxs]
    bys = betay[idxs]
    r2s = (bxs - x)*(bxs - x) + (bys - y)*(bys - y)
    # print("len(r2s) ",len(r2s))
    # idxs_r2_lessthan_1 = where( r2s <= 1 )[0] # wrong, not r2s <= 1
    idxs_r2_lessthan_1 = where( r2s <= rho2 )[0]
    r2s_inside_source = r2s[idxs_r2_lessthan_1]
    # print("len(r2s_inside_source) ",len(r2s_inside_source))
    return r2s_inside_source

@jit(nopython=True)
def limb_linear(r2, coe1, coe2):
    # when call this function, r should normalize to source raduis rho 
    
    # return (coe1 + coe2 * sqrt(1-r2))


    # return 2*sqrt(r2)*(coe1 + coe2 * sqrt(1-r2))
    return coe1 + coe2 * sqrt(1-r2)


class MicroLens(object):
    def __init__(self, mlens, xlens, ylens):
        self.mlens = array(mlens)
        self.xs = array(xlens)
        self.ys = array(ylens)
        self.betax, self.betay = None, None
        self.mag = None

    def gotstr(self, mlens, xs, ys, xlim, ylim, raynum, strlen = 30):
        tempstr = "MLENS{}XLENS{}YLENS{}XLIM{}YLIM{}RAYNUM{}".format(self.mlens, self.xs, self.ys, xlim, ylim, raynum)
        tempstr = "".join(x for x in tempstr if x.isalnum())
        lenstr = len(tempstr)
        if strlen > lenstr:
            self.betathetaxy_string = tempstr
            # return tempstr
        else:
            step = lenstr // strlen
            step = 1 if step == 0 else step
            self.betathetaxy_string = tempstr[::step]
            # return tempstr[::step]

    def ray_generator(self, xlim, ylim, raynum, datatype = np.double):
        """
        generate rays on lens plane at range xlim and ylim,
        with number of rays at unit thetaE specified by raynum
        """

        # self.betathetaxy_string = self.gotstr(self.mlens, self.xs, self.ys, xlim, ylim, raynum)
        self.raynum_perarea = raynum * raynum

        self.thetax_filename = datapath + "thetax" + self.betathetaxy_string + ".npy"
        self.thetay_filename = datapath + "thetay" + self.betathetaxy_string + ".npy"
        
        if os.path.exists(self.thetax_filename):
            print("Loading thetax thetay ...")
            print(self.thetax_filename)
            self.thetax = load(self.thetax_filename)
            self.thetay = load(self.thetay_filename)
        else:
            xlen = xlim[1] - xlim[0]
            ylen = ylim[1] - ylim[0]
            xyarea = xlen * ylen
            
            # totalraynum = raynum_perarea * xyarea
            print("Generating thetax thetay ...")
            x = np.linspace(xlim[0],xlim[1], int(xlen*raynum) ).astype(datatype)#np.single, np.double
            y = np.linspace(ylim[0],ylim[1], int(ylen*raynum) ).astype(datatype)

            # self.thetax, self.thetay =  genxy(x, y)# shape: (3,)

            X,Y = np.meshgrid(x,y)
            self.thetax, self.thetay = X.reshape(1,-1)[0], Y.reshape(1,-1)[0]

            print("Saving thetax thetay ...")
            save(self.thetax_filename, self.thetax)
            save(self.thetay_filename, self.thetay)


    # def ray_loader(self, xlim, ylim, raynum):
    #     xlen = xlim[1] - xlim[0]
    #     ylen = ylim[1] - ylim[0]
    #     xyarea = xlen * ylen
    #     self.raynum_perarea = raynum * raynum
    #     # self.thetax, self.thetay = 

    # def ray_saver(self):
    #     # "".join(x for x in s if x.isalnum())

    #     pass

    def inverse_ray_shooting(self):

        self.betax_filename = datapath + "betax" + self.betathetaxy_string + ".npy"
        self.betay_filename = datapath + "betay" + self.betathetaxy_string + ".npy"

        if os.path.exists(self.betax_filename):
            print("Loading betax betay ...")
            print(self.betax_filename)
            self.betax = load(self.betax_filename)
            self.betay = load(self.betay_filename)
        else:
            print("Generating betax betay ...")

            # self.betax, self.betay = self.thetax.copy(), self.thetay.copy()

            self.betax, self.betay = inverse_ray_shooting_use_numba(self.thetax, self.thetay, self.mlens, self.xs, self.ys)

            # for i in range(len(self.xs)):
            #     Ri2s = (self.thetax-self.xs[i])**2 + (self.thetay-self.ys[i])**2
            #     self.betax -= self.mlens[i]*(self.thetax - self.xs[i])/Ri2s
            #     self.betay -= self.mlens[i]*(self.thetay - self.ys[i])/Ri2s

            print("Saving betax betay ...")
            save(self.betax_filename, self.betax)
            save(self.betay_filename, self.betay)


    def limb_linear_coefficients(self, u = None, Gamma = None):
        if not u:
            if not Gamma:
                print("At least given u or Gamma parameter of the limb-darkening profile")
            else:
                u = gamma_to_u(Gamma)
            # now we have u
        self.limb_linear_coe1 = (3 - 3*u)/(3 - u)
        self.limb_linear_coe2 = 3*u/(3 - u)

    # def limb_linear(self, r2):
    #     '''
    #     f(r) = 0, if r > 1
    #     f(r) = 3/(3-a)*( 1 - a + a * sqrt(1-r^2) ), if r<=1
    #          = self.limb_linear_coe1 + self.limb_linear_coe2 * sqrt(1-r^2)

    #     r2 = r^2

    #     to speed up the code, r2 need to be clean and let no r2>1
    #     '''
    #     return self.limb_linear_coe1 + self.limb_linear_coe2 * np.sqrt(1-r2) 


    def lightcurve(self, xscenters, yscenters, rho, uselimb = False, u = None, Gamma = None):

        '''
        magnification at a certain source position
        '''
        rho2 = rho*rho


        # 下面的方式有解析的方法
        # r2s_inside_source = findidxs(self.thetax, self.thetay , 0, 0, rho, rho2)
        # avg_raynum_in_source = len(r2s_inside_source)
        # print("ray num if no lens: numerical = {}, analytical = {}".format(avg_raynum_in_source, M_PI*rho2 * self.raynum_perarea))

        # 解析的方法
        avg_raynum_in_source = M_PI*rho2 * self.raynum_perarea


        # find out the index of those rays which inside the square around source center

        print("Calculating ray shooting lightcurve")

        if not uselimb:
            mus = []
            muslimb = []
            with tqdm(total=len(xscenters)) as pbar:
                for x,y in zip(xscenters, yscenters):
                    # leftx = x - rho
                    # lefty = x + rho #// use abs
                    # https://stackoverflow.com/questions/13869173/numpy-find-index-of-the-elements-within-range
                    

                    # idxs = np.where(np.logical_and( np.abs(self.betax - x) <= rho, np.abs(self.betay - y) <= rho))[0]
                    # bxs = self.betax[idxs]
                    # bys = self.betay[idxs]
                    # r2s = (bxs - x)*(bxs - x) + (bys - y)*(bys - y)
                    # idxs_r2_lessthan_1 = np.where( r2s <= 1 )[0]


                    r2s_inside_source = findidxs(self.betax, self.betay , x, y, rho, rho2)
                    raynumifLens = len(r2s_inside_source)

                    # print("idxs_r2_lessthan_1 :",idxs_r2_lessthan_1)

                    # print("x,y = ({:.2e},{:.2e}); {} / {} = {}".format(x, y, raynumifLens, raynumifNolens, raynumifLens/raynumifNolens ))

                    mus.append( raynumifLens/avg_raynum_in_source )
                    pbar.update(1)
            muslimb = mus
        else:
            self.limb_linear_coefficients(u = u, Gamma = Gamma)
            # 下面的方式有无解析的方法？

            # limbres = limb_linear(r2s_inside_source/rho2, self.limb_linear_coe1, self.limb_linear_coe2)
            # raynumifNolens_limb = npsum(limbres)
            # print("ray num if no lens with limb darkening:, ", raynumifNolens_limb) # 3600 ???

            mus = []
            muslimb = []
            with tqdm(total=len(xscenters)) as pbar:
                for x,y in zip(xscenters, yscenters):


                    # r2s_inside_source = findidxs(self.thetax, self.thetay , x, y, rho)
                    # raynumifLens = len(r2s_inside_source)
                    # print("ray num if no lens:, ", raynumifLens) # 3600 ???                    
                    
                    r2s_inside_source = findidxs(self.betax, self.betay , x, y, rho, rho2)

                    raynumifLens = len(r2s_inside_source)
                    mus.append( raynumifLens/avg_raynum_in_source )  


                    limbres = limb_linear(r2s_inside_source/rho2, self.limb_linear_coe1, self.limb_linear_coe2)
                    raynumifLens_limb = npsum(limbres)

                    muslimb.append( raynumifLens_limb/avg_raynum_in_source ) 
                    pbar.update(1)

        return mus, muslimb




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
    def gamma_to_u(self, gamma):
        return 3. * gamma / (2. + gamma)
    def u_to_gamma(self, u):
        return (2. * u) / (3. - u)




if __name__ == '__main__':

    # s = 0.8
    # q = 0.1
    # mlens = [1.0/(1.0+q), q/(1.0 + q)]
    # zlens = [-q * s / (1.0 + q), 0, s / (1.0 + q), 0]

    TRIL = TripleLensing.TripleLensing()

    # lens masses: m1, m2, m3
    mlens = [0.968738798957637, 0.028093425169771, 0.003167775872591]
    # #lens positions: x1, y1, x2, y2, x3, y3
    zlens = [-0.039343051506317, 0, 1.356656948493683, 0, 0.638936196010800, -0.950873946634155]

    xs = [zlens[i] for i in range(0, len(zlens)) if i%2==0 ]
    ys = [zlens[i] for i in range(0, len(zlens)) if i%2==1 ]
    # print("xs {}, ys {}".format(xs, ys))
    
    testlens = MicroLens(mlens, xs, ys)

    xlim = (-1.5, 1.5)
    ylim = (-1.5, 1.5)
    raynum = 5000

    testlens.gotstr(mlens, xs, ys, xlim, ylim, raynum, strlen = 300)


    # triple lens event parameters
    t0 = 7494.153;
    u0 = 0.021;
    tE = 74.62;
    s2 = 1.396;
    q2 = 0.029;
    alpha = 2.948; #//rad
    s3 = 1.168;
    q3 = 3.27e-3;
    psi = 5.332; #//rad
    rs = 0.01#0.22e-3 #0.005#
    salpha = np.sin(alpha)
    calpha = np.cos(alpha)
    params = [t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs]
    # source position
    ts = np.linspace(7470, 7510, 200) #
    tn = (ts - t0) / tE;
    y1s = u0 * salpha + tn * calpha;
    y2s = u0 * calpha - tn * salpha;

    Gamma = 0.51#0.51 #0.51
    limbu = testlens.gamma_to_u(Gamma)

    saveFig_and_npz = False

    lightcurve_raysht_path = datapath + "lightcurve_raysht_rs{:.1e}".format(rs) + testlens.betathetaxy_string


    # test ray shooting
    # testlens.ray_generator(xlim, ylim, raynum)
    # testlens.inverse_ray_shooting()
    # musnolimb, muslimb = testlens.lightcurve(y1s[0:5], y2s[0:5], rs, uselimb = 0, u = limbu)


    load_rayshtlkv = True
    if load_rayshtlkv:
        print("Loading raysht lightcurve ...")
        # npzfile = load(lightcurve_raysht_path + ".npz" )


        # local 5000
        # npzfile = load("./data/rayshtdata/lightcurve_raysht_MLENS09687388002809343000316778XLENS00393430513566569506389362YLENS00095087395XLIM1515YLIM1515RAYNUM5000.npz" ) # rs 1e-2

        #from dpc
        # npzfile = load("./data/rayshtdata/fromdpc/lightcurve_raysht_rs1.0e-02MLENS09687388002809343000316778XLENS00393430513566569506389362YLENS00095087395XLIM1515YLIM1515RAYNUM10000.npz" ) # rs 1e-2
        # npzfile = load( "./data/rayshtdata/fromdpc/lightcurve_raysht_rs1.0e-03MLENS09687388002809343000316778XLENS00393430513566569506389362YLENS00095087395XLIM1515YLIM1515RAYNUM10000.npz" ) # rs 1e-3

        # from dpcg
        npzfile = load( "./data/rayshtdata/fromdpc/lightcurve_raysht_rs1.0e-02MLENS09687388002809343000316778XLENS00393430513566569506389362YLENS00095087395XLIM1515YLIM1515RAYNUM40000.npz" ) # rs 1e-2
        # npzfile = load( "./data/rayshtdata/fromdpc/lightcurve_raysht_rs1.0e-03MLENS09687388002809343000316778XLENS00393430513566569506389362YLENS00095087395XLIM1515YLIM1515RAYNUM40000.npz" ) # rs 1e-3

        musnolimb = npzfile["musnolimb"]
        muslimb = npzfile["muslimb"]

        # np.savez(lightcurve_raysht_path, ts=ts[start:end], musnolimb=musnolimb, muslimb = muslimb)
    else:
        testlens.ray_generator(xlim, ylim, raynum)
        testlens.inverse_ray_shooting()
        musnolimb, muslimb = testlens.lightcurve(y1s[start:end], y2s[start:end], rs, uselimb = 1, u = limbu)

        np.savez(lightcurve_raysht_path, ts=ts[start:end], musnolimb=musnolimb, muslimb = muslimb)

    # input("finished saving {}>>>".format(lightcurve_raysht_path))

    tempstr = "MLENS{}XLENS{}YLENS{}XLIM{}YLIM{}".format(testlens.mlens, testlens.xs, testlens.ys, xlim, ylim)
    tempstr = "".join(x for x in tempstr if x.isalnum())
    RelTolLimb = 1e-4
    AbsTolLimb = 1e-5
    start, end = 0, -1
    lightcurve_triple_path = datapath + "lightcurve_triple_RelTolLimb1e-4_rs{:.1e}".format(rs) + tempstr

    load_trilkv = 1
    if load_trilkv:
        print("Loading triple lightcurve ...")
        # npzfile = load(lightcurve_triple_path + ".npz" )
        # npzfile = load("./data/rayshtdata/lightcurve_triple_RelTolLimb1e-4_rs1.0e-02MLENS09687388002809343000316778XLENS00393430513566569506389362YLENS00095087395XLIM1515YLIM1515_bakNpts200_secnum90.npz")
        npzfile = load("./data/response_figure/data/newtest_lkv_limb20201103.npz")
        # npzfile = load("./data/response_figure/data/newtest_lkv_limb20201103temp.npz")


        bimusnolimb = npzfile["bimusnolimb"]
        bimuslimb = npzfile["bimuslimb"]   

        # np.savez(lightcurve_triple_path, ts=ts[start:end], bimuslimb=bimuslimb, bimusnolimb = bimusnolimb) # temp     
    # computing light curve
    else:
        print("generating light curve by triple lensing without limb darkening...")
        # mus = TRIL.TriLightCurve(params, y1s, y2s)
        t0 = time.time()
        bimusnolimb = TRIL.TriLightCurve(mlens, zlens, y1s[start:end], y2s[start:end], rs, 180, 2, 1e-10, 1e-3)
        t1 = time.time()
        dt1 = t1 - t0
        print("\t\t dt1 = {} ...".format(dt1))

        print("generating light curve by triple lensing with limb darkening...")
        t0 = time.time()
        bimuslimb = TRIL.TriLightCurveLimb(mlens, zlens, y1s[start:end], y2s[start:end], rs, 180, 2, 1e-10, 1e-3, RelTolLimb, AbsTolLimb, limbu)
        t1 = time.time()
        dt2 = t1 - t0

        bimuslimb = array(bimuslimb)
        bimusnolimb = array(bimusnolimb)
        np.savez("./data/response_figure/data/newtest_lkv_limb20201103temp", ts=ts[start:end], bimuslimb=bimuslimb, bimusnolimb = bimusnolimb)
        # np.savez(lightcurve_triple_path, ts=ts[start:end], bimuslimb=bimuslimb, bimusnolimb = bimusnolimb)

    # from MulensModel.binarylens import BinaryLens
    # bilens = BinaryLens(mass_1=1/(1.0+q), mass_2=q/(1.0 + q), separation=s)
    # bimuslimb = getbi_light_curve(bilens, y1s[start:end], y2s[start:end], rs, u_limb_darkening = limbu, accuracy = 1e-4)
    # bimusnolimb = getbi_light_curve(bilens, y1s[start:end], y2s[start:end], rs, u_limb_darkening = 0, accuracy = 1e-4)

        print("\t\t\t time no limb: {}, time has limb {}".format(dt1, dt2))


    print("rayshooting no limb[]: ", musnolimb[:10])
    print("rayshooting limb[]: ", muslimb[:10])

    print("vbbl no limb[]: ", bimusnolimb[:10])
    print("vbbl limb[]: ", bimuslimb[:10])




    # main, gs = pltlkv(ts[start:end], musnolimb, label = "ray shooting no limb")

    fig = plt.subplots(figsize=(16,8), dpi=100)
    gs = gridspec.GridSpec(3,1,height_ratios=[5,1,1])
    plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0, wspace = 0)
    main = plt.subplot(gs[0])


    main.plot(ts[start:end], np.log10(musnolimb), color="b",linewidth = 1, label = "No limb-darkening, rayshooting")
    
    main.plot(ts[start:end], np.log10( bimusnolimb ), "-.b",linewidth = 1, label = "No limb-darkening, our code")
    # main.scatter( ts[start:end], np.log10(musnolimb) ,marker = ".", s=30, facecolors='none', edgecolors='b')

    main.plot(ts[start:end], np.log10( muslimb ), "r",linewidth = 1, label = "With limb-darkening, rayshooting")
    main.plot(ts[start:end], np.log10( bimuslimb ), "-.r",linewidth = 1, label = "With limb-darkening, our code")
    # main.scatter( ts[start:end], np.log10(muslimb) ,marker = ".", s=30, facecolors='none', edgecolors='r')


    main.set_ylabel(r"log($\mu$)", fontdict=font)
    main.set_xlabel('HJD - 2450000', fontdict=font)
    main.tick_params(axis='both', labelsize = legend_tick_size, direction = "in")  
    # main.set_xticks([])
    plt.setp( main.get_xticklabels(), visible=False)


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
        $\rho$ = {:.1e}
        Gamma = {}
        """.format(t0, u0, tE, s2, q2, s3, q3, alpha, psi, rs, Gamma)

    main.text(0.0, 0.3, msg, transform=main.transAxes, fontdict = font)

    # plt.legend(loc = 1)
    leg = main.legend(prop = font2, frameon=False, loc="upper right") #, ncol=2


    nolimbcolor = "b"
    ax_traj = plt.subplot(gs[1],sharex=main)
    ax_traj.axhline(0,color=nolimbcolor)
    # ax_traj.set_xlabel('HJD - 2450000', fontdict=font, color="k")
    ax_traj.set_ylabel('Residual',fontdict=font, color="k")
    ax_traj.xaxis.set_tick_params(labelsize=17)
    ax_traj.yaxis.set_tick_params(labelsize=17)
    ax_traj.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax_traj.plot( ts[start:end], (bimusnolimb-musnolimb)/musnolimb, c=nolimbcolor, linewidth = 2 , label = "No limb-darkening")
    ax_traj.scatter( ts[start:end], (bimusnolimb-musnolimb)/musnolimb ,marker = ".", s=10, facecolors='none', edgecolors=nolimbcolor)
    ax_traj.legend(loc = 1)
    # ax_traj.set_ylim(-5e-5, 5e-5)
    # ax_traj.set_ylim(-1e-4, 1e-4)
    # ax_traj.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
    # ax_traj.yaxis.major.formatter._useMathText = True

    def formatnum(x, pos):
        if x != 0:
            return '$%.0f$x$10^{-5}$' % (x/1e-5)
        else:
            return 0
    formatter = FuncFormatter(formatnum)
    ax_traj.yaxis.set_major_formatter(formatter)

    # ax_traj.set_xticks([])
    plt.setp( ax_traj.get_xticklabels(), visible=False)


    limbcolor = "r"
    ax_traj = plt.subplot(gs[2],sharex=main)
    ax_traj.axhline(0,color=limbcolor)
    ax_traj.set_xlabel('HJD - 2450000', fontdict=font, color="k")
    ax_traj.set_ylabel('Residual',fontdict=font, color="k")
    ax_traj.xaxis.set_tick_params(labelsize=17)
    ax_traj.yaxis.set_tick_params(labelsize=17)
    ax_traj.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax_traj.plot( ts[start:end], (bimuslimb-muslimb)/muslimb, c=limbcolor, linewidth = 2, label = "With limb-darkening")
    ax_traj.scatter( ts[start:end], (bimuslimb-muslimb)/muslimb ,marker = ".", s=10, facecolors='none', edgecolors=limbcolor)
    ax_traj.legend(loc = 1)

    ax_traj.set_ylim(-1.5e-4, 1.5e-4)
    # ax_traj.set_xticks([])

    # ax_traj.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
    # ax_traj.yaxis.major.formatter._useMathText = True
    def formatnum(x, pos):
        if x != 0:
            return '$%.0f$x$10^{-4}$' % (x/1e-4)
        else:
            return 0
    formatter = FuncFormatter(formatnum)
    ax_traj.yaxis.set_major_formatter(formatter)


    wh = "80%"
    inset_axes = iax(main, width=wh, # width = 30% of parent_bbox
                            height=wh, # height : 1 inch
                            loc=1,
                            bbox_to_anchor=(0.25,0,0.5,0.5), bbox_transform=main.transAxes,
                            ) #1 top right, 

    inset_axes.plot(ts[start:end], np.log10(musnolimb), color="b",linewidth = 1)
    
    inset_axes.plot(ts[start:end], np.log10( bimusnolimb ), "-.b",linewidth = 2)
    # inset_axes.scatter( ts[start:end], np.log10(musnolimb) ,marker = ".", s=30)

    inset_axes.plot(ts[start:end], np.log10( muslimb ), "r",linewidth = 1)
    inset_axes.plot(ts[start:end], np.log10( bimuslimb ), "-.r",linewidth = 1)
    # inset_axes.scatter(ts[start:end], np.log10(muslimb) ,marker = ".", s=30)
    inset_axes.xaxis.set_tick_params(labelsize=14)
    inset_axes.yaxis.set_tick_params(labelsize=14)
    inset_axes.tick_params(axis='both', labelsize = legend_tick_size, direction = "in")  

    # inset_axes.plot(criticalx, criticaly, '--', color='red', markersize=1)
    # inset_axes.plot(causticsx, causticsy, '-', color='red', markersize=1)
    # inset_axes.plot(XS, YS, 'k')
    # inset_axes.plot(imgXS, imgYS, '.', color=cl, markersize=1)
    # # inset_axes.plot(imgXS, imgYS, '.', color='magenta', markersize=1)
    # inset_axes.plot(falseimgXS, falseimgYS, '.', color='b', markersize=1)
    # inset_axes.plot([xy[0]for xy in z], [xy[1]for xy in z], '+', color='k', markersize=15)
    # inset_axes.tick_params(axis='both', labelsize = 12, direction="in")

    inset_axes.set_xlim(7495,7498)
    inset_axes.set_ylim(1.2,1.9)
    showx = [7495, 7496, 7497, 7498]
    inset_axes.set_xticks(showx)
    inset_axes.set_xticklabels(showx, minor=False)

    main.set_xlim(7470,7510)


    # np.savez(lightcurve_path, ts=ts[start:end], musnolimb=musnolimb, muslimb = muslimb, bimuslimb = bimuslimb, bimusnolimb = bimusnolimb)

    # plt.savefig("./data/response_figure/newtest_lkv_limb20201103"  + ".png", dpi=300)
    if saveFig_and_npz:
        plt.savefig(datapath + testlens.betathetaxy_string  + ".png", dpi=300)
    plt.show()