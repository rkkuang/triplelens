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
    raynum = 1000

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
    rs = 0.05#0.22e-3 #0.005#
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

    map_raysht_path = datapath + "magmap_raysht_rs{:.1e}raynum{}".format(rs,raynum) + testlens.betathetaxy_string

    ImgSize = 128
    xlim0 = -0.2; xlim1 = 1; ylim0 = -0.6; ylim1 = 0.6;
    dx = (xlim1 - xlim0)/(ImgSize-1)
    dy = (ylim1 - ylim0)/(ImgSize-1)
    px = xlim0
    py = ylim0

    xsCenters = []
    ysCenters = []
    px = xlim0
    for i in range(ImgSize):
        py = ylim0
        for j in range(ImgSize):
            xsCenters.append(px)
            ysCenters.append(py)
            py += dy
        px += dx



    testlens.ray_generator(xlim, ylim, raynum)
    testlens.inverse_ray_shooting()
    musnolimb, muslimb = testlens.lightcurve(xsCenters, ysCenters, rs, uselimb = 0, u = limbu)

    muRayshoot = np.array(musnolimb)
    muRayshoot = muRayshoot.reshape(ImgSize, ImgSize)

    muRayshootlimb = np.array(muslimb)
    muRayshootlimb = muRayshootlimb.reshape(ImgSize, ImgSize)


    np.savez(map_raysht_path, mapnolimb=muRayshoot, maplimb = muRayshootlimb)