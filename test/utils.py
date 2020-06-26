import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import axes_grid1
import matplotlib
import matplotlib as mpl
import sys
from mpl_toolkits.axes_grid.inset_locator import inset_axes as iax

import TripleLensing
TRIL = TripleLensing.TripleLensing()

NLENS = 3
DEGREE = NLENS**2 + 1
import math
M_PI = math.pi
EPS = 1.0e-5


VERBOSE = False
verbose = False


font = {'family': 'Times New Roman',
    'color':  'k',
    'weight': 'normal',
    'size': 17,
    }
legend_tick_size = 17

font2 = {'family': 'Times New Roman',
    'weight': 'normal',
    'size': 20,
    }

def read_lens_system_triple(fileName):
    # ../data/lens_system_triple.dat
    # t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs
    f = open(fileName, "r")
    _ = f.readline()
    line1 = f.readline()
    f.close()
    lpara = line1.split(" ")
    lpara = (float(i) for i in lpara)
    spara = "t0,u0,tE,s2,q2,alpha,s3,q3,psi,rs,xsCenter,ysCenter".split(",")
    params = {}
    for s,l in zip(spara, lpara):
        params[s] = l
    return params

def readFile(fileName, column1, column2, expected_elem_each_row=4): #column number starts from 0
    x0=[]
    y0=[]
    f = open(fileName, "r")
    for line in f:
    #skip over empty lines or lines starting with spaces or spaces+#
        tempString=line.strip()
        if (tempString[0] =='#' or tempString[0] == '') or len(line.split())!=expected_elem_each_row:
            continue
        line = line.split()
        x = np.float(line[column1])
        y = np.float(line[column2])
        x0.append(x)
        y0.append(y)
    f.close()
    return x0,  y0

def read_cpplkv(fileName):
    #reading cpp generated light curve
    # four value, time in HJD, 2nd, 3rd are the coordinate (xs, ys) of source center ,the 4th is the corresponding magnification
    # fprintf(ftrilkv, "%.15f %.15f %.15f %.15f ", t_array[j], y1_array[j], y2_array[j], mag_array[j]);
    f = open(fileName, "r")
    line = f.readline()
    cols = line.split(" ")[:-1]
    f.close()
    cols = np.array([float(i) for i in cols])
    times = cols[::4]
    xs = cols[1::4]
    ys = cols[2::4]
    mags = cols[3::4]
    return times, xs, ys, mags

def read_cpplkv_adap(fileName):
    #reading cpp generated light curve
    # two rows, the first is the time HJD, the second is the corresponding magnification
    f = open(fileName, "r")
    line1 = f.readline()
    times = line1.split(" ")[:-1]
    line2 = f.readline()
    mags = line2.split(" ")[:-1]    
    f.close()
    times = np.array([float(i) for i in times])
    mags = np.array([float(i) for i in mags])
    return times, mags

def read_cppmap(fileName):
    #reading cpp generated magnification map
    # three value, the coordinate (xs, ys) of source center ,the 3rd is the corresponding magnification
    f = open(fileName, "r")
    line = f.readline()
    f.close()
    cols = line.split(" ")[:-1]
    cols = np.array([float(i) for i in cols])
    xs = cols[::3]
    ys = cols[1::3]
    mags = cols[2::3]
    return xs, ys , mags

def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    https://stackoverflow.com/questions/34017866/arrow-on-a-line-plot-with-matplotlib
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )

def draw_circle(ax, xsCenter, ysCenter, rs, color = "b"):
    xs = np.linspace(xsCenter-rs, xsCenter+rs, 100)
    ys = np.linspace(ysCenter-rs, ysCenter+rs, 100)
    r2 = rs**2
    X = []
    Y = []
    for x in xs:
        for y in ys:
            if (x- xsCenter)**2 + (y- ysCenter)**2 <= r2:
                X.append(x)
                Y.append(y)
    ax.scatter(X,Y, c = color)

def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
    """Add a vertical color bar to an image plot."""
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    return im.axes.figure.colorbar(im, cax=cax, **kwargs)


def fmt(x, pos):
    '''https://stackoverflow.com/questions/25983218/scientific-notation-colorbar-in-matplotlib'''
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    '''
    https://stackoverflow.com/questions/43324152/python-matplotlib-colorbar-scientific-notation-base
    '''
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
             self.format = r'$\mathdefault{%s}$' % self.format


# #https://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
# def add_subplot_axes(ax,rect,axisbg='w'):
#     fig = plt.gcf()
#     box = ax.get_position()
#     width = box.width
#     height = box.height
#     inax_position  = ax.transAxes.transform(rect[0:2])
#     transFigure = fig.transFigure.inverted()
#     infig_position = transFigure.transform(inax_position)    
#     x = infig_position[0]
#     y = infig_position[1]
#     width *= rect[2]
#     height *= rect[3]  # <= Typo was here
#     subax = fig.add_axes([x,y,width,height],facecolor=axisbg)
#     x_labelsize = subax.get_xticklabels()[0].get_size()
#     y_labelsize = subax.get_yticklabels()[0].get_size()
#     x_labelsize *= rect[2]**0.5
#     y_labelsize *= rect[3]**0.5
#     subax.xaxis.set_tick_params(labelsize=x_labelsize)
#     subax.yaxis.set_tick_params(labelsize=y_labelsize)
#     return subax







def readFile(fileName, column1, column2, expected_elem_each_row=4): #column number starts from 0
    x0=[]
    y0=[]
    f = open(fileName, "r")
    for line in f:
	#skip over empty lines or lines starting with spaces or spaces+#
        tempString=line.strip()
        if (tempString[0] =='#' or tempString[0] == '') or len(line.split())!=expected_elem_each_row:
            continue
        line = line.split()
        x = np.float(line[column1])
        y = np.float(line[column2])
        x0.append(x)
        y0.append(y)
    f.close()
    return x0,  y0

def FS_vs_VBBL():
    ImgSize = 100
    x_min = -0.75
    x_max = 0.75
    y_max = 0.75
    y_min = -0.75
    print("loading mu.dat ...")
    line2list = []
    f = open('./data/fmuFS.dat',"r")
    for i in range(ImgSize):
        line2 = f.readline()
        line2list += line2.split(" ")[:-1]
    f.close()
    muFS = np.array([float(i) for i in line2list])
    muFS = muFS.reshape((ImgSize, ImgSize))  

    line2list = []
    f = open('./data/fmuVBBL.dat',"r")
    for i in range(ImgSize):
        line2 = f.readline()
        line2list += line2.split(" ")[:-1]
    f.close()
    muVBBL = np.array([float(i) for i in line2list])
    muVBBL = muVBBL.reshape((ImgSize, ImgSize))

    line2list = []
    f = open('./data/fdtFS.dat',"r")
    for i in range(ImgSize):
        line2 = f.readline()
        line2list += line2.split(" ")[:-1]
    f.close()
    dtFS = np.array([float(i) for i in line2list])
    dtFS = dtFS.reshape((ImgSize, ImgSize))

    line2list = []
    f = open('./data/fdtVBBL.dat',"r")
    for i in range(ImgSize):
        line2 = f.readline()
        line2list += line2.split(" ")[:-1]
    f.close()
    dtVBBL = np.array([float(i) for i in line2list])
    dtVBBL = dtVBBL.reshape((ImgSize, ImgSize))

    cmap = "seismic"
    fig = plt.figure(figsize=(18,7))
    plt.subplot(231)
    plt.imshow(muVBBL,cmap = cmap, extent = [x_min , x_max, y_min , y_max])
    plt.title("Mu VBBL",fontdict = font)
    plt.colorbar()

    plt.subplot(232)
    # plt.imshow((muFS-muVBBL)/muVBBL,cmap = cmap)
    plt.imshow((muFS-muVBBL)/muVBBL,cmap = cmap, extent = [x_min , x_max, y_min , y_max])
    plt.title("Rel err",fontdict = font)

    plt.colorbar()

    plt.subplot(233)
    plt.imshow(dtFS/dtVBBL,cmap = cmap, extent = [x_min , x_max, y_min , y_max])
    # plt.title("Rel err by cpp",fontdict = font)
    plt.title("Time new/VBBL",fontdict = font)
    plt.colorbar()


    print(" total time VBBL: {}, new method: {}, total time ratio: {} ".format( np.sum(dtVBBL), np.sum(dtFS), np.sum(dtFS)/np.sum(dtVBBL) ))

    timeratio = (dtFS/dtVBBL).reshape(-1,1)
    # plt.figure(figsize=(8,8))
    plt.subplot(212)
    # plt.hist(timeratio, bins=range(min(timeratio), max(timeratio) + 10, 10))
    plt.hist(timeratio, bins=np.arange(min(timeratio), max(timeratio), 0.5))
    plt.show()





def plt_lightkv():
    # plt.figure()
    f = open("mu.dat", "r")
    # fprintf(fmu, "%f %f %d %f %f", xlim0, xlim1, ImgSize, k, b);
    head = f.readline().strip().split(" ")
    head = [float(i) for i in head]
    full_line = f.readline().strip().split(" ")
    full_line = [float(i) for i in full_line]
    pxs = np.linspace(head[0], head[1], int(head[2]))
    # plt.plot(pxs, full_line)
    # plt.show()
    f.close()
    return pxs, full_line, head

def plotcritcaus():
    print("Tracks function")
    # fig = plt.figure()
    # plt.clf()
    # ax = fig.add_subplot(111)      
    fig, ax = plt.subplots(figsize=(8,8))
    title="Critical lines (blue), caustics (red)\n source and lens (black), images (green) of triple lenses system\n yellow are the images of source centre"
    plt.suptitle(title)
    x, y = readFile("./data/caustics.dat", 0, 1, expected_elem_each_row=2)
    ax.plot(x, y, 'o', color='red', markersize=1)
    # plt.show()
    # input()

    x, y = readFile("./data/critical_curves.dat", 0, 1, expected_elem_each_row=2)
    ax.plot(x, y, 'o', color='blue', markersize=1)

    lensx, lensy = readFile("./data/lens_system.dat", 1, 2, expected_elem_each_row=3)
    lensm, _ = readFile("./data/lens_system.dat", 0, 2, expected_elem_each_row=3)
    # print()
    # print(lensx, lensy)
    for i in range(len(lensm)-1):
        plt.plot(lensx[1+i], lensy[1+i], 'o', color='k', markersize=5*lensm[i+1])
        plt.text(lensx[1+i], lensy[1+i],"lens{}".format(i+1))

    print("plot pureImgPoints.dat")

    # xc, yc = readFile("./data/pureImgPoints_centre.dat", 0, 1, expected_elem_each_row=2)
    # ax.plot(xc, yc, 'x', color='y', markersize=3)
    # print("src centre image number: ", len(xc))

    x, y = readFile("./data/pureImgPoints.dat", 0, 1, expected_elem_each_row=2)
    ax.plot(x, y, '.', color='green', markersize=3)

    f = open("./data/lens_system.dat", "r")
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
    x=xs+rs*np.cos(phi)
    y=ys+rs*np.sin(phi)

    ax.plot(x, y)
    plt.text(xs,ys,"src")
    plt.axis('equal')
    # plt.xlim(-1.5,1.5)
    # plt.ylim(-1.5,1.5)
    plt.show()

#python interface using TripleLensing module 


# def plot_critcaus_srcimgs_v3(mlens, z, xsCenter, ysCenter, rs, NLENS,nphi=2000, NPS=1000,secnum = 24, basenum = 50, scale = 10, pltfalseimg = False):
#     # non-uniform phis
#     critical , caustics = get_crit_caus(mlens, z, NLENS, NPS = NPS)
#     causticsx = np.array([xy[0]for xy in caustics])
#     causticsy = np.array([xy[1]for xy in caustics])
#     Phis = getphis_v3(mlens, z, xsCenter, ysCenter, rs, nphi, causticsx, causticsy,secnum=secnum, basenum = basenum, scale = scale)
#     Phis = Phis[0]
#     print("uniphi in image plane, len PHI:", len(Phis))
#     imgXS, imgYS, XS, YS, falseimgXS, falseimgYS = get_allimgs_v2(mlens, z, xsCenter, ysCenter, rs, NLENS, Phis)
#     fig, ax = plt.subplots(figsize=(6,6))

#     # critical , caustics = get_crit_caus(mlens, z, NLENS, NPS = NPS)
#     ax.plot([xy[0]for xy in critical], [xy[1]for xy in critical], 'o', color='blue', markersize=0.5)
#     ax.plot([xy[0]for xy in z], [xy[1]for xy in z], 'x', color='k', markersize=5)
#     # ax.plot([xy[0]for xy in caustics], [xy[1]for xy in caustics], 'o', color='red', markersize=1)
#     ax.plot(causticsx, causticsy, 'o', color='red', markersize=1)
#     ax.plot(XS, YS, '.', color="orange" , markersize=1)
#     ax.plot(imgXS, imgYS, '.', color="green", markersize=1)
#     ax.text(0.8,1,"rs:{:.0e}\nnphi:{}".format(rs,len(Phis)),fontdict = font)
#     for xy,m,i in zip(z, mlens, range(NLENS)):
#         ax.text(xy[0],xy[1],"m{}@{:.1e}".format(i,m,fontdict = font))

def plot_critcaus_srcimgs(mlens, zlens, xsCenter, ysCenter, rs,nphi=2000, NPS=4000,secnum = 360, basenum = 5, scale = 10, pltfalseimg = True, title = False, srctext = False, xy = (0.5, 0.9), inst = False, xylim = (-0.1,0.1,-0.1,0.1), wh = "32%"):
    # non-uniform phis
    z = [ [zlens[0], zlens[1]], [zlens[2], zlens[3]], [zlens[4], zlens[5]] ]
    nlens = len(mlens)
    fig, ax = plt.subplots(figsize=(8,8))
    # plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.8, left = 0.2, hspace = 0, wspace = 0)
    ax.tick_params(axis='both', labelsize = legend_tick_size, direction="in")

    critical , caustics = get_crit_caus(mlens, z, nlens, NPS = NPS)
    causticsx = np.array([xy[0]for xy in caustics])
    causticsy = np.array([xy[1]for xy in caustics])
    criticalx = [xy[0]for xy in critical]
    criticaly = [xy[1]for xy in critical]
    ax.plot(causticsx, causticsy, '-', color='red', markersize=1)
    ax.plot(criticalx, criticaly, '--', color='r', markersize=1)

    Phis = getphis_v3(mlens, z, xsCenter, ysCenter, rs, nphi, causticsx, causticsy,secnum=secnum, basenum = basenum, scale = scale)
    Phis = Phis[0]
    # print("uniphi in image plane, len PHI:", len(Phis))
    imgXS, imgYS, XS, YS, falseimgXS, falseimgYS = get_allimgs_v2(mlens, z, xsCenter, ysCenter, rs, nlens, Phis)

    ax.plot([xy[0]for xy in z], [xy[1]for xy in z], '+', color='k', markersize=5)
    ax.plot(XS, YS, '.', color="k" , markersize=1)
    ax.plot(imgXS, imgYS, '.', color='magenta', markersize=1)
    # ax.plot(imgXS, imgYS, '.', color='cyan', markersize=1)
    if pltfalseimg:
        ax.plot(falseimgXS, falseimgYS, '.', color='b', markersize=1)

    if srctext:
        for xy,m,i in zip(z, mlens, range(nlens)):
            ax.text(xy[0],xy[1],"m{}@{:.1e}".format(i+1,m,fontdict = font))
    ax.annotate('(${:.1e}$, ${:.1e}$)'.format(xsCenter,ysCenter), xy=xy, xycoords='axes fraction', fontsize=17,
                horizontalalignment='right', verticalalignment='bottom')
    plt.axis('equal')
    # plt.xlim(-1.1,1.6)
    # plt.ylim(-1.3,1.4)
    tit = """
    The three plus signs: the lens positions; The black circle: finite source
    The red solid and dashed curve: the caustics and critical curves.
    The cyan curves: true image trajectories; the blue curves: false image trajectories
    """
    if inst:
        inset_axes = iax(ax,
                            width=wh, # width = 30% of parent_bbox
                            height=wh, # height : 1 inch
                            loc=1) #1 top right, 
        inset_axes.plot(criticalx, criticaly, '--', color='red', markersize=1)
        inset_axes.plot(causticsx, causticsy, '-', color='red', markersize=1)
        inset_axes.plot(XS, YS, 'k')
        inset_axes.plot(imgXS, imgYS, '.', color='magenta', markersize=1)
        inset_axes.plot(falseimgXS, falseimgYS, '.', color='b', markersize=1)
        inset_axes.plot([xy[0]for xy in z], [xy[1]for xy in z], '+', color='k', markersize=5)
        inset_axes.tick_params(axis='both', labelsize = 12, direction="in")
        # inset_axes.set_xlim(-0.06,0.04)
        # inset_axes.set_ylim(-0.05,0.05)
        inset_axes.set_xlim(xylim[0],xylim[1])
        inset_axes.set_ylim(xylim[2],xylim[3])

    if title:
        plt.suptitle(tit, fontdict = font2)
    if 0:
        plt.savefig("./data/topo_{}_{}_rs{}.png".format(xsCenter,ysCenter,rs), dpi=300)

def get_crit_caus(mlens, z, NLENS, NPS = 200):
    zlens = [i[0] for i in z] + [i[1] for i in z]
    resxy = TRIL.outputCriticalTriple_list(mlens, zlens, NLENS, NPS)
    # critcaus_xy.sort()
    # print(critcaus_xy[-2:]) # NPS = 200, 1200, 1200

    #    // allxys:
    #    // num_critical_points
    #    // then critx_i, crity_i
    #    // numb_caustics_points
    #    //then causticsx_i, causticsy_i
    #    // [2, crx1,cry1,crx2,cry2, 1, cax1, cay1]
    # //idx [0, 1,   2,   3,   4,    5, 6,    7]

    # print("len resxy", len(resxy))
    # print(resxy[0], max(resxy[1:]))

    critical = []
    caustics = []
    numcrit = int(resxy[0])
    for i in range(numcrit):
        critical.append( [ resxy[2*i+1], resxy[2*i+2] ] )
    offset = 2*numcrit + 1
    numcaus = int( resxy[offset] )
    for i in range(numcaus):
        caustics.append([resxy[offset+2*i+1],resxy[offset+2*i+2]])
    return critical, caustics


def sol_len_equ_cpp(mlens, z, xsCenter, ysCenter, NLENS):
    zlens = [i[0] for i in z] + [i[1] for i in z]
    
    resxy = TRIL.solv_lens_equation(mlens, zlens, xsCenter, ysCenter, NLENS)
    # print(res)
    # res = [ [0]*2 ] * DEGREE # do not write like this, it will shallow copy
    res = [ [0, 0] for i in range(DEGREE)]
    # print(resxy)
    for i in range(DEGREE):
        # print(i, resxy[i], resxy[i+DEGREE])
        res[i][0] = resxy[i]
        res[i][1] = resxy[i+DEGREE]
        # print(res[i])
    # print(res)
    return res


def getphis_v3(mlens, z, xsCenter, ysCenter, rs, nphi, causticsx, causticsy,secnum = 24, basenum = 50, scale = 10, psf = [0.7,1,0.7]):
    # get phis non-uniformly 
    distype = None # away, almost, crossing
    dis = (causticsx - xsCenter)**2 + (causticsy - ysCenter)**2
    mindis = np.min(dis)
    mindixidx = np.argmin(dis)
    if mindis >= 4*rs**2:
        if VERBOSE or verbose:
            print("away from caustics")
        distype = "away"
        nphi = max(nphi, 32)
        mindis_ang = myatan(causticsx[mindixidx] - xsCenter, causticsy[mindixidx] - ysCenter)
        phi0 = mindis_ang+M_PI
        PHI = np.linspace(phi0,2.0*M_PI+phi0,nphi,endpoint = True)
        return PHI, distype, phi0, nphi
    elif mindis >= rs**2:
        if VERBOSE or verbose:
            print("almost caustic crossing")
        distype = "almost"
        # nphi = max(nphi, 256)
        nphi = max(nphi, 32)
        qnphi = int(nphi/4+0.5)
        scale = 1
        mindis_ang = myatan(causticsx[mindixidx] - xsCenter, causticsy[mindixidx] - ysCenter)
        PHI1 = np.linspace(mindis_ang-M_PI,mindis_ang - M_PI/3,scale*qnphi,endpoint = False)
        PHI2 = np.linspace(mindis_ang-M_PI/3,mindis_ang + M_PI/3,4*scale*qnphi,endpoint = False)
        # PHI3 = np.linspace(mindis_ang+M_PI/3,mindis_ang + M_PI,scale*qnphi,endpoint = False)
        PHI3 = np.linspace(mindis_ang+M_PI/3,mindis_ang + M_PI,scale*qnphi,endpoint = True) # must be true!!!!
        return np.concatenate([PHI1, PHI2, PHI3]), distype, mindis_ang, qnphi
    else:
        if VERBOSE or verbose:
            print("caustic crossing")
        # basenum = 50
        distype = "crossing"

        PHI = np.linspace(0,2.0*M_PI,secnum,endpoint = False)
        # print("PHI, ",PHI/np.pi*180)
        XS = xsCenter + rs*np.cos(PHI)
        YS = ysCenter + rs*np.sin(PHI)
        mus = []
        for xs, ys in zip(XS, YS):
            mus.append( muPoint(mlens, z, xs, ys, NLENS))
        npmus = np.array(mus)

        # maxmuidx = np.argmax(npmus)
        # maxmuidx = np.argmin(npmus)
        ratiomin = 1
        maxmuidx = (ratiomin * np.argmin(npmus) + (1-ratiomin)*np.argmax(npmus) ).astype(np.int)
        if VERBOSE:
            print("maxmuidx, ",maxmuidx)

        # psf = np.array([0.7, 1, 0.7])*0.2
        npmus = np.convolve(npmus/np.min(npmus),psf, 'same')        # 卷积

        npmus = npmus.astype(np.int)+1

        # print("len npmus: ",len(npmus), npmus)
        # input()

        # npmus = (npmus/np.min(npmus)).astype(np.int)
        # npmus = np.convolve(npmus,[1,1,1], 'same').astype(np.int)
        # print("len npmus: ",len(npmus), npmus)
        # input()

        # dphi = 2*M_PI/secnum
        # PHIfrom = PHI - dphi/2
        # PHIto = PHI + dphi/2
        # PHI1 = np.array([])
        # for i in range(1, secnum):
        #         PHI2 = np.linspace(PHIfrom[i],PHIto[i], npmus[i]*basenum,endpoint = False)
        #         PHI1 = np.concatenate([PHI1, PHI2])
        # PHI2 = np.linspace(PHIto[-1],PHIto[-1]+dphi, npmus[0]*basenum,endpoint = True)
        # PHI1 = np.concatenate([PHI1, PHI2])
        # # print("len PHI1: ",len(PHI1), PHI1/np.pi*180)
        # # input()
        # # PHI1.sort()
        # return PHI1, distype, npmus, secnum, PHI,  basenum

        secnumlist = list(range(maxmuidx, secnum))+list(range(maxmuidx))
        # print("secnumlist: ",secnumlist)
        offset = (list(range(secnum))<maxmuidx)*2*M_PI
        # print("offset, ", offset)


        dphi = 2*M_PI/secnum/2
        for i in secnumlist:
            # offset = 0
            # if i<maxmuidx:
            #     offset = 2*M_PI
            if i == secnumlist[0]:
                # PHI1 = np.linspace(offset[i]+2*M_PI*i/secnum,offset[i]+2*M_PI*(i+1)/secnum, npmus[i]*basenum,endpoint = False)
                PHI1 = np.linspace(offset[i]+PHI[i]-dphi,offset[i]+PHI[i]+dphi, npmus[i]*basenum,endpoint = False)
            elif i==secnumlist[-1]:
                PHI2 = np.linspace(offset[i]+PHI[i]-dphi,offset[i]+PHI[i]+dphi, npmus[i]*basenum,endpoint = True)
                # PHI2 = np.linspace(offset[i]+2*M_PI*i/secnum,offset[i]+2*M_PI*(i+1)/secnum, npmus[i]*basenum,endpoint = True)
                PHI1 = np.concatenate([PHI1, PHI2])
            else:
                # PHI2 = np.linspace(offset[i]+2*M_PI*i/secnum,offset[i]+2*M_PI*(i+1)/secnum, npmus[i]*basenum,endpoint = False)
                PHI2 = np.linspace(offset[i]+PHI[i]-dphi,offset[i]+PHI[i]+dphi, npmus[i]*basenum,endpoint = False)
                PHI1 = np.concatenate([PHI1, PHI2])
        if VERBOSE:
            print("len PHI1: ",len(PHI1))
        # for phi in PHI1:
        #     print(phi/np.pi*180)
        # input()
        PHI1 -= 2*M_PI
        #add on 2020.03.15
        if PHI1[0]+2*M_PI > PHI1[-1]:
            # PHI1 = np.concatenate([PHI1, np.array([PHI1[0]+2*M_PI])])
            PHI1 = np.concatenate([PHI1[:-1], np.linspace(PHI1[-1], PHI1[0]+2*M_PI, 4, endpoint=True )])
        return PHI1, distype, npmus, secnum, PHI, ratiomin, basenum

def get_allimgs_v2(mlens, z, xsCenter, ysCenter, rs, NLENS, Phis):
    # for non-uniform phis
    # allimgs = []
    XS = []
    YS = []
    imgXS = []
    imgYS = []
    falseimgXS = []
    falseimgYS = []
    # nphi = 1000
    # phi0 = 1.5
    # dphi = 2*math.pi/(nphi-1)
    # phi = phi0 - dphi
    for phi in Phis:
        # phi += dphi
        xs = xsCenter + rs*math.cos(phi)
        ys = ysCenter + rs*math.sin(phi)
        XS.append(xs)
        YS.append(ys)
        res = sol_len_equ_cpp(mlens, z, xs, ys, NLENS)
        for i in range(DEGREE):
            flag = trueSolution(mlens, z, xs, ys, res[i], cal_ang = False)
            if flag[0]:
                # allimgs.append(res[i])
                imgXS.append(res[i][0])
                imgYS.append(res[i][1])
            else:
                falseimgXS.append(res[i][0])
                falseimgYS.append(res[i][1])
    # for res in allimgs:
    #     imgXS.append(res[0])
    #     imgYS.append(res[1])   
    return imgXS, imgYS, XS, YS, falseimgXS, falseimgYS

def myatan(x,y):
    # print("inside myatan, x, y=",x,y)
    # let angle to be 0~2*M_PI
    if x >= 0 and y == 0:
        return 0
    if x == 0 and y > 0:
        return M_PI/2
    if x == 0 and y < 0:
        return 3*M_PI/2
    if y == 0 and x < 0:
        return M_PI
    ang = np.arctan(y/x)
    if ang > 0:
        if y > 0:
            return ang
        else:
            return M_PI + ang
    else:
        if y < 0:
            return 2*M_PI + ang
        else:
            return M_PI + ang

def trueSolution(mlens, zlens_list, xs, ys, z, cal_ang = True):
    # z = [x, y]
    zlens = []
    for i in range(len(zlens_list)):
        zlens.append( complex(zlens_list[i][0], zlens_list[i][1]) )
    z = complex(z[0], z[1])

    flag = 0
    Jxx = 1.0
    Jyy = 0.0
    Jxy = 0.0
    sum2 = 0.0
    sq = 0.0
    TINY = 1.0e-20
    # // mark imaginary soltion;
    lambda1, lambda2, thetaJ = 0,0,0
    mu = -1.0e10
    zs = complex(xs, ys)

    # dzs = zs - ( z - mlens[0] / conj(z - zlens[0]) - mlens[1] / conj(z - zlens[1]) - mlens[2] / conj(z - zlens[2]) )
    dzs = zs - z
    for i in range(len(mlens)):
        dzs += (mlens[i] / conj(z - zlens[i]))
    # // check the difference between zs(real source position) and the position computed from lens equation using solved image position z
    if abs(dzs) < EPS:
        flag = 1
        x = z.real
        y = z.imag
        for i in range(NLENS):
            dx = x - zlens[i].real
            dy = y - zlens[i].imag
            r2_1 = dx * dx + dy * dy + TINY#; // avoid zero in the denominator
            Jxx += mlens[i] * (dx * dx - dy * dy) / (r2_1 * r2_1)
            Jxy += 2.0 * mlens[i] * dx * dy / (r2_1 * r2_1)
        
        # //analytical results for the other components of the Jacobian
        Jyy = 2.0 - Jxx;
        Jyx = Jxy;
        mu = 1.0 / (Jxx * Jyy - Jxy * Jyx);
        # print("mu", mu)

        if cal_ang:
            sum2 = (Jxx + Jyy) / 2.0;
            sq = math.sqrt(sum2 * sum2 - 1.0 / mu);

            lambda1 = sum2 + sq;
            lambda2 = sum2 - sq;
            # lambda1 * lambda2 = 1/mu

            thetaJ = 0.5 * math.atan(2.0 * Jxy / (Jyy - Jxx + TINY))#// avoid zero in the denomintor
            if (thetaJ < 0.0):
                thetaJ += math.pi / 2.0
    # print([flag, mu, lambda1, lambda2, thetaJ])
    return [flag, mu, lambda1, lambda2, thetaJ]
def muPoint(mlens, z, xsCenter, ysCenter, NLENS):
    res = sol_len_equ_cpp(mlens, z, xsCenter, ysCenter, NLENS)
    mu = 0
    for i in range(DEGREE):
        flag = trueSolution(mlens, z, xsCenter, ysCenter, res[i], cal_ang = False)
        # mu = 0 #居然放在这里了，bug！
        if flag[0]:
            # print(flag)
            mu += abs(flag[1])
    return mu

def conj(z):
    return complex(z.real, -z.imag)
