import TripleLensing

TRIL = TripleLensing.TripleLensing()

#set up lens system
#lens masses: m1, m2, m3
mlens = [0.968738798957637, 0.028093425169771, 0.003167775872591]
#lens positions: x1, y1, x2, y2, x3, y3
zlens = [-0.039343051506317, 0, 1.356656948493683, 0, 0.638936196010800, -0.950873946634155]

#source center
xsCenter = -0.034747426672208
ysCenter = -0.026627816352184
#source radius
rs = 0.005

#compute the magnification:
mu = TRIL.tripleFS2python(mlens, zlens, xsCenter, ysCenter, rs)
print("finite source magnification: ", mu)


#plot lens system topology, critical curves and caustics
from utils import *
plot_critcaus_srcimgs(mlens, zlens, xsCenter, ysCenter, rs)


# test light curve generation:
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
rs = 0.22e-3;
salpha = np.sin(alpha)
calpha = np.cos(alpha)
params = [t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs]
# source position
ts = np.linspace(7470, 7510, 500)
tn = (ts - t0) / tE;
y1s = u0 * salpha + tn * calpha;
y2s = u0 * calpha - tn * salpha;
# computing light curve
print("generating light curve ...")
mus = TRIL.TriLightCurve(params, y1s, y2s)
pltlkv(ts, mus, params)
plt.show()


# import TripleLensing

# TRIL = TripleLensing.TripleLensing()

# #set up lens system
# #lens masses: m1, m2, m3
# mlens = [0.968738798957637, 0.028093425169771, 0.003167775872591]
# #lens positions: x1, y1, x2, y2, x3, y3
# zlens = [-0.039343051506317, 0, 1.356656948493683, 0, 0.638936196010800, -0.950873946634155]

# #source center
# xsCenter = 0
# ysCenter = 0
# #source radius
# rs = 0.5

# #compute the magnification:
# mu = TRIL.tripleFS2python(mlens, zlens, xsCenter, ysCenter, rs)
# print("finite source magnification: ", mu)


# #plot lens system topology, critical curves and caustics
# from utils import *
# plot_critcaus_srcimgs(mlens, zlens, xsCenter, ysCenter, rs, xy = (0.3, 0.9), inst = 1, xylim = (1.3,1.45,-0.075,0.075), wh = "22%")
# if 1:
#     plt.savefig("./data/topo_{}_{}_rs{}.png".format(xsCenter,ysCenter,rs), dpi=300)
# plt.show()