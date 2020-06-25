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
plt.show()