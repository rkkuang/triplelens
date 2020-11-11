import TripleLensing
from utils import *
import math

# 
#set up lens system
# #lens masses: m1, m2, m3
mlens = [0.968738798957637, 0.028093425169771, 0.003167775872591]
# #lens positions: x1, y1, x2, y2, x3, y3
zlens = [-0.039343051506317, 0, 1.356656948493683, 0, 0.638936196010800, -0.950873946634155]


# #lens masses: m1, m2, m3
# mlens = [0.968738798957637, 0.028093425169771]
# #lens positions: x1, y1, x2, y2, x3, y3
# zlens = [-0.039343051506317, 0, 1.356656948493683, 0]



TRIL = TripleLensing.TripleLensing()
# TRIL.reset(mlens, zlens)

#source center
xsCenter = -0.034747426672208
ysCenter = -0.026627816352184
#source radius
rs = 0.005




##### test lens equation solving


# rs = 0.0002200000000000000078323
# xsCenter0 = 0.1452752184254198497548316
# ysCenter0 = 0.0070811503507445944238796

rs = 0.0002200000000000000078323
xsCenter0 = 0.1452752184254198497548316
ysCenter0 = 0.0070811503507445944238796
mlens = [0.968738798957637, 0.028093425169771, 0.003167775872591]
zlens = [-0.039343051506317, 0, 1.356656948493683, 0, 0.6389361960108, -0.950873946634155]

# Zlens[i] = complex(zlens[i], zlens[i + NLENS]);

z = [ [zlens[0], zlens[1]], [zlens[2], zlens[3]], [zlens[4], zlens[5]] ]

res = sol_len_equ_cpp(mlens, z, xsCenter, ysCenter, len(mlens))
nimgs = 0
for i in range(len(res)):
    dzs = checkLensEqu(mlens, z, xsCenter, ysCenter, res[i])
    # print("i = {}, solution = {}, dzs = {:.5e}".format(i, res[i], dzs))
    if dzs < 1e-5:
        nimgs += 1
print("number of solutions which dzs < 1e-5: ", nimgs)



Npnt = 300
phi = math.pi/2
dphi = math.pi*2/(Npnt)


# rs = 1e-6
# xsCenter0 = 0.1452752184254198497548316
# ysCenter0 = 0.0071611503507445944238796

rs = 1e-7 # 1e-7就不行了,本应出现 6 个 解，但是全部只有 4 个
xsCenter0 = 0.1452752184254198497548316
ysCenter0 = 0.0071621503507445944238796

for j in range(1, Npnt+1):
    phi = dphi*(j-1)
    xsCenter = xsCenter0 + rs*math.cos(phi)
    ysCenter = ysCenter0 + rs*math.sin(phi)
    res = sol_len_equ_cpp(mlens, z, xsCenter, ysCenter, len(mlens))
    nimgs = 0
    for i in range(len(res)):
        dzs = checkLensEqu(mlens, z, xsCenter, ysCenter, res[i])
        print("\t\ti = {}, solution = {}, dzs = {:.5e}".format(i, res[i], dzs))
        if dzs < 1e-5:
            nimgs += 1
    print("phi = {:.3f}, number of solutions which dzs < 1e-5: {}".format(phi*180/math.pi, nimgs))


plot_critcaus_srcimgs(mlens, zlens, xsCenter0, ysCenter0, rs)
# 要么是 6 个解 dzs < 1e-5, 要么是 4 个
plt.show()



input("check lens equation solving first <<<<<<")