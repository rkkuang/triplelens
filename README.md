# Microlensing light-curve computation for finite source star lensed by triple lens system

This method is based on the image track connecting method, which is different from traditional inverse ray shooting routine.

### cpp software

After compile with "make" command, the test program can be run as:

```shell
./bin/testtriple
```

Here we show how the closed image tracks are connected:

<img align="right" src="./doc/connected_track_eg.gif" width="350" height="350"><img align="left" src="./doc/connected_track_eg2.gif" width="350" height="350"> 















```


















```

### python interface

first run:

```shell
sudo python3 setup.py install
```

to build the python module "TripleLensing"

then, run the following:

```python
#test/testpymodule.py
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
```

the output will be:

```shell
finite source magnification:  16.98207968801614
```

We also offer an interface to visualize the lens system configuration and critical curves and caustics:

```python
#simply call:
from utils import *
plot_critcaus_srcimgs(mlens, zlens, xsCenter, ysCenter, rs)
plt.show()
```

The resultant image is as follows:

![](./doc/critcaus.png)

We can also generate a light curve:

```python
from utils import *
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
ts = np.linspace(7470, 7510, 1000)
tn = (ts - t0) / tE;
y1s = u0 * salpha + tn * calpha;
y2s = u0 * calpha - tn * salpha;
# computing light curve
print("generating light curve ...")
mus = TRIL.TriLightCurve(params, y1s, y2s)
pltlkv(ts, mus, params)
plt.show()
```

Here we show the resultant light curve:

![](./doc/lkv.png)

