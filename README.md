# Microlensing light-curve computation for finite source star lensed by triple lens system

This method is based on the image track connecting method, which is different from traditional inverse ray shooting routine.

### cpp software

After compile with "make" command, the test program can be run as:

```
./bin/testtriple
```

Here we show how the closed image tracks are connected:

<img align="right" src="./doc/connected_track_eg.gif" width="350" height="350"><img align="left" src="./doc/connected_track_eg2.gif" width="350" height="350"> 















```










```

### python interface

first run:

```
sudo python3 setup.py install
```

to build the python module "TripleLensing"

then, run the following:

```
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

```
finite source magnification:  16.98207968801614
```


