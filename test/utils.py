import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import axes_grid1
import matplotlib

VERBOSE = False

font = {'family' : 'serif',
        # 'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
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

def tracks():
    print("Tracks function")
    # fig = plt.figure()
    # plt.clf()
    # ax = fig.add_subplot(111)      
    fig, ax = plt.subplots(figsize=(12,7))
    title="Critical lines (blue), caustics (red)\n source and lens (black), images (green) of triple lenses system\n yellow are the images of source centre"
    plt.suptitle(title)
    x, y = readFile("./data/caustics.dat", 0, 1, expected_elem_each_row=2)
    ax.plot(x, y, 'o', color='red', markersize=1)
    # plt.show()
    # input()

    x, y = readFile("./data/critical_curves.dat", 0, 1, expected_elem_each_row=2)
    ax.plot(x, y, 'o', color='blue', markersize=1)

    # x, y = readFile("images.dat", 0, 1)
    # # plt.plot(x, y, color='green', markersize=1)
    # plt.plot(x, y, '+', color='green', markersize=1 )



    # f = open("images.dat", "r")
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

    lensx, lensy = readFile("./data/lens_system.dat", 1, 2, expected_elem_each_row=3)
    lensm, _ = readFile("./data/lens_system.dat", 0, 2, expected_elem_each_row=3)
    # print(lensx, lensy)
    for i in range(3):
        plt.plot(lensx[1+i], lensy[1+i], 'o', color='k', markersize=5*lensm[i+1])
        plt.text(lensx[1+i], lensy[1+i],"lens{}".format(i+1))


    try:
        f = open("./data/images.dat", "r")
        #total image tracks
        _ = f.readline()
        full_line = f.readline()
        # print("full_line:",full_line)
        line=full_line.split()
        totalTracks=np.int(line[0])

        print("total number of tracks: ", totalTracks)

        #read image each image track
        ntrack=0

        while (ntrack<totalTracks):
            full_line = f.readline()
            # print("full_line:",full_line)
            line=full_line.split()
            number_of_lines=np.int(line[0])
            #  fileImage = fopen("images.dat", "w");
            # fprintf(fileImage, "%f %f %f\n", xsCenter, ysCenter, rs);
            # fprintf(fileImage, "%d\n", imageTracks->length);

            # for(_curve *c=imageTracks->first;c;c=c->next) { //curves
            #   fprintf(fileImage, "%d\n", c->length);
            #   for(_point *p=c->first;p;p=p->next) { //point
            #     fprintf(fileImage, "%f %f %f %f \n", p->x1, p->x2, p->phi, p->mu);
            #   }
            # }       
            i=0
            x0=[]
            y0=[]
            while (i<number_of_lines):
                full_line = f.readline()
                # print("full_line:",full_line)
                line=full_line.split()    
                
                x = np.float(line[0])
                y = np.float(line[1])
                
                x0.append(x)
                y0.append(y)
                i += 1
                
            # plt.plot(x0, y0, color='magenta')
            ax.plot(x0, y0, '+', color='green', markersize=3)
            ntrack +=1
        f.close()
    except:
        print("plot pureImgPoints.dat")

        xc, yc = readFile("./data/pureImgPoints_centre.dat", 0, 1, expected_elem_each_row=2)
        ax.plot(xc, yc, 'x', color='y', markersize=3)
        print("src centre image number: ", len(xc))

        x, y = readFile("./data/pureImgPoints.dat", 0, 1, expected_elem_each_row=2)
        ax.plot(x, y, '+', color='green', markersize=3)

        plt.figure()
        x = np.array(x)
        y = np.array(y)
        Clusters = cluster(x,y)
        dimback = 0
        idx = 0
        for cl in Clusters:
            cl = list(set(cl))
            print("number of points in track_{}: {}".format(idx,len(cl)))
            dimback += len(cl)
            tempx = x[np.array(cl)]
            tempy = y[np.array(cl)]
            plt.scatter(tempx, tempy,s=0.2)
            plt.text(tempx[0],tempy[0],"track_{}".format(idx))
            idx+=1
        print("dimback: ",dimback)
        plt.suptitle("Clustered image tracks")



    plt.grid(False)
    plt.xlabel('x')
    plt.ylabel('y')

    scale = 1e6
    Totalarea = areaFunc(x , y, Clusters, xc, yc, scale = 1e6)

    Mag = Totalarea/(np.pi* (rs*scale)**2 )
    print("Mag: {}".format(Mag))


    try:
        pxs, full_line, head = plt_lightkv()
        ax.plot(pxs, head[3]*pxs+head[4], 'k-')

        # ax2 = fig.add_axes([.54, .44, .4, .4], aspect=1)
        # ax2.plot(pxs, full_line, 'k-')
        # ax2.set_xlim(0,1000)
        # # ax2.set_ylim(-1,1)

        # plt.figure()
        ax2 = fig.add_axes([.58, .5, .4, .4], aspect=1,frameon=False)
        ax2.plot(pxs, full_line, 'k-')
        # ax2.set_xlim(0,1000)
        # ax2.set_ylim(-1,1)
        ax2.set_aspect(1./ax2.get_data_ratio())

        
        
        
    except:
        pass  
    
    plt.show()

def cluster(x,y):

    thresholdvalue = 0.01

    

    x = np.array(x) 
    y = np.array(y) 
    dim = x.shape[0] #3676
    print("dim", dim)

    dismatrix = np.ones((dim,dim))*1e2
    for i in range(dim):
        for j in range(i):
                dismatrix[i,j] = ( (x[i]-x[j])**2 + (y[i]-y[j])**2 )**0.5
    for j in range(dim):
        dismatrix[0:j, j] = dismatrix[ j, 0:j ]

    dismatrix_bak = np.copy(dismatrix)

    Clusters = [] # save index which belongs to a cluster
    na_mark = np.ones(dim).astype(np.int)

    Clusters.append([])
    clusternum = 0
    processidx = 0
    na_mark[0] = 0
    Clusters[0].append(0)
    threshold_max = []

    checked = [0]

    thresholdvalue = np.mean(dismatrix.min(axis=0))*3
    if VERBOSE:
        print("thresholdvalue", thresholdvalue)

    threshold = [thresholdvalue]

    while np.sum(na_mark)>0:
        mindis = dismatrix[:,processidx].min()
        closestidx = np.where(dismatrix[:,processidx]==mindis)[0][0]
        checked.append(processidx)
        if closestidx in checked:
            na_mark[closestidx] = 0
            dismatrix[closestidx, processidx] = 1e2
            dismatrix[processidx,closestidx] = 1e2
            continue

        if  na_mark[closestidx] == 0:
            dismatrix[closestidx, processidx] = 1e2
            dismatrix[processidx,closestidx] = 1e2

            threshold.append(mindis)
        else:
            if ( (mindis <= max(threshold)) or ( min(dismatrix[np.array(Clusters[clusternum]),closestidx]) < max(threshold) ) ):
                dismatrix[closestidx, processidx] = 1e2
                dismatrix[processidx,closestidx] = 1e2
                Clusters[clusternum].append(closestidx)
                threshold.append(mindis)
            else:
                # print("starting a new image track")
                Clusters.append([])
                clusternum += 1
                Clusters[clusternum].append(closestidx)
                threshold_max.append(  2 * max(threshold))
                threshold = [thresholdvalue]
        processidx = closestidx
        na_mark[closestidx] = 0

    lenCluster = [len(i) for i in Clusters]
    lensortidx = sorted(range(len(lenCluster)), key=lambda k: lenCluster[k])
    threshold_max.append(  2 * max(threshold))
    # if len(threshold_max)<len(lenCluster):
        # threshold_max.append(threshold_max[-1])
    if VERBOSE:
        print("threshold_max: ",threshold_max)
        print("len(Clusters): before loop ", len(Clusters))
    smallclusteridx = None
    mark = 0
    currentlen = len(lensortidx)
    i=0
    step = len(threshold_max)
    while i < step:
        smallclusteridx = lensortidx[mark]
        if VERBOSE:
            print('\nstep:{}, lenCluster:{}'.format(i, lenCluster))
            print("lensortidx, smallclusteridx: ",lensortidx, smallclusteridx)
        if mark == len(lensortidx)-1:
            break

        for pnt in Clusters[smallclusteridx]:
            mindis = 1e2
            mindisindx = None
            for largeclusteridx in lensortidx[mark:]:
                if largeclusteridx != smallclusteridx:
                    tempmin = min(dismatrix_bak[np.array(Clusters[largeclusteridx]), pnt ])
                    if mindis > tempmin:
                        mindis = tempmin
                        mindisindx = largeclusteridx
            # if a image track contains too little points, e.g. just 1, it must be merged to other track
            # if mindis < threshold_max[mindisindx] or lenCluster[smallclusteridx]<=10:
            if mindis < threshold_max[mindisindx] or ((lenCluster[smallclusteridx]/lenCluster[lensortidx[-1]]<=0.6) and len(lenCluster)>4):
                largeclusteridx = mindisindx
                Clusters[largeclusteridx]+=Clusters[smallclusteridx]
                if VERBOSE:
                    print("merging smallcluster{} to largecluster{}: ", smallclusteridx, largeclusteridx)
                del Clusters[smallclusteridx]
                del threshold_max[smallclusteridx]
                # del lenCluster[smallclusteridx]
                lenCluster = [len(i) for i in Clusters]
                lensortidx = sorted(range(len(lenCluster)), key=lambda k: lenCluster[k])
                break
        
        if currentlen == len(lensortidx):
            # no change for this image, go to the next one
            mark += 1
            if VERBOSE:
                print("mark+=1, mark=", mark)
        else:
            currentlen = len(lenCluster)
        i+=1
    
    #drop tracks which only contains little point

    # NLENS = 3
    # print(lenCluster)
    # if not (len(Clusters) - NLENS - 1)%2 == 0:
    #     Clusters[lensortidx[1]]+=Clusters[lensortidx[0]]
    #     del Clusters[lensortidx[0]]
    return Clusters



def areaFunc(x , y, Clusters, xc, yc, scale=1e6):
    # x, y: source edge points' images
    # Cluster: individual image Tracks
    # xc, yc: source center point's images
    lenCluster = len(Clusters)
    lenCentreImgs = len(xc)

    # compute center of each track

    trackcenters = []
    xsyslist = []
    for track in Clusters:
        # print(track)
        # print(type(track))
        # input()
        # track = list(set(track))
        xsyslist.append([ x[np.array(track)], y[np.array(track)] ] )
        trackcenters.append([  np.mean(xsyslist[-1][0]) , np.mean(xsyslist[-1][1]) ])

    xcyc_which_track = []
    for cx,cy in zip(xc,yc):
        tempdis = []
        for trackcxy in trackcenters:
            tempdis.append( (cx-trackcxy[0])**2+(cy-trackcxy[1])**2 )
        sortdis = sorted(range(len(tempdis)), key=lambda k: tempdis[k])
        xcyc_which_track.append(sortdis[0])
        print("center({},{}), imageTracks_{}".format(cx,cy,sortdis[0]))
    # which means [xc[i], yc[i]] corresponding to xcyc_which_track[i] track
    # and if you want to plot the center of track i, just try! to plot xc[i],yc[i]

    if lenCluster == lenCentreImgs:
        print("lenCluster == lenCentreImgs")
    else:
        print("lenCluster != lenCentreImgs")

    Totalarea = 0

    for i in range( lenCentreImgs ):
        trackidx = xcyc_which_track[i]
        data = xsyslist[trackidx]
        try:
            elliparas = fitellpise(data, plot=True, scale = scale)
            if not np.isnan(elliparas[1]):
                area = np.pi*elliparas[1]*elliparas[2]
            else:
                area = 0
            # print("area: ",area)
            Totalarea += area
        except:
            # we need to use other way to compute the area rather than ellipse fitting
            print("we need to use other way to compute the area rather than ellipse fitting for track_{}".format(trackidx))

        # Totalarea += area

    return Totalarea







def fitellpise(data, plot=False, scale = 1e6):
    '''
        Args
        ----
        data (list:list:float): list of two lists containing the x and y data of the
            ellipse. of the form [[x1, x2, ..., xi],[y1, y2, ..., yi]]
        Returns
        ------
        coef (list): list of the coefficients describing an ellipse
           [a,b,c,d,f,g] corresponding to ax**2+2bxy+cy**2+2dx+2fy+g
    '''
    # import ellipses as el
    # # import numpy as np
    # # import matplotlib.pyplot as plt
    # from matplotlib.patches import Ellipse

    # data = el.make_test_ellipse()

    lsqe = el.LSqEllipse()
    # print("type(data): ", type(data))
    # print(len(data))
    
    data = [scale*i.astype(np.float64) for i in data]

    # data = [list(i) for i in data]
    # print("\n\n\ndata[0]: ", data[0])
    # print("\n\n\ndata[1]: ", data[1])
    lsqe.fit(data)

    center, width, height, phi = lsqe.parameters()

    center = [i/scale for i in center]
    width /= scale
    height /= scale
    # phi /= scale

    # print(center, width, height, phi)

    data = [i/scale for i in data]

    if plot:        
        ax = plt.gca()
        plt.scatter(data[0], data[1],c='r',s=0.5)
        ellipse = Ellipse(xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi),
                       edgecolor='b', fc='None', lw=0.5, label='Fit', zorder = 2)
        ax.add_patch(ellipse)
    return lsqe.parameters()


def all_tracks():
    fig = plt.figure()
    plt.clf()
    plt.subplot(111)      
    x, y = readFile("caustics.dat", 0, 1, expected_elem_each_row=2)
    plt.plot(x, y, 'o', color='red', markersize=1)
    # plt.subplot(132)
    x, y = readFile("critical_curves.dat", 0, 1, expected_elem_each_row=2)
    plt.plot(x, y, 'o', color='blue', markersize=1)
    # plt.subplot(133)
    # x, y = readFile("allimages.dat", 0, 1)
    x, y = readFile("images.dat", 0, 1)
    plt.plot(x, y, '+', color='black', markersize=1 )

    #   # plot lenses:
  # fileImage1 = fopen("lens_system.dat", "w");
  # fprintf(fileImage1, "%f %f %f\n", xsCenter, ysCenter, rs);
  # // fprintf(fileImage1, "%f %f %f\n", mlens[0], mlens[1], mlens[2]);
  # // fprintf(fileImage1, "%f %f %f\n", (zlens[0].re), (zlens[1].re), zlens[2].re);
  # // fprintf(fileImage1, "%f %f %f\n", (zlens[0].im), (zlens[1].im), (zlens[2].im));
  # for (int i=0; i<3; i++){
  # fprintf(fileImage1, "%f %f %f\n", mlens[i], zlens[i].re, zlens[i].im);
  # }
    lesnx, lensy = readFile("lens_system.dat", 1, 2, expected_elem_each_row=3)
    print(lesnx, lensy)
    plt.plot(lesnx[1:], lensy[1:], 'o', color='k', markersize=3)


#    x, y = readFile("test2.dat", 0, 1)
#    ax.plot(x, y, '+', color='cyan', markersize=1 )
   
    plt.grid(False)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()


def cluster_bak(x,y):
    thresholdvalue = 0.02
    threshold = [thresholdvalue]

    x = np.array(x) 
    y = np.array(y) 
    dim = x.shape[0] #3676
    print("dim", dim)

    dismatrix = np.ones((dim,dim))*1e2
    for i in range(dim):
        for j in range(i):
                dismatrix[i,j] = ( (x[i]-x[j])**2 + (y[i]-y[j])**2 )**0.5
    for j in range(dim):
        dismatrix[0:j, j] = dismatrix[ j, 0:j ]

    dismatrix_bak = np.copy(dismatrix)
    # dismatrix = np.ones((dim,dim))*1e2
    # for i in range(dim):
    #     for j in range(i):
    #         dismatrix[i,j] = ( (x[i]-x[j])**2 + (y[i]-y[j])**2 )**0.5
    # plt.figure()
    # plt.imshow(dismatrix)
    # plt.colorbar()
    Clusters = [] # save index which belongs to a cluster
    # na = np.linspace(0, dim, dim, endpoint=False).astype(np.int)
    # na_cut = np.linspace(0, dim, dim, endpoint=False).astype(np.int)

    # na = [i for i in range(dim)]
    na_mark = np.ones(dim).astype(np.int)

    # markmindis = dismatrix[:,0].min()
    # markmindis = 10e3

    Clusters.append([])
    clusternum = 0
    processidx = 0
    na_mark[0] = 0
    Clusters[0].append(0)
    threshold_max = []

    checked = [0]
    while np.sum(na_mark)>0:
        # print("sum na_mark:", np.sum(na_mark))
        mindis = dismatrix[:,processidx].min()
        closestidx = np.where(dismatrix[:,processidx]==mindis)[0][0]
        checked.append(processidx)
        if closestidx in checked:
            na_mark[closestidx] = 0
            dismatrix[closestidx, processidx] = 1e2
            dismatrix[processidx,closestidx] = 1e2
            continue

        if  na_mark[closestidx] == 0:
            # Clusters[clusternum].append(processidx)
            dismatrix[closestidx, processidx] = 1e2
            dismatrix[processidx,closestidx] = 1e2
            # processidx = closestidx

            threshold.append(mindis)
        else:
            if ( (mindis < max(threshold)) or ( min(dismatrix[np.array(Clusters[clusternum]),closestidx]) < max(threshold) ) ):
                dismatrix[closestidx, processidx] = 1e2
                dismatrix[processidx,closestidx] = 1e2
                Clusters[clusternum].append(closestidx)
                # processidx = closestidx
                threshold.append(mindis)
            else:
                # print("starting a new track, mindis = ",mindis)
                # input()
                # lenCluster.append(len(Clusters[clusternum]))
                Clusters.append([])
                clusternum += 1
                Clusters[clusternum].append(closestidx)
                # print("max(threshold): ",max(threshold))
                threshold_max.append(  2 * max(threshold))
                threshold = [thresholdvalue]
                # processidx = closestidx

        processidx = closestidx
        na_mark[closestidx] = 0

    lenCluster = [len(i) for i in Clusters]
    lensortidx = sorted(range(len(lenCluster)), key=lambda k: lenCluster[k])
    if len(threshold_max)<len(lenCluster):
        threshold_max.append(threshold_max[-1])
    print("threshold_max: ",threshold_max)
    print("len(Clusters): before loop ", len(Clusters))
    # print("Clusters[-1]",Clusters[-1])
    smallclusteridx = None
    mark = 0
    currentlen = len(lensortidx)
    # for i in range(len(threshold_max)):
    i=0
    step = len(threshold_max)
    while i < step:
    # while min(lenCluster)<lenCluster[lensortidx[cut]]:
        # for smallclusteridx in lensortidx[:-4]:
        print('\nstep:{}, lenCluster:{}'.format(i, lenCluster))

        smallclusteridx = lensortidx[mark]
        print("lensortidx, smallclusteridx: ",lensortidx, smallclusteridx)
        if mark == len(lensortidx)-1:
            break

        # mergeflag = False
        for pnt in Clusters[smallclusteridx]:
            mindis = 1e2
            # for largeclusteridx in lensortidx[cut:]:
            mindisindx = None
            for largeclusteridx in lensortidx[mark:]:
                if largeclusteridx != smallclusteridx:
                    tempmin = min(dismatrix_bak[np.array(Clusters[largeclusteridx]), pnt ])
                    if mindis > tempmin:
                        mindis = tempmin
                        mindisindx = largeclusteridx
            if mindis < threshold_max[mindisindx]:
                largeclusteridx = mindisindx
                # mergeflag = True
                Clusters[largeclusteridx]+=Clusters[smallclusteridx]
                print("merging smallcluster{} to largecluster{}: ", smallclusteridx, largeclusteridx)
                # print("len(Clusters): before del ", len(Clusters))
                del Clusters[smallclusteridx]
                # print("len(Clusters): after del ", len(Clusters))
                del threshold_max[smallclusteridx]
                # lenCluster[largeclusteridx]+=lenCluster[smallclusteridx]
                del lenCluster[smallclusteridx]
                lenCluster = [len(i) for i in Clusters]
                lensortidx = sorted(range(len(lenCluster)), key=lambda k: lenCluster[k])
                break
        
        if currentlen == len(lensortidx):
            # no change for this image, go to the next one
            mark += 1
            print("mark+=1, mark=", mark)
        else:
            currentlen = len(lenCluster)

        i+=1
        # Clusters.append([])
        # clusternum = 0
        # Clusters[0].append(0)
        # for j in range(dim-1):
        #     mindis = dismatrix[:,j].min()
        #     print("mindis: ", mindis)
        #     closestidx = np.where(dismatrix[:,j]==mindis)[0][0]
        #     # if (mindis < threshold*markmindis) and (closestidx not in Clusters[clusternum]):
        #     if ( (mindis < max(threshold)) or ( min(dismatrix[np.array(Clusters[clusternum]),closestidx]) < max(threshold) ) ):
        #         if closestidx not in Clusters[clusternum] :
        #             Clusters[clusternum].append(closestidx)
        #             threshold.append(mindis)
        #     else:
        #         print("starting a new track, mindis = ",mindis)
        #         input()
        #         Clusters.append([])
        #         clusternum += 1
        #         Clusters[clusternum].append(closestidx)
        #         # markmindis = 10e3
        #         threshold = [thresholdvalue]
        #         # threshold.append(mindis)
        # print("clusternum: ",clusternum)
        # return Clusters


        # na_cut = list(na_cut)
        # clusternum = 0
        # Clusters.append([])
        # markmindis = dismatrix[:,0].min()
        # while(len(na_cut)>0):
        #     processidx = na_cut[0]
        #     del na_cut[0]
        #     Clusters[clusternum].append(processidx)
        #     mindis = dismatrix[:,processidx].min()
        #     closestidx = np.where(dismatrix[:,processidx]==mindis)[0][0]
        #     if mindis < 2*markmindis:
        #         Clusters[clusternum].append(closestidx)
        #         del na_cut[]
    return Clusters

  
def ogle():
    fig = plt.figure()
    plt.clf()
    ax = fig.add_subplot(111)      

    import numpy as np
    import pandas as pd
    event='ogle-2008-blg-270.dat'
    event='ogle-2012-blg-0207.dat'    
    event='ogle-2012-blg-0442.dat'    

    df = pd.read_table(event, delim_whitespace=True, header=None)     


    x=df[0]
    x=x-2450000
    y=df[1]
    yerr=df[2]

    xmin=np.min(x)
    xmax=np.max(x)

    ymin=np.min(y)
    ymax=np.max(y)

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymax+0.1*(ymax-ymin), ymin-0.1*(ymax-ymin)])
    errorbar(x, y, yerr, marker='o', mfc='red', mec='green', ms=4, mew=4)