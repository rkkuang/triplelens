# implement a version that use array, rather than linked list, to represent the source/image boundaries.
# then transfer this part of python code to C++ version
# which is helpful when calculating the magnifications with GPU

# procedures
# now, you only need the lens equation solver, and assuming you already got an "array" of PHIs around the source limb --> as the initial state for the python code that we need to implement here.

# i.e., you just need to implement the following functions with array data structure
# double TripleLensing::TripleMag(double xsCenter, double ysCenter, double rs) 
# and _sols *TripleLensing::outputTracks_v2_savehalf(double xsCenter, double ysCenter, double rs, _linkedarray *PHI, _sols **prevstore)

# in first version, you can use exactly the same version as cpp
# later, you can optimize like use the same source index, rather than source position? --> how to track the source index when # of points doubled?

from numpy import cos
from numpy import sin
from utils import sol_len_equ_cpp, trueSolution, get_crit_caus
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

VERBOSE = 1
EPS = 1e-5

class pyTriple:
    def __init__(self, mlens, zlens):
        self.mlens = mlens
        self.nlens = len(mlens)
        self.zlens = zlens
        self.DEGREE = self.nlens ** 2 + 1
        self.correct_parity = -self.nlens + 1

    def get_closest_idx(self, preProv_x, preProv_y, Prov_x, Prov_y):
        # Hungarian Algorithm ?

        # for each point i (Prov_x[i], Prov_y[i]), find the index in preProv_x, preProv_y so that the point j is closest to i
        has_removed = np.zeros(self.DEGREE)
        res_idx = np.zeros(self.DEGREE).astype(int)
        #return np.arange(0, 10).astype(int)

        for i in range(self.DEGREE):
            currx = preProv_x[i]
            curry = preProv_y[i]
            mindis = 1e100
            minidx = 0
            for j in range(self.DEGREE):
                if has_removed[j] == 0:
                    dis = (currx - Prov_x[j])**2 + (curry - Prov_y[j])**2
                    if dis<mindis:
                        mindis = dis
                        minidx = j
            res_idx[i] = minidx
            has_removed[minidx] = 1
            #print(minidx, mindis)
        #print(res_idx)
        #input(">>>")
        return res_idx


    def saveimage(self, allSolutions_flag, allSolutions_mu, allSolutions_absdzs, curr_nimages, curr_total_parity, j):
        # find one true solutions that has been classified as false solution (due to the simple criterion absdzs > EPS)
        ifsave = False
        saveindex = 0
        # print("1st \t\t nimages = %d, and parity = %d which are wrong\n"%(curr_nimages, curr_total_parity))
        # for i in range(self.DEGREE):
        #     print(i, "flag, mu, absdzs = %d, %f, %e"%(allSolutions_flag[i,j], allSolutions_mu[i,j], allSolutions_absdzs[i,j]))
        

        # find the index list which mu != -1e10 and flad = 0
        # waitinglist = 
        # directly find the index
        saveindex = 0
        minabsdzs = 1e10
        for i in range(self.DEGREE):
            if allSolutions_flag[i, j] == 0 and allSolutions_mu[i, j] != -1e10:
                if allSolutions_absdzs[i, j] < minabsdzs:
                    minabsdzs = allSolutions_absdzs[i, j]
                    saveindex = i
        # after save the solution, whether the nimages and parity are correct
        if (curr_nimages + 1 - self.nlens)%2 == 1 and (curr_total_parity+self.one_parity(allSolutions_mu[saveindex,j])) == self.correct_parity:
            if VERBOSE: print("we can save the solution %d, %f, %e"%(saveindex, allSolutions_mu[saveindex,j], allSolutions_absdzs[saveindex,j]))
            ifsave = True

        # input(">>>>>")

        return ifsave, saveindex

    def one_parity(self, mu):
        if mu == -1e10:
            return 0
        return 1 if mu > 0 else -1


    def outputTracks_v2_savehalf(self, xsCenter, ysCenter, rs, PHI, prevstore = None, prevstore_x = None, prevstore_y = None, prevstore_srcx = None, prevstore_srcy = None, prevstore_mu = None, prevstore_flag = None, prevstore_absdzs = None, mindphi = None):
        # prevstore is the solution arrays from previous step
        nphi = len(PHI)

        if not mindphi:
            mindphi = np.min(PHI[1:] - PHI[:-1]) # positive
        
        mindsource = rs*2*np.sin((0.5*mindphi)) # minimum distance between two sampled points at the source boundary
        if VERBOSE: print("mindsource = %e"%mindsource) # 1.097354e-08, mindsource = 1.047547e-04 mindsource ** 2 = 1.097354e-08
        mindsource = mindsource*mindsource
        if VERBOSE: print("mindsource ** 2 = %e"%mindsource)

        if VERBOSE: print(nphi)

        allSolutions_x = np.zeros((self.DEGREE, nphi)) # 1 row = 1 mixed image track
        allSolutions_y = np.zeros((self.DEGREE, nphi))
        allSolutions_flag = np.zeros((self.DEGREE, nphi)).astype(int)
        allSolutions_srcx = np.zeros((self.DEGREE, nphi)) 
        allSolutions_srcy = np.zeros((self.DEGREE, nphi)) 
        allSolutions_absdzs = np.zeros((self.DEGREE, nphi)) 
        allSolutions_mu = np.zeros((self.DEGREE, nphi)) 

        Prov_x = np.zeros(self.DEGREE) # store the current image point of the i-th image track
        Prov_y = np.zeros(self.DEGREE)
        Prov_flag = np.zeros(self.DEGREE).astype(int)
        Prov_srcx = np.zeros(self.DEGREE)
        Prov_srcy = np.zeros(self.DEGREE)
        Prov_absdzs = np.zeros(self.DEGREE)
        Prov_mu = np.zeros(self.DEGREE)

        preProv_x = np.zeros(self.DEGREE) # store the last image point of the i-th image track
        preProv_y = np.zeros(self.DEGREE)

        attach_idx = np.zeros(self.DEGREE)

        posflagnums = np.zeros(self.DEGREE) # how many true image points in this image Track

        #onebyone_prev = False
        curr_total_parity = 0 # total parity
        curr_nimages = 0 # how many true solutions for a certain source position

        if not prevstore: # equavalent to cpp code: if (ftime)
            # solve lens equations for the first time 
            for j in range(nphi):
                phi = PHI[j]
                xs = xsCenter + rs * cos(phi)
                ys = ysCenter + rs * sin(phi)

                res = sol_len_equ_cpp(self.mlens, self.zlens, xs, ys, self.nlens, self.DEGREE)

                if j == 0:
                    curr_total_parity, curr_nimages = 0, 0
                    for i in range(self.DEGREE):
                        flag = trueSolution(self.mlens, self.zlens, xs, ys, res[i], cal_ang = False, NLENS = self.nlens) # [flag, mu, absdzs, xxx]
                        # print(res[i], allSolutions_x.shape )
                        allSolutions_x[i, j] = res[i][0]
                        allSolutions_y[i, j] = res[i][1]
                        allSolutions_flag[i, j] = flag[0]
                        allSolutions_mu[i, j] = flag[1]
                        allSolutions_srcx[i, j] = xs
                        allSolutions_srcy[i, j] = ys
                        allSolutions_absdzs[i, j] = flag[2]
                        posflagnums[i] += flag[0]

                        if flag[0]: curr_total_parity += self.one_parity(flag[1])
                        curr_nimages += flag[0]
                    if (curr_nimages - self.nlens)%2 != 1 or curr_total_parity != self.correct_parity:
                        # there might be some true images being classified as false images
                        ifsave, saveindex = self.saveimage(allSolutions_flag, allSolutions_mu, allSolutions_absdzs, curr_nimages, curr_total_parity, j)
                        if ifsave:
                            allSolutions_flag[saveindex, j] = 1
                            posflagnums[saveindex] += 1

                else:
                    # for j > 0, you need to attached the closest solution to the existing allSolutions
                    curr_total_parity, curr_nimages = 0, 0
                    for i in range(self.DEGREE):
                        flag = trueSolution(self.mlens, self.zlens, xs, ys, res[i], cal_ang = False, NLENS = self.nlens) #

                        Prov_x[i] = res[i][0]
                        Prov_y[i] = res[i][1]
                        Prov_flag[i] = flag[0]
                        Prov_mu[i] = flag[1]
                        Prov_srcx[i] = xs
                        Prov_srcy[i] = ys
                        Prov_absdzs[i] = flag[2]

                        preProv_x[i] = allSolutions_x[i, j-1] # the last point in the i-th image track
                        preProv_y[i] = allSolutions_y[i, j-1]

                        if flag[0]: curr_total_parity += self.one_parity(flag[1])
                        curr_nimages += flag[0]

                    attach_idx = self.get_closest_idx(preProv_x, preProv_y, Prov_x, Prov_y)

                    for i in range(self.DEGREE):
                        allSolutions_x[i, j] = Prov_x[attach_idx[i]]
                        allSolutions_y[i, j] = Prov_y[attach_idx[i]]
                        allSolutions_flag[i, j] = Prov_flag[attach_idx[i]]
                        allSolutions_mu[i, j] = Prov_mu[attach_idx[i]]
                        allSolutions_srcx[i, j] = Prov_srcx[attach_idx[i]]
                        allSolutions_srcy[i, j] = Prov_srcy[attach_idx[i]]
                        allSolutions_absdzs[i, j] = Prov_absdzs[attach_idx[i]]   
                        posflagnums[i] += Prov_flag[attach_idx[i]]

                    if (curr_nimages - self.nlens)%2 != 1 or curr_total_parity != self.correct_parity:
                        # there might be some true images being classified as false images
                        ifsave, saveindex = self.saveimage(allSolutions_flag, allSolutions_mu, allSolutions_absdzs, curr_nimages, curr_total_parity, j)
                        if ifsave:
                            allSolutions_flag[saveindex, j] = 1
                            posflagnums[saveindex] += 1

        else: # use previously calculated result
            scan_prev_j = 0
            for j in range(0, nphi, 2):
                for i in range(self.DEGREE):
                    allSolutions_x[i, j] = prevstore_x[i, scan_prev_j]
                    allSolutions_y[i, j] = prevstore_y[i, scan_prev_j]
                    allSolutions_flag[i, j] = prevstore_flag[i, scan_prev_j]
                    allSolutions_mu[i, j] = prevstore_mu[i, scan_prev_j]
                    allSolutions_srcx[i, j] = prevstore_srcx[i, scan_prev_j]
                    allSolutions_srcy[i, j] = prevstore_srcy[i, scan_prev_j]
                    allSolutions_absdzs[i, j] = prevstore_absdzs[i, scan_prev_j]

                    posflagnums[i] += prevstore_flag[i, scan_prev_j]
                scan_prev_j += 1

            for j in range(1, nphi, 2):
                phi = PHI[j]
                xs = xsCenter + rs * cos(phi)
                ys = ysCenter + rs * sin(phi)                 
                res = sol_len_equ_cpp(self.mlens, self.zlens, xs, ys, self.nlens, self.DEGREE)

                # for j > 0, you need to attached the closest solution to the existing allSolutions
                curr_total_parity, curr_nimages = 0, 0
                for i in range(self.DEGREE):
                    flag = trueSolution(self.mlens, self.zlens, xs, ys, res[i], cal_ang = False, NLENS = self.nlens)

                    Prov_x[i] = res[i][0]
                    Prov_y[i] = res[i][1]
                    Prov_flag[i] = flag[0]
                    Prov_mu[i] = flag[1]
                    Prov_srcx[i] = xs
                    Prov_srcy[i] = ys
                    Prov_absdzs[i] = flag[2]

                    preProv_x[i] = allSolutions_x[i, j-1] # the last point in the i-th image track
                    preProv_y[i] = allSolutions_y[i, j-1]

                    if flag[0]: curr_total_parity += self.one_parity(flag[1])
                    curr_nimages += flag[0]

                attach_idx = self.get_closest_idx(preProv_x, preProv_y, Prov_x, Prov_y)

                # attach correct index to each of the tracks in allSolutions
                for i in range(self.DEGREE):
                    allSolutions_x[i, j] = Prov_x[attach_idx[i]]
                    allSolutions_y[i, j] = Prov_y[attach_idx[i]]
                    allSolutions_flag[i, j] = Prov_flag[attach_idx[i]]
                    allSolutions_mu[i, j] = Prov_mu[attach_idx[i]]
                    allSolutions_srcx[i, j] = Prov_srcx[attach_idx[i]]
                    allSolutions_srcy[i, j] = Prov_srcy[attach_idx[i]]
                    allSolutions_absdzs[i, j] = Prov_absdzs[attach_idx[i]]   

                    posflagnums[i] += Prov_flag[attach_idx[i]]

                if (curr_nimages - self.nlens)%2 != 1 or curr_total_parity != correct_parity:
                    # there might be some true images being classified as false images
                    ifsave, saveindex = self.saveimage(allSolutions_flag, allSolutions_mu, allSolutions_absdzs, curr_nimages, curr_total_parity, j)
                    if ifsave:
                        allSolutions_flag[saveindex, j] = 1
                        posflagnums[saveindex] += 1

        # return allSolutions_x, allSolutions_y, allSolutions_srcx, allSolutions_srcy, allSolutions_mu, allSolutions_flag, allSolutions_absdzs, None, None

        # now, select and connect true segments
        # initialize array to store the information of true image segments:
        # i.e., the i-th track, the index of head, the index of tail, (the length of the segments, this is not necessary)
        ntrue_segments = 0
        # true_segments_info = np.zeros((10 * self.DEGREE, 3 )).astype(int)
        true_segments_info = np.zeros((self.DEGREE * self.DEGREE, 3 )).astype(int)
        # each row contains where this segments is, a, b, c,
        npure_close_segments = 0

        imgTrack_type = np.zeros(self.DEGREE) # 0: pure false, 1: pure true, 2: mixed
        for i in range(self.DEGREE):
            # iter over each image Track, which may contains only true images, false images, or mixed
            if posflagnums[i] <= 1 or ( (1.0 * posflagnums[i] / nphi < 0.5) and abs(allSolutions_mu[i,0]) < EPS  and abs(allSolutions_mu[i,nphi-1]) < EPS ):
                if VERBOSE: print(i, "pure false")
                pass # pure false image Track
            elif posflagnums[i] == nphi: #(1.0 * posflagnums[i] / nphi) > 0.99:
                imgTrack_type[i] = 1
                if VERBOSE: print(i, "pure true")
                true_segments_info[ntrue_segments][0] = i
                true_segments_info[ntrue_segments][1] = 0
                true_segments_info[ntrue_segments][2] = nphi - 1 # index 0, 1, ..., nphi - 1 are true image points
                ntrue_segments += 1
                npure_close_segments += 1
            else:
                imgTrack_type[i] = 2
        # now for each of the type == 2 image Tracks, select out all true image segments
        for i in range(self.DEGREE):
            if imgTrack_type[i] == 2:
                # select true image tracks
                if VERBOSE: print(i, "now you only need to proceed this posflagnum %d, posflagnum/length %f, first->flag %d, mu %f, x1 x2 = %f, %f"%(posflagnums[i], 1.0 * posflagnums[i] / nphi, allSolutions_flag[i][0], allSolutions_mu[i][0], allSolutions_x[i][0], allSolutions_y[i][0]))
                # iter over this mixed track, record the indexes
                scan_true_idx = 0
                while scan_true_idx < nphi:
                    if allSolutions_flag[i][scan_true_idx]:
                        true_segments_info[ntrue_segments][0] = i
                        true_segments_info[ntrue_segments][1] = scan_true_idx
                        while scan_true_idx < nphi and allSolutions_flag[i][scan_true_idx] == 1:
                            scan_true_idx += 1   
                        true_segments_info[ntrue_segments][2] = scan_true_idx - 1
                        ntrue_segments += 1                        
                    else:
                        while scan_true_idx < nphi and allSolutions_flag[i][scan_true_idx] == 0:
                            scan_true_idx += 1

        # now, connect segments into closed track
        # if a segment has length == npoints, then pass, do not need to handle this, we just need to add this segment to the closed_image_info
        # otherwise, you need to connect the current segment
        # by the way, you need to record how many segments left, to be connected.

        # initialize array to store the information of true-close-image boundaries
        # there are at most self.DEGREE true-close-image boundaries
        nclosed_image = 0
        # closed_image_info = np.zeros(( min(self.DEGREE, ntrue_segments), 2 + 3 * (ntrue_segments - npure_close_segments + 1) )) # at most self.DEGREE closed-images, #### npure_close_segments is wrong number, not all segments with length nphi is closed
        closed_image_info = np.zeros(( min(self.DEGREE, ntrue_segments), 2 + 3 * (ntrue_segments + 1) )) # at most self.DEGREE closed-images, 
        # each row contains how this closed-image is built from segments, i.e.,
        # (parity, from_n_segments, a, b, c, ...)
        # a is the allSolutions index where the segments belongs to
        # b, c is the head_index, tail_index (if b>c, then means this segment is actually reversed)
        nfinal_closed_image = 0
        already_done_segments = np.zeros(ntrue_segments)
        open_seg_leftover = ntrue_segments

        if VERBOSE: print("true_segments_info before sort by length = ", true_segments_info[:ntrue_segments,:])
        # do we need to order by the length of the segments in ntrue_segments?
        # bubble sort? // pros: when create new image track, you always start with the longest one among the remaining segments
        tmp_segments_length = np.zeros(ntrue_segments)
        for i in range(ntrue_segments): tmp_segments_length[i] = (true_segments_info[i,2] - true_segments_info[i, 1]) + 1
        sorted_lenidx = np.argsort(tmp_segments_length)[::-1]
        true_segments_info = true_segments_info[sorted_lenidx]


        if VERBOSE: print("ntrue_segments = ", ntrue_segments)
        if VERBOSE: print(true_segments_info[:ntrue_segments,:])
        if VERBOSE: # print info of head tail at each segments
            for iv in range(ntrue_segments):
                print(">>>", iv, "true_segments_info: ", true_segments_info[iv, :], "x, y, mu, flag, absdzs")
                jv, hidv, tidv = true_segments_info[iv, :]
                print("head ", allSolutions_x[jv, hidv], allSolutions_y[jv, hidv], allSolutions_mu[jv, hidv], allSolutions_flag[jv, hidv], allSolutions_absdzs[jv, hidv])
                print("tail ", allSolutions_x[jv, tidv], allSolutions_y[jv, tidv], allSolutions_mu[jv, tidv], allSolutions_flag[jv, tidv], allSolutions_absdzs[jv, tidv])


        ### first find all closed tracks
        if_creat_new = True
        for i in range(ntrue_segments):
            if already_done_segments[i] or open_seg_leftover <= 0:
                continue
            j, hid, tid = true_segments_info[i, :]
            head1, tail1 = [j, hid], [j, tid]

            if hid == 0 and tid == nphi - 1 and self.head_tail_close(head1, tail1, allSolutions_x, allSolutions_y)<EPS**2:
                closed_image_info[nfinal_closed_image][0] = 1 if allSolutions_mu[j, hid] > 0 else -1 # final->first->mu > 0 ? 1 : -1;
                closed_image_info[nfinal_closed_image][1] = 1 # this closed image is originate from 1 segment
                closed_image_info[nfinal_closed_image][2] = j # which allSolutions_x index this image belongs to
                closed_image_info[nfinal_closed_image][3] = hid #
                closed_image_info[nfinal_closed_image][4] = tid #

                already_done_segments[i] = 1
                open_seg_leftover -= 1
                nfinal_closed_image += 1

                if VERBOSE: print(i, ">>>>>> initialize new segments, start from seg ", i, "open_seg_leftover = ", open_seg_leftover, "hid, tid = ", hid, tid)
                # to check whether a track is closed or not
                if VERBOSE: print(i, "close track", "head", allSolutions_x[head1[0], head1[1]], allSolutions_y[head1[0], head1[1]], allSolutions_x[tail1[0], tail1[1]], allSolutions_y[tail1[0], tail1[1]],)

        if_creat_new = True
        while open_seg_leftover > 0:
            for i in range(ntrue_segments):
                if already_done_segments[i] or open_seg_leftover <= 0:
                    continue
                j, hid, tid = true_segments_info[i, :]

                if if_creat_new:
                    head1, tail1 = [j, hid], [j, tid]

                    if hid == tid and open_seg_leftover > 1:
                        if VERBOSE: print("<<<<<< segment %d has only 1 data point, we do not start from this"%(i))
                        continue


                    closed_image_info[nfinal_closed_image][0] = 1 if allSolutions_mu[j, hid] > 0 else -1 # final->first->mu > 0 ? 1 : -1;
                    closed_image_info[nfinal_closed_image][1] = 1 # this closed image is originate from 1 segment
                    closed_image_info[nfinal_closed_image][2] = j # which allSolutions_x index this image belongs to
                    closed_image_info[nfinal_closed_image][3] = hid #
                    closed_image_info[nfinal_closed_image][4] = tid #

                    already_done_segments[i] = 1
                    open_seg_leftover -= 1

                    if if_creat_new and VERBOSE: print(i, ">>>>>> initialize new segments, start from seg ", i, "open_seg_leftover = ", open_seg_leftover, "hid, tid = ", hid, tid)
                    if_creat_new = False
                
                # if hid == 0 and tid == nphi - 1 and self.head_tail_close(head1, tail1, allSolutions_x, allSolutions_y)<EPS**2:
                if self.head_tail_close(head1, tail1, allSolutions_x, allSolutions_y)<EPS**2:
                    # to check whether a track is closed or not
                    if VERBOSE: print(i, "close track", "head", allSolutions_x[head1[0], head1[1]], allSolutions_y[head1[0], head1[1]], allSolutions_x[tail1[0], tail1[1]], allSolutions_y[tail1[0], tail1[1]],)
                    nfinal_closed_image += 1
                    if_creat_new = True
                    continue
                else:
                    ifcontinue = False
                    # now, test whether we can connect segment i with "other segments"
                    for i2 in range(ntrue_segments):
                        if open_seg_leftover > 0 and (not already_done_segments[i2]):
                            # # judge whether two segments (i, i2) can be connected together
                            j2, hid2, tid2 = true_segments_info[i2, :]
                            head2, tail2 = [j2, hid2], [j2, tid2]
                            # call a function to test whether these two segments can be connected
                            itype = self.if_two_segments_connect(head1, tail1, head2, tail2, allSolutions_x, allSolutions_y, EPS**2)
                            if itype > 0:
                                # we can connect seg_i with seg_i2
                                existing_seg_n = closed_image_info[nfinal_closed_image][1]
                                offset = int(3*existing_seg_n)
                                if itype == 3: # tail connect with head
                                    closed_image_info[nfinal_closed_image][2 + offset] = j2
                                    closed_image_info[nfinal_closed_image][3 + offset] = hid2
                                    closed_image_info[nfinal_closed_image][4 + offset] = tid2
                                    tail1 = tail2
                                elif itype == 4: # tail connect with tail
                                    closed_image_info[nfinal_closed_image][2 + offset] = j2
                                    closed_image_info[nfinal_closed_image][3 + offset] = tid2
                                    closed_image_info[nfinal_closed_image][4 + offset] = hid2
                                    tail1 = head2
                                elif itype == 2: # head connect with tail
                                    # move prev segments behind
                                    closed_image_info[nfinal_closed_image][5: 5+offset] = closed_image_info[nfinal_closed_image][2: 2+offset]
                                    closed_image_info[nfinal_closed_image][2] = j2
                                    closed_image_info[nfinal_closed_image][3] = hid2
                                    closed_image_info[nfinal_closed_image][4] = tid2
                                    head1 = head2
                                elif itype == 1: # head connect with head
                                    # move prev segments behind
                                    closed_image_info[nfinal_closed_image][5: 5+offset] = closed_image_info[nfinal_closed_image][2: 2+offset]
                                    closed_image_info[nfinal_closed_image][2] = j2
                                    closed_image_info[nfinal_closed_image][3] = tid2
                                    closed_image_info[nfinal_closed_image][4] = hid2
                                    head1 = tail2
                                closed_image_info[nfinal_closed_image][1] += 1
                                already_done_segments[i2] = 1
                                open_seg_leftover -= 1
                                #continue # you have connect seg2 with seg1, now proceed to the next segment
                                # check whether current is already closed
                                if self.head_tail_close(head1, tail1, allSolutions_x, allSolutions_y)<EPS**2 or open_seg_leftover<=0:
                                    nfinal_closed_image += 1
                                    if_creat_new = True
                                else:
                                    if VERBOSE: print("break 322, open_seg_leftover, nfinal_closed_image, if_creat_new = ", open_seg_leftover, nfinal_closed_image, if_creat_new)
                                    ifcontinue = True # connect once, then try whether we can connect again
                                    break
                    if ifcontinue: continue # if connected above, then continue, otherwise try jump
                    if VERBOSE: print("330 continue")
                    # for i2 in range(ntrue_segments):
                    #     if open_seg_leftover > 0 and (not already_done_segments[i2]):
                    # now, test whether we can jump from segment i to "other segments"
                    # different from connect, jump is more trivial, need to find the best place to jump
                    canwejump = 0
                    bestjumptype = 0
                    bestjump_i2 = 0
                    bestjump_fac_mu = 1e10 # find the minimum how_close, shoule be very close to 2.0
                    bestjump_dis = 1e10
                    for i3 in range(ntrue_segments):
                        if open_seg_leftover > 0 and (not already_done_segments[i3]):
                            # # judge whether two segments (i, i2) can be connected together
                            j2, hid2, tid2 = true_segments_info[i3, :]
                            head2, tail2 = [j2, hid2], [j2, tid2]
                            # call a function to test whether these two segments can be connected
                            
                            itype, how_close, jumpdis = self.if_two_segments_jump(head1, tail1, head2, tail2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu)

                            if bestjump_fac_mu > how_close and bestjump_dis >= jumpdis:
                                bestjump_fac_mu = how_close
                                bestjump_dis = jumpdis
                                bestjump_i2 = i3
                                bestjumptype = itype
                            if itype > 0 and jumpdis <= mindsource: # jumpdis has an upper limit: the mimimum distances between two sample points in the source limb, dphi
                                canwejump = 1
                                # if bestjump_fac_mu > how_close and bestjump_dis >= jumpdis:
                                #     bestjump_fac_mu = how_close
                                #     bestjump_dis = jumpdis
                                #     bestjump_i2 = i3
                                #     bestjumptype = itype
                            else:
                                if VERBOSE:
                                    print("440, we can not jump, the bestjump_fac_mu = %f, bestjump_dis = %e"%(bestjump_fac_mu, bestjump_dis))
                                    print("\talready_done_segments = ", already_done_segments)
                                    print("\tcurrent head: x, y, mu, flag, srcx, srcy = ", allSolutions_x[head1[0],head1[1]], allSolutions_y[head1[0],head1[1]],allSolutions_mu[head1[0],head1[1]],allSolutions_flag[head1[0],head1[1]],allSolutions_srcx[head1[0],head1[1]], allSolutions_srcy[head1[0],head1[1]])
                                    print("\tcurrent tail: x, y, mu, flag, srcx, srcy = ", allSolutions_x[tail1[0],tail1[1]], allSolutions_y[tail1[0],tail1[1]],allSolutions_mu[tail1[0],tail1[1]],allSolutions_flag[tail1[0],tail1[1]],allSolutions_srcx[tail1[0],tail1[1]], allSolutions_srcy[tail1[0],tail1[1]])                                   
                                    print("\ttmp seg %d head: x, y, mu, flag, srcx, srcy = "%i3, allSolutions_x[head2[0],head2[1]], allSolutions_y[head2[0],head2[1]],allSolutions_mu[head2[0],head2[1]],allSolutions_flag[head2[0],head2[1]],allSolutions_srcx[head2[0],head2[1]], allSolutions_srcy[head2[0],head2[1]])
                                    print("\ttmp seg %d tail: x, y, mu, flag, srcx, srcy = "%i3, allSolutions_x[tail2[0],tail2[1]], allSolutions_y[tail2[0],tail2[1]],allSolutions_mu[tail2[0],tail2[1]],allSolutions_flag[tail2[0],tail2[1]], allSolutions_srcx[tail2[0],tail2[1]], allSolutions_srcy[tail2[0],tail2[1]])
                    if canwejump:
                        #print(i, "best jump to ", bestjump_i2, "bestjump_fac_mu = ", bestjump_fac_mu)
                        
                        # can find some place to jump, 
                        j2, hid2, tid2 = true_segments_info[bestjump_i2, :]
                        head2, tail2 = [j2, hid2], [j2, tid2]             

                        if VERBOSE: print(">>> best jump to ", bestjump_i2, "bestjump_fac_mu = ", bestjump_fac_mu, 'type = ', bestjumptype, allSolutions_x[head2[0], head2[1]], allSolutions_y[head2[0], head2[1]], allSolutions_x[tail2[0], tail2[1]], allSolutions_y[tail2[0], tail2[1]], "bestjump_dis = ", bestjump_dis)           

                        # we can connect seg_i with seg_i2
                        existing_seg_n = closed_image_info[nfinal_closed_image][1]
                        offset = int(3*existing_seg_n)
                        if bestjumptype == 3: # tail connect with head
                            closed_image_info[nfinal_closed_image][2 + offset] = j2
                            closed_image_info[nfinal_closed_image][3 + offset] = hid2
                            closed_image_info[nfinal_closed_image][4 + offset] = tid2
                            tail1 = tail2
                        elif bestjumptype == 4: # tail connect with tail
                            closed_image_info[nfinal_closed_image][2 + offset] = j2
                            closed_image_info[nfinal_closed_image][3 + offset] = tid2
                            closed_image_info[nfinal_closed_image][4 + offset] = hid2
                            tail1 = head2
                        elif bestjumptype == 2: # head connect with tail
                            # move prev segments behind
                            closed_image_info[nfinal_closed_image][5: 5+offset] = closed_image_info[nfinal_closed_image][2: 2+offset]
                            closed_image_info[nfinal_closed_image][2] = j2
                            closed_image_info[nfinal_closed_image][3] = hid2
                            closed_image_info[nfinal_closed_image][4] = tid2
                            head1 = head2
                        elif bestjumptype == 1: # head connect with head
                            # move prev segments behind
                            closed_image_info[nfinal_closed_image][5: 5+offset] = closed_image_info[nfinal_closed_image][2: 2+offset]
                            closed_image_info[nfinal_closed_image][2] = j2
                            closed_image_info[nfinal_closed_image][3] = tid2
                            closed_image_info[nfinal_closed_image][4] = hid2
                            head1 = tail2
                        closed_image_info[nfinal_closed_image][1] += 1
                        already_done_segments[bestjump_i2] = 1
                        open_seg_leftover -= 1

                        # check whether current is already closed
                        if self.head_tail_close(head1, tail1, allSolutions_x, allSolutions_y)<EPS**2 or open_seg_leftover<=0:
                            nfinal_closed_image += 1
                            if_creat_new = True
                        else:
                            if VERBOSE: print("break 398, open_seg_leftover, nfinal_closed_image", open_seg_leftover, nfinal_closed_image)
                            if VERBOSE: print(already_done_segments)
                            continue # after jump, not close, and open_seg_leftover > 0

                    else:
                        # we can not connect, and we can not jump, we need to test whether current
                        # we need to test whether it is closed, other wise, we regard this as a close image anyway
                        # new test whether head1, tail1 is close, and whether there is no open segments left
                        #if self.head_tail_close(head1, tail1, allSolutions_x, allSolutions_y, EPS**2):
                        if VERBOSE: print("410, open_seg_leftover", open_seg_leftover)
                        if open_seg_leftover > 0:
                            nfinal_closed_image += 1
                            if_creat_new = True

        if VERBOSE: print("nfinal_closed_image = ", nfinal_closed_image)
        if VERBOSE: print(closed_image_info[:nfinal_closed_image,:].astype(int))

        return allSolutions_x, allSolutions_y, allSolutions_srcx, allSolutions_srcy, allSolutions_mu, allSolutions_flag, allSolutions_absdzs, true_segments_info[:ntrue_segments,:].astype(int), closed_image_info[:nfinal_closed_image,:].astype(int)

    # def select_connect_segments(self, allSolutions_x, allSolutions_y, allSolutions_srcx, allSolutions_srcy, allSolutions_mu, allSolutions_flag, allSolutions_absdzs):


    def head_tail_close(self, head1, tail1, allSolutions_x, allSolutions_y):
        x1 = allSolutions_x[head1[0], head1[1]]
        y1 = allSolutions_y[head1[0], head1[1]]
        x2 = allSolutions_x[tail1[0], tail1[1]]
        y2 = allSolutions_y[tail1[0], tail1[1]]
        #print("x1, y1, x2, y2", x1, y1, x2, y2)
        dis = self.if_dis_close(x1, y1, x2, y2)
        return dis

    def if_dis_close(self, x1, y1, x2, y2):
        dis = (x1 - x2)**2 + (y1 - y2)**2
        return dis

    def if_head_tail_jump(self, head1, tail1, allSolutions_srcx, allSolutions_srcy, allSolutions_mu):
        srcxprev = allSolutions_srcx[head1[0], head1[1]]
        srcyprev = allSolutions_srcy[head1[0], head1[1]]
        srcxpost = allSolutions_srcx[tail1[0], tail1[1]]
        srcypost = allSolutions_srcy[tail1[0], tail1[1]]
        muprev = allSolutions_mu[head1[0], head1[1]]
        mupost = allSolutions_mu[tail1[0], tail1[1]]
        ifjump, how_close, dis = self.if_jump( srcxprev, srcyprev, srcxpost, srcypost, muprev, mupost)
        return ifjump, how_close, dis

    def if_jump(self, srcxprev, srcyprev, srcxpost, srcypost, muprev, mupost):
        how_close = ( abs(muprev / mupost) + abs(mupost/muprev) )
        dis = self.if_dis_close(srcxprev, srcyprev, srcxpost, srcypost)
        
        ifjump = False
        if muprev * mupost < 0.0:
            if dis <= 1e-15:
                ifjump = True
            else:
                ifjump = how_close < 2.5

        return ifjump, how_close, dis

    def if_two_segments_connect(self, head1, tail1, head2, tail2, allSolutions_x, allSolutions_y, EPS2):
        # head1 = [head1_j, head1_idx]

        #print(head1, tail1, head2, tail2)

        # connect seg1 and seg2
        # itype 0 = can not connect
        # itype 1 = head head --> how to handle this? --> you need to insert seg2 into the front of seg1, and reverse seg2
        # itype 2 = head tail --> how to handle this? --> you need to insert seg2 into the front of seg1, do not need to reverse seg2
        # itype 3 = tail head --> how to handle this? --> the easiest scenerio, just add seg2 behind seg1
        # itype 4 = tail tail --> how to handle this? --> add seg2 behind seg1, and reverse seg2

        itype = 0
        # we prefer connect, rather than jump, so, first test image connectivity; if we cannot connect, then try whether we can "jump"
        # besides, you need to have a seperate function for jump

        # we prefer to find the easiest scenario
        if self.head_tail_close( tail1, head2, allSolutions_x, allSolutions_y)<EPS2:
            itype = 3 # tail connect with head
            if VERBOSE: print("tail connect with head", itype, tail1, head2, allSolutions_x[tail1[0], tail1[1]], allSolutions_y[tail1[0], tail1[1]] )

        if self.head_tail_close( tail1, tail2, allSolutions_x, allSolutions_y)<EPS2:
            itype = 4 # tail connect with tail
            if VERBOSE: print("tail connect with tail", itype, tail1, tail2, allSolutions_x[tail1[0], tail1[1]], allSolutions_y[tail1[0], tail1[1]])

        if self.head_tail_close( head1, tail2, allSolutions_x, allSolutions_y)<EPS2:
            itype = 2 # head connect with tail
            if VERBOSE: print("head connect with tail", itype, head1, tail2, allSolutions_x[tail2[0], tail2[1]], allSolutions_y[tail2[0], tail2[1]])

        if self.head_tail_close( head1, head2, allSolutions_x, allSolutions_y)<EPS2:
            itype = 1 # head connect with head
            if VERBOSE: print("head connect with head", itype, head1, head2, allSolutions_x[head1[0], head1[1]], allSolutions_y[head1[0], head1[1]])
        return itype


    def if_two_segments_jump(self, head1, tail1, head2, tail2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu):
        itype = 0
        # we prefer connect, rather than jump, so, first test image connectivity; if we cannot connect, then try whether we can "jump"
        # besides, you need to have a seperate function for jump

        # we prefer to find the easiest scenario
        besttype = 0
        bestjump_fac_mu = 1e10
        bestjump_dis= 1e10

        ifjump, how_close, dis = self.if_head_tail_jump( tail1, head2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu)
        if bestjump_fac_mu > how_close and bestjump_dis >= dis:
            bestjump_fac_mu = how_close
            bestjump_dis = dis
        if ifjump:     besttype = 3 # tail connect with head

        ifjump, how_close, dis = self.if_head_tail_jump( tail1, tail2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu)
        if bestjump_fac_mu > how_close and bestjump_dis >= dis:
            bestjump_fac_mu = how_close
            bestjump_dis = dis
        if ifjump:     besttype = 4 # tail connect with tail
            

        ifjump, how_close, dis =  self.if_head_tail_jump( head1, tail2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu)
        if bestjump_fac_mu > how_close and bestjump_dis >= dis:
            bestjump_fac_mu = how_close
            bestjump_dis = dis
        if ifjump:     besttype = 2 # head connect with tail
        

        ifjump, how_close, dis =  self.if_head_tail_jump( head1, head2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu)
        if bestjump_fac_mu > how_close and bestjump_dis >= dis:
            bestjump_fac_mu = how_close
            bestjump_dis = dis
        if ifjump:     besttype = 1 # head connect with head

        # ifjump, how_close, dis = self.if_head_tail_jump( tail1, head2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu)
        # if ifjump and bestjump_fac_mu > how_close and bestjump_dis >= dis:
        #     bestjump_fac_mu = how_close
        #     bestjump_dis = dis
        #     besttype = 3 # tail connect with head

        # ifjump, how_close, dis = self.if_head_tail_jump( tail1, tail2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu)
        # if ifjump and bestjump_fac_mu > how_close and bestjump_dis >= dis:
        #     bestjump_fac_mu = how_close
        #     bestjump_dis = dis
        #     besttype = 4 # tail connect with tail
            

        # ifjump, how_close, dis =  self.if_head_tail_jump( head1, tail2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu)
        # if ifjump and bestjump_fac_mu > how_close and bestjump_dis >= dis:
        #     bestjump_fac_mu = how_close
        #     bestjump_dis = dis
        #     besttype = 2 # head connect with tail
        

        # ifjump, how_close, dis =  self.if_head_tail_jump( head1, head2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu)
        # if ifjump and bestjump_fac_mu > how_close and bestjump_dis >= dis:
        #     bestjump_fac_mu = how_close
        #     bestjump_dis = dis
        #     besttype = 1 # head connect with head

        if besttype == 3 and VERBOSE: print("tail jump to head", besttype, tail1, head2, bestjump_fac_mu, bestjump_dis)
        if besttype == 4 and VERBOSE: print("tail jump to tail", besttype, tail1, tail2, bestjump_fac_mu,bestjump_dis)
        if besttype == 2 and VERBOSE: print("head jump to tail", besttype, head1, tail2, bestjump_fac_mu,bestjump_dis)
        if besttype == 1 and VERBOSE: print("head jump to head", besttype, head1, head2, bestjump_fac_mu,bestjump_dis)

        return besttype, bestjump_fac_mu, bestjump_dis

    def show_connected_tracks_static(self, allSolutions_x, allSolutions_y, true_segments_info, step = 100,xlim=(-1.5,1.5),ylim=(-1.5,1.5), txt = 1, Narrows_each_track = 15, colors = [], txtstr = "H{}", head_center = "head", onlytrue = False, mus = None, srcx = None, srcy = None, showfalse = False):
        # onlytrue = 1, plot only true image segments

        #print(true_segments_info.shape)

        if onlytrue == False:
            nsegments, nphi, = allSolutions_x.shape
            print("nsegments, nphi", nsegments, nphi)
            true_segments_info = np.zeros((nsegments, 3)).astype(int)
            true_segments_info[:,0] = np.arange(nsegments).astype(int)
            true_segments_info[:,2] = nphi

        #ntrack = len(allSolutions_x[:,0])
        nsegments = true_segments_info.shape[0]
        headx, heady, tailx, taily = [],[],[],[]
        for jj in range(nsegments):
            j, hid, tid = true_segments_info[jj, 0], true_segments_info[jj, 1], true_segments_info[jj, 2]
            if onlytrue == False: tid -= 1
            if VERBOSE: print(j, "track length = ", tid - hid + 1, hid, tid)
            headx.append( allSolutions_x[j, hid] )    
            heady.append( allSolutions_y[j, hid] )
            tailx.append( allSolutions_x[j, tid] )
            taily.append( allSolutions_y[j, tid] )

        fig, ax = plt.subplots(figsize=(8,8))
        plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.15, hspace = 0, wspace = 0)
        ax.tick_params(axis='both', labelsize = 17, direction="in")
        
        critical , caustics = get_crit_caus(self.mlens, self.zlens, self.nlens, NPS = 3000)
        ax.plot([xy[0]for xy in caustics], [xy[1]for xy in caustics], 'o', color='red', markersize=2)
        ax.plot([xy[0]for xy in critical], [xy[1]for xy in critical], 'o', color='k', markersize=1)
        ax.plot([xy[0]for xy in self.zlens], [xy[1]for xy in self.zlens], 'x', color='k', markersize=5)
        PHIS = np.linspace(0, 2*np.pi, len(allSolutions_x[0,:]))
        XS = xsCenter + rs * cos(PHIS)
        YS = ysCenter + rs * sin(PHIS)
        ax.plot(XS, YS, '.', color="k" , markersize=1)

        ncolor = len(colors)
        #print("ncolor = ", ncolor)
        scal = 5e-3

        if showfalse: # show false images
            ax.plot(allSolutions_x, allSolutions_y, '.', color='gray', markersize=0.5)

        for hx,hy,tx,ty,jj in zip(headx, heady, tailx, taily, range(len(headx))):
            # npts = len(allSolutions_x[j,:]) # number of points in that track
            j, hid, tid = true_segments_info[jj, 0], true_segments_info[jj, 1], true_segments_info[jj, 2]
            npts = tid - hid + 1
            if txt:
                if head_center == "head": 
                    ax.text(hx-1.1*scal,hy-3*scal,"H{}_{:.2f}_{:.2e}_{:.2e}".format(jj,mus[j,hid],srcx[j,hid],srcy[j,hid]), color=colors[jj%ncolor], fontsize = 7)
                    ax.text(tx+1.1*scal,ty+3*scal,"T{}_{:.2f}_{:.2e}_{:.2e}".format(jj,mus[j,tid],srcx[j,tid],srcy[j,tid]), color=colors[jj%ncolor], fontsize = 7)
            ax.plot(hx, hy, '*', color=colors[jj%ncolor], markersize=20)
            ax.plot(tx, ty, 'o', color=colors[jj%ncolor], markersize=20, fillstyle = "none")
            line = ax.plot(allSolutions_x[j, hid:tid+1], allSolutions_y[j, hid:tid+1], '.', markersize=2, color=colors[jj%ncolor])[0]
            if 1:
                try:
                    for i in range(hid , hid+npts, int((npts-5)/Narrows_each_track) ):
                        # print("i:", i)
                        try:
                            x0 = allSolutions_x[j, i]
                            y0 = allSolutions_y[j, i]
                            x1 = allSolutions_x[j, i+3]
                            y1 = allSolutions_y[j, i+3]
                        except:
                            continue
                        # plt.arrow(x0, y0 , x1,  y1 - y0, shape='full', lw=0, length_includes_head=True, head_width=.05, color = colors[j%ncolor])
                        line.axes.annotate('',
                            xytext=(x0, y0),
                            xy=(x1, y1),
                            # arrowprops=dict(arrowstyle="->", color=colors[j%ncolor]),
                            arrowprops=dict(arrowstyle="-|>", color=colors[jj%ncolor]),
                            # https://matplotlib.org/3.3.2/tutorials/text/annotations.html
                            size=10
                        )
                except:
                    pass
                    print("might be wrong at 595", j, hid, tid, npts)

        ax.set_xlabel(r"$x_1/ \theta_E $", fontsize = 17, fontname='Times New Roman')
        ax.set_ylabel(r"$x_2/ \theta_E $", fontsize = 17,fontname='Times New Roman')

        ax.set_xlim(xlim[0],xlim[1])
        ax.set_ylim(ylim[0],ylim[1])


    def show_closed(self, allSolutions_x, allSolutions_y, closed_image_info, step = 100,xlim=(-1.5,1.5),ylim=(-1.5,1.5), txt = 1, Narrows_each_track = 15, colors = [], txtstr = "H{}", head_center = "head", ax = None, xs = 0, ys = 0, showfalse = False):

        nsegments = closed_image_info.shape[0]
        headx, heady, tailx, taily = [],[],[],[]
        for jj in range(nsegments):
            parity, nsegs = closed_image_info[jj, 0], closed_image_info[jj, 1]

            hj, hid = closed_image_info[jj, 2], closed_image_info[jj, 3]
            tj, tid = closed_image_info[jj, int(2+3*nsegs-3)], closed_image_info[jj, int(2+3*nsegs-1)]
            if VERBOSE: print(jj, "hj, hid, tj, tid", hj, hid, tj, tid)
            #print(j, "track length = ", tid - hid + 1)
            headx.append( allSolutions_x[hj, hid] )    
            heady.append( allSolutions_y[hj, hid] )
            tailx.append( allSolutions_x[tj, tid] )
            taily.append( allSolutions_y[tj, tid] )

        if not ax:
            fig, ax = plt.subplots(figsize=(8,8))
            plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.15, hspace = 0, wspace = 0)
            ax.tick_params(axis='both', labelsize = 17, direction="in")
        
        critical , caustics = get_crit_caus(self.mlens, self.zlens, self.nlens, NPS = 3000)
        ax.plot([xy[0]for xy in caustics], [xy[1]for xy in caustics], '-', color='red', markersize=2) # ax.plot(causticsx, causticsy, '-', color='red', markersize=1)
        ax.plot([xy[0]for xy in critical], [xy[1]for xy in critical], '--', color='red', markersize=1) # ax.plot(criticalx, criticaly, '--', color='r', markersize=1)
        ax.plot([xy[0]for xy in self.zlens], [xy[1]for xy in self.zlens], 'x', color='k', markersize=5)
        PHIS = np.linspace(0, 2*np.pi, len(allSolutions_x[0,:]))
        XS = xs + rs * cos(PHIS)
        YS = ys + rs * sin(PHIS)
        ax.plot(XS, YS, '.', color="k" , markersize=1)

        ncolor = len(colors)

        if showfalse: # show false images
            ax.plot(allSolutions_x, allSolutions_y, '.', color='gray', markersize=0.5)

        if txt: ax.annotate('(${}$, ${}$)'.format(xs,ys), xy=(0.3, 0.9), xycoords='axes fraction', fontsize=17,horizontalalignment='right', verticalalignment='bottom')

        scal = 5e-3

        for hx,hy,tx,ty,jj in zip(headx, heady, tailx, taily, range(len(headx))):
            # npts = len(allSolutions_x[j,:]) # number of points in that track

            if txt:
                    if head_center == "head": 
                        ax.text(hx-1.1*scal,hy-3*scal,txtstr.format(jj), color=colors[jj%ncolor], fontsize = 17)
                        ax.text(tx+1.1*scal,ty+3*scal,"T{}".format(jj), color=colors[jj%ncolor], fontsize = 17)
            ax.plot(hx, hy, '*', color=colors[jj%ncolor], markersize=20)
            ax.plot(tx, ty, 'o', color=colors[jj%ncolor], markersize=20, fillstyle = "none")

            parity, nsegs = closed_image_info[jj, 0], closed_image_info[jj, 1]
            NPT = 0
            for k in range(nsegs):
                offset = k*3
                j, hid, tid = closed_image_info[jj, 2+offset: 5+offset]
                
                npts = abs(tid - hid) + 1
                NPT += npts

                if hid < tid:
                    isgn = 1
                    line = ax.plot(allSolutions_x[j, hid:tid+1], allSolutions_y[j, hid:tid+1], '.', markersize=2, color=colors[jj%ncolor])[0]
                else:
                    isgn = -1
                    line = ax.plot(allSolutions_x[j, hid:tid:-1], allSolutions_y[j, hid:tid:-1], '.', markersize=2, color=colors[jj%ncolor])[0]

                if npts > 50:
                    Narrows_each_track0 = Narrows_each_track
                else:
                    Narrows_each_track0 = 2

                if 1:
                    try:
                        for i in range(hid , tid, int(isgn*(npts-5)/Narrows_each_track0) ):
                            # print("i:", i)
                            try:
                                x0 = allSolutions_x[j, i]
                                y0 = allSolutions_y[j, i]
                                x1 = allSolutions_x[j, int(i+isgn*3)]
                                y1 = allSolutions_y[j, int(i+isgn*3)]
                            except:
                                continue
                            if isgn>0 and int(i+isgn*3)>tid:
                                continue
                            if isgn<0 and int(i+isgn*3)<hid:
                                continue
                            # plt.arrow(x0, y0 , x1,  y1 - y0, shape='full', lw=0, length_includes_head=True, head_width=.05, color = colors[j%ncolor])
                            line.axes.annotate('',
                                xytext=(x0, y0),
                                xy=(x1, y1),
                                # arrowprops=dict(arrowstyle="->", color=colors[j%ncolor]),
                                arrowprops=dict(arrowstyle="-|>", color=colors[jj%ncolor]),
                                # https://matplotlib.org/3.3.2/tutorials/text/annotations.html
                                size=10
                            )
                    except:
                        print("might be wrong at source position", xs, ys)
            if VERBOSE: print(jj, "track total length = ", NPT)
        ax.set_xlabel(r"$x_1/ \theta_E $", fontsize = 17, fontname='Times New Roman')
        ax.set_ylabel(r"$x_2/ \theta_E $", fontsize = 17,fontname='Times New Roman')

        ax.set_xlim(xlim[0],xlim[1])
        ax.set_ylim(ylim[0],ylim[1])

    def areaFunc(self, allSolutions_x, allSolutions_y, closed_image_info, rs):
        nclosed_image = closed_image_info.shape[0]
        area_array = np.zeros(nclosed_image)
        total_area = 0

        for jj in range(nclosed_image):
            parity, nsegs = closed_image_info[jj, 0], closed_image_info[jj, 1]
            NPT = 0
            area = 0
            for k in range(nsegs):
                offset = k*3
                j, hid, tid = closed_image_info[jj, 2+offset: 5+offset]
                npts = abs(tid - hid) + 1
                if npts < 2:
                    print(jj, k, "this segment has only %d point"%npts)
                    continue # just in case 
                NPT += npts
                if hid < tid:
                    isgn = 1
                else:
                    isgn = -1                    
                scanpnt = hid
                x1, y1 = allSolutions_x[j, scanpnt], allSolutions_y[j, scanpnt]
                scanpnt += isgn
                while scanpnt != tid+isgn:
                    x2, y2 = allSolutions_x[j, scanpnt], allSolutions_y[j, scanpnt]
                    area += 0.5 * (y1 + y2) * (x2 - x1)
                    x1, y1 = x2, y2
                    scanpnt += isgn
            if nsegs > 1: # if nsegs > 1, you need to add the area between two segments
                for k in range(nsegs - 1):
                    offset = k*3
                    j, hid, tid = closed_image_info[jj, 2+offset: 5+offset]
                    x1, y1 = allSolutions_x[j, tid], allSolutions_y[j, tid]
                    offset = k*3 + 3
                    j, hid, tid = closed_image_info[jj, 2+offset: 5+offset]
                    x2, y2 = allSolutions_x[j, hid], allSolutions_y[j, hid]
                    area += 0.5 * (y1 + y2) * (x2 - x1)
            # the tail and head
            hj, hid = closed_image_info[jj, 2], closed_image_info[jj, 3]
            tj, tid = closed_image_info[jj, int(2+3*nsegs-3)], closed_image_info[jj, int(2+3*nsegs-1)]            
            x1, y1 = allSolutions_x[tj, tid], allSolutions_y[tj, tid]
            x2, y2 = allSolutions_x[hj, hid], allSolutions_y[hj, hid]
            area += 0.5 * (y1 + y2) * (x2 - x1)

            if VERBOSE: print(jj, "track total length = ", NPT, "parity = ", parity, "area = ", area)
            area *= parity

            total_area += area

            #area_array[jj] = area
        return abs(total_area/ (np.pi * rs * rs) )

if __name__ == "__main__":

    #set up lens system
    #fractional lens masses: m1, m2, m3
    mlens = [0.968738798957637, 0.028093425169771, 0.003167775872591]
    #lens positions: x1, y1, x2, y2, x3, y3
    zlens = [-0.039343051506317, 0, 1.356656948493683, 0, 0.638936196010800, -0.950873946634155]

    #set up lens system
    #fractional lens masses: m1, m2, m3
    mlens = [0.968738798957637, 0.020093425169771, 0.003167775872591, 0.008]
    #lens positions: x1, y1, x2, y2, x3, y3
    zlens = [-0.039343051506317, 0, 1.356656948493683, 0, 0.638936196010800, -0.950873946634155, -0.638936196010800, 0.950873946634155,]

    #set up lens system
    #fractional lens masses: m1, m2, m3
    mlens = [0.968738798957637, 0.028093425169771]
    #lens positions: x1, y1, x2, y2, x3, y3
    zlens = [-0.039343051506317, 0, 1.356656948493683, 0]


    #source center
    xsCenter = -0.034747426672208
    ysCenter = -0.026627816352184

    #source radius
    rs = 0.005
    rs = 0.1

    zlens0 = []
    for i in range(len(mlens)):
        zlens0.append([zlens[2*i],zlens[2*i+1]])
    zlens = zlens0
    print(zlens)

    # colors = cm.rainbow(np.linspace(0,1,10))  #["blue", "green", "black", "red", "orange", "salmon", "lime"] * 2
    #["blue", "green", "black", "red", "orange", "salmon", "lime"] * 2

    pyTRIL = pyTriple(mlens, zlens)

    print("pyTRIL.nlens = ", pyTRIL.nlens)

    # self, xsCenter, ysCenter, rs, PHI, prevstore = None

    PHI = np.linspace(0, 2*np.pi, 3000)
    mindphi = np.min(PHI[1:] - PHI[:-1]) # positive

    xsCenter = -0.12
    xsCenter = 0.6777777777777777777777777777
    ysCenter = 0


    #source center
    xsCenter = -0.034747426672208
    ysCenter = -0.026627816352184
    #source radius
    rs = 0.005

    rs = 0.1
    xsCenter =  -0.1171717171717172
    xsCenter = -0.13131313131313133
    ysCenter = 0
    # xsCenter = -0.4
    rs = 0.05
    # xsCenter = -0.004040404040404066
    # xsCenter = 0.010101010101010055
    # xsCenter = -0.018181818181818243
    xsCenter = -0.02131313131313133


    if 1: # show static
        for iteri in range(1):
            if VERBOSE: print("iteri = ", iteri)
            if iteri == 0:
                prevstore = False

                allSolutions_x, allSolutions_y, allSolutions_srcx, allSolutions_srcy, allSolutions_mu, allSolutions_flag, allSolutions_absdzs, true_segments_info, closed_image_info = pyTRIL.outputTracks_v2_savehalf(xsCenter, ysCenter, rs, PHI, prevstore, mindphi = mindphi)
            else:
                prevstore = True

                nphi = len(PHI)
                midPHI = 0.5*(PHI[:-1] + PHI[1:])
                newPHI = np.zeros( nphi + nphi - 1 )
                newPHI[::2] = PHI
                newPHI[1::2] = midPHI
                PHI = newPHI
                allSolutions_x, allSolutions_y, allSolutions_srcx, allSolutions_srcy, allSolutions_mu, allSolutions_flag, allSolutions_absdzs, true_segments_info, closed_image_info = pyTRIL.outputTracks_v2_savehalf(xsCenter, ysCenter, rs, PHI, prevstore, prevstore_x = allSolutions_x, prevstore_y = allSolutions_y, prevstore_srcx = allSolutions_srcx, prevstore_srcy = allSolutions_srcy, prevstore_mu = allSolutions_mu, prevstore_flag = allSolutions_flag, prevstore_absdzs = allSolutions_absdzs)
        if 1:
            colors = cm.rainbow(np.linspace(0,1,allSolutions_x.shape[0]))
            colors = cm.rainbow(np.linspace(0,1,true_segments_info.shape[0]))
            pyTRIL.show_connected_tracks_static(allSolutions_x, allSolutions_y, true_segments_info, colors = colors, xlim = (-1.4, 2), ylim = (-1.7, 1.7), onlytrue = True, mus = allSolutions_mu, srcx = allSolutions_srcx, srcy = allSolutions_srcy, showfalse = True)
            plt.show()
        else: # show closed image boundaries
            colors = cm.seismic(np.linspace(0,1,closed_image_info.shape[0]))
            colors = cm.rainbow(np.linspace(0,1,closed_image_info.shape[0]))
            pyTRIL.show_closed(allSolutions_x, allSolutions_y, closed_image_info,txt = 1,showfalse = True, colors = colors, xlim = (-1.4, 2), ylim = (-1.7, 1.7), xs = xsCenter, ys = ysCenter)
            mu = pyTRIL.areaFunc(allSolutions_x, allSolutions_y, closed_image_info, rs)
            print("mu = ", mu)

            plt.show()


    # show movie
    if 0:

        import numpy as np
        import matplotlib
        # matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.animation as manimation
        from matplotlib import gridspec
        from tqdm import tqdm

        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='Movie Test', artist='Matplotlib',
                        comment='Movie support!')
        writer = FFMpegWriter(fps=4, metadata=metadata)

        prevstore = False

        ysCenter = 0

        # fig, axs = plt.subplots(figsize=(16,8), nrows=1, ncols=2)
        # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.15, hspace = 0, wspace = 0)
        # ax = axs[0]
        # ax2 = axs[1]

        fig = plt.figure(figsize=(15,7))
        ax = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0, wspace = 0.15)

        ax.tick_params(axis='both', labelsize = 17, direction="in")
        ax2.tick_params(axis='both', labelsize = 17, direction="in")

        if 0:
            for xsCenter in np.linspace(-0.2, 0.9, 10):

                allSolutions_x, allSolutions_y, allSolutions_srcx, allSolutions_srcy, allSolutions_mu, allSolutions_flag, allSolutions_absdzs, true_segments_info, closed_image_info = pyTRIL.outputTracks_v2_savehalf(xsCenter, ysCenter, rs, PHI, prevstore)

                colors = cm.rainbow(np.linspace(0,1,closed_image_info.shape[0]))
                pyTRIL.show_closed(allSolutions_x, allSolutions_y, closed_image_info, colors = colors, xlim = (-1.2, 1.8), ax = ax)
                plt.show()

        # top_lkv: /Users/anything/THU/astro/softwares/gravlens/triplelens/test/gentopos_movie2.py


        rs = 0.1
        rs = 0.05
        xlim = -0.4, 1
        ylim = 0.2, 1.65
        xsCenters = np.linspace(xlim[0], xlim[1], 100) # NPHI = 3000 # first try with 10 points, to locate the peak logmu, so that you can set_ylim

        
        testsuffix = 'rs_%s'%rs
        # testsuffix = '4lens'
        testsuffix = 'tmp'
        testsuffix = '2lens'
        mufilename = "../data/pymu_22Nov01_%s.txt"%testsuffix

        if 0:
            if os.path.exists(mufilename):
                mus = np.loadtxt(mufilename)
            else:

                mus = []
                for idx in range(len(xsCenters)):
                    try:
                        allSolutions_x, allSolutions_y, allSolutions_srcx, allSolutions_srcy, allSolutions_mu, allSolutions_flag, allSolutions_absdzs, true_segments_info, closed_image_info = pyTRIL.outputTracks_v2_savehalf(xsCenters[idx], ysCenter, rs, PHI, prevstore)
                        mu = pyTRIL.areaFunc(allSolutions_x, allSolutions_y, closed_image_info, rs)
                    except:
                        print(idx, "wrong, xsCenter = ", xsCenters[idx]) # 20 wrong, xsCenter =  -0.1171717171717172
                        mu = 1
                    mus.append(mu)
                np.savetxt(mufilename, mus)
            mus = np.log10(mus)


            ax2.plot(xsCenters, np.log10(mus))
            plt.show()
            input(">>>>")

        mus = np.zeros(len(xsCenters))
        def update_figure(idx):
            # print("xsCenters[idx] = ", xsCenters[idx])
            ax.clear()
            ax2.clear()

            allSolutions_x, allSolutions_y, allSolutions_srcx, allSolutions_srcy, allSolutions_mu, allSolutions_flag, allSolutions_absdzs, true_segments_info, closed_image_info = pyTRIL.outputTracks_v2_savehalf(xsCenters[idx], ysCenter, rs, PHI, prevstore)

            colors = cm.rainbow(np.linspace(0,1,closed_image_info.shape[0]))
            pyTRIL.show_closed(allSolutions_x, allSolutions_y, closed_image_info, showfalse = True, colors = colors, txt = 0, xlim = (-1.4, 2), ylim = (-1.7, 1.7), xs = xsCenters[idx], ys=ysCenter, ax = ax)
            mu = pyTRIL.areaFunc(allSolutions_x, allSolutions_y, closed_image_info, rs)
            mus[idx] = np.log10(mu)

            #ax2.plot(xsCenters, mus, '.', color = "gray")
            ax2.plot(xsCenters[:idx], mus[:idx], '-', color = "red")
            ax2.set_xlabel(r"$x/ \theta_E $", fontsize = 17, fontname='Times New Roman')
            ax2.set_ylabel(r"log$(\mu)$", fontsize = 17,fontname='Times New Roman')
        
            # ax.set_xlim(-1.35,1.85)
            # ax.set_ylim(-1.6,1.6)

            ax2.set_xlim(xlim[0], xlim[1])
            ax2.set_ylim(ylim[0], ylim[1])

        # https://matplotlib.org/1.4.1/examples/animation/moviewriter.html
        with writer.saving(fig, "check_img_22Nov01_%s.mp4"%testsuffix, dpi=300):
            with tqdm(total=len(xsCenters)) as pbar:
                for j in range(len(xsCenters)):
                    update_figure(j)
                    writer.grab_frame()
                    pbar.update(1)

        np.savetxt(mufilename, mus)



### debug:

        # print(i, i2, "can not attach")

        # xprev, yprev, srcxprev, srcyprev, muprev = allSolutions_x[j,hid], allSolutions_y[j,hid], allSolutions_srcx[j,hid], allSolutions_srcy[j,hid], allSolutions_mu[j,hid]
        # print(i, "head", xprev, yprev, srcxprev, srcyprev, muprev)
        # xprev, yprev, srcxprev, srcyprev, muprev = allSolutions_x[j,tid], allSolutions_y[j,tid], allSolutions_srcx[j,tid], allSolutions_srcy[j,tid], allSolutions_mu[j,tid]
        # print(i, "tail", xprev, yprev, srcxprev, srcyprev, muprev)

        # xpost, ypost, srcxpost, srcypost, mupost = allSolutions_x[j2,hid2], allSolutions_y[j2,hid2], allSolutions_srcx[j2,hid2], allSolutions_srcy[j2,hid2], allSolutions_mu[j,hid2] **** --> allSolutions_mu[******j2******,hid2]
        # print(i2, "head", xpost, ypost, srcxpost, srcypost, mupost)
        # xpost, ypost, srcxpost, srcypost, mupost = allSolutions_x[j2,tid2], allSolutions_y[j2,tid2], allSolutions_srcx[j2,tid2], allSolutions_srcy[j2,tid2], allSolutions_mu[j,tid2] **** --> allSolutions_mu[******j2******,hid2]
        # print(i2, "tail", xpost, ypost, srcxpost, srcypost, mupost)

        # 4 6 can not attach
        # 4 head 0.19764838415833735 0.9440370274343007 -0.029747426672208 -0.034747426672208 -33.00279379640709
        # 4 tail 0.458212046074449 0.8396768164456602 -0.030895595065191295 -0.03155941894461588 -1283.8566291520983
        # 6 head 0.7206507298965913 0.6120205177529883 -0.029747426672208 -0.034747426672208 -33.00279379640709
        # 6 tail 0.47148464005208884 0.8317471393445778 -0.030895595065191295 -0.03155941894461588 -1283.8566291520983
        # f1 = trueSolution(pyTRIL.mlens, pyTRIL.zlens, -0.030895595065191295, -0.03155941894461588, [0.458212046074449, 0.8396768164456602], cal_ang = False) 
        # f2 = trueSolution(pyTRIL.mlens, pyTRIL.zlens, -0.030895595065191295, -0.03155941894461588, [0.47148464005208884, 0.8317471393445778], cal_ang = False) 
        # print("f1", f1)
        # print("f2", f2)
        # f1 [1, -1283.8566291520983, 1.0775462328294101e-12, 0, 0, 0]
        # f2 [1, 1276.8817477735402, 1.2045701799888764e-12, 0, 0, 0] --- 6 tail has wrong magnification




















