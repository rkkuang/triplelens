
    def if_two_segments_connect(self, head1, tail1, head2, tail2, allSolutions_x, allSolutions_y, allSolutions_srcx, allSolutions_srcy, allSolutions_mu):
        # head1 = [head1_j, head1_idx]

        # connect seg1 and seg2
        # itype 0 = can not connect
        # itype 1 = head head --> how to handle this? --> you need to insert seg2 into the front of seg1, and reverse seg2
        # itype 2 = head tail --> how to handle this? --> you need to insert seg2 into the front of seg1, do not need to reverse seg2
        # itype 3 = tail head --> how to handle this? --> the easiest scenerio, just add seg2 behind seg1
        # itype 4 = tail tail --> how to handle this? --> add seg2 behind seg1, and reverse seg2

        h1_j, h1_idx = head1
        t1_j, t1_idx = tail1
        h2_j, h2_idx = head2
        t2_j, t2_idx = tail2

        itype = 0
        # we prefer connect, rather than jump, so, first test image connectivity; if we cannot connect, then try whether we can "jump"
        # besides, you need to have a seperate function for jump

        # we prefer to find the easiest scenario
        xprev, yprev, srcxprev, srcyprev, muprev = allSolutions_x[t1_j, t1_idx], allSolutions_y[t1_j, t1_idx], allSolutions_srcx[t1_j, t1_idx], allSolutions_srcy[t1_j, t1_idx], allSolutions_mu[t1_j, t1_idx]

        xpost, ypost, srcxpost, srcypost, mupost = allSolutions_x[h2_j, h2_idx], allSolutions_y[h2_j, h2_idx], allSolutions_srcx[h2_j, h2_idx], allSolutions_srcy[h2_j, h2_idx], allSolutions_mu[h2_j, h2_idx]
        if self.if_dis_close( xprev, yprev, xpost, ypost, EPS**2):
            itype = 3 # tail connect with head
            print("tail connect with head", itype, xprev, yprev, xpost, ypost)
        # elif self.if_jump(srcxprev, srcyprev, srcxpost, srcypost, muprev, mupost, EPS**2):
        #     itype = 3 
        #     print(i, i2, "tail jump to head", srcxprev, srcyprev, srcxpost, srcypost, muprev, mupost)

        xpost, ypost, srcxpost, srcypost, mupost = allSolutions_x[t2_j, t2_idx], allSolutions_y[t2_j, t2_idx], allSolutions_srcx[t2_j, t2_idx], allSolutions_srcy[t2_j, t2_idx], allSolutions_mu[t2_j, t2_idx]
        if self.if_dis_close( xprev, yprev, xpost, ypost, EPS**2):
            itype = 4 # tail connect with tail
            print("tail connect with tail", itype, xprev, yprev, xpost, ypost)
        # elif self.if_jump(srcxprev, srcyprev, srcxpost, srcypost, muprev, mupost, EPS**2):
        #     itype = 4 # tail connect with tail
        #     print(i, i2, "tail jump to tail", srcxprev, srcyprev, srcxpost, srcypost, muprev, mupost)

        xprev, yprev, srcxprev, srcyprev, muprev = allSolutions_x[h1_j, h1_idx], allSolutions_y[h1_j, h1_idx], allSolutions_srcx[h1_j, h1_idx], allSolutions_srcy[h1_j, h1_idx], allSolutions_mu[h1_j, h1_idx]
        xpost, ypost, srcxpost, srcypost, mupost = allSolutions_x[t2_j, t2_idx], allSolutions_y[t2_j, t2_idx], allSolutions_srcx[t2_j, t2_idx], allSolutions_srcy[t2_j, t2_idx], allSolutions_mu[t2_j, t2_idx]
        if self.if_dis_close( xprev, yprev, xpost, ypost, EPS**2):
            itype = 2 # head connect with tail
            print("head connect with tail", itype, xprev, yprev, xpost, ypost)
        # elif self.if_jump(srcxprev, srcyprev, srcxpost, srcypost, muprev, mupost, EPS**2):
        #     itype = 2 # head connect with tail
        #     print(i, i2, "head jump to tail", srcxprev, srcyprev, srcxpost, srcypost, muprev, mupost)

        xpost, ypost, srcxpost, srcypost, mupost = allSolutions_x[h2_j, h2_idx], allSolutions_y[h2_j, h2_idx], allSolutions_srcx[h2_j, h2_idx], allSolutions_srcy[h2_j, h2_idx], allSolutions_mu[h2_j, h2_idx]
        if self.if_dis_close( xprev, yprev, xpost, ypost, EPS**2):
            print("head connect with head", itype, xprev, yprev, xpost, ypost)
            itype = 1 # head connect with head
        # elif self.if_jump(srcxprev, srcyprev, srcxpost, srcypost, muprev, mupost, EPS**2):
        #     itype = 1 # head connect with head
        #     print(i, i2, "head jump to head", srcxprev, srcyprev, srcxpost, srcypost, muprev, mupost)


        if itype == 0:
            pass
            print("can not attach")

            xprev, yprev, srcxprev, srcyprev, muprev = allSolutions_x[h1_j, h1_idx], allSolutions_y[h1_j, h1_idx], allSolutions_srcx[h1_j, h1_idx], allSolutions_srcy[h1_j, h1_idx], allSolutions_mu[h1_j, h1_idx]
            print("head", xprev, yprev, srcxprev, srcyprev, muprev)
            xprev, yprev, srcxprev, srcyprev, muprev = allSolutions_x[t1_j, t1_idx], allSolutions_y[t1_j, t1_idx], allSolutions_srcx[t1_j, t1_idx], allSolutions_srcy[t1_j, t1_idx], allSolutions_mu[t1_j, t1_idx]
            print("tail", xprev, yprev, srcxprev, srcyprev, muprev)

            xpost, ypost, srcxpost, srcypost, mupost = allSolutions_x[h2_j, h2_idx], allSolutions_y[h2_j, h2_idx], allSolutions_srcx[h2_j, h2_idx], allSolutions_srcy[h2_j, h2_idx], allSolutions_mu[h2_j, h2_idx]
            print("head", xpost, ypost, srcxpost, srcypost, mupost)
            xpost, ypost, srcxpost, srcypost, mupost = allSolutions_x[t2_j, t2_idx], allSolutions_y[t2_j, t2_idx], allSolutions_srcx[t2_j, t2_idx], allSolutions_srcy[t2_j, t2_idx], allSolutions_mu[t2_j, t2_idx]
            print("tail", xpost, ypost, srcxpost, srcypost, mupost)


        return itype

