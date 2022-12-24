// #include<stdio.h>
#include<memory>
#include<math.h>
#include "VBBinaryLensingLibrary.h"

// output an image
#define TRUE_IMAGE 1
#define FALSE_IMAGE -1
// #define VERBOSE
#define _PRINT_ERRORS2

// #define verbose
// #define combine_verbose
#define combine_verbose_v 0

// #define CALCULATE_ANGLE

// #define NLENS 3 // NLENS and DEGREE should not be fixed ?
// #define DEGREE (NLENS*NLENS+1)
#define EPS 1.0e-5 // segment close threshold, 1.0e-5 is ok
#define EPS2 1.0e-10 //(EPS*EPS)

#define SOLEPS 1.0e-5 // true or false solution of lens equation solving, 1.0e-5 is a bit too strict at some cases
// #define SOLEPS1e2 1e-3 // threshold used to judge whether we use Newton method to polish the root using the original lens equation instead of the polynomial
// #define SOLEPS1e3 1e-2 // threshold used to judge whether we should taken the result of Newton polish method 


// 2021.10.07 to do: perhaps the timing when NEWTON POLISH used should be changed?
#define POLISH_USE_NEWTON // used for small mass ratio
#define SOLEPS1e2 1e-4 // if we find one solution has absdsz > SOLEPS1e2, we use NEWTON polish
#define SOLEPS1e3 1e-4// be cautious, newton polish may lead to a dilemma that two different solutions converge to the same position, this will lead to bugs during the track-connect process



#define MISS_SOL_THRESHOLD 1

#define JUMP_SEARCH_DEPTH 5 // search depth for the situation that we can not jump

#define one_24 (double)1.0/24


// #define parabcorr

#define EPS_CLOSE 1.0e-10   // how accurate the root should be
#define MAXIT_NEWTON 50     // maximum steps in the Newton-Raphson method


#define MAXARRLEN 5120 // 512000 = 1000 * 2**9 
#define NLENSFIX 3
#define DEGREEFIX (NLENSFIX*NLENSFIX+1)
#define CORRECT_PARITY (1-NLENSFIX)

class _linkedarray {
public:
    int length;
    Node *first;
    Node *last;
    _linkedarray(void);
    ~_linkedarray(void);

    void concatenate(_linkedarray *);
    void append(double);
    void linspace(double phi0, double phiend, int nphi,  int endpoint);
    void print();
    void extend();
    void append_before(Node *place, double x);
    void append_behind(Node *place, double x);
    // Node *insert(double th);
};


class TripleLensing {
    int finalnphi, nimages, ftime, FSflag, nsolution, trackcnt, parity, degrees, adanp0 = 2, MAXNPS = 2000, secnum_priv, basenum_priv, clpairparity1, clpairparity2, special_flag = 0;
    double TINY, ph1, ph2, ph3, adaerrTol = 0.0005, timerrTol = 1e-3, tempdis, mindis, M_PI2 = 2 * M_PI, absdzs;
    // absdzs is used to store the abs(dzs) in trueSolution()
    double r2_1, r2_2, x, y, dx_db, dy_db, dx_db2, dy_db2 , Jxx, Jyy, Jxy, rho2, areaSource, phi0, x1, x2, y1, y2, phi1, phi2, dphi, dphi3, ds1, ds2, subArea, pos_area, rep, imp, x_xj, y_yj, xy_j2, relerr_priv; // trueSolution
    complex zs, dzs, dzsdz, zsc, z1, z2, z3, z1bar, z2bar, z3bar, zsbar; // zc[NLENS]
    complex tempzc, tempzc2, tempzc3, tempzc4, J1c, dz1, dz2;
    // complex p[NLENS + 1][NLENS], q[NLENS + 1][NLENS], p_const[NLENS + 1][NLENS];
    // complex temp[NLENS + 1], temp_const1[NLENS + 1][NLENS + 1], temp_const2[NLENS + 1][NLENS + 1][NLENS + 1], temp_const22[NLENS + 1];
    // complex ctemp[DEGREE + 1];
    // complex qtemp[DEGREE + 1], qtemp2[NLENS + 1];
    // complex ptemp[DEGREE + 1];
    _point *pscan;
    double therr;
    int NPS;
    double *mlens, *Zlens;
public:
    int nphi, secnum, maxmuidx,  pntnum, CQ, finalNPS;
    unsigned short int NLENS, DEGREE, ifjump, flag, basenum, distype, ifFinite, area_quality = 1;
    // area_quality = 1 means this magnification is ok, otherwise (mainly caused by insufficient samples around the source edge)
    // if area_quality == 0, mean parity in areaFunc might be wrong
    // if area_quality == 2, means return from a looser threshold
    // if area_quality == 3, means return muPS

    // complex p[NLENS + 1][NLENS], q[NLENS + 1][NLENS], p_const[NLENS + 1][NLENS], temp[NLENS + 1], temp_const1[NLENS + 1][NLENS + 1], temp_const2[NLENS + 1][NLENS + 1][NLENS + 1], temp_const22[NLENS + 1], ctemp[DEGREE + 1], qtemp[DEGREE + 1], qtemp2[NLENS + 1], ptemp[DEGREE + 1], zr[DEGREE], coefficients[DEGREE + 1];
    complex **p, **q, **p_const, *temp, **temp_const1, ***temp_const2, *temp_const22, *ctemp, *qtemp, *qtemp2, *ptemp, *zr, *coefficients, *zc;


    double quad_err, quaderr_Tol, area, mu0, mu, relerr_mag, lambda1, lambda2, thetaJ, eacherr, muTotal, xs, ys, phi, SD, MD, CD, muPS, ds, dJ2, RelTolLimb = 1e-3, AbsTolLimb = 1e-4;
    
    complex *zlens;
    // complex zr[DEGREE];
    // complex coefficients[DEGREE + 1];
    complex J1, J2, J3, dJ, firstterm, secondterm, thirdterm, J1c2, dJ5, dy, dz;


    TripleLensing();
    ~TripleLensing();
    TripleLensing(double mlens[], complex zlens[]);
    TripleLensing(double mlens[], double Zlens[]);

    void setnlens(int nlens);
    void reset(double mlens[], complex zlens[]);
    void reset3(double mlens[], complex zlens[]);
    void reset2(double mlens[],  double Zlens[]);
    double TriplePS(double xs, double ys);
    double gould(double xsCenter, double ysCenter, double rs, double Gamma, double *A2rho2, double *A4rho4);
    double gould(double xsCenter, double ysCenter, double rs, double Gamma);


    double tripleFS2python(double mlens[], double zlens[], double xsCenter, double ysCenter, double rs);//core
    double tripleFS2python(double xsCenter, double ysCenter, double rs);//core

    void solv_lens_equation(double *zrxy, double *mlens, double *zlens, double xs, double ys, int nlens);
    void outputCriticalTriple_list(double resxy[], double mlens[], double Zlens[], int nlens, int NPS);
    void outputCriticalBinary_list(double resxy[], double s, double q, int NPS);
    void TriLightCurve(double *pr, double *mags, double *y1s, double *y2s, int np);

    // double tripleFS_v2_savehalf_quadtest(double xsCenter, double ysCenter, double rs);
    double TripleMag(double xsCenter, double ysCenter, double rs);
    double TripleMagFull(double xsCenter, double ysCenter, double rs);
    double TripleMagDark(double xsCenter, double ysCenter, double RSv, double a1);

    double tripleQuatrapoleTest(double xs, double ys, double rs);
    _linkedarray *getphis_v3(double xsCenter, double ysCenter, double rs);
    _sols *outputTracks_v2_savehalf(double xsCenter, double ysCenter, double rs, _linkedarray *PHI, _sols **prevstore);

    void bisec(_curve *lightkv, _point *scan2, _point *scan3, double errorTol_critical, double tE_inv, double u0, double salpha, double calpha, double t0, double rs);
    double angcos(_point *p1, _point *p2, _point *p3, _point *p4);

    void TripleLkvAdap(double rs, _curve *lightkv, double t_start, double t_end, double alpha, double tE, double t0 , double u0);


    int trueSolution(double xs, double ys, complex z, double *mu);
    int trueSolution_qtest(double xs, double ys, complex z);
    int trueSolution_nomu(double xs, double ys, complex z, double true_solution_threshold);
    int trueSolution_withmu(double xs, double ys, complex z, double true_solution_threshold, double *mu);

    void outImgPoints(double xsCenter, double ysCenter, double rs, int nphi);
    _sols *outputTracks(double xsCenter, double ysCenter, double rs, int nphi, double phi0);

    void outputImagesTriple(double xsCenter, double ysCenter, int nphi, double rs);


    double areaFunc(_sols *track, double rs, int finalnphi, double muPS, double muPrev, int *quality);
    double areaFunc_parab(_sols *track, double rs);

    void polynomialCoefficients(double xs, double ys, complex c[]);

    void tripleFS2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, int secnum, int basenum, double quaderr_Tol, double relerr_mag, double mags[], int Np) ; //core
    void tripleGould2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, double Gamma, double mags[], int Np);
    void tripleFSLimb2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, int secnum, int basenum, double quaderr_Tol, double relerr_mag, double mags[], int Np, double limba1, double RelTolLimb, double AbsTolLimb) ; //core
    void triple_num_real_sol2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double true_solution_threshold , double numbers_mups[], int Np);
    void findCloseImages(double mlens[], complex zlens[], double xs, double ys, complex *z, bool *imageFound);
    void newtonStep(double mlens[], complex zlens[], double xs, double ys, complex z, complex *dz, double *mu);
    void outsys(double mlens[], complex zlens[], double t0, double u0, double tE, double s2, double q2, double alpha, double s3, double q3, double psi, double rs, double xsCenter, double ysCenter);


    double ARRPHI[MAXARRLEN], allSolutions_x[DEGREEFIX][MAXARRLEN], allSolutions_y[DEGREEFIX][MAXARRLEN], allSolutions_srcx[DEGREEFIX][MAXARRLEN], allSolutions_srcy[DEGREEFIX][MAXARRLEN], allSolutions_mu[DEGREEFIX][MAXARRLEN], allSolutions_absdzs[DEGREEFIX][MAXARRLEN], Prov_x[DEGREEFIX], Prov_y[DEGREEFIX], Prov_srcx[DEGREEFIX], Prov_srcy[DEGREEFIX], Prov_absdzs[DEGREEFIX], Prov_mu[DEGREEFIX], preProv_x[DEGREEFIX], preProv_y[DEGREEFIX], connectEPS, bestconnect_dis, connectdis, bestjump_fac_mu, bestjump_dis, how_close, jumpdis, srcxprev, srcyprev, srcxpost, srcypost, muprev, mupost, currx, curry;
    bool allSolutions_flag[DEGREEFIX][MAXARRLEN], has_removed[DEGREEFIX],already_done_segments[DEGREEFIX * DEGREEFIX], if_creat_new, canweconnect, ifcontinue, canwejump, ifjumpbool;
    unsigned short int Prov_flag[DEGREEFIX], attach_idx[DEGREEFIX], res_idx[DEGREEFIX], nclosed_image, ntrue_segments, open_seg_leftover, npure_close_segments, nfinal_closed_image, imgTrack_type[DEGREEFIX],true_segments_info[DEGREEFIX * DEGREEFIX][3], *temp_true_segments_info[DEGREEFIX*DEGREEFIX], bestconnect_type, bestconnect_i2, i2, j2, existing_seg_n, i3, i4;
    int posflagnums[DEGREEFIX], tmp_segments_length[DEGREEFIX*DEGREEFIX], hid, tid, hid2, tid2, head1[2], tail1[2], head2[2], tail2[2], itype, scan_true_idx, continue_left, sorted_lenidx[DEGREEFIX*DEGREEFIX], bestjumptype, bestjump_i2, closed_image_info[DEGREEFIX][2+3*(DEGREEFIX*DEGREEFIX+1)];

    // closed_image_info[DEGREEFIX][2+3*(DEGREEFIX*DEGREEFIX+1)] must be int type, because may save -1 parity
    // closed_image_info = np.zeros(( min(DEGREE, ntrue_segments), 2 + 3 * (ntrue_segments + 1) ));
    // already_done_segments = np.zeros(ntrue_segments)

    double arrTripleMag(double xsCenter, double ysCenter, double rs);
    void get_arrphi(double xsCenter, double ysCenter, double rs, int *retnphi);
    void arrlinspace(double *arr, double phi0, double phiend, int insertidx, int nphi, int endpoint);
    void arroutputTracks(double xsCenter, double ysCenter, double rs, bool prevstore, int nphi, double mindphi);//allSolutions_x, allSolutions_y, allSolutions_srcx, allSolutions_srcy, allSolutions_mu, allSolutions_flag, allSolutions_absdzs
    //return true_segments_info[:ntrue_segments,:].astype(int), closed_image_info[:nfinal_closed_image,:].astype(int)

    void extend_arrphi(int *retnphi);
    void saveimage(int curr_nimages, int curr_total_parity, int j, bool *ifsave, int *saveindex);
    void get_closest_idx(double preProv_x[], double preProv_y[], double Prov_x[], double Prov_y[], unsigned short int res_idx[]);
    void if_two_segments_connect(int head1[], int tail1[],int head2[], int tail2[], double EPS_2, int *itype, double *dis);

    void if_head_tail_jump(int head1[], int tail1[], bool *ifjumpbool, double *how_close, double *dis);
    void if_two_segments_jump(int head1[], int tail1[],int head2[], int tail2[], int *itype, double *bestjump_fac_mu, double *bestjump_dis);
    double head_tail_close(int head1[], int tail1[]);
    void arrareaFunc(double *area);
};



void addPoint(_curve *final, _point *p);
_curve *newSegment(_sols *imageTracks);


_curve *jump_over_caustics(_curve *final, _sols *imageTracks, int *head, int *tail, _curve *markedCurve, bool checkfinalhead, int *ifcontinueFromjump,  _sols * connected, int *ifkeep_jumping);
_curve *jump_over_caustics_deeper(_curve *final, _sols *imageTracks, int *head, int *tail, _curve *markedCurve, bool checkfinalhead, int *ifcontinueFromjump,  _sols * connected, int *ifkeep_jumping);
_curve *connect_head_or_tail(int *ifcontinue, _curve *final , _sols *imageTracks, _sols *connected, bool checkfinalhead, int *ifkeep_connect);
_point *pnt_reach_head_or_tail(_point * p, _curve * Prov, _sols * imageTracks, int *reach);

void printOneTrack(_curve *curve);
void printAllTracks(_sols *track);
void saveTracks(_sols *track, int nphi);
void saveTracks_before(_sols *track, int nphi);
void outputCriticalTriple(double mlens[], complex zlens[], int nlens, int npnt = 3000);
void outputCaustics_for_crossSection(double mlens[], complex zlens[], int nlens, int npnt);

void outputCriticalTriple_list(double allxys[], double mlens[], double zlens[], int nlens, int NPS);

void get_crit_caus(double mlens[], complex zlens[], int nlens, int NPS, double *criticalx, double *criticaly, double *causticsx, double *causticsy, int *numcritical);

void multiply(complex a[], int na, complex b[], int nb, complex c[]);
void multiply_z(complex c[], complex a, int n);
// void multiply_z_v2(complex c[][NLENS + 1], complex a, int n, int firstdim);
void multiply_z_v2(complex **c, complex a, int n, int firstdim);
void multiply_zSquare(complex c[], complex a, int n);


int one_parity(double mu);
void myargsort(int arr[], int index[], int n, int mode);
void myswap(int *a, int *b);
void mydbswap(double *a, double *b);


