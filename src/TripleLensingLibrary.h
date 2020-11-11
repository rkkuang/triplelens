// #include<stdio.h>
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

#define NLENS 3
#define DEGREE (NLENS*NLENS+1)
#define EPS 1.0e-5
#define one_24 (double)1.0/24
// #define parabcorr


#define EPS_CLOSE 1.0e-10   // how accurate the root should be
#define MAXIT_NEWTON 50     // maximum steps in the Newton-Raphson method

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
    int finalnphi, nimages, ftime, FSflag, nsolution, trackcnt, parity, degrees, adanp0 = 2, MAXNPS = 2000, adanp = 1, secnum_priv, clpairparity1, clpairparity2;
    double TINY, ph1, ph2, ph3, adaerrTol = 0.0005, timerrTol = 1e-5, tempdis, mindis, M_PI2 = 2 * M_PI;
    double r2_1, r2_2, x, y, dx_db, dy_db, dx_db2, dy_db2 , Jxx, Jyy, Jxy, Jyx, rho2, areaSource, phi0, x1, x2, y1, y2, phi1, phi2, dphi, dphi3, ds1, ds2, subArea, pos_area, rep, imp, x_xj, y_yj, xy_j2, m13, m23, m33, m12, m22, m32, m1, m2, m3, relerr_priv; // trueSolution
    complex zs, dzs, dzsdz, zsc, zc[NLENS], z1, z2, z3, z1bar, z2bar, z3bar, zsbar, z13, z23, z33, z12, z22, z32, zsmod, z1barz2barzs, z1barz3barzs, z2barz3barzs, z1barzs, z2barzs, z3barzs, z1barzsmod, z2barzsmod, z3barzsmod, z1barzsbar, z2barzsbar, z3barzsbar, z1barz2barzsmod, z1barz3barzsmod, z2barz3barzsmod, z1barz2barzsbar, z1barz3barzsbar, z2barz3barzsbar;
    complex tempzc, tempzc2, tempzc3, tempzc4, J1c, dz1, dz2;
    complex p[NLENS + 1][NLENS], q[NLENS + 1][NLENS], p_const[NLENS + 1][NLENS];
    complex temp[NLENS + 1], temp_const1[NLENS + 1][NLENS + 1], temp_const2[NLENS + 1][NLENS + 1][NLENS + 1], temp_const22[NLENS + 1];
    complex ctemp[DEGREE + 1];
    complex qtemp[DEGREE + 1], qtemp2[NLENS + 1];
    complex ptemp[DEGREE + 1];
    _point *pscan;
    double corrquad2, corrquad, therr;
    int NPS;
    double *mlens, *Zlens;
public:
    int nphi, secnum, basenum, distype, maxmuidx, flag,  pntnum, CQ, finalNPS, ifFinite;
    double quad_err, quaderr_Tol, area, mu0, mu, relerr_mag, lambda1, lambda2, thetaJ, eacherr, muTotal, xs, ys, phi, SD, MD, CD, muPS, ds, dJ2, RelTolLimb = 0, AbsTolLimb = 0;
    unsigned short int ifjump;
    complex *zlens;
    complex zr[DEGREE];
    complex coefficients[DEGREE + 1];
    complex J1, J2, J3, dJ, firstterm, secondterm, thirdterm, J1c2, dJ5, dy, dz;


    TripleLensing();
    TripleLensing(double mlens[], complex zlens[]);
    TripleLensing(double mlens[], double Zlens[]);
    // double TripleMag(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs);
    // double TriplePS(double mlens[], complex zlens[], double xsCenter, double ysCenter);
    void reset(double mlens[], complex zlens[]);
    void reset3(double mlens[], complex zlens[]);
    void reset2(double mlens[],  double Zlens[]);
    double TriplePS(double xs, double ys);
    double gould(double xsCenter, double ysCenter, double rs, double Gamma);


    double tripleFS2python(double mlens[], double zlens[], double xsCenter, double ysCenter, double rs);//core
    double tripleFS2python(double xsCenter, double ysCenter, double rs);//core

    void solv_lens_equation(double *zrxy, double *mlens, double *zlens, double xs, double ys, int nlens);
    void outputCriticalTriple_list(double resxy[], double mlens[], double Zlens[], int nlens, int NPS);
    void TriLightCurve(double *pr, double *mags, double *y1s, double *y2s, int np);

    // double tripleFS_v2_savehalf_quadtest(double xsCenter, double ysCenter, double rs);
    double TripleMag(double xsCenter, double ysCenter, double rs);
    double TripleMagFull(double xsCenter, double ysCenter, double rs);
    double TripleMagDark(double xsCenter, double ysCenter, double RSv, double a1);

    double tripleQuatrapoleTest(double xs, double ys, double rs);
    _linkedarray *getphis_v3(double xsCenter, double ysCenter, double rs);
    _sols *outputTracks_v2_savehalf(double xsCenter, double ysCenter, double rs, _linkedarray *PHI, _sols **prevstore);

    // void bisec(_curve *lightkv, _point *scan2, _point *scan3, double errTol, double errorTol_critical, double tE_inv, double u0, double salpha, double calpha, double t0, double mlens [], double zlens [], double rs, int nphi, int secnum, int basenum, double quaderr_Tol, int MAXNPS, double timerrTol);

    void bisec(_curve *lightkv, _point *scan2, _point *scan3, double errorTol_critical, double tE_inv, double u0, double salpha, double calpha, double t0, double rs);


    // void tri_ada_sample_1D( double mlens[], double zlens[], double rs, int nphi, int secnum, int basenum, double quaderr_Tol, _curve *lightkv, int np0, int MAXNPS, double errTol, double t_start, double t_end, double alpha, double tE, double t0 , double u0, int np, double timerrTol);

    void TripleLkvAdap(double rs, _curve *lightkv, double t_start, double t_end, double alpha, double tE, double t0 , double u0);


    int trueSolution(double xs, double ys, complex z);
    int trueSolution_qtest(double xs, double ys, complex z);

    void outImgPoints(double xsCenter, double ysCenter, double rs, int nphi);
    _sols *outputTracks(double xsCenter, double ysCenter, double rs, int nphi, double phi0);

    void outputImagesTriple(double xsCenter, double ysCenter, int nphi, double rs);


    double areaFunc(_sols *track, double rs);
    double areaFunc_parab(_sols *track, double rs);

    // void polynomialCoefficients(double mlens[], complex zlen[], double xs, double ys, complex c[], int nlens, int degree);
    void polynomialCoefficients(double xs, double ys, complex c[]);
    // void TripleCoefs(double xs, double ys, complex *coefs);
    // void BinaryCoefs(double xs, double ys, complex *coefs);

    void tripleFS2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, int secnum, int basenum, double quaderr_Tol, double relerr_mag, double mags[], int Np) ; //core
    void tripleFSLimb2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, int secnum, int basenum, double quaderr_Tol, double relerr_mag, double mags[], int Np, double limba1, double RelTolLimb, double AbsTolLimb) ; //core


};



void addPoint(_curve *final, _point *p);
_curve *newSegment(_sols *imageTracks);



// _linkedarray *getphis_v3(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int secnum, int basenum, int *distype, int *maxmuidx);

// _linkedarray *thetas2phis(_thetas * thetas);


_curve *jump_over_caustics(_curve *final, _sols *imageTracks, int *head, int *tail, _curve *markedCurve, bool checkfinalhead = false);
_curve *connect_head_or_tail(int *ifcontinue, _curve *final , _sols *imageTracks, _sols *connected, bool checkfinalhead = false);
_point *pnt_reach_head_or_tail(_point * p, _curve * Prov, _sols * imageTracks, int *reach);

void printOneTrack(_curve *curve);
void printAllTracks(_sols *track);
void saveTracks(_sols *track);
void saveTracks_before(_sols *track);
void outputCriticalTriple(double mlens[], complex zlens[], int nlens, int npnt = 3000);
void outputCaustics_for_crossSection(double mlens[], complex zlens[], int nlens, int npnt);

void outputCriticalTriple_list(double allxys[], double mlens[], double zlens[], int nlens, int NPS);
void outsys(double mlens[], complex zlens[], double t0, double u0, double tE, double s2, double q2, double alpha, double s3, double q3, double psi, double rs, double xsCenter, double ysCenter);

// void outputImagesTriple(double mlens[], complex zlens[], double xsCenter, double ysCenter, int nphi, double rs);

// int trueSolution(double mlens[], complex zlens[], double xs, double ys, complex z, double *mu, double *lambda1, double *lambda2, double *thetaJ, int nlens, complex *J1, complex *J2, complex *dJ, complex *J3);

// _sols *outputTracks(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, double phi0);

// _sols *outputTracks_v2_savehalf(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, _linkedarray *PHI, _sols **prevstore, int firsttime, unsigned short int *ifjump);

// void lensEquation(double mlens[], complex zlens[], complex z, complex *zs, double *mu);

// double areaFunc_v4(_sols *track, double *E1, double *E2, double *E3, double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, _linkedarray *PHI, double errTol, int NPS);

// double areaFunc(_sols *track, double rs);

// void outImgPoints(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi);

// double triplePS(double mlens[], complex zlens[], double xs, double ys, int *nimages) ;

// double tripleQuatrapoleTest(double mlens[], complex zlens[], double xs, double ys, int *nimages, double rs, double *quaderr);

// double tripleFS2python(double mlens[], double zlens[], double xsCenter, double ysCenter, double rs, int nphi, int secnum, int basenum, double quaderr_Tol = 1e-2);//core


void get_crit_caus(double mlens[], complex zlens[], int nlens, int NPS, double *criticalx, double *criticaly, double *causticsx, double *causticsy, int *numcritical);

// complex _LL_Nlens(complex y, complex z, double mlens[], complex zlens[]);


void multiply(complex a[], int na, complex b[], int nb, complex c[]);
void multiply_z(complex c[], complex a, int n);
void multiply_z_v2(complex c[][NLENS + 1], complex a, int n, int firstdim);
void multiply_zSquare(complex c[], complex a, int n);
// void polynomialCoefficients(double mlens[], complex zlen[], double xs, double ys, complex c[], int nlens, int degree);

// void tri_ada_sample_1D( double mlens[], double zlens[], double rs, int nphi, int secnum, int basenum, double quaderr_Tol, _curve *lightkv, int np0, int MAXNPS, double errTol, double t_start, double t_end, double alpha, double tE, double t0 , double u0, int np, double timerrTol);
// void bisec(_curve *lightkv, _point *scan2, _point *scan3, double errTol, double errorTol_critical, double tE_inv, double u0, double salpha, double calpha, double t0, double mlens [], double zlens [], double rs, int nphi, int secnum, int basenum, double quaderr_Tol, int MAXNPS, double timerrTol);