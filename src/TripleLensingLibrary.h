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
#define parabcorr

#define EPS_CLOSE 1.0e-10   // how accurate the root should be
#define MAXIT_NEWTON 50     // maximum steps in the Newton-Raphson method

class TripleLensing {
    int finalnphi, nimages;
public:
    int nphi, secnum, basenum;
    double quad_err, quaderr_Tol;
    TripleLensing();
    double TripleMag(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs);
    double TriplePS(double mlens[], complex zlens[], double xsCenter, double ysCenter);
    void TripleLkvAdap(double mlens[], double zlens[], double rs, _curve *lightkv, int np0, int MAXNPS, double errTol, double t_start, double t_end, double alpha, double tE, double t0 , double u0, int np, double timerrTol);
    // double Zlens[2 * NLENS] = { zlens[0].re, zlens[0].im, zlens[1].re, zlens[1].im, zlens[2].re, zlens[2].im };
    double tripleFS2python(double mlens[], double zlens[], double xsCenter, double ysCenter, double rs);//core
};

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

void addPoint(_curve *final, _point *p);
_curve *newSegment(_sols *imageTracks);



_linkedarray *getphisNonfirst_caustic_crossing(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int secnum, int basenum, int maxmuidx);
_linkedarray *getphis_v3(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int secnum, int basenum, int *distype, int *maxmuidx);
_linkedarray *thetas2phis(_thetas * thetas);


_curve *jump_over_caustics(_curve *final, _sols *imageTracks, int *head, int *tail, _curve *markedCurve);
_curve *connect_head_or_tail(int *ifcontinue, _curve *final , _sols *imageTracks, _sols *connected);
_point *pnt_reach_head_or_tail(_point * p, _curve * Prov, _sols * imageTracks, int *reach);
void printOneTrack(_curve *curve);
void printAllTracks(_sols *track);
void outputCriticalTriple(double mlens[], complex zlens[], int nlens);
void outputImagesTriple(double mlens[], complex zlens[], double xsCenter, double ysCenter, int nphi, double rs);

int trueSolution(double mlens[], complex zlens[], double xs, double ys, complex z, double *mu, double *lambda1, double *lambda2, double *thetaJ, int nlens, complex *J1, complex *J2, complex *dJ, complex *J3);
_sols *outputTracks(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, double phi0);
_sols *outputTracks_v2_savehalf(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, _linkedarray *PHI, _sols **prevstore, int firsttime, unsigned short int *ifjump);
// _sols *outputTracks_mod(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs,int nphi, double phi0);
void lensEquation(double mlens[], complex zlens[], complex z, complex *zs, double *mu);
double areaFunc_v4(_sols *track, double *E1, double *E2, double *E3, double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, _linkedarray *PHI, double errTol, int NPS);
double areaFunc(_sols *track, double rs, int paracorr = 0);

void outImgPoints(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi);

double triplePS(double mlens[], complex zlens[], double xs, double ys, int *nimages) ;
double tripleQuatrapoleTest(double mlens[], complex zlens[], double xs, double ys, int *nimages, double rs, double *quaderr);
double tripleFS_v2_savehalf_quadtest(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int *finalnphi, int secnum, int basenum, double *quad_err, double quaderr_Tol = 1e-2);
double tripleFS2python(double mlens[], double zlens[], double xsCenter, double ysCenter, double rs, int nphi, int secnum, int basenum, double quaderr_Tol = 1e-2);//core
void get_crit_caus(double mlens[], complex zlens[], int nlens, int NPS, double *criticalx, double *criticaly, double *causticsx, double *causticsy, int *numcritical);
void outputCriticalTriple_list(double allxys[], double mlens[], double zlens[], int nlens, int NPS);
complex _LL_Nlens(complex y, complex z, double mlens[], complex zlens[]);
void tri_ada_sample_1D( double mlens[], double zlens[], double rs, int nphi, int secnum, int basenum, double quaderr_Tol, _curve *lightkv, int np0, int MAXNPS, double errTol, double t_start, double t_end, double alpha, double tE, double t0 , double u0, int np, double timerrTol);
void bisec(_curve *lightkv, _point *scan2, _point *scan3, double errTol, double errorTol_critical, double tE_inv, double u0, double salpha, double calpha, double t0, double mlens [], double zlens [], double rs, int nphi, int secnum, int basenum, double quaderr_Tol, int MAXNPS, double timerrTol);