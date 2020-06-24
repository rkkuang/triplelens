// #include<stdio.h>
#include "VBBinaryLensingLibrary.h"

// using namespace nc;
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

// void multiply(complex a[], int na, complex b[], int nb, complex c[]);
// void multiply_z(complex c[], complex a, int n);
// void multiply_zSquare(complex c[], complex a, int n);
// void polynomialCoefficients(double mlens[], complex zlen[], double xs, double ys, complex c[]);




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
_linkedarray *getphis(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, double causticsx[], double causticsy[], int secnum, int basenum, int *distype, int *maxmuidx, int numcritical);
_linkedarray *getphis_v2(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int secnum, int basenum, int *distype, int *maxmuidx);
_linkedarray *getphis_v3(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int secnum, int basenum, int *distype, int *maxmuidx);
_linkedarray *judge_distance(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, double causticsx[], double causticsy[], int secnum, int basenum, int *distype, int *maxmuidx, int numcritical);
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
_sols *outputTracks_v2(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, _linkedarray *PHI);
_sols *outputTracks_v2_savehalf(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, _linkedarray *PHI, _sols **prevstore,int firsttime, unsigned short int *ifjump);
// _sols *outputTracks_mod(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs,int nphi, double phi0);
void lensEquation(double mlens[], complex zlens[], complex z, complex *zs, double *mu);
double areaFunc_v2(_sols *track, double *E1, double *E2, double *E3, double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, _linkedarray *PHI);
double areaFunc_v3(_sols *track, double *E1, double *E2, double *E3, double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, _linkedarray *PHI);
double areaFunc_v4(_sols *track, double *E1, double *E2, double *E3, double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, _linkedarray *PHI, double errTol, int NPS);
double areaFunc(_sols *track, double rs, int paracorr=0);

void outpureImgPoints(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi);

double gould(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, double Gamma);
double gould_v2(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, double Gamma, int secd_or4th, double *correction, double *mups);
void findCloseImages(double mlens[], complex zlens[], double xs, double ys, complex z[]);
double triplePS(double mlens[], complex zlens[], double xs, double ys, int *nimages) ;
double tripleQuatrapoleTest(double mlens[], complex zlens[], double xs, double ys, int *nimages, double rs, double *quaderr);
double gould(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, double Gamma);
double tripleFS(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, double phi0, int *finalnphi);
double tripleFS_v2(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int *finalnphi, double *causticsx, double *causticsy, int secnum, int basenum, int numcritical);
double tripleFS_v2_gould(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int *finalnphi, double *causticsx, double *causticsy, int secnum, int basenum, int numcritical);
double tripleFS_v2_gould_error(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int *finalnphi, double *causticsx, double *causticsy, int secnum, int basenum, int numcritical);
double tripleFS_combine(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int *finalnphi, int nphi0, double phi0);
double tripleFS_v2_savehalf(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int *finalnphi, double *causticsx, double *causticsy, int secnum, int basenum, int numcritical);
double tripleFS_v2_savehalf_gould(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int *finalnphi, int secnum, int basenum, double *quad_err, double quaderr_Tol = 1e-2);
double tripleFS2python(double mlens[], double zlens[], double xsCenter, double ysCenter, double rs, int nphi, int secnum, int basenum, double quaderr_Tol = 1e-2);//core
double tripleFS_v3(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int *finalnphi, int secnum, int basenum, int numcritical);

double myatan(double x, double y);
void get_crit_caus(double mlens[], complex zlens[], int nlens, int NPS, double *criticalx, double *criticaly, double *causticsx, double *causticsy, int *numcritical);
void outputCriticalTriple_list(double allxys[], double mlens[], double zlens[], int nlens, int NPS);

_curve *NewImagesNlens(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs , _theta *, complex  *coefficients);
complex _LL_Nlens(complex y, complex z, double mlens[], complex zlens[]);
void OrderImagesNlens(_sols *Sols, _curve *Newpts);
double Mag_Nlens(double mlens[], complex zlens[], double xsCenter, double ysCenter, double RSv, _sols **Images);
double Mag_Nlensv2(double mlens[], complex zlens[], double xsCenter, double ysCenter, double RSv, _sols **Images);
void tri_ada_sample_np_20200606( double mlens[], double zlens[], double rs, int nphi, int secnum, int basenum, double quaderr_Tol, _curve *lightkv, int np0, int MAXNPS, double errTol, double t_start, double t_end, double alpha, double tE, double t0 , double u0, int np, double timerrTol);
void bisec(_curve *lightkv,_point *scan2,_point *scan3,double errTol,double errorTol_critical, double tE_inv, double u0, double salpha, double calpha, double t0,double mlens [], double zlens [],double rs,int nphi,int secnum,int basenum,double quaderr_Tol, int MAXNPS, double timerrTol);