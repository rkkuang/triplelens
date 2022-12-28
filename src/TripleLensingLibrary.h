// #include<stdio.h>
#include<memory>
#include<math.h>

#define EPS 1.0e-5 // segment close threshold, 1.0e-5 is ok
#define EPS2 1.0e-10 //(EPS*EPS)
#define SOLEPS 1.0e-5 // true or false solution of lens equation solving, 1.0e-5 is a bit too strict at some cases

#define EPS_CLOSE 1.0e-10   // how accurate the root should be
#define MAXIT_NEWTON 50     // maximum steps in the Newton-Raphson method

#define MAXARRLEN 5120 // 512000 = 1000 * 2**9 
#define NLENSFIX 3
#define DEGREEFIX (NLENSFIX*NLENSFIX+1)
#define CORRECT_PARITY (1-NLENSFIX)

#define MR 8
#define MT 10
#define MAXIT (MT*MR)
#define MAXM 30
#define FRAC_ERR 2.0e-17 //Fractional Error for double precision, 2.0e-15


class complex {
public:
    double re;
    double im;
    complex(double, double);
    complex(double);
    complex(void);
};

class TripleLensing {
    int finalnphi, nimages, ftime, FSflag, nsolution, trackcnt, parity, degrees, adanp0 = 2, MAXNPS = 2000, secnum_priv, basenum_priv, clpairparity1, clpairparity2, special_flag = 0;
    double TINY, ph1, ph2, ph3, adaerrTol = 0.0005, timerrTol = 1e-3, tempdis, mindis, M_PI2 = 2 * M_PI, absdzs;
    // absdzs is used to store the abs(dzs) in trueSolution()
    double r2_1, r2_2, x, y, dx_db, dy_db, dx_db2, dy_db2 , Jxx, Jyy, Jxy, rho2, areaSource, phi0, x1, x2, y1, y2, phi1, phi2, dphi, dphi3, ds1, ds2, subArea, pos_area, rep, imp, x_xj, y_yj, xy_j2, relerr_priv; // trueSolution
    complex zs, dzs, dzsdz, zsc, z1, z2, z3, z1bar, z2bar, z3bar, zsbar; // zc[NLENS]
    complex tempzc, tempzc2, tempzc3, tempzc4, J1c, dz1, dz2;
    double therr;
    int NPS;
    double *mlens, *Zlens;


    void cmplx_laguerre(complex *, int, complex *, int &, bool &);
    void cmplx_newton_spec(complex *, int, complex *, int &, bool &);
    void solve_quadratic_eq(complex &, complex &, complex *);
    void solve_cubic_eq(complex &, complex &, complex &, complex *);

public:
    int nphi, secnum, maxmuidx,  pntnum, CQ, finalNPS;
    unsigned short int NLENS, DEGREE, ifjump, flag, basenum, distype, ifFinite, area_quality = 1;
    complex **p, **q, **p_const, *temp, **temp_const1, ***temp_const2, *temp_const22, *ctemp, *qtemp, *qtemp2, *ptemp, *zr, *coefficients, *zc;

    double quad_err, quaderr_Tol, area, mu0, mu, relerr_mag, lambda1, lambda2, thetaJ, eacherr, muTotal, xs, ys, phi, SD, MD, CD, muPS, ds, dJ2, RelTolLimb = 1e-3, AbsTolLimb = 1e-4;
    
    complex *zlens;
    complex J1, J2, J3, dJ, firstterm, secondterm, thirdterm, J1c2, dJ5, dy, dz;


    TripleLensing();
    ~TripleLensing();
    TripleLensing(double mlens[], complex zlens[]);
    TripleLensing(double mlens[], double Zlens[]);

    void setnlens(int nlens);
    void reset3(double mlens[], complex zlens[]);
    void reset2(double mlens[],  double Zlens[]);
    double TriplePS(double xs, double ys);
    double gould(double xsCenter, double ysCenter, double rs, double Gamma, double *A2rho2, double *A4rho4);
    double gould(double xsCenter, double ysCenter, double rs, double Gamma);


    double tripleFS2python(double mlens[], double zlens[], double xsCenter, double ysCenter, double rs);//core
    double tripleFS2python(double xsCenter, double ysCenter, double rs);//core

    void solv_lens_equation(double *zrxy, double *mlens, double *zlens, double xs, double ys, int nlens);

    double tripleQuatrapoleTest(double xs, double ys, double rs);


    int trueSolution(double xs, double ys, complex z, double *mu);
    int trueSolution_qtest(double xs, double ys, complex z);
    int trueSolution_withmu(double xs, double ys, complex z, double true_solution_threshold, double *mu);

    void outImgPoints(double xsCenter, double ysCenter, double rs, int nphi);
    void polynomialCoefficients(double xs, double ys, complex c[]);

    void tripleFS2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, int secnum, int basenum, double quaderr_Tol, double relerr_mag, double mags[], int Np) ; //core
    void tripleGould2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, double Gamma, double mags[], int Np);
    void triple_num_real_sol2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double true_solution_threshold , double numbers_mups[], int Np);
    void outsys(double mlens[], complex zlens[], double t0, double u0, double tE, double s2, double q2, double alpha, double s3, double q3, double psi, double rs, double xsCenter, double ysCenter);


    // ########################################################################################
    // variables and functions related to magnification calculations with array structure, begin
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
    // variables and functions related to magnification calculations with array structure, end
    // ########################################################################################

    // Skowron & Gould root calculator
    void cmplx_roots_gen(complex *, complex *, int, bool, bool);
    void cmplx_laguerre2newton(complex *, int, complex *, int &, bool &, int);
};


void multiply(complex a[], int na, complex b[], int nb, complex c[]);
void multiply_z(complex c[], complex a, int n);
// void multiply_z_v2(complex c[][NLENS + 1], complex a, int n, int firstdim);
void multiply_z_v2(complex **c, complex a, int n, int firstdim);
void multiply_zSquare(complex c[], complex a, int n);


int one_parity(double mu);
void myargsort(int arr[], int index[], int n, int mode);
void myswap(int *a, int *b);
void mydbswap(double *a, double *b);

double abs(complex);
double absint(int);
double absf(double);
complex conj(complex);
complex sqrt(complex);
double real(complex);
double imag(complex);
complex expcmplx(complex);
complex cbrt(complex);
complex operator+(complex, complex);
complex operator-(complex, complex);
complex operator*(complex, complex);
complex operator/(complex, complex);
complex operator+(complex, double);
complex operator-(complex, double);
complex operator*(complex, double);
complex operator/(complex, double);
complex operator+(double, complex);
complex operator-(double, complex);
complex operator*(double, complex);
complex operator/(double, complex);
complex operator+(int, complex);
complex operator-(int, complex);
complex operator*(int, complex);
complex operator/(int, complex);
complex operator+(complex, int);
complex operator-(complex, int);
complex operator*(complex, int);
complex operator/(complex, int);
complex operator-(complex);
bool operator==(complex, complex);
bool operator!=(complex, complex);