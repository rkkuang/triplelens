#include "VBBinaryLensingLibrary.h"
#include "TripleLensingLibrary.h"
#include <stdio.h>
#include <math.h>
#include <time.h>

VBBinaryLensing VBBL;

TripleLensing::TripleLensing() {
    nphi = 32;
    secnum = 45;
    basenum = 1;
    quad_err = 0;
    quaderr_Tol = 1e-4;
    CQ = 1;
    relerr_mag = 1e-3;
    TINY = 1.0e-20;
}

TripleLensing::TripleLensing(double mlens[], complex zlens[]) {
    nphi = 32;
    secnum = 45;
    basenum = 1;
    quad_err = 0;
    quaderr_Tol = 1e-4;
    CQ = 1;
    relerr_mag = 1e-3;
    TINY = 1.0e-20;


    this -> mlens = mlens;
    this -> zlens = zlens;

    //used in polynomialCoefficients
    for (int i = 0; i < NLENS; i++) {
        zc[i] = conj(zlens[i]);
    }

    temp_const1[0][0] = 1.0;
    temp_const2[0][0][0] = 1.0;
    for (int i = 0; i < NLENS; i++) {

        for (int j = 0; j < NLENS; j++) {
            multiply_z_v2(temp_const1, zlens[j], j, i);
        }

        for (int i2 = 0; i2 <= NLENS; i2++) temp_const1[i + 1][i2] = temp_const1[i][i2];
        temp_const1[i + 1][0] = 1.0;


        //  //numerator, m[i] * (z-z1)*(z-z2) ... (z-zn)
        for (int j = 0; j <= NLENS; j++) {
            p_const[j][i] = mlens[i] * temp_const1[i][j];
        }

        for (int j = 0; j < NLENS; j++) {
            /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */

            degrees = 0;
            for (int k = 0; k < NLENS; k++) {
                if (k == j) continue;
                multiply_z_v2(temp_const2[i], zlens[k], degrees, j);
                degrees++;
            }

            for (int i2 = 0; i2 <= NLENS; i2++) temp_const2[i][j + 1][i2] = temp_const2[i][j][i2];
            temp_const2[i][j + 1][0] = 1.0;

        }

        for (int i2 = 0; i2 <= NLENS; i2++) {
            temp_const2[i + 1][0][i2] = temp_const2[i][NLENS][i2];
        }
        temp_const2[i + 1][0][0] = 1.0;
    }


}

TripleLensing::TripleLensing(double mlens[], double Zlens[]) {
    nphi = 32;
    secnum = 45;
    basenum = 1;
    quad_err = 0;
    quaderr_Tol = 1e-4;
    CQ = 1;
    relerr_mag = 1e-3;
    TINY = 1.0e-20;


    complex zlens[NLENS];
    for (int i = 0; i < NLENS ; i++) {
        zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
    }


    for (int i = 0; i < NLENS; i++) {
        zc[i] = conj(zlens[i]);
    }

    temp_const1[0][0] = 1.0;
    temp_const2[0][0][0] = 1.0;
    for (int i = 0; i < NLENS; i++) {

        for (int j = 0; j < NLENS; j++) {
            multiply_z_v2(temp_const1, zlens[j], j, i);
        }

        for (int i2 = 0; i2 <= NLENS; i2++) temp_const1[i + 1][i2] = temp_const1[i][i2];
        temp_const1[i + 1][0] = 1.0;


        //  //numerator, m[i] * (z-z1)*(z-z2) ... (z-zn)
        for (int j = 0; j <= NLENS; j++) {
            p_const[j][i] = mlens[i] * temp_const1[i][j];
        }

        for (int j = 0; j < NLENS; j++) {
            /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */

            degrees = 0;
            for (int k = 0; k < NLENS; k++) {
                if (k == j) continue;
                multiply_z_v2(temp_const2[i], zlens[k], degrees, j);
                degrees++;
            }

            for (int i2 = 0; i2 <= NLENS; i2++) temp_const2[i][j + 1][i2] = temp_const2[i][j][i2];
            temp_const2[i][j + 1][0] = 1.0;

        }

        for (int i2 = 0; i2 <= NLENS; i2++) {
            temp_const2[i + 1][0][i2] = temp_const2[i][NLENS][i2];
        }
        temp_const2[i + 1][0][0] = 1.0;
    }


    this -> mlens = mlens;
    this -> zlens = zlens;
    this -> Zlens = Zlens;
}


void TripleLensing::reset2(double mlens[],  double Zlens[]) {

    complex zlens[NLENS];
    for (int i = 0; i < NLENS ; i++) {
        zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
    }


    for (int i = 0; i < NLENS; i++) {
        zc[i] = conj(zlens[i]);
    }

    temp_const1[0][0] = 1.0;
    temp_const2[0][0][0] = 1.0;
    for (int i = 0; i < NLENS; i++) {
        for (int j = 0; j < NLENS; j++) {
            multiply_z_v2(temp_const1, zlens[j], j, i);
        }

        for (int i2 = 0; i2 <= NLENS; i2++) temp_const1[i + 1][i2] = temp_const1[i][i2];
        temp_const1[i + 1][0] = 1.0;

        //  //numerator, m[i] * (z-z1)*(z-z2) ... (z-zn)
        for (int j = 0; j <= NLENS; j++) {
            p_const[j][i] = mlens[i] * temp_const1[i][j];
        }

        for (int j = 0; j < NLENS; j++) {
            /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */

            degrees = 0;
            for (int k = 0; k < NLENS; k++) {
                if (k == j) continue;
                multiply_z_v2(temp_const2[i], zlens[k], degrees, j);
                degrees++;
            }

            for (int i2 = 0; i2 <= NLENS; i2++) temp_const2[i][j + 1][i2] = temp_const2[i][j][i2];
            temp_const2[i][j + 1][0] = 1.0;

        }

        for (int i2 = 0; i2 <= NLENS; i2++) {
            temp_const2[i + 1][0][i2] = temp_const2[i][NLENS][i2];
        }
        temp_const2[i + 1][0][0] = 1.0;
    }
    this -> mlens = mlens;
    this -> zlens = zlens;
    this -> Zlens = Zlens;
}


void TripleLensing::reset(double mlens[], complex zlens[]) {
    nphi = 32;
    secnum = 45;
    basenum = 1;
    quad_err = 0;
    quaderr_Tol = 1e-4;
    CQ = 1;
    relerr_mag = 1e-3;
    TINY = 1.0e-20;


    this -> mlens = mlens;
    this -> zlens = zlens;

    //used in polynomialCoefficients
    for (int i = 0; i < NLENS; i++) {
        zc[i] = conj(zlens[i]);
    }

    temp_const1[0][0] = 1.0;
    temp_const2[0][0][0] = 1.0;
    for (int i = 0; i < NLENS; i++) {

        for (int j = 0; j < NLENS; j++) {
            multiply_z_v2(temp_const1, zlens[j], j, i);
        }

        for (int i2 = 0; i2 <= NLENS; i2++) temp_const1[i + 1][i2] = temp_const1[i][i2];
        temp_const1[i + 1][0] = 1.0;


        //  //numerator, m[i] * (z-z1)*(z-z2) ... (z-zn)
        for (int j = 0; j <= NLENS; j++) {
            p_const[j][i] = mlens[i] * temp_const1[i][j];
        }

        for (int j = 0; j < NLENS; j++) {
            /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */

            degrees = 0;
            for (int k = 0; k < NLENS; k++) {
                if (k == j) continue;
                multiply_z_v2(temp_const2[i], zlens[k], degrees, j);
                degrees++;
            }

            for (int i2 = 0; i2 <= NLENS; i2++) temp_const2[i][j + 1][i2] = temp_const2[i][j][i2];
            temp_const2[i][j + 1][0] = 1.0;

        }

        for (int i2 = 0; i2 <= NLENS; i2++) {
            temp_const2[i + 1][0][i2] = temp_const2[i][NLENS][i2];
        }
        temp_const2[i + 1][0][0] = 1.0;
    }
}


void TripleLensing::reset3(double mlens[], complex zlens[]) {
    // just update mlens and zlens
    this -> mlens = mlens;
    this -> zlens = zlens;

    //used in polynomialCoefficients
    for (int i = 0; i < NLENS; i++) {
        zc[i] = conj(zlens[i]);
    }


    temp_const1[0][0] = 1.0;
    temp_const2[0][0][0] = 1.0;
    for (int i = 0; i < NLENS; i++) {


        for (int j = 0; j < NLENS; j++) {
            multiply_z_v2(temp_const1, zlens[j], j, i);
        }

        for (int i2 = 0; i2 <= NLENS; i2++) temp_const1[i + 1][i2] = temp_const1[i][i2];
        temp_const1[i + 1][0] = 1.0;


        //  //numerator, m[i] * (z-z1)*(z-z2) ... (z-zn)
        for (int j = 0; j <= NLENS; j++) {
            p_const[j][i] = mlens[i] * temp_const1[i][j];
        }

        for (int j = 0; j < NLENS; j++) {
            /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */

            degrees = 0;
            for (int k = 0; k < NLENS; k++) {
                if (k == j) continue;
                multiply_z_v2(temp_const2[i], zlens[k], degrees, j);
                degrees++;
            }

            for (int i2 = 0; i2 <= NLENS; i2++) temp_const2[i][j + 1][i2] = temp_const2[i][j][i2];
            temp_const2[i][j + 1][0] = 1.0;

        }

        for (int i2 = 0; i2 <= NLENS; i2++) {
            temp_const2[i + 1][0][i2] = temp_const2[i][NLENS][i2];
        }
        temp_const2[i + 1][0][0] = 1.0;
    }

}

void TripleLensing::tripleFSLimb2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, int secnum, int basenum, double quaderr_Tol, double relerr_mag, double mags[], int Np, double limba1, double RelTolLimb, double AbsTolLimb)
{
    complex zlens[NLENS];
    for (int i = 0; i < NLENS ; i++) {
        zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
    }

    reset3(mlens, zlens);

    this->secnum = secnum;
    this->quaderr_Tol = quaderr_Tol;
    this->basenum = basenum;
    this->relerr_mag = relerr_mag;
    this->RelTolLimb = RelTolLimb;
    this->AbsTolLimb = AbsTolLimb;

    for (int i = 0; i < Np; i++) {
        mags[i] = TripleMagDark(xsCenters[i], ysCenters[i], rs, limba1);
    }

}


#define EPS_CLOSE 1.0e-10 // how accurate the root should be
#define MAXIT_NEWTON 50   // maximum steps in the Newton-Raphson method

// Given the lens information, the source position and a guess of the image position z, use the Newton-Raphson method to find
// the next step to be taken, *dz, and its magnification *mu
//
void newtonStep(double mlens[], complex zlens[], double xs, double ys, complex z, complex *dz, double *mu)
{
  complex zs;
  complex dzs;
  complex dzsdz;

  int i;
  // int flag = 0;     //

  double r2_1;//, r2_2;
  double x, y;
  double dx, dy;
  double Jxx = 1.0, Jyy = 0.0, Jxy = 0.0, Jyx; // Jxx + Jyy = 2.0, Jxy = Jyx
  // double sum2 = 0.0, sq = 0.0;
  static double TINY = 1.0e-20;
  double dxs, dys;

  zs = complex(xs, ys);
  dzs = zs - ( z - mlens[0] / conj(z - zlens[0]) - mlens[1] / conj(z - zlens[1]) - mlens[2] / conj(z - zlens[2]) );

  dxs = real(dzs);
  dys = imag(dzs);


  // calculate the Jacobian 2x2 matrix (Jxx, Jxy, Jyx, Jyy)
  x = real(z);
  y = imag(z);

  for (i = 0; i < NLENS; i++) {
    dx = x - real(zlens[i]);
    dy = y - imag(zlens[i]);

    r2_1 = dx * dx + dy * dy + TINY; // avoid zero in the denominator

    Jxx += mlens[i] * (dx * dx - dy * dy) / (r2_1 * r2_1);
    Jxy += 2.0 * mlens[i] * dx * dy / (r2_1 * r2_1);
  }

  //
  //analytical results for the other components of the Jacobian
  //
  Jyy = 2.0 - Jxx;
  Jyx = Jxy;
  *mu = 1.0 / (Jxx * Jyy - Jxy * Jyx);

  // now calculate the change we need to make using Newton's iteration
  dx = (*mu) * ( Jyy * dxs - Jxy * dys);
  dy = (*mu) * (-Jxy * dxs + Jxx * dys);

  *dz = complex(dx, dy);
}


//     find the images close to a initial guess
void findCloseImages(double mlens[], complex zlens[], double xs, double ys, complex *z, bool *imageFound)
{
  int i;
  complex dz;
  double mu;

    *imageFound = false;


    for (i = 0; i < 50; i++) {
      newtonStep(mlens, zlens, xs, ys, *z, &dz, &mu);

      *z = *z + dz;

      if (abs(dz) < EPS_CLOSE) {
        *imageFound = true;
        break;      // we are done in the iteration
      }
    }

#ifdef VERBOSE
    if (*imageFound) {
      fprintf(stderr, "image found at position: %le %le with magnification %le\n", real(*z), imag(*z), mu);
    } else {
      fprintf(stderr, "failed to find image \n");
    }
#endif

}


#undef EPS_CLOSE
#undef MAXIT_NEWTON


void TripleLensing::tripleFS2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, int secnum, int basenum, double quaderr_Tol, double relerr_mag, double mags[], int Np)
{
    complex zlens[NLENS];
    for (int i = 0; i < NLENS ; i++) {
        zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
    }

    reset3(mlens, zlens);

    this->secnum = secnum;
    this->quaderr_Tol = quaderr_Tol;
    this->basenum = basenum;
    this->relerr_mag = relerr_mag;

    // for (int i = 0; i < Np; i++) {
    //     mags[i] = TripleMag(xsCenters[i], ysCenters[i], rs);
    // }

    if (rs>1e-10){
        for (int i = 0; i < Np; i++) {
            mags[i] = TripleMag(xsCenters[i], ysCenters[i], rs);
        }
    }else{
        for (int i = 0; i < Np; i++) {
            mags[i] = TriplePS(xsCenters[i], ysCenters[i]);
        }        
    }

}

// Gould approximation calculate interface to Python
void TripleLensing::tripleGould2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, double Gamma, double mags[], int Np)
{
    double A2rho2, A4rho4;
    int Np2 = 2*Np;
    complex zlens[NLENS];
    for (int i = 0; i < NLENS ; i++) {
        zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
    }
    reset3(mlens, zlens);
    for (int i = 0; i < Np; i++) {
        mags[i] = gould(xsCenters[i], ysCenters[i], rs, Gamma, &A2rho2, &A4rho4);
        mags[Np+i] = A2rho2;
        mags[Np2+i] = A4rho4;
    }
}


// return the number of real solutions to python
void TripleLensing::triple_num_real_sol2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double true_solution_threshold , double numbers_mups[], int Np)
{
    double mu = 0., muTotal = 0.;
    int flag = 0;
    complex zlens[NLENS];

    complex zr[DEGREE];
    complex coefficients[DEGREE + 1];
    int true_imgnum_cnt = 0;

    for (int i = 0; i < NLENS ; i++) {
        zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
    }

    reset3(mlens, zlens);

    for (int i = 0; i < Np; i++) {
        muTotal = 0.;

        polynomialCoefficients(xsCenters[i], ysCenters[i], coefficients);
        VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

        true_imgnum_cnt = 0;
        for (int j = 0; j < DEGREE; j++) {
            // true_imgnum_cnt += trueSolution_nomu(xsCenters[i], ysCenters[i], zr[j], true_solution_threshold);
            flag =  trueSolution_withmu(xsCenters[i], ysCenters[i], zr[j], true_solution_threshold, &mu);
            true_imgnum_cnt += flag;
            muTotal += ( flag ? abs(mu) : 0  );

        }
        numbers_mups[i] = true_imgnum_cnt;
        numbers_mups[i + Np] = muTotal;
    }

}



// double TripleLensing::tripleFS2python(double mlens[], double Zlens[], double xsCenter, double ysCenter, double rs) {
double TripleLensing::tripleFS2python(double xsCenter, double ysCenter, double rs) {
    mu = TripleMag(xsCenter, ysCenter, rs);
    return mu;
}

double TripleLensing::tripleFS2python(double mlens[], double zlens[], double xsCenter, double ysCenter, double rs) {
    reset2(mlens, zlens);
    mu = TripleMag(xsCenter, ysCenter, rs);
    return mu;
}

void TripleLensing::solv_lens_equation(double zrxy[], double mlens[], double zlens[], double xs, double ys, int nlens) {
    complex zr[DEGREE];
    complex coefficients[DEGREE + 1];

    complex Zlens[NLENS];
    int i;
    for (i = 0; i < NLENS; i++) {
        Zlens[i] = complex(zlens[i], zlens[i + NLENS]);
    }
    reset3(mlens, Zlens);
    polynomialCoefficients(xs, ys, coefficients);
    VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

    for (i = 0; i < DEGREE; i++) {
        zrxy[i] = zr[i].re;
        zrxy[i + DEGREE] = zr[i].im;
    }
}


void TripleLensing::TriLightCurve(double *pr, double *mags, double *y1s, double *y2s, int np) {
    // pr: t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs
    double s2 = pr[3], q2 = pr[4], s3 = pr[6], q3 = pr[7], psi = pr[8], rs = pr[9];

    // parameter conversion, from s2, q2,, s3, q3, psi to mlens[i], zlens[i]
    double inv1andq2 = 1 / (1 + q2);
    double inv1andq2q3 = 1 / (1 + q2 + q3);
    mlens[0] = inv1andq2q3;
    mlens[1] = q2 * inv1andq2q3;
    mlens[2] = q3 * inv1andq2q3;
    zlens[0] = complex(- q2 * inv1andq2 * s2 , 0);
    zlens[1] = complex( inv1andq2 * s2 , 0);
    zlens[2] = complex(zlens[0].re + s3 * cos(psi), s3 * sin(psi));

    if (rs>1e-10){
        for (int i = 0; i < np; i++) {
            mags[i] = TripleMag(y1s[i], y2s[i], rs);
        }
    }else{
        for (int i = 0; i < np; i++) {
            mags[i] = TriplePS(y1s[i], y2s[i]);
        }        
    }

}


//output in two files the triple critical curves and caustics
void TripleLensing::outputCriticalTriple_list(double allxys[], double mlens[], double zlens[], int nlens, int NPS)
// void outputCriticalTriple_list(double allxys[], double mlens[], double zlens[], int NLENS, int NPS)
{
    VBBinaryLensing VBBL;
    VBBL.outputCriticalTriple_list(allxys, mlens, zlens, nlens, NPS);
}

void TripleLensing::outputCriticalBinary_list(double resxy[], double s, double q, int NPS){
    VBBinaryLensing VBBL;
    VBBL.outputCriticalBinary_list(resxy, s, q, NPS);    
}


void outsys(double mlens[], complex zlens[], double t0, double u0, double tE, double s2, double q2, double alpha, double s3, double q3, double psi, double rs, double xsCenter, double ysCenter) {

    // output triple lensing event parameters
    FILE *file, *fileImage1;
    file = fopen("data/lens_system_triple.dat", "w");
    fileImage1 = fopen("data/lens_system.dat", "w");
    fprintf(file, "t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs, xsCenter, ysCenter\n");
    fprintf(file, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs, xsCenter, ysCenter);
    fprintf(fileImage1, "%.15f %.15f %.15f\n", xsCenter, ysCenter, rs);
    for (int i = 0; i < NLENS; i++) {
        fprintf(fileImage1, "%.15f %.15f %.15f\n", mlens[i], zlens[i].re, zlens[i].im);
    }
    fclose(file);
    fclose(fileImage1);
}

double TripleLensing::angcos(_point *p1, _point *p2, _point *p3, _point *p4) {
    // calculate the cos between the vector p1 --> p2, and vector p3 --> p4;
    double x1, y1, x2, y2;
    x1 = p2->x1 - p1->x1;
    y1 = p2->x2 - p1->x2;

    x2 = p4->x1 - p3->x1;
    y2 = p4->x2 - p3->x2;

    return (x1 * x2 + y1 * y2) / ( sqrt(  (x1 * x1 + y1 * y1) * ( x2 * x2 + y2 * y2 ) ) );

}

// tripleFS_v2_savehalf_quadtest(mlens, zlens, xsCenter, ysCenter, rs, nphi, &finalnphi, secnum, basenum, &quad_err, quaderr_Tol);
double TripleLensing::TripleMag(double xsCenter, double ysCenter, double rs) {
    // VBBinaryLensing VBBL;

    finalNPS = 0;
    rho2 = rs * rs;
    areaSource = M_PI * rho2;
    ftime = 1;
    relerr_priv = this->relerr_mag;
    secnum_priv = this->secnum;
    quad_err = this->quaderr_Tol;
    basenum_priv = this->basenum;
    finalnphi = 0;
    int area_quality_local = 1;

    double _curr_relerr_priv = relerr_priv;

    // 2021.06.08
    for (int jj = 0; jj < DEGREE + 1; jj++) {
        zr[jj] = complex(0, 0);
    }

    muPS = tripleQuatrapoleTest(xsCenter, ysCenter, rs); // this is not robust, didn't check the number of true solution is 4, 6, 8, or 10

// #ifdef VERBOSE
//     char arr[1024];
//     fprintf(stderr, "please input >>> ");
//     scanf("%c%*c", &arr[0]);
// #endif


//     if (CQ * quad_err <= quaderr_Tol) {
//         muPS = gould(xsCenter, ysCenter, rs, 0);
// #ifdef VERBOSE
//         fprintf(stderr, "quad_err %f, using gould approximation, mu gould = %f\n", quad_err, muPS);
// #endif
//         ifFinite = 0;
//         return muPS;
//     } 

    if ( CQ * quad_err <= 1e-1 * quaderr_Tol) {
#ifdef VERBOSE
        fprintf(stderr, "quad_err %f,using point source magnification, muPS = %f\n", quad_err, muPS);
#endif
        ifFinite = 0;
        return muPS;
    }
    else if (CQ * quad_err <= quaderr_Tol) {
        muPS = gould(xsCenter, ysCenter, rs, 0);
#ifdef VERBOSE
        fprintf(stderr, "quad_err %f, using gould approximation, mu gould = %f\n", quad_err, muPS);
#endif
        ifFinite = 0;
        return muPS;
    } else {
        ifFinite = 1;
#ifdef VERBOSE
        fprintf(stderr, "quad_err %f,using caustics crossing computation\n", quad_err);
#endif
        mu0 = muPS;
        secnum_priv = fmax(fmin((int)(2 * secnum_priv * ( abs(log10(quad_err * 1e4)) * mu0 * mu0) / abs(log10(rs))) , 45), secnum_priv);
        _curr_relerr_priv = fmax( _curr_relerr_priv , _curr_relerr_priv / fmax( quad_err, mu0 * mu0 ) * abs(log10(rs))  );


#ifdef VERBOSE
        printf("point magnification: %f, relerr_priv = %.3e\n", mu0, _curr_relerr_priv);
#endif

        _sols *imageTracks;//  = new _sols;
        _sols *prevstore;// = new _sols;
        _linkedarray *phis;// = new _linkedarray;
        phis = getphis_v3(  xsCenter,  ysCenter,  rs);

//         for (int i = 0; i < 10; i++) {
        for (int i = 0; i < 7; i++) { // 2022.06.05
            if (i > 0) { // 2021.06.08
                delete imageTracks;
                //imageTracks = new _sols;
#ifdef VERBOSE
                fprintf(stderr, "prevstore->length in tripleFS_v2_savehalf %d, phis->length %d\n", prevstore->length, phis->length);
#endif
            }


            imageTracks = outputTracks_v2_savehalf(xsCenter,  ysCenter,  rs, phis, &prevstore);

            finalnphi = phis->length;
#ifdef VERBOSE
            saveTracks(imageTracks, finalnphi);// temp
            fprintf(stderr, "\t\t\t imageTracks saved. \n");
#endif


            if (!imageTracks) {
#ifdef VERBOSE
                fprintf(stderr, "\n\nin tripleFS, i= %d, imageTracks = Null (be careful), nphi = %d, xsCenter=%f, ysCenter = %f\n\n", i, finalnphi, xsCenter, ysCenter);
#endif
                if (distype == 2) {
                    delete phis;
                    //phis = new _linkedarray;
                    maxmuidx = (int)(absint(maxmuidx + 1)) % secnum;
                    basenum_priv *= 2;
                    phis = getphis_v3( xsCenter,  ysCenter, ys);
                } else {
                    phis->extend();
                }
                i -= 1;
                ftime = 1;
                continue;
            }
            ftime = 0;

#ifdef parabcorr
            area = areaFunc_parab(imageTracks, rs);
#else
            area = areaFunc(imageTracks, rs, finalnphi, muPS, mu0, &area_quality_local);
#endif
            mu = area / areaSource;

            area_quality = area_quality_local;
            if (area_quality_local == 0) {
                _curr_relerr_priv *= 2;
            }

#ifdef VERBOSE
            fprintf(stderr, "in tripleFS, i= %d, mu0= %f, mu= %f, nphi = %d, xsCenter=%f, ysCenter = %f, nimages = %d, abs(mu - mu0) / mu = %.3e, errTol = %.3e, \n", i, mu0, mu, finalnphi, xsCenter, ysCenter, imageTracks->length, abs(mu - mu0) / mu, _curr_relerr_priv);
#endif
#ifdef verbose
            fprintf(stderr, "in tripleFS, i= %d, mu0= %f, mu= %f, nphi = %d, xsCenter=%f, ysCenter = %f, nimages = %d\n", i, mu0, mu, finalnphi, xsCenter, ysCenter, imageTracks->length);
#endif
            if (abs(mu - mu0) / mu < _curr_relerr_priv) {
                break;
            } else if (abs(mu - mu0) / mu < _curr_relerr_priv * i  ||  phis->length > 5e5 ) { // add "|  imageTracks->length > 1e6" on 2022.06.04, when logq is too small <~ -6, the memory usage increases due to he amount of sampling points on the source boundary
#ifdef VERBOSE
                fprintf(stderr, "return mu = %f, because mu0 = %f, mu = %f, relerr_priv*i = %.3e this might be a hard case, area_quality_local = %d\n", mu, mu0, mu, _curr_relerr_priv * i, area_quality_local);

#endif
                delete imageTracks;
                delete phis;
                delete prevstore;

                area_quality = 2; // return because of looser threshold, quality might be low
                return (mu);
//             } else if ( (i > 2 && (0.5 * ( abs(mu / mu0) + abs(mu0 / mu) > 2.1)) ) ) {
            } else if ( (i > 2 && ( 0.5*(abs(mu/mu0) + abs(mu0/mu)) > 2.1 ) ) ) { // modified on 2022.06.04
#ifdef VERBOSE
                fprintf(stderr, "return muPS = %f, because mu0 = %f, mu = %f are too arbitrary, this might be a hard case, area_quality_local=%d\n", muPS, mu0, mu, area_quality_local);
#endif
                delete imageTracks;
                delete phis;
                delete prevstore;

                area_quality = 3; // return point source magnification
                return (muPS);
            }

            mu0 = mu;
            phis->extend();

        }
#ifdef VERBOSE
        printf("mu = %f\n", mu);
        // saveTracks(imageTracks);
#endif

        // show image track animation here:

        finalNPS = finalnphi;
        nimages = imageTracks->length;
        delete imageTracks;
        delete phis;
        delete prevstore;
        return (mu);
    }

}


//
//   implement Gould (2008) approximation, perhaps we can use static tables for phi to get rid of cos and sin
//
double TripleLensing::gould(double xsCenter, double ysCenter, double rs, double Gamma)
{
    int i;
    double phi, dphi;
    double muRho[8];
    double muRhoHalf[4];
    double AhalfPlus, APlus, Across, A0;
    double A2rho2, A4rho4;// gould, 10%,
    double xs, ys;
    // int nimages;
    double mag;
    A0 = TriplePS(xsCenter, ysCenter);

    // 8 points at full radius
    dphi = M_PI * 0.25;
    for (i = 0; i < 8; i++) {
        phi = dphi * (i - 1);
        xs = xsCenter + rs * cos(phi);
        ys = ysCenter + rs * sin(phi);
        muRho[i] = TriplePS(xs, ys);
    }

    // four points at half radius
    dphi = M_PI * 0.5;
    for (i = 0; i < 4; i++) {
        phi = dphi * (i - 1);
        xs = xsCenter + 0.5 * rs * cos(phi);
        ys = ysCenter + 0.5 * rs * sin(phi);
        muRhoHalf[i] = TriplePS(xs, ys);
    }

    AhalfPlus = (muRhoHalf[0] + muRhoHalf[1] + muRhoHalf[2] + muRhoHalf[3]) * 0.25 - A0;
    APlus  = (muRho[0] + muRho[2] + muRho[4] + muRho[6]) * 0.25 - A0;
    Across = (muRho[1] + muRho[3] + muRho[5] + muRho[7]) * 0.25 - A0;

    A2rho2 = (16.0 * AhalfPlus - APlus) / 3.0;
    A4rho4 = (APlus + Across) / 2.0 - A2rho2;
#ifdef VERBOSE
    printf("I am in gould, second and fourth order corrections: %f %f\n", A2rho2 / 2.0, A4rho4 / 3.0);
#endif

    mag = A0 + A2rho2 / 2.0 * (1.0 - 0.2 * Gamma)  + A4rho4 / 3.0 * (1.0 - 11.0 / 35.0 * Gamma);
    return (mag);
}

double TripleLensing::gould(double xsCenter, double ysCenter, double rs, double Gamma, double *A2rho2, double *A4rho4)
{
    int i;
    double phi, dphi;
    double muRho[8];
    double muRhoHalf[4];
    double AhalfPlus, APlus, Across, A0;
    // double A2rho2, A4rho4;// gould, 10%,
    double xs, ys;
    // int nimages;
    double mag;
    A0 = TriplePS(xsCenter, ysCenter);

    // 8 points at full radius
    dphi = M_PI * 0.25;
    for (i = 0; i < 8; i++) {
        phi = dphi * (i - 1);
        xs = xsCenter + rs * cos(phi);
        ys = ysCenter + rs * sin(phi);
        muRho[i] = TriplePS(xs, ys);
    }

    // four points at half radius
    dphi = M_PI * 0.5;
    for (i = 0; i < 4; i++) {
        phi = dphi * (i - 1);
        xs = xsCenter + 0.5 * rs * cos(phi);
        ys = ysCenter + 0.5 * rs * sin(phi);
        muRhoHalf[i] = TriplePS(xs, ys);
    }

    AhalfPlus = (muRhoHalf[0] + muRhoHalf[1] + muRhoHalf[2] + muRhoHalf[3]) * 0.25 - A0;
    APlus  = (muRho[0] + muRho[2] + muRho[4] + muRho[6]) * 0.25 - A0;
    Across = (muRho[1] + muRho[3] + muRho[5] + muRho[7]) * 0.25 - A0;

    *A2rho2 = (16.0 * AhalfPlus - APlus) / 3.0;
    *A4rho4 = (APlus + Across) / 2.0 - *A2rho2;
#ifdef VERBOSE
    printf("I am in gould, second and fourth order corrections: %f %f\n", A2rho2 / 2.0, A4rho4 / 3.0);
#endif

    mag = A0 + *A2rho2 / 2.0 * (1.0 - 0.2 * Gamma)  + *A4rho4 / 3.0 * (1.0 - 11.0 / 35.0 * Gamma);
    return (mag);
}


double TripleLensing::TripleMagDark(double xsCenter, double ysCenter, double RSv, double a1) {
    double Mag = -1.0, Magold = 0., Tolv = this->AbsTolLimb;
    double tc, lb, rb, lc, rc, cb = 0, cc, r2, cr2, scr2;
    int c = 0, flag;
    double currerr, maxerr;
    annulus *first, *scan, *scan2;
    int nannuli, nannold, totNPS = 1, minannuli = 1;
    // _sols *Images;

    double y_1 = xsCenter;
    double y_2 = ysCenter;
    while ((Mag < 0.9) && (c < 3)) {

        first = new annulus;
        first->bin = 0.;
        first->cum = 0.;// cumulative function F(r)
        first->Mag = TriplePS(y_1, y_2);
        first->nim = nimages;//Images->length; // nimages is a global variable in TripleLensing

        first->f = 3 / (3 - a1);// r = 0, so using point source mag computation above, Images is a _sols object, its length = nim -- number of image tracks
        first->err = 0;
        first->prev = 0;

        first->next = new annulus;
        scan = first->next;
        scan->prev = first;
        scan->next = 0;
        scan->bin = 1.;
        scan->cum = 1.;

        // scan->Mag = TripleMagFull(xsCenter, ysCenter, RSv);
        scan->Mag = TripleMag(xsCenter, ysCenter, RSv);

        totNPS += finalNPS;//NPS;
        scan->nim = nimages; //Images->length;
        // delete Images;
        scan->f = first->f * (1 - a1); // r = 1, so using FS mag computation above
        if (scan->nim == scan->prev->nim) {
            scan->err = fabs((scan->Mag - scan->prev->Mag) * (scan->prev->f - scan->f) / 4);
        }
        else {
            scan->err = fabs((scan->Mag) * (scan->prev->f - scan->f) / 4);
        }

        Magold = scan->Mag;
        Mag = scan->Mag;

        currerr = scan->err;
        flag = 0;

        nannuli = nannold = 1;
        while (((flag < nannold + 5) && (currerr > Tolv) && (currerr > this->RelTolLimb * Mag)) || (nannuli < minannuli)) {
            // update/find the annuli with the maximum error
            maxerr = 0;
            for (scan2 = first->next; scan2; scan2 = scan2->next) {
                if (scan2->err > maxerr) {
                    maxerr = scan2->err;
                    scan = scan2;
                }
            }
            nannuli++;
            Magold = Mag;
            Mag -= (scan->Mag * scan->bin * scan->bin - scan->prev->Mag * scan->prev->bin * scan->prev->bin) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
            currerr -= scan->err;
            rb = scan->bin; // current bin rb, right bin
            rc = scan->cum; // right cumulative function
            lb = scan->prev->bin; // prev bin lb, left bin
            lc = scan->prev->cum; // left cumulative function
            tc = (lc + rc) / 2; // Mid Cumulative Function as target
            do {
                cb = rb + (tc - rc) * (rb - lb) / (rc - lc); // center bin, bin actually is the corresponding radius in range (0,1)
                r2 = cb * cb;
                cr2 = 1 - r2;
                scr2 = sqrt(cr2);
                cc = (3 * r2 * (1 - a1) - 2 * a1 * (scr2 * cr2 - 1)) / (3 - a1);
                // cumulative function F(r) from equ(45) at new r = cb, integral of 2rf(r)
                if (cc > tc) {
                    rb = cb;
                    rc = cc;
                }
                else {
                    lb = cb;
                    lc = cc;
                }
            } while (fabs(cc - tc) > 1.e-5); // loop until equal partition the cumulative function

            if (RSv > 1e-3 && RSv * cb < 1e-4) {
                fprintf(stderr, "\t\t\t return because rho = %.25e\n", RSv * cb);
                return Magold;
            }

            scan->prev->next = new annulus;
            scan->prev->next->prev = scan->prev;
            scan->prev = scan->prev->next;
            scan->prev->next = scan;
            scan->prev->bin = cb;
            scan->prev->cum = cc;
            scan->prev->f = first->f * (1 - a1 * (1 - scr2));//Bozza 2010 42

            // scan->prev->Mag = TripleMagFull(xsCenter, ysCenter, RSv * cb);
            scan->prev->Mag = TripleMag(xsCenter, ysCenter, RSv * cb);

            totNPS += finalNPS;//NPS;

            scan->prev->nim = nimages;//Images->length;
            if (scan->prev->prev->nim == scan->prev->nim) {
                scan->prev->err = fabs((scan->prev->Mag - scan->prev->prev->Mag) * (scan->prev->prev->f - scan->prev->f) * (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin) / 4);
            }
            else {
                scan->prev->err = fabs((scan->prev->bin * scan->prev->bin * scan->prev->Mag - scan->prev->prev->bin * scan->prev->prev->bin * scan->prev->prev->Mag) * (scan->prev->prev->f - scan->prev->f) / 4);
            }
            if (scan->nim == scan->prev->nim) {
                scan->err = fabs((scan->Mag - scan->prev->Mag) * (scan->prev->f - scan->f) * (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin) / 4);
            }
            else {
                scan->err = fabs((scan->bin * scan->bin * scan->Mag - scan->prev->bin * scan->prev->bin * scan->prev->Mag) * (scan->prev->f - scan->f) / 4);
            }
            rb = (scan->Mag + scan->prev->prev->Mag - 2 * scan->prev->Mag);
            scan->prev->err += fabs(rb * (scan->prev->prev->f - scan->prev->f) * (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin));
            scan->err += fabs(rb * (scan->prev->f - scan->f) * (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin));
#ifdef _PRINT_ERRORS_DARK
            printf("\n%d", nimages);
#endif

            Mag += (scan->bin * scan->bin * scan->Mag - cb * cb * scan->prev->Mag) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
            Mag += (cb * cb * scan->prev->Mag - scan->prev->prev->bin * scan->prev->prev->bin * scan->prev->prev->Mag) * (scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin);
            currerr += scan->err + scan->prev->err;

            if (fabs(Magold - Mag) * 2 < Tolv) {
                flag++;
            }
            else {
                flag = 0;
                nannold = nannuli;
            }

        }

        Tolv /= 10;
        c++;
    }
    NPS = totNPS;
    therr = currerr;
    finalNPS = NPS;
    return Mag;
}



double TripleLensing::TripleMagFull(double xsCenter, double ysCenter, double rs) {
    finalNPS = 0;
    rho2 = rs * rs;
    areaSource = M_PI * rho2;
    ftime = 1;
    relerr_priv = relerr_mag;
    secnum_priv = secnum;
    int area_quality;

    muPS = TriplePS(xsCenter, ysCenter);
    ifFinite = 1;
    mu0 = muPS;
#ifdef VERBOSE
    printf("point magnification: %f, relerr_priv = %.3e\n", mu0, relerr_priv);
#endif
    _sols *imageTracks  = new _sols;
    _sols *prevstore  = new _sols;
    _linkedarray *phis = new _linkedarray;

    phis = getphis_v3(xsCenter,  ysCenter,  rs);
    for (int i = 0; i < 10; i++) {
        delete imageTracks;
        imageTracks = new _sols;
#ifdef VERBOSE
        fprintf(stderr, "prevstore->length in tripleFS_v2_savehalf %d\n", prevstore->length);
#endif
        imageTracks = outputTracks_v2_savehalf(xsCenter,  ysCenter,  rs, phis, &prevstore);

        finalnphi = phis->length;
        if (!imageTracks) {
#ifdef VERBOSE
            fprintf(stderr, "\n\nin TripleMagFull, i= %d, imageTracks = Null, nphi = %d, xsCenter=%f, ysCenter = %f, rs = %f\n\n", i, finalnphi, xsCenter, ysCenter, rs);
#endif
            if (distype == 2) {
                delete phis;
                phis = new _linkedarray;
                maxmuidx = (int)(absint(maxmuidx + 1)) % secnum;
                basenum *= 2;
                phis = getphis_v3( xsCenter,  ysCenter, ys);
            } else {
                phis->extend();
            }
            i -= 1;
            ftime = 1;
            continue;
        }
        ftime = 0;

#ifdef parabcorr
        area = areaFunc_parab(imageTracks, rs);
#else
        area = areaFunc(imageTracks, rs, finalnphi, muPS, mu0, &area_quality);
#endif

        mu = area / areaSource;
#ifdef VERBOSE
        fprintf(stderr, "in TripleMagFull, i= %d, mu0= %f, mu= %f, nphi = %d, xsCenter=%f, ysCenter = %f, nimages = %d, rs = %f\n", i, mu0, mu, finalnphi, xsCenter, ysCenter, imageTracks->length, rs);
#endif
#ifdef verbose
        fprintf(stderr, "in TripleMagFull, i= %d, mu0= %f, mu= %f, nphi = %d, xsCenter=%f, ysCenter = %f, nimages = %d, rs = %f\n", i, mu0, mu, finalnphi, xsCenter, ysCenter, imageTracks->length, rs);

#endif
        if (abs(mu - mu0) / mu < relerr_priv) break;
        mu0 = mu;
        phis->extend();
    }
#ifdef VERBOSE
    printf("mu = %f\n", mu);
    saveTracks(imageTracks, finalnphi);
#endif


    finalNPS = finalnphi;
    nimages = imageTracks->length;
    delete imageTracks;
    delete phis;
    delete prevstore;
    return (mu);

}


// tripleQuatrapoleTest Bozza 2018 equation 39 with c_Q = 1, add on
double TripleLensing::tripleQuatrapoleTest(double xs, double ys, double rs) {

    polynomialCoefficients(xs, ys, coefficients);

    VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

    muTotal = 0.0;
    eacherr = 0.0;
    // rho2 = rs * rs;

    nimages = 0;
    quad_err = 0.0;//quad_err
    for (int i = 0; i < DEGREE; i++) {
        flag = trueSolution_qtest(xs, ys, zr[i]);
        nimages = nimages + flag;
        if (flag == 1) {
            muTotal = muTotal + abs(mu);
            // fprintf(stderr, "muTotal = %f \n", muTotal);
            // fprintf(stderr, "mu, J1.re, J1.im, J2.re, J2.im, dJ.re, dJ.im, J3.re, J3.im : %f, %f, %f, %f, %f, %f, %f, %f, %f\n", mu, J1.re, J1.im, J2.re, J2.im, dJ.re, dJ.im, J3.re, J3.im  );
            eacherr = 0.0;
            J1c = conj(J1);
            J1c2 = J1c * J1c;
            firstterm = J1c2 * J1c * J2 * J2;
            dJ2 = dJ.re * dJ.re;
            secondterm = (3.0 - 3.0 * dJ.re + dJ2 / 2.0) * (J2 * conj(J2));
            thirdterm = dJ.re * J1c2 * J3;

            eacherr += abs( (6.0 * firstterm - 2.0 * secondterm + 2.0 * thirdterm).re ); //muQI
            eacherr += abs( 6.0 * firstterm.im );
            // quad_err += eacherr * rho2 / abs(dJ2 * dJ2 * dJ.re);
            quad_err += eacherr / (abs(dJ2 * dJ2 * dJ.re) + TINY); // 2021.06.08, add +TINY
        }
    }
    quad_err *= rho2; // so if rho = 0, we would use point source magnification
    return (muTotal);
}


void addPoint(_curve *final, _point *p)
{
    final->append(p->x1, p->x2);
    final->last->phi = p->phi;
    final->last->mu = p->mu;
    final->last->zs = p->zs;
    final->last->flag = p->flag ;
    final->last->thetaJ = p->thetaJ;

    final->last->ds = p->ds;
    final->last->dz = p->dz;
    final->last->absdzs = p->absdzs;
    final->last->closepairparity = p->closepairparity;

}

_curve *newSegment(_sols *imageTracks)
{
    _curve *final;
    final = new _curve;

    for (_point *p = (imageTracks->first)->first; p; p = p->next) {
        final->append(p->x1, p->x2);
        final->last->phi = p->phi;
        final->last->mu = p->mu;
        final->last->zs = p->zs;
        final->last->flag = p->flag ;
        final->last->thetaJ = p->thetaJ;
        final->last->absdzs = p->absdzs;

        final->last->ds = p->ds;
        final->last->dz = p->dz;
        final->last->closepairparity = p->closepairparity;
    }

    return (final);
}

#define _newSegment \
if(imageTracks->length > 1){\
    final = imageTracks->first;\
    imageTracks->first = imageTracks->first->next;\
    imageTracks->first->prev = 0;\
    imageTracks->length --;\
}else if (imageTracks->length == 1){\
    final = imageTracks->first;\
    imageTracks->first = 0;\
    imageTracks->length --;\
}


double TripleLensing::TriplePS(double xs, double ys) {
    double mu;
    double mulist[DEGREE];
    int total_parity = 0, flaglist[DEGREE], absdzslist[DEGREE];

    polynomialCoefficients(xs, ys, coefficients);
    VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);
    muTotal = 0.0;
    nimages = 0;
    for (int i = 0; i < DEGREE; i++) {
        flag = trueSolution(xs, ys, zr[i], &mu);
        mulist[i] = mu;
        flaglist[i] = flag;
        absdzslist[i] = absdzs;

        nimages += flag;
        muTotal += ( flag ? abs(mu) : 0  );
        total_parity += (flag ? ( mu > 0 ? 1 : -1 ) : 0 );
    }

    /* remove on 2021.11.13
    int k = 0;
    int missing_sol_idx = 0;
    int worst_true_sol_idx = 0;
    double missing_sol_absdzs = 1e3, missing_sol_mu;
    double worst_tru_sol_absdzs = -1e1, worst_tru_sol_mu;

// correction according to the number of true images and total parity
    if ( (nimages - NLENS) % 2 != 1 || total_parity != -2) {
        for (k = 0 ; k < DEGREE; k++) {
            if ( flaglist[k] == 1 &&  absdzslist[k] > worst_tru_sol_absdzs) {
                worst_tru_sol_absdzs =  absdzslist[k];
                worst_true_sol_idx = k;
                worst_tru_sol_mu =  mulist[k];
            }
            if ( flaglist[k] == 0 &&  absdzslist[k] < missing_sol_absdzs) {
                missing_sol_absdzs =  absdzslist[k];
                missing_sol_mu = mulist[k];
                missing_sol_idx = k;
            }
        }

        // if total parity is wrong, we try to remove solutions
        if ((total_parity - ( worst_tru_sol_mu > 0 ? 1 : -1  )  ) == -2 && worst_tru_sol_absdzs > SOLEPS * 0.5) {
            flaglist[worst_true_sol_idx] = -1;
            nimages -= 1;
            muTotal -= abs( worst_tru_sol_mu );
        }
        else if ( (total_parity +  ( missing_sol_mu > 0 ? 1 : -1  )  ) == -2 && missing_sol_absdzs < 2 * SOLEPS) {
            flaglist[missing_sol_idx] = 1;
            nimages += 1;
            muTotal += abs(missing_sol_mu);
        } else {
        }
    }
    */


    return (muTotal);
}



// out put pure image points around the edge of the source
// void TripleLensing::outImgPoints(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi)
void TripleLensing::outImgPoints(double xsCenter, double ysCenter, double rs, int nphi)
{

    fprintf(stderr, "Generating image positions corresponding to source boundary ...");
    /* now find the closed curves for all real and imaginary roots for a circular source */
    double mu;
    nsolution = 0;
    FSflag = 0;

    phi0 = 1.5;
    dphi = 2.0 * M_PI / (nphi - 1);
    phi = -dphi + phi0;

    //output centre image:
    FILE *fileImage0;
    fileImage0 = fopen("./data/pureImgPoints_centre.dat", "w");
    xs = xsCenter;
    ys = ysCenter;
    // polynomialCoefficients(mlens, zlens, xs, ys, coefficients, NLENS, DEGREE);
    polynomialCoefficients(xs, ys, coefficients);
    VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);
    for (int i = 0; i < DEGREE; i++) {
        // flag = trueSolution(mlens, zlens, xs, ys, zr[i], &mu, &lambda1, &lambda2, &thetaJ, NLENS, &J1, &J2, &dJ, &J3);
        flag = trueSolution(xs, ys, zr[i], &mu);
        if (flag == 1) {        // true solution
            fprintf(fileImage0, "%f %f\n", zr[i].re, zr[i].im);
        }
    }
    fclose(fileImage0);

    FILE *fileImage1, *fileImage2;
    fileImage1 = fopen("data/pureImgPoints.dat", "w");
    fileImage2 = fopen("data/pureImgPoints_falseimg.dat", "w");

    for (int j = 0; j < nphi; j++) {
        phi += dphi;
        xs = xsCenter + rs * cos(phi);
        ys = ysCenter + rs * sin(phi);
        // polynomialCoefficients(mlens, zlens, xs, ys, coefficients, NLENS, DEGREE);
        polynomialCoefficients(xs, ys, coefficients);
        VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);
        nimages = 0;
        for (int i = 0; i < DEGREE; i++) {
            // flag = trueSolution(mlens, zlens, xs, ys, zr[i], &mu, &lambda1, &lambda2, &thetaJ, NLENS, &J1, &J2, &dJ, &J3);
            flag = trueSolution( xs, ys, zr[i], &mu);
            nimages += flag;
            if (flag == 1) {        // true solution
                fprintf(fileImage1, "%f %f\n", zr[i].re, zr[i].im);
            } else {
                fprintf(fileImage2, "%f %f\n", zr[i].re, zr[i].im);
            }
        }
    }
    fclose(fileImage1);
    fclose(fileImage2);
    fprintf(stderr, " ... done.\n");
}

//output in two files the triple critical curves and caustics 2020.09.08
void outputCaustics_for_crossSection(double mlens[], complex zlens[], int nlens, int npnt)
{

    _sols *track;

    track = VBBL.outCritTriple(mlens, zlens, npnt, nlens);

    FILE *fp;

    fp = fopen("data/_allTracks.dat", "w");
    for (_curve *c = track->first; c; c = c->next) { //curves
        fprintf(fp, "%d ", c->length);
    }
    fprintf(fp, "\n");

    for (_curve *c = track->first; c; c = c->next) { //curves
        for (_point *p = c->first; p; p = p->next) { //point
            fprintf(fp, "%f %f %f %f ", p->ds, p->dJ, p->x1, p->x2);// caustics x, caustics y, critical x, critical y
        }
    }
    fclose(fp);


}


// #define _c_parity \
// if (c->length>1){\
//             if (c->first->mu > 0 && c->first->next->mu > 0) {\
//                 c->parity = 1;\
//             }\
//             else if(c->first->mu < 0 && c->first->next->mu < 0) {\
//                 c->parity = -1;\
//             }else if (c->first->mu > 0 && c->first->next->mu < 0 && c->first->next->next->mu < 0){\
//                 c->parity = -1;\
//             }else if (c->first->mu < 0 && c->first->next->mu > 0 && c->first->next->next->mu > 0){\
//                 c->parity = 1;\
//             }\
// }else{\
//     c->parity = (c->first->mu > 0) ? 1:-1;\
// }

// revised on 2021.06.08
#define _c_parity \
if (c->length>2){\
            if (c->first->mu > 0 && c->first->next->mu > 0) {\
                c->parity = 1;\
            }\
            else if(c->first->mu < 0 && c->first->next->mu < 0) {\
                c->parity = -1;\
            }else if ( ( (c->first->mu > 0 ? 1:-1) + (c->first->next->mu > 0 ? 1:-1) + (c->first->next->next->mu > 0 ? 1:-1))  > 0){\
                c->parity = 1;\
                special_flag = 1;\
            }else{\
                c->parity = -1;\
                special_flag = 1;\
            }\
}else{\
    c->parity = (c->first->mu > 0) ? 1:-1;\
}



#define _final_parity \
if (final->parity == 0){\
final->parity = (final->first->mu > 0 ? 1:-1);}\

#define _Prov_parity \
Prov->parity = (Prov->first->mu > 0 ? 1:-1);\

#define _Prov2_parity \
if(Prov2->length>1){\
if (Prov2->first->mu > 0 && Prov2->first->next->mu > 0) {\
    Prov2->parity = 1;\
}\
else if (Prov2->first->mu < 0 && Prov2->first->next->mu < 0) {\
    Prov2->parity = -1;\
}else if (Prov2->first->mu > 0 && Prov2->first->next->mu < 0 && Prov2->first->next->next->mu < 0) {\
    Prov2->parity = -1;\
}else if (Prov2->first->mu < 0 && Prov2->first->next->mu > 0 && Prov2->first->next->next->mu > 0) {\
    Prov2->parity = 1;\
}\
}else{\
    Prov2->parity = (Prov2->first->mu > 0) ? 1:-1;\
}


// _sols *outputTracks_v2_savehalf
_sols *TripleLensing::outputTracks_v2_savehalf(double xsCenter, double ysCenter, double rs, _linkedarray *PHI, _sols **prevstore) {
    nsolution = 0;
    int j, k, missing_sol_idx, worst_true_sol_idx, total_parity;
    double missing_sol_absdzs, missing_sol_mu , worst_tru_sol_absdzs, worst_tru_sol_mu;
    // double mu, lambda1, lambda2, thetaJ;
    double mu;

    bool imageFound;
    complex tmpz;

    _sols *allSolutions, *imageTracks;
    allSolutions = new _sols;
    imageTracks = new _sols;

    //auto allSolutions = make_unique<_sols>();
    //auto imageTracks = make_unique<_sols>();

    // unique_ptr<_sols> allSolutions(new _sols);
    // unique_ptr<_sols> imageTracks(new _sols);




    _curve *Prov;// = new _curve; // 21 Sep 11
    _curve *Prov2 = new _curve;
    _curve *tempProv;// = new _curve;

    _point *pisso, *missing_sol, *pscan1, *pscan2; // *worst_true_sol
    _point *scanpisso;
    // double SD, MD, CD;

    // connect all images, including fake images, initialize ten _curve object to save solutions from polynomial solving
    
    for (int i = 0; i < DEGREE; i++) {
        Prov = new _curve;
        allSolutions->append(Prov);
    }

    pntnum = PHI->length;

    if (ftime) {
        *prevstore = NULL;
        *prevstore = new _sols;
        Node *scan;
        scan = PHI->first;

        phi = scan->value;
        xs = xsCenter + rs * cos(phi);
        ys = ysCenter + rs * sin(phi);

        polynomialCoefficients(xs, ys, coefficients);
        VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

        // first point in all solution tracks
        tempProv = new _curve;
        (*prevstore)->append(tempProv);
        // tempProv = nullptr;
        // delete tempProv;
        Prov = allSolutions->first;
        // the head of all tracks
        nimages = 0;
        total_parity = 0;
        for (int i = 0; i < DEGREE; i++) {
            flag = trueSolution(xs, ys, zr[i], &mu);

            // polish using Newton method, 2021.09.18
            #ifdef POLISH_USE_NEWTON
            if (flag == 0 && absdzs < SOLEPS1e2){  // we only polish when it is necessary
             tmpz = complex(zr[i].re, zr[i].im);
             findCloseImages(this->mlens, this->zlens, xs, ys, &tmpz, &imageFound);
             if (imageFound && abs(tmpz - zr[i])< SOLEPS1e3) {
                zr[i] = tmpz;
                flag = trueSolution(xs, ys, zr[i], &mu);
                }
            }
            #endif


            nimages += flag;

            Prov->append(zr[i].re, zr[i].im);
            Prov->last->phi = phi;
            Prov->last->mu = mu;
            Prov->last->zs = complex(xs, ys);
            Prov->last->thetaJ = thetaJ;
            Prov->last->absdzs = absdzs;

            (*prevstore)->last->append(zr[i].re, zr[i].im);// (*prevstore) = multiple_curves, (*prevstore)->last = one curve, (*prevstore)->last->last = one point
            (*prevstore)->last->last->phi = phi;
            (*prevstore)->last->last->mu = mu;
            (*prevstore)->last->last->zs = complex(xs, ys);
            (*prevstore)->last->last->thetaJ = thetaJ;
            (*prevstore)->last->last->absdzs = absdzs;

            if (flag == 1) {      // true solution
                total_parity +=  (mu > 0 ? 1 : -1);
                Prov->last->flag = 1;
                (*prevstore)->last->last->flag = 1 ;

                // add on 2021.06.10
                Prov->posflagnum += 1;

#ifdef parabcorr
                dy = complex(-sin(phi), cos(phi)) * rs;
                // dz = (dy - conj(J1) * conj(dy)) / dJ.re;
                dz = (dy - J1 * conj(dy)) / dJ.re;
                ds = (imag(dy * dz * dz * J2) + rho2) / dJ.re;
                Prov->last->ds = ds;
                Prov->last->dz = dz;
                (*prevstore)->last->last->ds = ds;
                (*prevstore)->last->last->dz = dz;
#endif

            } else {    // false solution
                (*prevstore)->last->last->flag = -1 ;
                Prov->last->flag = -1 ;
            }
            Prov = Prov->next; // this might be wrong
        }

        // whether nimages = 4, 6, 8, or 10;
        if ( (nimages - NLENS) % 2 != 1 || total_parity != -2) { // check the number of solutions is correct or not
#ifdef VERBOSE
            fprintf(stderr, "1st \t\t nimages = %d, and parity = %d which are wrong\n\n", nimages, total_parity);
#endif
            // you need to find out one more solution, and change their flag to 1
            k = 0;
            missing_sol_idx = 0;
            worst_true_sol_idx = 0;
            missing_sol_absdzs = 1e3;
            worst_tru_sol_absdzs = -1e1;
            for (_curve *_tmp = allSolutions->first; _tmp; _tmp = _tmp->next) {
#ifdef VERBOSE
                fprintf(stderr, "k = %d, missing_sol->absdzs %.5e, missing_sol->mu %.5e\n", k, _tmp->first->absdzs,  _tmp->first->mu);
#endif
                if ( _tmp->first->flag == 1 &&  _tmp->first->absdzs > worst_tru_sol_absdzs) {
                    worst_tru_sol_absdzs =  _tmp->first->absdzs;
                    worst_true_sol_idx = k;
                    worst_tru_sol_mu =  _tmp->first->mu;
                }
                if ( _tmp->first->flag == -1 &&  _tmp->first->absdzs < missing_sol_absdzs) {
                    missing_sol_absdzs =  _tmp->first->absdzs;
                    missing_sol_mu = _tmp->first->mu;
                    missing_sol_idx = k;
                }
                k = k + 1;
            }
#ifdef VERBOSE
            fprintf(stderr , "missing_sol_idx = %d, missing_sol_mu = %f\n\n", missing_sol_idx, missing_sol_mu);
            fprintf(stderr , "worst_true_sol_idx = %d, worst_tru_sol_mu = %f\n\n", worst_true_sol_idx, worst_tru_sol_mu);
#endif

            if ((total_parity - ( worst_tru_sol_mu > 0 ? 1 : -1  )  ) == -2  && worst_tru_sol_absdzs > SOLEPS * 0.5) {
                // pscan1 = Prov->first;
                _curve *_tmp = allSolutions->first;
                pscan2 = (*prevstore)->last->first;
                k = 0;
                while (k != worst_true_sol_idx ) {
                    // pscan1 = pscan1->next;
                    _tmp = _tmp->next;
                    pscan2 = pscan2->next;
                    k += 1;
                }
                // 定位到了 missing solution，将其 flag 变为 1
                // pscan1->flag = 0;
                _tmp->first->flag = -1; // 0 or -1? -- -1!
                pscan2->flag = -1;
                nimages -= 1;
                _tmp->posflagnum -= 1;

            }
            else if ( (total_parity +  ( missing_sol_mu > 0 ? 1 : -1  )  ) == -2 && missing_sol_absdzs < 2 * SOLEPS) {
                // 不是最优的实现方法
                // pscan1 = Prov->first;
                _curve *_tmp = allSolutions->first;
                pscan2 = (*prevstore)->last->first;
                k = 0;
                while (k != missing_sol_idx) {
                    _tmp = _tmp->next;
                    pscan2 = pscan2->next;
                    k += 1;
                }
                // 定位到了 missing solution，将其 flag 变为 1
                // pscan1->flag = 1;
                _tmp->first->flag = 1;
                pscan2->flag = 1;
                nimages += 1;

                _tmp->posflagnum += 1;

            } else {
#ifdef VERBOSE
                fprintf(stderr, "\t\t\t ****** 1278, 1st, can not add or delete solution, \n");
#endif
            }
        }

        scan->nimages = nimages;
        scan = scan->next;

        for (j = 1; j < pntnum; j++) {
            phi = scan->value;
            xs = xsCenter + rs * cos(phi);
            ys = ysCenter + rs * sin(phi);

            // polynomialCoefficients(mlens, zlens, xs, ys, coefficients, NLENS, DEGREE);
            polynomialCoefficients(xs, ys, coefficients);
            VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

            delete Prov2;
            Prov2 = new _curve; // use to store the ten solutions in current phi
            tempProv = new _curve;
            (*prevstore)->append(tempProv);// prevstore 在这只是用于保存以前解出来的解, 后续 phi->extend 之后不用再重新解这些 lens equations; 每次的 10 个解作为一条curve 保存

        // tempProv = nullptr;
        // delete tempProv;

            nimages = 0;
            total_parity = 0;
            for (int i = 0; i < DEGREE; i++) {
                // flag = trueSolution(mlens, zlens, xs, ys, zr[i], &mu, &lambda1, &lambda2, &thetaJ, NLENS, &J1, &J2, &dJ, &J3);
                flag = trueSolution(xs, ys, zr[i], &mu);

            // polish using Newton method, 2021.09.18
            #ifdef POLISH_USE_NEWTON
            if (flag == 0 && absdzs < SOLEPS1e2){  // we only polish when it is necessary
             tmpz = complex(zr[i].re, zr[i].im);
             findCloseImages(this->mlens, this->zlens, xs, ys, &tmpz, &imageFound);
             if (imageFound && abs(tmpz - zr[i])< SOLEPS1e3) {
                zr[i] = tmpz;
                flag = trueSolution(xs, ys, zr[i], &mu);
                }
            }
            #endif

                nimages += flag;

                Prov2->append(zr[i].re, zr[i].im);
                Prov2->last->phi = phi;
                Prov2->last->mu = mu;
                Prov2->last->zs = complex(xs, ys);
                Prov2->last->thetaJ = thetaJ;
                Prov2->last->absdzs = absdzs;

                (*prevstore)->last->append(zr[i].re, zr[i].im);
                (*prevstore)->last->last->phi = phi;
                (*prevstore)->last->last->mu = mu;
                (*prevstore)->last->last->zs = complex(xs, ys);
                (*prevstore)->last->last->thetaJ = thetaJ;
                (*prevstore)->last->last->absdzs = absdzs;


                if (flag == 1) {      // true solution
                    total_parity += ( mu > 0 ? 1 : -1 );

                    Prov2->last->flag = +1 ;
                    (*prevstore)->last->last->flag = +1 ;

#ifdef parabcorr
                    dy = complex(-sin(phi), cos(phi)) * rs;
                    // dz = (dy - conj(J1) * conj(dy)) / dJ.re;
                    dz = (dy - J1 * conj(dy)) / dJ.re;
                    ds = (imag(dy * dz * dz * J2) + rho2) / dJ.re;
                    Prov2->last->ds = ds;
                    Prov2->last->dz = dz;
                    (*prevstore)->last->last->ds = ds;
                    (*prevstore)->last->last->dz = dz;
#endif

                } else {
                    Prov2->last->flag = -1;
                    (*prevstore)->last->last->flag = -1;
                }
            }


            // whether nimages = 4, 6, 8, or 10;
            if ( (nimages - NLENS) % 2 != 1 || total_parity != -2) { // check the number of solutions is correct or not
#ifdef VERBOSE
                fprintf(stderr, "\n2nd \t\t nimages = %d, and parity = %d which are wrong\n\n", nimages, total_parity);
                // you need to find out one more solution, and change their flag to 1
#endif

                k = 0;
                missing_sol_idx = 0;
                worst_true_sol_idx = 0;

                missing_sol_absdzs = 1e3;
                worst_tru_sol_absdzs = -1e1;
                for (missing_sol = Prov2->first; missing_sol; missing_sol = missing_sol->next) {
#ifdef VERBOSE
                    fprintf(stderr, "idx = %d,  absdzs = %.5e, mu = %.17g, (x,y) = %.17g, %.17g\n", k, missing_sol->absdzs, missing_sol->mu, missing_sol->x1, missing_sol->x2);
#endif
                    if (missing_sol->flag == 1 && missing_sol->absdzs > worst_tru_sol_absdzs) {
                        worst_tru_sol_absdzs = missing_sol->absdzs;
                        worst_true_sol_idx = k;
                        worst_tru_sol_mu = missing_sol->mu;
                    }

                    // find missing solution
                    if (missing_sol->flag == -1 && missing_sol->absdzs < missing_sol_absdzs) {
                        missing_sol_absdzs = missing_sol->absdzs;
                        missing_sol_mu = missing_sol->mu;
                        missing_sol_idx = k;
                    }
                    k = k + 1;
                }
#ifdef VERBOSE
                fprintf(stderr , "missing_sol_idx = %d, missing_sol_mu = %f\n\n", missing_sol_idx, missing_sol_mu);
                fprintf(stderr , "worst_true_sol_idx = %d, worst_tru_sol_mu = %f\n\n", worst_true_sol_idx, worst_tru_sol_mu);
#endif

                if ((total_parity - ( worst_tru_sol_mu > 0 ? 1 : -1  )  ) == -2 && worst_tru_sol_absdzs > SOLEPS * 0.5) {
#ifdef VERBOSE
                    fprintf(stderr, "2nd \t\t remove a false positive\n");
#endif
                    pscan1 = Prov2->first;
                    pscan2 = (*prevstore)->last->first;
                    k = 0;
                    while (k != worst_true_sol_idx ) {
                        pscan1 = pscan1->next;
                        pscan2 = pscan2->next;
                        k += 1;
                    }
                    pscan1->flag = -1;
                    pscan2->flag = -1; // not 0
                    nimages -= 1;
                }
                else if ( (total_parity  +  ( missing_sol_mu > 0 ? 1 : -1  )  ) == -2 && missing_sol_absdzs < 2 * SOLEPS) {
                    // if (  missing_sol_absdzs / worst_tru_sol_absdzs < MISS_SOL_THRESHOLD ) {
                    // 不是最优的实现方法
#ifdef VERBOSE
                    fprintf(stderr, "2nd \t\t add a missing solution\n");
#endif
                    pscan1 = Prov2->first;
                    pscan2 =  (*prevstore)->last->first;
                    k = 0;
                    while (k != missing_sol_idx) {
                        pscan1 = pscan1->next;
                        pscan2 = pscan2->next;
                        k += 1;
                    }
                    // 定位到了 missing solution，将其 flag 变为 1
                    pscan1->flag = 1;
                    pscan2->flag = 1;
                    nimages += 1;
                } else {

#ifdef VERBOSE
                    fprintf(stderr, "2nd \t\t nimages = %d, and parity = %d which are wrong, srcxy = (%.17g, %.15g)\n\n", nimages, total_parity, xs, ys);
                    k = 0;
                    for (missing_sol = Prov2->first; missing_sol; missing_sol = missing_sol->next) {
                        fprintf(stderr, "idx = %d,  absdzs = %.5e, mu = %.5e\n", k, missing_sol->absdzs, missing_sol->mu);
                    }
                    fprintf(stderr, "\t\t\t ****** 1411, 2nd, can not add or delete solution, \n");
#endif
                }
            }


            scan->nimages = nimages;
            scan = scan->next;

            // attach the current solutions (both true and false) positions to last ones with the closest distance
            for (Prov = allSolutions->first; Prov;  Prov = Prov->next) {
                Prov2->closest(Prov->last, &pisso);
                Prov2->drop(pisso);
                Prov->append(pisso);
                if (pisso->flag == 1) {Prov->posflagnum += 1;}
            }
                delete Prov2;
                Prov2 = nullptr; // 2021.09.26
        }
#ifdef VERBOSE
        fprintf(stderr, "*prevstore->length %d\n", (*prevstore)->length);
#endif

    } else {
        //firsttime = 0 // when it is possible, use prerviously obtained solutions
#ifdef VERBOSE
        fprintf(stderr, "not first time\n");
#endif
        _sols *tempcurrstore;
        tempcurrstore = new _sols;

        Node *scan;
        scan = PHI->first;
        _curve *scanstorecurve;
#ifdef VERBOSE
        fprintf(stderr, "*prevstore->length in not firsttime %d\n", (*prevstore)->length);
#endif

        scanstorecurve = (*prevstore)->first;
#ifdef VERBOSE
        fprintf(stderr, "scanstorecurve->length = %d\n", scanstorecurve->length);
#endif

        phi = scan->value;
        xs = xsCenter + rs * cos(phi);
        ys = ysCenter + rs * sin(phi);


        // first point in all solution tracks
        tempProv = new _curve;
        tempcurrstore->append(tempProv);
        // tempProv = nullptr;
        // delete tempProv;

        //Prov = new _curve;
        Prov = allSolutions->first;
        nimages = 0;
        total_parity = 0;
        scanpisso = scanstorecurve->first;

        for (int i = 0; i < DEGREE; i++) {
            flag = (int)((scanpisso->flag + 1) / 2.0);
            nimages += flag;

            Prov->append(scanpisso->x1, scanpisso->x2);
            Prov->last->phi = scanpisso->phi;
            Prov->last->mu = scanpisso->mu;
            Prov->last->zs = scanpisso->zs; //complex(xs, ys);
            Prov->last->thetaJ = scanpisso->thetaJ;
            Prov->last->flag = scanpisso->flag;
            Prov->last->absdzs = scanpisso->absdzs;

            Prov->last->ds = scanpisso->ds;
            Prov->last->dz = scanpisso->dz;

            if (scanpisso->flag) {Prov->posflagnum = 1;}

            tempcurrstore->last->append(scanpisso->x1, scanpisso->x2);
            tempcurrstore->last->last->phi = scanpisso->phi;
            tempcurrstore->last->last->mu = scanpisso->mu;
            tempcurrstore->last->last->zs = scanpisso->zs; //complex(xs, ys);
            tempcurrstore->last->last->thetaJ = scanpisso->thetaJ;
            tempcurrstore->last->last->flag = scanpisso->flag;
            tempcurrstore->last->last->absdzs = scanpisso->absdzs;

            tempcurrstore->last->last->ds = scanpisso->ds;
            tempcurrstore->last->last->dz = scanpisso->dz;

            Prov = Prov->next;
            scanpisso = scanpisso->next;
        }
        scan->nimages = nimages;
        scan = scan->next;

        for (j = 1; j < pntnum; j++) {
            phi = scan->value;
            xs = xsCenter + rs * cos(phi);
            ys = ysCenter + rs * sin(phi);

            if ( (j % 2) == 1 ) {
                // polynomialCoefficients(mlens, zlens, xs, ys, coefficients, NLENS, DEGREE);
                polynomialCoefficients(xs, ys, coefficients);
                VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);
                delete Prov2;
                Prov2 = new _curve;
                tempProv = new _curve;
                tempcurrstore->append(tempProv);
        // tempProv = nullptr;
        // delete tempProv;

                nimages = 0;
                total_parity = 0;
                for (int i = 0; i < DEGREE; i++) {
                    // flag = trueSolution(mlens, zlens, xs, ys, zr[i], &mu, &lambda1, &lambda2, &thetaJ, NLENS, &J1, &J2, &dJ, &J3);
                    flag = trueSolution(xs, ys, zr[i], &mu);

            // polish using Newton method, 2021.09.18
            #ifdef POLISH_USE_NEWTON
            if (flag == 0 && absdzs < SOLEPS1e2){  // we only polish when it is necessary
             tmpz = complex(zr[i].re, zr[i].im);
             findCloseImages(this->mlens, this->zlens, xs, ys, &tmpz, &imageFound);
             if (imageFound && abs(tmpz - zr[i])< SOLEPS1e3) {
                zr[i] = tmpz;
                flag = trueSolution(xs, ys, zr[i], &mu);
                }
            }
            #endif

                    nimages += flag;
                    Prov2->append(zr[i].re, zr[i].im);
                    Prov2->last->phi = phi;

                    Prov2->last->mu = mu;
                    Prov2->last->zs = complex(xs, ys);
                    Prov2->last->thetaJ = thetaJ;
                    Prov2->last->absdzs = absdzs;

                    tempcurrstore->last->append(zr[i].re, zr[i].im);
                    tempcurrstore->last->last->phi = phi;
                    tempcurrstore->last->last->mu = mu;
                    tempcurrstore->last->last->zs = complex(xs, ys);
                    tempcurrstore->last->last->thetaJ = thetaJ;
                    tempcurrstore->last->last->absdzs = absdzs;

                    if (flag == 1) {      // true solution
                        total_parity += (mu > 0 ? 1 : -1);
                        Prov2->last->flag = +1 ;
                        tempcurrstore->last->last->flag = +1 ;
#ifdef parabcorr
                        dy = complex(-sin(phi), cos(phi)) * rs;
                        // dz = (dy - conj(J1) * conj(dy)) / dJ.re;
                        dz = (dy - J1 * conj(dy)) / dJ.re;
                        ds = (imag(dy * dz * dz * J2) + rho2) / dJ.re;
                        Prov2->last->ds = ds;
                        Prov2->last->dz = dz;
                        tempcurrstore->last->last->ds = ds;
                        tempcurrstore->last->last->dz = dz;
#endif

                    } else {
                        Prov2->last->flag = -1;
                        tempcurrstore->last->last->flag = -1;
                    }
                }

                // whether nimages = 4, 6, 8, or 10;
                if ( (nimages - NLENS) % 2 != 1 || total_parity != -2) { // check the number of solutions is correct or not
#ifdef VERBOSE
                    fprintf(stderr, "3rd \t\t nimages = %d, and parity = %d which are wrong\n\n", nimages, total_parity);
#endif
                    // you need to find out one more solution, and change their flag to 1

                    k = 0;
                    missing_sol_idx = 0;
                    worst_true_sol_idx = 0;

                    missing_sol_absdzs = 1e3;
                    worst_tru_sol_absdzs = -1e1;
                    for (missing_sol = Prov2->first; missing_sol; missing_sol = missing_sol->next) {
                        if (missing_sol->flag == 1 && missing_sol->absdzs > worst_tru_sol_absdzs) {
                            worst_tru_sol_absdzs = missing_sol->absdzs;
                            worst_true_sol_idx = k;
                            worst_tru_sol_mu = missing_sol->mu;
                        }

                        if (missing_sol->flag == -1 && missing_sol->absdzs < missing_sol_absdzs) {
                            missing_sol_absdzs = missing_sol->absdzs;
                            missing_sol_idx = k;
                            missing_sol_mu = missing_sol->mu;
                        }
                        k = k + 1;
                    }

                    if ((total_parity - ( worst_tru_sol_mu > 0 ? 1 : -1  )  ) == -2 && worst_tru_sol_absdzs > SOLEPS * 0.5) {
                        pscan1 = Prov2->first;
                        pscan2 = tempcurrstore->last->first;
                        k = 0;
                        while (k != worst_true_sol_idx ) {
                            pscan1 = pscan1->next;
                            pscan2 = pscan2->next;
                            k += 1;
                        }
                        // 定位到了 missing solution，将其 flag 变为 1
                        pscan1->flag = -1;
                        pscan2->flag = -1;
                        nimages -= 1;
                    }


                    else if ( (total_parity +  ( missing_sol_mu > 0 ? 1 : -1  )  ) == -2 && missing_sol_absdzs < 2 * SOLEPS) {
                        // 不是最优的实现方法
                        pscan1 = Prov2->first;
                        pscan2 = tempcurrstore->last->first;
                        k = 0;
                        while (k != missing_sol_idx) {
                            pscan1 = pscan1->next;
                            pscan2 = pscan2->next;
                            k += 1;
                        }
                        // 定位到了 missing solution，将其 flag 变为 1
                        pscan1->flag = 1;
                        pscan2->flag = 1;
                        nimages += 1;
                    } else {
#ifdef VERBOSE
                        fprintf(stderr, "\t\t\t ****** 1613, 3rd, can not add or delete solution, \n");
#endif
                    }
                }

                scan->nimages = nimages;
                scan = scan->next;

                // attach the current solutions positions to last ones with the closest distance
                for (Prov = allSolutions->first; Prov;  Prov = Prov->next) {
                    Prov2->closest(Prov->last, &pisso);
                    Prov2->drop(pisso);
                    Prov->append(pisso);
                    if (pisso->flag == 1) {Prov->posflagnum += 1;}
                }
                delete Prov2;
                Prov2 = nullptr; // 2021.09.26

            } else {
                delete Prov2;
                Prov2 = new _curve;
                tempProv = new _curve;
                tempcurrstore->append(tempProv);
                // tempProv = nullptr;
                // delete tempProv;

                nimages = 0;
                scanstorecurve = scanstorecurve->next;
                scanpisso = scanstorecurve->first;
                for (int i = 0; i < DEGREE; i++) {
                    flag = (int)((scanpisso->flag + 1) / 2);
                    nimages += flag;
                    Prov2->append(scanpisso->x1, scanpisso->x2);
                    Prov2->last->phi = scanpisso->phi;
                    Prov2->last->mu = scanpisso->mu;
                    Prov2->last->zs = scanpisso->zs; //complex(xs, ys);
                    Prov2->last->thetaJ = scanpisso->thetaJ;
                    Prov2->last->flag = scanpisso->flag;
                    Prov2->last->absdzs = scanpisso->absdzs;

                    Prov2->last->ds = scanpisso->ds;
                    Prov2->last->dz = scanpisso->dz;

                    tempcurrstore->last->append(scanpisso->x1, scanpisso->x2);
                    tempcurrstore->last->last->phi = scanpisso->phi;
                    tempcurrstore->last->last->mu = scanpisso->mu;
                    tempcurrstore->last->last->zs = scanpisso->zs; //complex(xs, ys);
                    tempcurrstore->last->last->thetaJ = scanpisso->thetaJ;
                    tempcurrstore->last->last->flag = scanpisso->flag;
                    tempcurrstore->last->last->absdzs = scanpisso->absdzs;

                    tempcurrstore->last->last->ds = scanpisso->ds;
                    tempcurrstore->last->last->dz = scanpisso->dz;

                    scanpisso = scanpisso->next;
                }
                scan->nimages = nimages;
                scan = scan->next;

                // attach the current solutions positions to last ones with the closest distance
                for (Prov = allSolutions->first; Prov;  Prov = Prov->next) {
                    Prov2->closest(Prov->last, &pisso);
                    Prov2->drop(pisso);
                    Prov->append(pisso);
                    if (pisso->flag == 1) {Prov->posflagnum += 1;}
                }
                delete Prov2;
                Prov2 = nullptr; // 2021.09.26
            }

        }
        delete (*prevstore);//200711
        //(*prevstore) = new _sols;;//200711

        (*prevstore) = tempcurrstore;
        tempcurrstore = nullptr;
        delete tempcurrstore;
    }

#ifdef VERBOSE
    printAllTracks(allSolutions);
    printf("allSolutions->length %d\n", allSolutions->length);
    k = 0;
    for (_curve *tmpsc = allSolutions->first; tmpsc; tmpsc = tmpsc->next) {
        fprintf(stderr, "\t\t %d, track->length %d, posflagnum %d, percentage %f  , \n", k, tmpsc->length, tmpsc->posflagnum, 100.0 * tmpsc->posflagnum / tmpsc->length);
        k = k + 1;
    }
#endif




    //allocate memory for potentially true image tracks, we make a factor four more segments just in case
    // imageTracks = new _sols;

    //2021.06.10, before proceed, we directly drop out pure false image curves;
    // this must be slightly faster
    // Prov = allSolutions->first;
    // while (Prov) {
    //     Prov2 = Prov->next;
    //     // if (1.0*Prov->posflagnum/Prov->length < 0.01) {
    //     if (Prov->posflagnum == 0) {
    //         allSolutions->drop(Prov);
    //         delete Prov;
    //         Prov = Prov2;
    //     } else {Prov = Prov->next;}
    // }

    //_curve *tmp1802;//

    Prov = allSolutions->first;
    while (Prov) {
        if (Prov->posflagnum <= 1 || ( (1.0 * Prov->posflagnum / Prov->length < 0.5) && abs(Prov->first->mu) < EPS  && abs(Prov->last->mu) < EPS )  ) {
            // the second condition is for a tiny image, where the flag can not be correctly judged
            // what's the area of this tiny image, compared to the source area?

            Prov2 = Prov->next;
            allSolutions->drop(Prov);
            delete Prov;
            Prov = Prov2;
        } else if (1.0 * Prov->posflagnum / Prov->length > 0.99) {
            // append Prov to imageTracks
            Prov2 = Prov->next;

            // tmp1802 = Prov;
            // allSolutions->drop(Prov);
            // imageTracks->append(tmp1802);

            allSolutions->drop(Prov);
            imageTracks->append(Prov);

            Prov = Prov2;
        }
        else {Prov = Prov->next;}
    }

#ifdef VERBOSE
    for (_curve *ttmp = imageTracks->first; ttmp; ttmp = ttmp->next) {
        fprintf(stderr, "1811 imageTracks ttmp length = %d, posflagnum = %d \n", ttmp->length, ttmp->posflagnum);
    }
    for (_curve *ttmp = allSolutions->first; ttmp; ttmp = ttmp->next) {
        fprintf(stderr, "1815 allSolutions ttmp length = %d, posflagnum = %d \n", ttmp->length, ttmp->posflagnum);
    }
#endif

    if (allSolutions->length == 0)  {
        for (Prov = imageTracks->first; Prov; Prov = Prov->next) {
            _Prov_parity
        }
        Prov2 = NULL;
        Prov = NULL;
        delete Prov2;
        delete Prov;
        delete allSolutions;
        return (imageTracks);
    }



    //
    // Solution tracks may contain only true images, some only false solutions, while some may be mixed
    //
    complex z, zs;
    int mixedTracks = 0; // number of mixed tracks
    _curve *tempcurve;


    int trueImages = 0, falseImages = 0, total_scaned, pos_scaned, scan_percentage;
    int firstTime = 1;
    int disfirstTime = 1, flagPrevious;
    int previousImage = 0;
    // int Parity;  // within the same curve, we require the parity to be the same
    // double ratio;
    double muPrevious , disPrevious;//mu
    scan_percentage = fmin(20, 0.05 * pntnum); // 10% of toal length



//     //  ************************************** begin version 1 *************************************** //
//     // Prov2 = imageTracks->first;
//     tempcurve = new _curve;
//     imageTracks->append(tempcurve);
//     // use a new _curve to store the remaining segments
//     Prov2 = imageTracks->last;
//     Prov = allSolutions->first;
//     int fistpls = 1;
//     double helpangcos = 0.0;
//     while (Prov) {
// #ifdef VERBOSE
//         fprintf(stderr, "now you only need to proceed this posflagnum %d, posflagnum/length %f, first->flag %d, mu %f \n", Prov->posflagnum, 1.0 * Prov->posflagnum / Prov->length, Prov->first->flag, Prov->first->mu);
// #endif
//         trueImages = 0;
//         falseImages = 0;
//         firstTime = 1;
//         disfirstTime = 1;
//         previousImage = 0;
//         for (_point *p = Prov->first; p; p = p->next) {
//             if (p->flag == 1) { // true solutions, should return identical source radius
//                 trueImages++;
//                 if (firstTime) {  // we enter the track the first time
//                     firstTime = 0;
//                     flagPrevious = p->flag;
//                     muPrevious = p->mu;
//                 } else {
//                     // parity has changed, start a new curve
//                     // actually there should be a furthur more test: test that whether p and p->next is close enough, sometimes they will be very far away from each other, in this case we should break the current Prov and start a new one.
//                     // now this situation is done by adding a testing afterwards, see the "imageTracks_checked" below
//                     // if (flagPrevious * p->flag < 0.0 && muPrevious * p->mu < 0.0)  {
//                     if ( muPrevious * p->mu < 0.0)  {
//                         // add one more condition about the angle, when crossing the caustics
//                         helpangcos = 1;
//                         if (p->next && p->prev) {
//                             if (p->next->next && p->prev->prev) {
//                                 helpangcos = angcos(p->prev->prev, p->prev, p->next, p->next->next);
//                             }
//                         }
//                         fprintf(stderr, "\t and p->prev->x1 x2 = %f, %f, parity = %d, mu = %f\n", p->prev->x1, p->prev->x2, p->prev->flag, p->prev->mu);
//                         fprintf(stderr, "\t break segment, helpangcos = %f, p->x1 x2 = %f, %f, parity = %d, mu = %f\n", helpangcos, p->x1, p->x2, p->flag, p->mu);
//                         fprintf(stderr, "\t and p->next->x1 x2 = %f, %f, parity = %d, mu = %f\n", p->next->x1, p->next->x2, p->next->flag, p->next->mu);

//                         // maybe we can also use the sign information of helpangcos ?
//                         if (  abs(helpangcos) < 0.95) {
//                             printf("\t --> start a new segment\n");
//                             if (!Prov2->next) {
//                                 tempcurve = new _curve;
//                                 imageTracks->append(tempcurve);
//                             }
//                             Prov2 = Prov2->next;
//                             flagPrevious = p->flag;
//                             p->mu = muPrevious > 0 ? abs(p->mu):-abs(p->mu); //
//                             // muPrevious = p->mu;
//                             disfirstTime = 1;
//                         }
//                     }
//                     // begin ************************
//                     if (disfirstTime && Prov2->length > 1) {
//                         disPrevious = *(Prov2->last) - *(Prov2->last->prev); // sqrt(dis)
//                         disfirstTime = 0;
//                     }
//                     if (Prov2->length > 1 && ((*(p) - * (Prov2->last)) >=  1e8 * disPrevious) ) {
//                         if (!Prov2->next) {
//                             tempcurve = new _curve;
//                             imageTracks->append(tempcurve);
//                         }
//                         Prov2 = Prov2->next;
//                         flagPrevious = p->flag;
//                         muPrevious = p->mu;
//                         disfirstTime = 1;
//                     }
//                     // end ************************
//                 }

//                 addPoint(Prov2, p); // add point p to the track
//                 //2021.06.11 you need to update:
//                 flagPrevious = p->flag;
//                 muPrevious = p->mu;
//                 previousImage = TRUE_IMAGE;
//             } else {
//                 falseImages++;
//                 //first time that we enter this
//                 if (firstTime) {
//                     previousImage = FALSE_IMAGE;
//                     firstTime = 0;
//                     flagPrevious = p->flag; // 2021.06.10
//                     muPrevious = p->mu;
//                 }
//                 // already a true segment existing, we cut the segment, and start a new one
//                 if (previousImage == TRUE_IMAGE) { // previously we had true images, now we are encountering on false solutions
//                     if (!Prov2->next) {
//                         tempcurve = new _curve;
//                         imageTracks->append(tempcurve);
//                         disfirstTime = 1;
//                     }
//                     Prov2 = Prov2->next;  // we start a new image track
//                     previousImage = FALSE_IMAGE;
//                 }
//             }
//         }
//         if (trueImages) {
//             if (!Prov2->next) {
//                 tempcurve = new _curve;
//                 imageTracks->append(tempcurve);
//                 disfirstTime = 1;
//             }
//             Prov2 = Prov2->next;
//         }
//         if (trueImages && falseImages) mixedTracks += 1; // there are tracks that are mixed
//         fprintf(stderr, "can you here 1911 \n");
//         Prov = Prov->next;
//     }
//     // finish spliting the segments.

//     //  ****************************************** end version 1 ************************************************* //


    //rewrite the split procedure, 2021.06.11
//     //  ****************************************** begin version 2 ************************************************* //
//     _point *p2;
//     // Prov2 = imageTracks->first;
//     tempcurve = new _curve;
//     imageTracks->append(tempcurve);
//     // use a new _curve to store the remaining segments
//     Prov2 = imageTracks->last;
//     Prov = allSolutions->first;
//     // for (Prov = allSolutions->first; Prov;  Prov = Prov->next) {
//     int fistpls = 1, enterfirstTime = 1 , need_new_seg, currentTrueOrFalse; // when currentTrueOrFalse = 1, we add the point p,  otherwise we just ignore p
//     double helpangcos = 0.0;
//     _point *runner, *anchor;
//     int can_percentage;
//     while (Prov) {
//         mixedTracks = 1; // if you can go inside this loop, means you have _curve in allSolutions that contains both +1 and -1 parity points
// #ifdef VERBOSE
//         fprintf(stderr, "now you only need to proceed this posflagnum %d, posflagnum/length %f \n", Prov->posflagnum, 1.0 * Prov->posflagnum / Prov->length);
// #endif
//         firstTime = 1;
//         need_new_seg = 0;
//         disfirstTime = 1;
//         previousImage = 0;
//         currentTrueOrFalse = 0; // currently, we are scan at a point with +1 flag or -1 flag
//         p2 = Prov->first;
//         // fprintf(stderr, "can you here 2076, p2->flag = %d?\n", p2->flag);
//         while (p2->next && p2->flag != 1) {
//             p2 = p2->next;
//         }
//         while (p2) {
//             if (need_new_seg == 1 && firstTime == 0) {
//                 if (!Prov2->next) {
//                     tempcurve = new _curve;
//                     imageTracks->append(tempcurve);
//                 }
//                 Prov2 = Prov2->next;
//             }
//             if (firstTime) {  // we enter the track the first time
//                 firstTime = 0;
//                 flagPrevious = p2->flag;
//                 muPrevious = p2->mu;
//                 currentTrueOrFalse = 1;//p2->flag > 0 ? 1 : 0; // this must be 1
//             } else {
//                 currentTrueOrFalse = 1;
//                 // parity has changed, start a new curve
//                 // actually there should be a furthur more test: test that whether p and p->next is close enough, sometimes they will be very far away from each other, in this case we should break the current Prov and start a new one.
//                 // now this situation is done by adding a testing afterwards, see the "imageTracks_checked" below
//                 if (flagPrevious * p2->flag < 0.0 || muPrevious * p2->mu < 0.0)  { // use or, because sometimes +1 flag has neg mu!!!
//                     // if previously we have a true segment, and now p satisfy the condition to break

//                     // we need to decide whether this is a real place to break
//                     // we test by check whether most of the following points have  parity -1
//                     total_scaned = 0; // total number of scanned point
//                     pos_scaned = 0; // total number of positive flag in it;
//                     anchor = p2;

//                     while (p2->next) {
//                         if (total_scaned < scan_percentage) {
//                             total_scaned += 1;
//                             pos_scaned += ( p2->flag == 1 ? 1 : 0 );
//                             p2 = p2->next;
//                         } else {break;}
//                     }

//                     if (!(p2 ->next )) {
//                         #ifdef VERBOSE
//                         fprintf(stderr, "2124 break from here\n");
//                         #endif
//                         break;
//                     }

//                     #ifdef VERBOSE
//                     fprintf(stderr, "\t total_scaned = %d, pos_scaned = %d, percentage = %f\n", total_scaned, pos_scaned,  1.0 * pos_scaned / total_scaned * 100);
//                     #endif
//                     if (pos_scaned == total_scaned){
//                         // we can add all these points to Prov2
//                         while(anchor != p2->prev){ // leave p2 to add by the below code
//                             addPoint(Prov2, anchor);
//                             need_new_seg = 0;
//                             anchor = anchor->next;
//                         }

//                         // p2 = p2->next; //
//                         flagPrevious = 1;
//                         muPrevious = p2->mu;
//                         currentTrueOrFalse = 1;
//                         #ifdef VERBOSE
//                         fprintf(stderr, "2149 here\n");
//                         #endif

//                         // continue;
//                     }else if (pos_scaned == 0){
//                         // we can discard them totally
//                         #ifdef VERBOSE
//                         fprintf(stderr, "2154 here\n");
//                         #endif
//                         currentTrueOrFalse = 0;
//                         // p2 = p2->next; // this will do in the following if    currentTrueOrFalse = 0;
//                         // do not need to update flagPrevious and muPrevious
//                         // continue;// do not continue, you can do better, you can scan from p2 to check whether the whole segments are -1 flag
//                         // muPrevious, flagPrevious remains the same?
//                         flagPrevious = -1; // must be 1 -- bug !!!!
//                         muPrevious = p2->mu;
//                     }else if (1.0 * pos_scaned / total_scaned > 0.8) {
//                         #ifdef VERBOSE
//                         fprintf(stderr, "2163 here\n");
//                         #endif
//                         // this point is safe
//                         currentTrueOrFalse = 1;
//                         p2 = anchor; // this go first
//                         muPrevious = p2->mu;
//                         flagPrevious = p2->flag; // assuming we have +1 flag
//                     } else { // we check from this point until it reaches another place of a true segments
//                         // we should break at here
//                         #ifdef VERBOSE
//                         fprintf(stderr, "2171 here\n");
//                         #endif
//                         currentTrueOrFalse = 0;
//                         p2 = anchor;
//                         muPrevious = p2->mu;
//                         flagPrevious = p2->flag;
//                     }
// //                     // add one more condition about the angle, when crossing the caustics
// //                     helpangcos = 1;
// //                     if (p2->next && p2->prev) {
// //                         if (p2->next->next && p2->prev->prev) {
// //                             helpangcos = angcos(p2->prev->prev, p2->prev, p2->next, p2->next->next);
// //                         }
// //                     }
// // #ifdef VERBOSE
// //                     fprintf(stderr, "\t helpangcos = %f\n", helpangcos);
// //                     fprintf(stderr, "\t and p->prev->x1 x2 = %f, %f, parity = %d, mu = %f, phi = %f \n", p2->prev->x1, p2->prev->x2, p2->prev->flag, p2->prev->mu, p2->prev->phi * 180 / M_PI);
// //                     fprintf(stderr, "\t break segment, p->x1 x2 = %f, %f, parity = %d, mu = %f, phi = %f \n", p2->x1, p2->x2, p2->flag, p2->mu, p2->phi * 180 / M_PI);
// //                     fprintf(stderr, "\t and p->next->x1 x2 = %f, %f, parity = %d, mu = %f, phi = %f \n", p2->next->x1, p2->next->x2, p2->next->flag, p2->next->mu, p2->next->phi * 180 / M_PI);
// // #endif

// //                     if (p2->prev->flag == 1 && p2->flag == -1 && p2->next->flag == 1) {
// //                         fprintf(stderr, "\t save one at 2128, \n");
// //                         currentTrueOrFalse = 1;
// //                         // this point can connect, although it has a flag -1;
// //                         flagPrevious = 1;
// //                         muPrevious = p2->next->mu; // use p2->next rather than p2
// //                     } else if (p2->prev->flag == -1 && p2->flag == 1 && p2->next->flag == -1 && helpangcos > 0.8) {
// //                         fprintf(stderr, "\t save one at 2128, \n");
// //                         currentTrueOrFalse = 1;
// //                         // this point can connect, although it has a flag -1;
// //                         flagPrevious = -1;
// //                         muPrevious = p2->next->mu; // use p2->next rather than p2
// //                     }
// //                     // maybe we can also use the sign information of helpangcos ?
// //                     else if (  abs(helpangcos) < 0.95 && p2->next->flag != 1) { // sometime its is just one flag@-1 point between two flad@1 point
// //                         // if p2->flag == -1, this means that now it is very dangerous

// //                         // if (1) {
// // // #ifdef VERBOSE
// // //                         fprintf(stderr, "\t and p->prev->x1 x2 = %f, %f, parity = %d, mu = %f\n", p2->prev->x1, p2->prev->x2, p2->prev->flag, p2->prev->mu);
// // //                         fprintf(stderr, "\tbreak segment, helpangcos = %f, p->x1 x2 = %f, %f, parity = %d, mu = %f\n", helpangcos, p2->x1, p2->x2, p2->flag, p2->mu);
// // //                         fprintf(stderr, "\t and p->next->x1 x2 = %f, %f, parity = %d, mu = %f\n", p2->next->x1, p2->next->x2, p2->next->flag, p2->next->mu);
// // // #endif
// //                         if (!Prov2->next) {
// //                             tempcurve = new _curve;
// //                             imageTracks->append(tempcurve);
// //                         }
// //                         Prov2 = Prov2->next;
// //                         flagPrevious = p2->flag;
// //                         muPrevious = p2->mu;
// //                         disfirstTime = 1;
// //                         currentTrueOrFalse = 0; // now we are at a false segment
// //                     }
//                 } else {
//                     // begin ************************
//                     if (disfirstTime && Prov2->length > 1) {
//                         disPrevious = *(Prov2->last) - *(Prov2->last->prev); // sqrt(dis)
//                         disfirstTime = 0;
//                     }
//                     if (Prov2->length > 1 && ((*(p2) - * (Prov2->last)) >=  1e8 * disPrevious) ) { // if the next point is too far, we also break it
//                         fprintf(stderr, "\t special case 2158\n");
//                         if (!Prov2->next) {
//                             tempcurve = new _curve;
//                             imageTracks->append(tempcurve);
//                         }
//                         Prov2 = Prov2->next;
//                         flagPrevious = p2->flag;
//                         muPrevious = p2->mu;
//                         disfirstTime = 1;
//                         currentTrueOrFalse = p2->flag > 0 ? 1 : 0;
//                     }
//                     // end ************************
//                 }
//             }
//             if (currentTrueOrFalse) {
//                 addPoint(Prov2, p2);
//                 need_new_seg = 0;
//                 p2 = p2->next;
//                 continue;
//             } // add point p to the track only when we are at a true segment
//             else {
//                 need_new_seg = 1; // some point can not be add to current segment, so you need a new Prov2 to store the follwing segments.
//                 // now since you have enconter a "fasle" point, you can scan and check where you can jump to

//                 // fprintf(stderr, "2221, p->flag %d, flagPrevious %d, p2->next->flag %d\n", p2->flag, flagPrevious, p2->next->flag);


//                 while (p2->next) { // not while (p2), because otherwise you will end in inquiring p2->flag, i.e. 0->flag, which is wrong
//                     // if (p2->flag == muPrevious){ // wrong , not muPrevious!
//                     // fprintf(stderr, "p2->flag %d, p2->mu %fm flagPrevious = %d\n", p2->flag, p2->mu, flagPrevious);
//                     p2 = p2->next;
//                     if (p2->flag != flagPrevious) {
//                         #ifdef VERBOSE
//                         fprintf(stderr, "2251, p->flag %d, flagPrevious %d\n", p2->flag, flagPrevious);
//                         #endif
//                         flagPrevious = p2->flag;
//                         muPrevious = p2->mu;
//                         break;
//                     }
//                 }
//                 // flagPrevious = 1; // p2->flag;
//                 // muPrevious = 1; //p2->mu;
//                 // continue;
//             }

//             if (!(p2 ->next )) {
//                 #ifdef VERBOSE
//                 fprintf(stderr, "2241 break from here\n");
//                 #endif
//                 break;
//             }
//             // }
//             // else {
//             //     falseImages++;
//             //     //first time that we enter this
//             //     if (firstTime) {
//             //         previousImage = FALSE_IMAGE;
//             //         firstTime = 0;
//             //         flagPrevious = p->flag; // 2021.06.10
//             //         muPrevious = p->mu;
//             //     }
//             //     // already a true segment existing, we cut the segment, and start a new one
//             //     if (previousImage == TRUE_IMAGE) { // previously we had true images, now we are encountering on false solutions
//             //         if (!Prov2->next) {
//             //             tempcurve = new _curve;
//             //             imageTracks->append(tempcurve);
//             //             disfirstTime = 1;
//             //         }
//             //         Prov2 = Prov2->next;  // we start a new image track
//             //         previousImage = FALSE_IMAGE;
//             //     }
//             // }
//         }
//         // if (trueImages) {
//         //     if (!Prov2->next) {
//         //         tempcurve = new _curve;
//         //         imageTracks->append(tempcurve);
//         //         disfirstTime = 1;
//         //     }
//         //     Prov2 = Prov2->next;
//         // }

//         // if (trueImages && falseImages) mixedTracks += 1; // there are tracks that are mixed

//         Prov = Prov->next;
//         // a Prov is checked, start a new Prov2 to store
//         if (!Prov2->next) {
//             tempcurve = new _curve;
//             imageTracks->append(tempcurve);
//         }
//         Prov2 = Prov2->next;
//     }
// finish spliting the segments.
//     //  ****************************************** end version 2 ************************************************* //



    //  ************************************** begin version 3 *************************************** //
#ifdef VERBOSE
    Prov = allSolutions->first;
    k = 0;
    while (Prov) {
        fprintf(stderr, ">>>>>><<<<<<<>>>>>>> %d, first->flag, mu = %d, %f, last->flag, mu = %d, %f\n", k,  Prov->first->flag, Prov->first->mu, Prov->last->flag, Prov->last->mu);
        fprintf(stderr, ">>>>>><<<<<<<>>>>>>> %d, first->x1,x2 %f, %f, last->x1, x2 = %f, %f\n", k,  Prov->first->x1, Prov->first->x2, Prov->last->x1, Prov->last->x2);
        k += 1;
        Prov = Prov->next;
    }
#endif

    tempcurve = new _curve;
    imageTracks->append(tempcurve);
    // tempcurve = nullptr;
    // delete tempcurve;

    // use a new _curve to store the remaining segments
    Prov2 = imageTracks->last;
    Prov = allSolutions->first;
    _point *anchor, *p2;
    int  have_add = 0, currentTrueOrFalse = 1; //fistpls = 1,
    double helpangcos = 0.0;
    while (Prov) {
#ifdef VERBOSE
        fprintf(stderr, "now you only need to proceed this posflagnum %d, posflagnum/length %f, first->flag %d, mu %f, x1 x2 = %f, %f \n", Prov->posflagnum, 1.0 * Prov->posflagnum / Prov->length, Prov->first->flag, Prov->first->mu, Prov->first->x1, Prov->first->x2);
#endif
        trueImages = 0;
        falseImages = 0;
        firstTime = 1;
        disfirstTime = 1;
        previousImage = 0;
        p2 = Prov->first;
        while (p2) {
            if (firstTime) {  // we enter the track the first time
                firstTime = 0;
                flagPrevious = p2->flag;
                muPrevious = p2->mu;
                if (flagPrevious == 1) {
                    addPoint(Prov2, p2); // add point p to the track
                    previousImage = TRUE_IMAGE;
                    have_add = 1;
                    if (p2->next) {
                        muPrevious = p2->mu;
                        p2 = p2->next;
                        continue;
                    } else {
                        break;
                    }
                } else {
                    if (1.0 * Prov->posflagnum / Prov->length > 0.85) { // this is also not safe
                        addPoint(Prov2, p2); // add point p to the track
                        //2021.06.11 you need to update:
                        flagPrevious = 1;
                        previousImage = TRUE_IMAGE;
                        have_add = 1;
                        if (p2->next) {
                            muPrevious = p2->mu;
                            p2 = p2->next;
                            continue;
                        } else {
                            break;
                        }

                    } else {
                        firstTime = 1;
                        have_add = 0;
                        previousImage = FALSE_IMAGE;
                        while (p2->next && p2->next->flag == -1) {
                            p2 = p2->next;
                        }
                        muPrevious = p2->mu;
                    }
                }
            } else {
                currentTrueOrFalse = 1;
                // parity has changed, start a new curve
                // actually there should be a furthur more test: test that whether p and p2->next is close enough, sometimes they will be very far away from each other, in this case we should break the current Prov and start a new one.
                // now this situation is done by adding a testing afterwards, see the "imageTracks_checked" below
                if (flagPrevious * p2->flag < 0.0) {
                    if (muPrevious * p2->mu < 0.0) {
                        // add one more condition about the angle, when crossing the caustics
                        helpangcos = 1;
                        if (p2->next && p2->prev) {
                            if (p2->next->next && p2->prev->prev) {
                                helpangcos = angcos(p2->prev->prev, p2->prev, p2->next, p2->next->next);
#ifdef VERBOSE
                                fprintf(stderr, "\t 2433 and p2->prev->x1 x2 = %f, %f, parity = %d, mu = %f\n", p2->prev->x1, p2->prev->x2, p2->prev->flag, p2->prev->mu);
                                fprintf(stderr, "\t break segment, helpangcos = %f, p2->x1 x2 = %f, %f, parity = %d, mu = %f\n", helpangcos, p2->x1, p2->x2, p2->flag, p2->mu);
                                fprintf(stderr, "\t and p2->next->x1 x2 = %f, %f, parity = %d, mu = %f\n", p2->next->x1, p2->next->x2, p2->next->flag, p2->next->mu);
#endif

                                // maybe we can also use the sign information of helpangcos ?
                                if (  abs(helpangcos) < 0.95) {
#ifdef VERBOSE
                                    printf("\t 2441 --> start a new segment\n");
#endif
                                    if (Prov2->length > 0 && !Prov2->next) {
                                        // printf("\t 2535 --> start a new segment\n");
                                        tempcurve = new _curve;
                                        imageTracks->append(tempcurve);
                                        Prov2 = Prov2->next;
                                        //tempcurve = nullptr;
                                        //delete tempcurve;
                                    }

                                    // 在这做个额外的处理，后面大概率是false image，哪怕有也是只有几个了，可利用 Prov->posflagnum 的信息大概判断这是 大segment 还是 小 segment
                                    if (1.0 * Prov->posflagnum / Prov->length < 0.3 && flagPrevious == 1 && p2->flag == -1) {
                                        // you can be aggressive and drop more here
                                        // by scan more
                                        total_scaned = 0;
                                        pos_scaned = 0;

                                        anchor = p2;
                                        while (p2->next) {
                                            if (total_scaned < scan_percentage) {
                                                total_scaned += 1;
                                                pos_scaned += ( p2->flag == 1 ? 1 : 0 );
                                                p2 = p2->next;
                                            } else {break;}
                                        }
                                        if (!(p2 ->next )) {
                                            break;
                                        }
                                        if (1.0 * pos_scaned / total_scaned < 0.8) {

                                            while (p2->next && p2->next->flag == -1) {
                                                p2 = p2->next;
                                            }
                                            flagPrevious = -1; // now p2->next->flag = 1; but p2->flag still -1
                                            muPrevious = p2->mu;
                                        }
                                    }


                                    flagPrevious = p2->flag;
                                    muPrevious = p2->mu;
                                    disfirstTime = 1;
                                    firstTime = 1;
                                    if (p2->next) {
                                        p2 = p2->next;
                                        continue;
                                    } else {
                                        break;
                                    }
                                }

                            }
                        } else {
                            break;
                        }
                    } else {
                        // same flags, but opposite mu
                        // now flag from 1 --> -1, but mu is the same, we need to judge whether this point is a false negtive or true negtive
#ifdef VERBOSE
                        fprintf(stderr, "\t 2392 p2->flag = %d, p2->mu %f, p2->x1 x2 = %f, %f, p2->absdzs %.5e\n", p2->flag, p2->mu, p2->x1, p2->x2, p2->absdzs);
#endif
                        total_scaned = 0;
                        pos_scaned = 0;

                        anchor = p2;
                        while (p2->next) {
                            if (total_scaned < scan_percentage) {
                                total_scaned += 1;
                                pos_scaned += ( p2->flag == 1 ? 1 : 0 );
                                p2 = p2->next;
                            } else {break;}
                        }
                        if (!(p2 ->next )) {
                            break;
                        }
#ifdef VERBOSE
                        fprintf(stderr, "\t total_scaned = %d, pos_scaned = %d, percentage = %f\n", total_scaned, pos_scaned,  1.0 * pos_scaned / total_scaned * 100);
#endif
                        if (pos_scaned == total_scaned) {
                            // we can add all these points to Prov2
                            while (anchor != p2) { // leave p to add by the below code
                                addPoint(Prov2, anchor);
                                anchor = anchor->next;
                            }
                            while (p2->next && p2->next->flag == 1) {
                                addPoint(Prov2, p2);
                                p2 = p2->next;
                            }

                            if (p2 == Prov->last && p2->flag == 1) {
                                addPoint(Prov2, p2);
                            }
                            flagPrevious = 1;
                            // muPrevious = p2->prev->mu;
                            currentTrueOrFalse = 1;
                            if (p2->next) {
                                muPrevious = p2->mu;
                                p2 = p2->next;
                                continue;
                            } else {
                                break;
                            }
                        } else if (pos_scaned == 0) {
                            // we can discard them totally
                            // you can check whether from here, the segments is mainly positive or negative.
                            // scan from -1 to next point with flag +1
                            while (p2->next && p2->next->flag == -1) {
                                p2 = p2->next;
                            }


//                             if (p2 == Prov->last && p2->flag == -1) {
// #ifdef VERBOSE
//                                 fprintf(stderr, "2587 the last point has flag == -1 \n");
// #endif
//                             }

                            // now you should initiate a new Prov2:
                            if (Prov2->length > 0 &&  !Prov2->next) {
                                tempcurve = new _curve;
                                imageTracks->append(tempcurve);
                                disfirstTime = 1;
                                Prov2 = Prov2->next;
                                // tempcurve = nullptr;
                                // delete tempcurve;
                            }

                            currentTrueOrFalse = 0;
                            flagPrevious = -1;
                            muPrevious = p2->mu;
                            if (p2->next) {
                                p2 = p2->next;
                                continue;
                            } else {
                                break;
                            }
                        } else if (1.0 * pos_scaned / total_scaned > 0.8) {

                            // this point is safe
                            currentTrueOrFalse = 1;

                            addPoint(Prov2, anchor);

                            p2 = anchor; // this go first
                            muPrevious = p2->mu;
                            flagPrevious = p2->flag; // assuming we have +1 flag
                            if (p2->next) {
                                p2 = p2->next;
                                continue;
                            } else {
                                break;
                            }
                        } else { // we check from this point until it reaches another place of a true segments
                            // we should break at here

                            currentTrueOrFalse = 0;


                            // be aggresive
                            if (1.0 * Prov->posflagnum / Prov->length < 0.3 && flagPrevious == 1 && p2->flag == -1) {
                                while (p2->next && p2->next->flag == -1) {
                                    p2 = p2->next;
                                }
                                flagPrevious = -1; // now p2->next->flag = 1; but p2->flag still -1
                                muPrevious = p2->mu;
                            } else {
                                p2 = anchor;
                                muPrevious = p2->mu;
                                flagPrevious = p2->flag;
                            }

                            if (Prov2->length > 0 &&  !Prov2->next) {
                                tempcurve = new _curve;
                                imageTracks->append(tempcurve);
                                disfirstTime = 1;
                                Prov2 = Prov2->next;
                                // tempcurve = nullptr;
                                // delete tempcurve;
                            }

                            if (p2->next) {
                                p2 = p2->next;
                                continue;
                            } else {
                                break;
                            }
                        }
                    }
                }
                // begin ************************
                if (disfirstTime && Prov2->length > 1) {
                    disPrevious = *(Prov2->last) - *(Prov2->last->prev); // sqrt(dis)
                    disfirstTime = 0;
                }
                if (Prov2->length > 1 && ((*(p2) - * (Prov2->last)) >=  1e8 * disPrevious) ) {
                    if (!Prov2->next) {
                        tempcurve = new _curve;
                        imageTracks->append(tempcurve);
                        // tempcurve = nullptr;
                        // delete tempcurve;
                    }
                    Prov2 = Prov2->next;
                    flagPrevious = p2->flag;
                    muPrevious = p2->mu;
                    disfirstTime = 1;
                    p2 = p2->next;
                }
                // end ************************


                if (currentTrueOrFalse) {
                    addPoint(Prov2, p2); // add point p to the track
                    //2021.06.11 you need to update:
                    flagPrevious = 1; //p2->flag;
                    muPrevious = p2->mu;
                    previousImage = TRUE_IMAGE;
                    have_add = 1;
                }


            }

            if ( have_add ) {
                if (p2->next) {
                    p2 = p2->next;
                    continue;
                } else {
                    break;
                }
            } else {
                falseImages++;
                //first time that we enter this
                if (firstTime) {
                    previousImage = FALSE_IMAGE;
                    firstTime = 0;
                    flagPrevious = p2->flag; // 2021.06.10
                    muPrevious = p2->mu;
                }
                // already a true segment existing, we cut the segment, and start a new one
                if (previousImage == TRUE_IMAGE) { // previously we had true images, now we are encountering on false solutions
                    if (!Prov2->next) {
                        tempcurve = new _curve;
                        imageTracks->append(tempcurve);
                        disfirstTime = 1;
                        // tempcurve = nullptr;
                        // delete tempcurve;
                    }
                    Prov2 = Prov2->next;  // we start a new image track
                    previousImage = FALSE_IMAGE;
                }
                if (p2->next) {
                    p2 = p2->next;
                } else {
                    break;
                }
            }
        }
        Prov = Prov->next;

        if (Prov2->length > 0 &&  !Prov2->next) {
            tempcurve = new _curve;
            imageTracks->append(tempcurve);
            disfirstTime = 1;
            Prov2 = Prov2->next;
            // tempcurve = nullptr;
            // delete tempcurve;
        }


    }
    mixedTracks = 1;
    // finish spliting the segments.
    //  ****************************************** end version 3 ************************************************* //

    // tempcurve = NULL;
    // delete tempcurve;
    int badscaned = 0;

    // drop empty _curve
    Prov = imageTracks->first;
    while (Prov) {
        Prov2 = Prov->next;
        pos_scaned = 0;
        // count how many +1 flag points here
        for (_point *scan = Prov->first; scan; scan = scan->next) {
            pos_scaned += ( scan->flag == 1 ? 1 : 0 );
        }
        Prov->posflagnum = pos_scaned;
        total_scaned = Prov->length;

        if (!Prov->first || !Prov->last) { // add one more condition
            imageTracks->drop(Prov);
            delete Prov;
            Prov = Prov2;
        }
        else if (1.0 * pos_scaned / total_scaned < 0.05) {
            imageTracks->drop(Prov);
            delete Prov;
            Prov = Prov2;
        } else if ( total_scaned < 0.01 * pntnum ) {
            // we are dealing with a tiny segment, probably a bad one

            badscaned = 0;
            // count how many +1 flag points here
            for (_point *scan = Prov->first; scan; scan = scan->next) { // remove tiny bad segments
                badscaned += ( scan->absdzs > 1e-7 ? 1 : 0 );
            }

            if (badscaned > 0.9 * Prov->length) {
                imageTracks->drop(Prov);
                delete Prov;
                Prov = Prov2;
            } else {
                // 2021.08.06
                // this tiny segments is ok
                Prov = Prov->next;
            }
        }
        else {
            Prov = Prov->next;
        }
    }

#ifdef VERBOSE
    Prov = imageTracks->first;
    while (Prov) {

        fprintf(stderr, "---> Prov->length %d, first->absdzs %.5e\n", Prov->length, Prov->first->absdzs);
        Prov = Prov->next;
    }
#endif


//2021.06.13 // now we do pruning: remove false images from these segments
// 查找头，尾部 prung_depth 范围内有没有 absdzs ~ 0 的一对点
    int littlecnt = 0;
    int prung_depth = 0;
    double absdzsPrevious = 0;
    int prung_flag = 0;
    bool muPrevious2;
    Prov = imageTracks->first;
    while (Prov) {
// check head first
        // 如果 track 全部 flag 都 = 1，则不用处理, 这个条件还不够
        if (Prov->posflagnum == Prov->length && (Prov->first->absdzs < SOLEPS * 1e-2) && (Prov->last->absdzs < SOLEPS * 1e-2) ) {
            Prov = Prov->next;

#ifdef VERBOSE
            printf("2804 do not pruning\n");
#endif
            continue;
        }

#ifdef VERBOSE
        printf("2902 do pruning\n");
#endif
        // prung_depth = int(0.3 * Prov->length);

        // Prov->length - Prov->posflagnum 表示有多少个 flag = -1 的点

        prung_depth = fmax( int( 1.5 * (Prov->length - Prov->posflagnum) ), int(0.2 * Prov->length)  ) ; // 2021.09.01

        if (Prov->first->absdzs > SOLEPS * 1e-2) { // // 2021.09.01
            // we need to find out the separation between bad and good points
            littlecnt = 1; // counter for the number of bad images
            p2 = Prov->first;
            muPrevious2 = p2->mu > 0 ? true : false; // do not need update
            absdzsPrevious = p2->absdzs;
            prung_flag = 0;
            if (p2->next) {p2 = p2->next;}
            else break;
            // while (p2 && ( littlecnt < prung_depth )) {
            while (p2) { // 2021.09.01
                if (p2->absdzs > SOLEPS * 1e-2 && ( absdzsPrevious / p2->absdzs < 1e2 ) && !( muPrevious2 ^ ( p2->mu > 0 ? true : false ) ) ) {
                    // this means now p2 is normal
                    absdzsPrevious = p2->absdzs;
                    if (p2->next) {p2 = p2->next;}
                    else break;
                    littlecnt += 1;
                } else {
                    prung_flag = 1;
                    break;
                }
            }
#ifdef VERBOSE
            printf("2963 littlecnt = %d\n", littlecnt);
#endif
            // now p2 should be the new head
            // if (!( littlecnt == prung_depth || p2 == Prov->last ) && prung_flag) { // 2021.09.01 comment out
            if (prung_flag && littlecnt > 1 && (p2->mu)*muPrevious2 < 0) { // 2021.09.15 如果 littlecnt=1，则保留第一个点
                // 2021.09.12, we shold delete points between Prov->first and p2
                for (_point *tp = Prov->first->next; tp != p2->next; tp = tp->next) {
                    delete tp->prev;
                    // printf("delete one point 2938");
                }
                Prov->first = p2;
                p2->prev = 0;
                Prov->length -= littlecnt;
            }
        } else {
            // the head is good
        }
        if (Prov->last->absdzs > SOLEPS * 1e-2) {
            // we need to find out the separation between bad and good points
            littlecnt = 1; // counter for the number of bad images
            p2 = Prov->last;
            muPrevious2 = p2->mu > 0 ? true : false; // do not need update
            absdzsPrevious = p2->absdzs;
            prung_flag = 0;
            if (p2->prev) {p2 = p2->prev;}
            else break;
            // while (p2 && ( littlecnt < prung_depth)) {
            while (p2) { // 2021.09.01
                if (p2->absdzs > SOLEPS * 1e-2 && ( absdzsPrevious / p2->absdzs < 1e2 ) && !( muPrevious2 ^ ( p2->mu > 0 ? true : false ) )) {
                    // this means now p2 is normal
                    absdzsPrevious = p2->absdzs;
                    if (p2->prev) {p2 = p2->prev;}
                    else break;
                    littlecnt += 1;
                } else {
                    prung_flag = 1;
                    break;
                }
            }
#ifdef VERBOSE
            printf("3003 littlecnt = %d\n", littlecnt);
#endif

            // now p2 should be the new head
            // if (!( littlecnt == prung_depth || p2 == Prov->first ) && prung_flag) {
            if (prung_flag && littlecnt > 1 && (p2->mu)*muPrevious2 < 0 ) { // 2021.09.01, 2021.09.15
                // 2021.09.12, we shold delete points between p2 and Prov->last
                for (_point *tp = Prov->last->prev; tp != p2->prev; tp = tp->prev) {
                    // delete tp; // wrong
                    delete tp->next; // 2021.09.15
                }
                Prov->last = p2;
                p2->next = 0;
                Prov->length -= littlecnt;
            }
        } else {
            // the head is good
        }
        Prov = Prov->next;
    }

// show head and tail information of each track in imageTracks
#ifdef VERBOSE
    int ctk = 0, ctk2 = 0, maxprint = 10;
    for (_curve *scan = imageTracks->first; scan; scan = scan->next) {
        printf("\n\n-----> track %d, length %d\n", ctk, scan->length);
        printf("\n\n-----> head\n");
        ctk2 = 0;
        for (_point *p = scan->first; p; p = p->next) {
            if (ctk2 < maxprint) {
                printf("p->flag %d, p->mu %f, p->absdzs = %.5e\n", p->flag, p->mu , p->absdzs );

                ctk2 += 1;
            } else {
                break;
            }
        }
        printf("\n\n-----> tail\n");
        ctk2 = 0;
        for (_point *p = scan->last; p; p = p->prev) {
            if (ctk2 < maxprint) {
                printf("p->flag %d, p->mu %f, p->absdzs = %.5e\n", p->flag, p->mu , p->absdzs);
                ctk2 += 1;
            } else {
                break;
            }
        }
        ctk += 1;
    }
#endif



#ifdef VERBOSE
    saveTracks_before(imageTracks, pntnum);
#endif

#ifdef VERBOSE
    fprintf(stderr, "\nimageTracks->length %d before linking to closed continuous track\n", imageTracks->length);//for debug
#endif


    // now the following is actually no longer needed,
    // The easy case: every curve is self-closed without caustics-crossing, we are done
    if (mixedTracks == 0)  {
        if ( (imageTracks->length - (NLENS + 1)) % 2 == 0) { // double check the number of tracks should be 4 + 2*n in this case
            for (Prov = imageTracks->first; Prov; Prov = Prov->next) {
                _Prov_parity
            }
            Prov2 = NULL;
            Prov = NULL;
            delete Prov2;
            delete Prov;
            delete allSolutions;
            return (imageTracks);
        } else {
#ifdef VERBOSE
            fprintf(stderr, "image tracks should be %d, but is %d\n", NLENS + 1, imageTracks->length);
#endif
            for (Prov = imageTracks->first; Prov; Prov = Prov->next) {
                _Prov_parity
            }
            Prov2 = NULL;
            Prov = NULL;
            delete Prov2;
            delete Prov;
            delete allSolutions;
            return imageTracks;
        }
    }

    //
    // The hard case, we have caustic crossing and mixed tracks
    //
#ifdef VERBOSE
    if (mixedTracks) printf("we have mixed tracks\n");
#endif

    //
    // We first connect all self-closed, non-caustic crossing image tracks
    //
    _sols *connected;   // connected image tracks
    // _curve *tempfinal;

    connected = new _sols;

    //
    //now attach self-closed, non-caustic crossing image tracks
    //
    int selfClosed = 0;
    Prov = imageTracks->first;
    tempProv = NULL;
    while (Prov) {
        SD = *(Prov->first) - *(Prov->last);
        // if (SD < EPS * EPS && Prov->length > 3) { // closed tracks
        if (SD < EPS * EPS && Prov->length > 3) { // add on 2021.06.10
            tempProv = Prov->next;
            imageTracks->drop(Prov);
            if (Prov->first->mu < 0 && abs(Prov->first->mu) < 1e-5 ) {
                // this is a tiny image
#ifdef VERBOSE
                fprintf(stderr, "A TINY image with length %d, final->mu = %.5e, last->mu = %.5e\n", Prov->length, Prov->first->mu, Prov->last->mu);
#endif
                delete Prov; // add on 2021.09.27
            } else {
                connected->append(Prov);
                connected->last->parity = connected->last->first->mu > 0 ? 1 : -1;
                selfClosed++;
#ifdef VERBOSE
                printf("entered, pointInTrack=%d imageTracks->length %d\n", Prov->length, imageTracks->length);
#endif

            }



            //
            // such image tracks must have number of points == number of points on the source circle
            //
            // special case may be there is only one single point around a narrow cusp, / we probably
            // need to increase the number of points to solve the problem
            //
#ifdef VERBOSE
            if (Prov->length != pntnum) {
                fprintf(stderr, "the number of points is %d, should be equal to pntnum=%d!\n", Prov->length, pntnum);
            }
#endif
            // drop the already connected curve from imageTracks
            Prov = tempProv;
        } else {
            Prov = Prov->next;
        }
    }
    tempProv = NULL;



    if (imageTracks->length == 0) {
#ifdef VERBOSE
        fprintf(stderr, "return beacause imageTracks.length==0, although mixedTracks = %d\n", mixedTracks);
#endif
        Prov2 = NULL;
        Prov = NULL;
        delete imageTracks;
        delete Prov;
        delete Prov2;
        delete allSolutions;
        return connected;
    } else {
#ifdef VERBOSE
        for (_curve *scan = imageTracks->first; scan; scan = scan->next) {
            fprintf(stderr, "\t\t scan length %d, scan->first->mu = %f\n", scan->length, scan->first->mu);
        }
#endif
    }


    // 2021.06.10
    // before proceed, add a sort procedure to make the longer track go first
    imageTracks->sort();


#ifdef VERBOSE
    k = 0;
    for (_curve *scan = imageTracks->first; scan; scan = scan->next) {
        fprintf(stderr, "\t\t %d, after sort scan length %d, scan->first->mu = %f\n", k, scan->length, scan->first->mu);
        k += 1;
    }

#endif

    //
    // now we deal with the final image track that are non-closed, presumbly associated with caustic crossing
    //
    _curve *final;

    _newSegment


    //
    //we first connect identical image positions and then try identical source positions
    //
    complex zhead, ztail;//dz
    int ifcontinue = 0;
    int ifkeep_connect = 0;

    int ifcontinue_Fromjump = 0;
    int ifkeep_jumping = 0;

    int one_more_life = 0;

    // use another way, less reverse
    int head = 0, tail = 0;
    _curve *markedCurve = nullptr;
    while ( imageTracks->length >= 1) {
        _final_parity // add on 2021.06.13
#ifdef VERBOSE
        // fprintf(stderr, "can not jump, reversing final and then check if connectWithHead or connectWithTail, or jump over caustics\n");
        fprintf(stderr, "\t right below while, imageTracks->length %d, final->length %d, final->parity = %d\n", imageTracks->length, final->length, final->parity);
#endif


        // be cautious, if current "final" has length 1, you can not tell its parity
        // because, when trying to connect or jump, the direction can go both way from just one point, and the parity might be messed up,

        //_final_parity
        ifkeep_connect = 1;
        while (ifkeep_connect) {
            final = connect_head_or_tail(&ifcontinue, final, imageTracks, connected, false, &ifkeep_connect);
        }
        if (ifcontinue) continue;
        else {
#ifdef VERBOSE
            // fprintf(stderr, "can not jump, reversing final and then check if connectWithHead or connectWithTail, or jump over caustics\n");
            fprintf(stderr, "1933, 1st connect null, check head\n");
#endif
            if (imageTracks->length > 0 ){ifkeep_connect = 1;}
            else{
                // we should stop here, 2021.09.14
                _final_parity
#ifdef VERBOSE
                fprintf(stderr, "\t >>>connected->append(final) place 3230, final->length = %d\n", final->length);
#endif
                connected->append(final);

                delete allSolutions;
                delete imageTracks;
                Prov2 = NULL;
                Prov = NULL;
                delete Prov;
                delete tempProv;
                delete Prov2;                                
                return connected;
            }
            while (ifkeep_connect) {
                final = connect_head_or_tail(&ifcontinue, final, imageTracks, connected, true, &ifkeep_connect);
            }
            if (ifcontinue) continue;
            else {

                // 2021.06.11, before jump, we need to check whether the remaining _curve has only 1 point
                if (imageTracks->first->length == 1) {
                    // because we have sort imageTracks by length, this mean all remaining curve has lengh = 1
                    // we just link them to final according to distance?

                    // what if some points can connect to previous connected "final", due to we ignore them?

                    // do not handle here, handle when you jump to that point,
                    //
                }
                ifkeep_jumping = 1;
                // final = jump_over_caustics(final, imageTracks, &head, &tail, markedCurve, false, &ifcontinue_Fromjump, connected, &ifkeep_jumping);
                while (ifkeep_jumping) {
                    final = jump_over_caustics(final, imageTracks, &head, &tail, markedCurve, false, &ifcontinue_Fromjump, connected, &ifkeep_jumping);
                    if (head || tail) {
                        one_more_life = 1; // 中间只要有一次跳了，就能再获得一次机会
                        // and we now should check connectivity first, rather than jump

// #ifdef VERBOSE
//                         printf("--> 2901, final->last->x1,x2 = %f, %f\n", final->last->x1, final->last->x2);
//                         for (_curve *tt = imageTracks->first; tt; tt = tt->next) {
//                             printf("--> 2901, tt length = %d\n", tt->length);
//                             printf("--> 2901, tt head->x1,x2 = %f, %f, dzs = %.5e\n", tt->first->x1, tt->first->x2, abs( complex(tt->first->x1, tt->first->x2) - complex(final->last->x1, final->last->x2) ));
//                             printf("--> 2901, tt tail->x1,x2 = %f, %f, dzs = %.5e\n", tt->last->x1, tt->last->x2, abs( complex(tt->last->x1, tt->last->x2) - complex(final->last->x1, final->last->x2)) );
//                         }
// #endif

                        final = connect_head_or_tail(&ifcontinue, final, imageTracks, connected, false, &ifkeep_connect);
                        if (ifcontinue) continue;
                    }
                }
                if (ifcontinue_Fromjump) continue;
                else {      // we failed to jump at head, check whether we can jump from tail
                    // we are not yet self closed, we need  to jump over caustics, we seek a point with the same source position but opposite parity, and smallest magnification ratio
                    ifkeep_jumping = 1;
                    while (ifkeep_jumping) {
                        final = jump_over_caustics(final, imageTracks, &head, &tail, markedCurve, true, &ifcontinue_Fromjump, connected, &ifkeep_jumping);
                        if (head || tail) {
                            one_more_life = 1; // 中间只要有一次跳了，就能再获得一次机会

                            final = connect_head_or_tail(&ifcontinue, final, imageTracks, connected, true, &ifkeep_connect);
                            if (ifcontinue) continue;

                        }
                    }
                    if (ifcontinue_Fromjump) continue;
                    else {
                        // the most difficult part, we can neither connect or jump

                        // 2021.06.10, jumper by searching deeper, it might be that they can connect when we search a little bit deeper
                        if (final->length > JUMP_SEARCH_DEPTH ) {
                            // ifkeep_jumping = 1;

                            // do not use while here, this is not a channel that we can use unlimitedly, only try when the previous steps are fail

                            // while (ifkeep_jumping) {
                            final = jump_over_caustics_deeper(final, imageTracks, &head, &tail, markedCurve, false, &ifcontinue_Fromjump, connected, &ifkeep_jumping);
                            if (head || tail) {
                                one_more_life = 1; //
                            }
                            // }
                            if (ifcontinue_Fromjump) continue;
                            else {
                                if (final->length > JUMP_SEARCH_DEPTH ) { // it is possible that after jump_over_caustics_deeper, final is even shorter
                                    // ifkeep_jumping = 1;
                                    // while (ifkeep_jumping) {
                                    final = jump_over_caustics_deeper(final, imageTracks, &head, &tail, markedCurve, true, &ifcontinue_Fromjump, connected, &ifkeep_jumping);
                                    if (head || tail) {
                                        one_more_life = 1; // 中间只要有一次跳了，就能再获得一次机会
                                    }
                                    // }
                                    if (ifcontinue_Fromjump) continue;
                                    else {
                                        // I have noting to do, go to error correction part.
                                    }
                                }
                            }
                        }


                        //ERROR correction begin. We forgot to check the case when the track is already closed when
                        //the first and last points share the same source position.
                        dzs = final->first->zs - final->last->zs;
                        if (abs(dzs) < EPS) {
                            if (final->length >= PHI->length && imageTracks->length == 0) { // this is wrong?
                                _final_parity
                                connected->append(final);


                                delete allSolutions;
                                delete imageTracks;
                                Prov2 = NULL;
                                Prov = NULL;
                                delete Prov;
                                delete tempProv;
                                delete Prov2;
                                return connected;
                            }
                        }

                        //
                        if (one_more_life) {

// #ifdef VERBOSE
//                             fprintf(stderr, "\t\t I have one more chance \n");

//                             fprintf(stderr, "final->length = %d, final->first place = %f, %f\n", final->length, final->first->x1, final->first->x2);
//                             fprintf(stderr, "                     final->last place = %f, %f\n", final->last->x1, final->last->x2);

//                             for (_curve *ttmp = imageTracks->first; ttmp; ttmp = ttmp->next) {
//                                 fprintf(stderr, "curve remains length %d first place = %f, %f\n", ttmp->length, ttmp->first->x1, ttmp->last->x2);
//                                 fprintf(stderr, "                         last place = %f, %f\n", ttmp->last->x1, ttmp->last->x2);
//                             }
// #endif


                            one_more_life = 0;
                            continue;
                        }

                        // one possible scenario,

                        if ( ( abs(final->first->phi - final->last->phi) < EPS ||  abs( abs(final->first->phi - final->last->phi) - 2.0 * M_PI) < EPS) ) {

                            // we check whether the slope changes are smooth using four previous points
                            int i = 0;
                            double slope[2], slope_if_jumped;
                            double scatter;
#ifdef VERBOSE
                            fprintf(stderr, "number of points so far %d\n", final->length);
#endif
                            if (final->length < 3) { // final->length might even be 1

                                // just ingnore current "final"
                                _newSegment
                                _final_parity
                            } else {
                                for (_point *p = final->last;  i <= 1; i++, p = p->prev) {
                                    slope[i] = (p->prev->x2 - p->x2) / (p->prev->x1 -  p->x1);
#ifdef VERBOSE
                                    printf("slope: %f %f %d %f\n", p->x1, p->x2, i, slope[i]);
                                    printf("slope=%f\n", tan(p->thetaJ));
#endif
                                }
                                scatter = abs( log10(abs(slope[1] / slope[0])) );
                                slope_if_jumped = (final->last->x2 - final->first->x2) / (final->last->x1 - final->first->x1);
                                if (abs(log10(abs(slope[0] / slope_if_jumped))) < 2.0 * scatter) {
                                    if (final->length > 10 && imageTracks->length == 0) {
                                        connected->append(final);
                                        if (imageTracks->length > 0) {
                                            _newSegment
                                            _final_parity
                                        }

                                    }

                                }
                            }
                        }

                        // final check, what if

#ifdef VERBOSE
                        fprintf(stderr, "we can do nothing with this segment, something went wrong, append final anyway\n");
                        fprintf(stderr, "\t >>>2028 connected->append(final) place 5, final->length = %d\n", final->length);
                        fprintf(stderr, "2029 due to may be small mass ratio (lack solution), some image tracks is not closed, just drop it, the area is small anyway, \n \t\t imageTracks.length = %d, connected.length = %d, final->length = %d\n", imageTracks->length, connected->length, final->length);
#endif
                        _final_parity
                        connected->append(final);

                        if (imageTracks->length > 0) {
                            connected->append(final);
                            _newSegment
                            _final_parity
                        }
                        continue;
                    }
                    //ERROR correction end
                }
            }
        }
    }


//check whether the first and last points are connected, if not, join the first and last points
    SD = *(final->first) - *(final->last);
#ifdef VERBOSE
    fprintf(stderr, "2047 imageTracks->length %d, connected->length %d\n", imageTracks->length, connected->length );
    printf("final->length = %d, joining the first and last point %f %f %f %f %f %f\n", final->length, final->first->x1, final->first->x2, final->first->phi,
           final->last->x1, final->last->x2, final->last->phi);
    fprintf(stderr, "SD = %.3e\n", SD);
#endif


    if (connected->last != final && final->length > 1) {
        _final_parity
#ifdef VERBOSE
        fprintf(stderr, "\t >>>connected->append(final) place 8, final->length = %d\n", final->length);
#endif
        connected->append(final); // attach the final curve
    }
    delete allSolutions;
    delete imageTracks;
    Prov2 = NULL;
    Prov = NULL;
    delete Prov;
    delete tempProv;
    delete Prov2;
    return connected;
}



#undef TRUE_IMAGE
#undef FALSE_IMAGE



_curve *jump_over_caustics_deeper(_curve * final, _sols * imageTracks, int *head, int *tail, _curve * markedCurve, bool checkfinalhead, int *ifcontinue_Fromjump, _sols * connected, int *ifkeep_jumping)
{
    *tail = 0;
    *head = 0;
    int cnt1, cnt2, mark_reduced_final, mark_reduced_sub;
    _point *pntscan, *pntscan_inner1, *pntscan_inner2, *markpnt_final, *markpnt_tobetrunc;
    int JUMP_SEARCH_DEPTH_FINAL = 5;
    int JUMP_SEARCH_DEPTH_OTHER = 5;

    double ratio, ratioMin = 1.0e10, muPrevious, disimg, mu;
    complex zs, dzs;
    complex zsPrevious, imgPrevious;
    if (checkfinalhead == false) {
        // 从 final 的 末尾开始找
        pntscan = final->last;
    } else {
        pntscan = final->first;
    }
    JUMP_SEARCH_DEPTH_FINAL = fmin( fmax(int(0.1 * final->length ), JUMP_SEARCH_DEPTH), 20 * JUMP_SEARCH_DEPTH );
    for (cnt1 = 0; cnt1 < JUMP_SEARCH_DEPTH_FINAL; cnt1++) {
        zsPrevious = pntscan->zs;
        muPrevious = pntscan->mu;
        imgPrevious = complex(pntscan->x1, pntscan->x2);

        for (_curve *Prov2 = imageTracks->first; Prov2; Prov2 = Prov2->next) {
            JUMP_SEARCH_DEPTH_OTHER =  fmin( fmax(int(0.1 * Prov2->length) ,  JUMP_SEARCH_DEPTH), 20 * JUMP_SEARCH_DEPTH );
            pntscan_inner1 = Prov2->first;
            pntscan_inner2 = Prov2->last;
            for (cnt2 = 0; cnt2 < JUMP_SEARCH_DEPTH_OTHER && (pntscan_inner1 || pntscan_inner2); cnt2++) {
                // we first test whether the last end point connects with the head of some curve
                zs = pntscan_inner1->zs;
                mu = pntscan_inner1->mu;

                dzs = zs - zsPrevious;

                disimg = abs(complex(pntscan_inner1->x1, pntscan_inner1->x2) - imgPrevious);

                // this condition is not safe enough, sometime it will link from just a normal place, so do not use this deeper search multiple times, just once (works after all previous trail fails), i.e. do not use while
                if (abs(dzs) < EPS && muPrevious * mu < 0.0 && abs(pntscan->phi - pntscan_inner1->phi) < EPS  && ( ( abs(muPrevious / mu) + abs(muPrevious / mu) ) < 2.5 ) ) { // connected with the head with opposite parity
                    // phi is exactly the same

                    ratio = abs( log10( abs(muPrevious / mu)) ) * disimg;
                    //  printf("head - ratio, ratioMin=%f %f\n", ratio, ratioMin);
                    if (ratio < ratioMin) {
                        ratioMin = ratio;
                        markedCurve = Prov2;
                        *head = 1;   // head has the minimum
                        *tail = 0;
                        markpnt_tobetrunc = pntscan_inner1;
                        mark_reduced_sub = cnt2;

                        markpnt_final = pntscan;
                        mark_reduced_final = cnt1; // mark how many points will be dropped from the "final"
                    }
                }

                // we now test whether the last end point connects with the tail of the next curve
                zs = pntscan_inner2->zs;
                mu = pntscan_inner2->mu;
                dzs = zs - zsPrevious;

                disimg = abs(complex(pntscan_inner2->x1, pntscan_inner2->x2) - imgPrevious);
                if (abs(dzs) < EPS && muPrevious * mu < 0.0 && abs( pntscan->phi - pntscan_inner2->phi) < EPS && ( ( abs(muPrevious / mu) + abs(muPrevious / mu) ) < 2.5 ) ) { // connected with the tail // remove on 2020.07.17
                    ratio = abs( log10( abs(muPrevious / mu)) ) * disimg;
                    if (ratio < ratioMin) {
                        ratioMin = ratio;
                        markedCurve = Prov2;
                        *tail = 1;   // tail has the minimum
                        *head = 0;
                        markpnt_tobetrunc = pntscan_inner2;
                        mark_reduced_sub = cnt2;

                        markpnt_final = pntscan;
                        mark_reduced_final = cnt1; // mark how many points will be dropped from the "final"

                    }
                }
                pntscan_inner1 = pntscan_inner1->next;
                pntscan_inner2 = pntscan_inner2->prev;
            }
        }


        if (checkfinalhead == false) {
            pntscan = pntscan->prev;
        } else {
            pntscan = pntscan->next;
        }
    }

#ifdef VERBOSE
    if (*head || *tail) {
        fprintf(stderr, "in jump_over_caustics_deeper, find a pair to jump from the tail of final\n");
        fprintf(stderr, "\t\t head = %d, tail = %d, checkfinalhead = %d \n", *head, *tail, checkfinalhead);

        fprintf(stderr, "\t\t mark_reduced_final %d, mark_reduced_sub %d\n", mark_reduced_final, mark_reduced_sub);
        fprintf(stderr, "\t\t at: final @ phi = %f, mu = %.5e, xy = %f, %f\n", markpnt_final->phi, markpnt_final->mu, markpnt_final->x1, markpnt_final->x2);
        fprintf(stderr, "\t\t at:   sub @ phi = %f, mu = %.5e, xy = %f, %f\n", markpnt_tobetrunc->phi, markpnt_tobetrunc->mu, markpnt_tobetrunc->x1, markpnt_tobetrunc->x2);
    }
#endif

    if (checkfinalhead == true) {
        // 用 final 的 head 去连


        if (*head == 1) {              // jumped over caustics to the head of a curve
#ifdef VERBOSE
            printf("2167 jumped to head of a curve: %f %f %f %f, phi=%.15f\n", markedCurve->first->x1, markedCurve->first->x2, muPrevious, markedCurve->first->mu, markedCurve->first->phi * 180 / M_PI);
            printf("jumped to head and going to tail at: %f %f\n", markedCurve->last->x1, markedCurve->last->x2);
            fprintf(stderr, "final->length before drop %d, connected->length = %d\n", final->length, connected->length);
            printf("jump_dis %.2e, final-head-tail-dis %.2e\n", abs(complex(markpnt_final->x1, markpnt_final->x2) -   complex(markpnt_tobetrunc->x1, markpnt_tobetrunc->x2)), abs(complex(final->first->x1, final->first->x2) -   complex(final->last->x1, final->last->x2)  ));
#endif

            // 2021.08.06 不仅要 head == 1, 而且要 连完之后，连接处的长度，不能比 原来 final 的首尾距离  长太多，如 10倍
            if ( abs(complex(markpnt_final->x1, markpnt_final->x2) -   complex(markpnt_tobetrunc->x1, markpnt_tobetrunc->x2)) <
                    10 *   abs(complex(final->first->x1, final->first->x2) -   complex(final->last->x1, final->last->x2)  ) ) {

                imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list
                // final head 连 head
                // _final_parity
                // 首先更新 parity and final:
                // final->parity = markpnt_final->mu > 0 ? 1 : -1; // wrong
                // or         -(markpnt_tobetrunc->mu > 0 ? 1:-1);
                final->first->next = markpnt_final->next;
                final->first = markpnt_final;
                markpnt_final->prev = 0;
                final->length -= mark_reduced_final;

                // 更新 marked:
                // markedCurve->first->next = markpnt_tobetrunc->next;

                // 2021.09.12 // 2021.09.15
                for (_point *tp = markedCurve->first->next; tp != markpnt_tobetrunc->next; tp = tp->next) {
                    delete tp->prev;
                }

                markedCurve->first = markpnt_tobetrunc;
                markpnt_tobetrunc->prev = 0;
                markedCurve->length -= mark_reduced_sub;


                markedCurve->reverse();
                final->first->closepairparity = (final->first->mu) > 0 ? 1 : -1;
                markedCurve->last->closepairparity = (markedCurve->last->mu) > 0 ? 1 : -1;
                markedCurve->last->next = final->first;
                final->first->prev = markedCurve->last;
                markedCurve->last = final->last;
                markedCurve->length += final->length;
                markedCurve->parity = final->parity;
                final = markedCurve;
            } else {

                *head = 0;
            }


        } else if (*tail == 1) { // jumped over caustics to the tail of a curve
#ifdef VERBOSE
            printf("2323 jumped to tail of a curve: %f %f %f %f, phi=%.15f\n", markedCurve->last->x1, markedCurve->last->x2, muPrevious, markedCurve->last->mu, markedCurve->last->phi * 180 / M_PI);
            printf("jumped to tail and going to head: %f %f\n", markedCurve->first->x1, markedCurve->first->x2);
            fprintf(stderr, "1984 imageTracks->length before drop %d, final->length %d, final->parity %d\n", imageTracks->length, final->length, final->parity);
            printf("jump_dis %.2e, final-head-tail-dis %.2e\n", abs(complex(markpnt_final->x1, markpnt_final->x2) -   complex(markpnt_tobetrunc->x1, markpnt_tobetrunc->x2)), abs(complex(final->first->x1, final->first->x2) -   complex(final->last->x1, final->last->x2)  ));
#endif

            if ( abs(complex(markpnt_final->x1, markpnt_final->x2) -   complex(markpnt_tobetrunc->x1, markpnt_tobetrunc->x2)) <
                    10 *   abs(complex(final->first->x1, final->first->x2) -   complex(final->last->x1, final->last->x2)  ) ) {

                final->first->closepairparity = final->first->mu > 0 ? 1 : -1;
                markedCurve->last->closepairparity = (markedCurve->last->mu) > 0 ? 1 : -1;

                imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list

                // final->parity = markpnt_final->mu > 0 ? 1 : -1; // wrong
                // or         -(markpnt_tobetrunc->mu > 0 ? 1:-1);
                final->first->next = markpnt_final->next;
                final->first = markpnt_final;
                markpnt_final->prev = 0;
                final->length -= mark_reduced_final;


                // markedCurve->last->prev = markpnt_tobetrunc->prev;
                // 2021.09.12, we shold delete points between Prov->last and p2
                for (_point *tp = markedCurve->last->prev; tp != markpnt_tobetrunc->prev; tp = tp->prev) {
                    delete tp->next;
                }
                markedCurve->last = markpnt_tobetrunc;
                markpnt_tobetrunc->next = 0;
                markedCurve->length -= mark_reduced_sub;


                markedCurve->last->next = final->first;
                final->first->prev = markedCurve->last;
                markedCurve->last = final->last;
                markedCurve->length += final->length;
                markedCurve->parity = final->parity;
                final = markedCurve;

            } else {
                *tail = 0;
            }

        }

    }

    if (checkfinalhead == false) {
        // here we jump over caustics from the tail of "final"

        if (*head == 1) {               // jumped over caustics to the head of a curve
#ifdef VERBOSE
            printf("2394 jumped to head of a curve: %f %f %f %f, final->last->mu %f, phi=%.15f\n", markedCurve->first->x1, markedCurve->first->x2, muPrevious, markedCurve->first->mu, final->last->mu, markedCurve->first->phi * 180 / M_PI);
            printf("jumped to head and going to tail at: %f %f\n", markedCurve->last->x1, markedCurve->last->x2);
            printf("jump_dis %.2e, final-head-tail-dis %.2e\n", abs(complex(markpnt_final->x1, markpnt_final->x2) -   complex(markpnt_tobetrunc->x1, markpnt_tobetrunc->x2)), abs(complex(final->first->x1, final->first->x2) -   complex(final->last->x1, final->last->x2)  ));
#endif

            if ( abs(complex(markpnt_final->x1, markpnt_final->x2) -   complex(markpnt_tobetrunc->x1, markpnt_tobetrunc->x2)) <
                    10 *   abs(complex(final->first->x1, final->first->x2) -   complex(final->last->x1, final->last->x2)  ) ) {


                imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list
                // drop point from tail
                final->last->prev = markpnt_final->prev;
                final->last = markpnt_final;
                markpnt_final->next = 0;
                final->length -= mark_reduced_final;

                // 更新 marked:
                // markedCurve->first->next = markpnt_tobetrunc->next;
                for (_point *tp = markedCurve->first->next; tp != markpnt_tobetrunc->next; tp = tp->next) {
                    delete tp->prev;
                }
                markedCurve->first = markpnt_tobetrunc;
                markpnt_tobetrunc->prev = 0;
                markedCurve->length -= mark_reduced_sub;

                final->last->closepairparity = (final->last->mu) > 0 ? 1 : -1;
                markedCurve->first->closepairparity = (markedCurve->first->mu) > 0 ? 1 : -1;

                // 2020.07.17, attach two segments directly instead of for loop above, to speed up the code.
                // _final_parity
                final->last->next = markedCurve->first;
                markedCurve->first->prev = final->last;
                final->last = markedCurve->last;
                final->length += markedCurve->length;

            } else {

                *head = 0;
            }


        } else if (*tail == 1) { // jumped over caustics to the tail of a curve
#ifdef VERBOSE
            printf("2441 jumped to tail of a curve: %f %f %f %f, final->last->mu %f, phi=%.15f\n", markedCurve->last->x1, markedCurve->last->x2, muPrevious, markedCurve->last->mu, final->last->mu, markedCurve->last->phi * 180 / M_PI);
            printf("jumped to tail and going to head: %f %f\n", markedCurve->first->x1, markedCurve->first->x2);
            printf("jump_dis %.2e, final-head-tail-dis %.2e\n", abs(complex(markpnt_final->x1, markpnt_final->x2) -   complex(markpnt_tobetrunc->x1, markpnt_tobetrunc->x2)), abs(complex(final->first->x1, final->first->x2) -   complex(final->last->x1, final->last->x2)  ));

#endif
            if ( abs(complex(markpnt_final->x1, markpnt_final->x2) -   complex(markpnt_tobetrunc->x1, markpnt_tobetrunc->x2)) <
                    10 * abs(complex(final->first->x1, final->first->x2) -   complex(final->last->x1, final->last->x2)  ) ) {


                imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list

                // jump from tail to tail;
                final->last->prev = markpnt_final->prev;
                final->last = markpnt_final;
                markpnt_final->next = 0;
                final->length -= mark_reduced_final;


                // drop points from tail of the markedCurve
                // markedCurve->last->prev = markpnt_tobetrunc->prev;
                for (_point *tp = markedCurve->last->prev; tp != markpnt_tobetrunc->prev; tp = tp->prev) {
                    delete tp->next;
                }
                markedCurve->last = markpnt_tobetrunc;
                markpnt_tobetrunc->next = 0;
                markedCurve->length -= mark_reduced_sub;



                final->last->closepairparity = (final->last->mu) > 0 ? 1 : -1;
                // _final_parity
                markedCurve->reverse();
                markedCurve->first->closepairparity = (markedCurve->first->mu) > 0 ? 1 : -1;
                final->last->next = markedCurve->first;
                markedCurve->first->prev = final->last;
                final->last = markedCurve->last;
                final->length += markedCurve->length;

            } else {

                *tail = 0;
            }
        }
    }


    if (*head || *tail) {

        *ifkeep_jumping = 0; // we have jumped once, let's check whether we can keep jumping
        return final;
    } else  {
        *ifkeep_jumping = 0;
        double SD;
        //    After jumping, we now check whether the tail and head of the current curve actually connect.
        SD = *(final->first) - *(final->last);

        if (SD < EPS * EPS) { // we are already connected

#ifdef VERBOSE
            fprintf(stderr, "jumped, head and tail connected!\n");
#endif
            _final_parity
#ifdef VERBOSE
            fprintf(stderr, "\t >>>connected->append(final) place 6, final->length = %d, imageTracks->length = %d\n", final->length, imageTracks->length);
#endif
            connected->append(final);

            //all segments connected, we are done
            if (imageTracks->length > 0) {
                _newSegment
                _final_parity
                *ifcontinue_Fromjump = 1;
            }
        }

    }
// this "final" is tough, needs to be further handle, so *ifcontinue_Fromjump should be 0, to let the code keep run, rather than "continue" the while loop
    *ifcontinue_Fromjump = 0;
    return final;
}


// head, tail, markedCurve = jump_over_caustics(final, imageTracks)
_curve *jump_over_caustics(_curve * final, _sols * imageTracks, int *head, int *tail, _curve * markedCurve, bool checkfinalhead, int *ifcontinue_Fromjump, _sols * connected, int *ifkeep_jumping)
{
    // int jumpOverCaustics = 0;
    *tail = 0;
    *head = 0;
    double ratio, ratioMin = 1.0e10, muPrevious, disimg, mu;
    complex zs, dzs;
    complex zsPrevious, imgPrevious;
    if (checkfinalhead == false) {
        zsPrevious = final->last->zs;
        muPrevious = final->last->mu;
        imgPrevious = complex(final->last->x1, final->last->x2);
#ifdef VERBOSE
        printf("2542, use final tail , I am in jump_over_caustics(), need to jump over caustics at: %f %f final->last->mu %f, phi=%.15f, zs: %f, %f, imageTracks-length = %d\n", final->last->x1, final->last->x2,  final->last->mu, final->last->phi * 180 / M_PI, zsPrevious.re, zsPrevious.im , imageTracks->length);
#endif
    } else {
        // check the head of final instead of tail of final. // 2020.07.17
        zsPrevious = final->first->zs;
        muPrevious = final->first->mu;
        imgPrevious = complex(final->first->x1, final->first->x2);
#ifdef VERBOSE
        printf("2550, use final head , I am in jump_over_caustics(), need to jump over caustics at first point of final: %f %f final->first->mu %f, phi=%.15f, zs: %f, %f, imageTracks-length = %d\n", final->first->x1, final->first->x2,  muPrevious, final->first->phi * 180 / M_PI, zsPrevious.re, zsPrevious.im, imageTracks->length);
#endif
    }

    for (_curve *Prov2 = imageTracks->first; Prov2; Prov2 = Prov2->next) {
        // we first test whether the last end point connects with the head of some curve
        zs = Prov2->first->zs;
        mu = Prov2->first->mu;
        dzs = zs - zsPrevious;

        disimg = abs(complex(Prov2->first->x1, Prov2->first->x2) - imgPrevious);
        // 这里 一刀切 abs(dzs) < EPS，不是最佳的策略，但不加也不行？
        // && ( ( abs(muPrevious / mu) + abs(muPrevious / mu) )<2.5 ) -- add on 2020.06.10

        // 2021.06.11, actually, if you can jump, you should have exactlly the same source positions.

        // if (abs(dzs) < EPS && muPrevious * mu < 0.0 && ( ( abs(muPrevious / mu) + abs(muPrevious / mu) ) < 2.5 ) ) { // connected with the head with opposite parity
        if (abs(dzs) < EPS && muPrevious * mu < 0.0 && ( ( abs(muPrevious / mu) + abs(mu/muPrevious) ) < 2.5 ) ) { // 2021.10.07
            // if (abs(dzs) < EPS) { // connected with the head with opposite parity // remove on 2020.07.17
            ratio = abs( log10( abs(muPrevious / mu)) ) * disimg;
            //  printf("head - ratio, ratioMin=%f %f\n", ratio, ratioMin);
            if (ratio < ratioMin) {
                ratioMin = ratio;
                markedCurve = Prov2;
                *head = 1;   // head has the minimum
                *tail = 0;
            }
        }

        // we now test whether the last end point connects with the tail of the next curve
        zs = Prov2->last->zs;
        mu = Prov2->last->mu;
        dzs = zs - zsPrevious;

        disimg = abs(complex(Prov2->last->x1, Prov2->last->x2) - imgPrevious);
        if (abs(dzs) < EPS && muPrevious * mu < 0.0 && ( ( abs(muPrevious / mu) + abs(mu/muPrevious) ) < 2.5 ) ) {
            ratio = abs( log10( abs(muPrevious / mu)) ) * disimg;
            //  printf("tail - ratio, ratioMin=%f %f\n", ratio, ratioMin);
            if (ratio < ratioMin) {
                ratioMin = ratio;
                markedCurve = Prov2;
                // fprintf(stderr, "head???\n");
                *tail = 1;   // tail has the minimum
                *head = 0;
            }
        }
    }


    if (checkfinalhead == true) {
        // 用 final 的 tail 去连
        if (*head == 1) {               // jumped over caustics to the head of a curve
#ifdef VERBOSE
            printf("2167 jumped to head of a curve: %f %f %f %f, phi=%.15f\n", markedCurve->first->x1, markedCurve->first->x2, muPrevious, markedCurve->first->mu, markedCurve->first->phi * 180 / M_PI);
            printf("jumped to head and going to tail at: %f %f\n", markedCurve->last->x1, markedCurve->last->x2);
            fprintf(stderr, "final->length before drop %d, connected->length = %d\n", final->length, connected->length);
#endif
            imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list

            _final_parity
            markedCurve->reverse();
            final->first->closepairparity = (final->first->mu) > 0 ? 1 : -1;
            markedCurve->last->closepairparity = (markedCurve->last->mu) > 0 ? 1 : -1;
            markedCurve->last->next = final->first;
            final->first->prev = markedCurve->last;
            markedCurve->last = final->last;
            markedCurve->length += final->length;
            markedCurve->parity = final->parity;
            final = markedCurve;
#ifdef VERBOSE
            fprintf(stderr, "3531 imageTracks->length after drop %d, final->length %d, final->parity %d\n", imageTracks->length, final->length, final->parity);
#endif
        } else if (*tail == 1) { // jumped over caustics to the tail of a curve
#ifdef VERBOSE
            printf("2196 jumped to tail of a curve: %f %f %f %f, phi=%.15f\n", markedCurve->last->x1, markedCurve->last->x2, muPrevious, markedCurve->last->mu, markedCurve->last->phi * 180 / M_PI);
            printf("jumped to tail and going to head: %f %f\n", markedCurve->first->x1, markedCurve->first->x2);
#endif
            final->first->closepairparity = final->first->mu > 0 ? 1 : -1;
            markedCurve->last->closepairparity = (markedCurve->last->mu) > 0 ? 1 : -1;
            imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list
            markedCurve->last->next = final->first;
            final->first->prev = markedCurve->last;
            markedCurve->last = final->last;
            markedCurve->length += final->length;
            markedCurve->parity = final->parity;
            final = markedCurve;

#ifdef VERBOSE
            fprintf(stderr, "3552 imageTracks->length after drop %d, final->length %d, final->parity %d\n", imageTracks->length, final->length, final->parity);
#endif
        }
    }

    if (checkfinalhead == false) {
        // here we jump over caustics from the tail of "final"
        if (*head == 1) {               // jumped over caustics to the head of a curve
#ifdef VERBOSE
            printf("2226 jumped to head of a curve: %f %f muPrev %f markedCurve->mu %f, phi=%.15f\n", markedCurve->first->x1, markedCurve->first->x2, muPrevious, markedCurve->first->mu, markedCurve->first->phi * 180 / M_PI);
            printf("\t\t abs(dzs) = %.5e \n", abs(  final->last->zs - markedCurve->first->zs  ) );
            printf("jumped to head and going to tail at: %f %f\n", markedCurve->last->x1, markedCurve->last->x2);

#endif

            final->last->closepairparity = (final->last->mu) > 0 ? 1 : -1;
            markedCurve->first->closepairparity = (markedCurve->first->mu) > 0 ? 1 : -1;
#ifdef VERBOSE
            fprintf(stderr, "markedCurve->length %d\n", markedCurve->length);
            fprintf(stderr, "final->length before drop %d, connected->length = %d\n", final->length, connected->length);
#endif
            imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list
            // 2020.07.17, attach two segments directly instead of for loop above, to speed up the code.
            // _final_parity
            final->last->next = markedCurve->first;
            markedCurve->first->prev = final->last;
            final->last = markedCurve->last;
            final->length += markedCurve->length;
#ifdef VERBOSE
            fprintf(stderr, "3585 imageTracks->length after drop %d, final->length %d, final->parity %d\n", imageTracks->length, final->length, final->parity);
#endif
        } else if (*tail == 1) { // jumped over caustics to the tail of a curve
#ifdef VERBOSE
            printf("2248 jumped to tail of a curve: %f %f muPrevious: %f markedCurve->last->mu: %f, final->last->mu %f, phi=%.15f\n", markedCurve->last->x1, markedCurve->last->x2, muPrevious, markedCurve->last->mu, final->last->mu, markedCurve->last->phi * 180 / M_PI);
            printf("jumped to tail and going to head: %f %f\n", markedCurve->first->x1, markedCurve->first->x2);
#endif

            final->last->closepairparity = (final->last->mu) > 0 ? 1 : -1;

            imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list

            // _final_parity
            markedCurve->reverse();
            markedCurve->first->closepairparity = (markedCurve->first->mu) > 0 ? 1 : -1;
            final->last->next = markedCurve->first;
            markedCurve->first->prev = final->last;
            final->last = markedCurve->last;
            final->length += markedCurve->length;
#ifdef VERBOSE
            fprintf(stderr, "3605 imageTracks->length after drop %d, final->length %d, final->parity %d\n", imageTracks->length, final->length, final->parity);
#endif
        }
    }


    if (*head || *tail) {
#ifdef VERBOSE
        fprintf(stderr, "2269 jumped, and I am tring to keep jumping!\n");
#endif
        *ifkeep_jumping = 1; // we have jumped once, let's check whether we can keep jumping
        return final;
    } else  {
        *ifkeep_jumping = 0;
        double SD;
        //    After jumping, we now check whether the tail and head of the current curve actually connect.
        SD = *(final->first) - *(final->last);

        if (SD < EPS * EPS) { // we are already connected

#ifdef VERBOSE
            fprintf(stderr, "jumped, head and tail connected!\n");
#endif
            _final_parity
#ifdef VERBOSE
            fprintf(stderr, "\t >>>connected->append(final) place 6, final->length = %d, imageTracks->length = %d\n", final->length, imageTracks->length);
#endif
            connected->append(final);

            //all segments connected, we are done
            if (imageTracks->length > 0) {
                _newSegment
                _final_parity
                *ifcontinue_Fromjump = 1;
            }
        }

    }
    // this "final" is tough, needs to be further handle, so *ifcontinue_Fromjump should be 0, to let the code keep run, rather than "continue" the while loop
    *ifcontinue_Fromjump = 0;
    return final;
}

_curve *connect_head_or_tail(int *ifcontinue, _curve * final, _sols * imageTracks, _sols * connected, bool checkfinalhead, int *ifkeep_connect)
{
    *ifcontinue = 0;
    int connectWithHead = 0, connectWithTail = 0;
    complex zsPrevious, z, zhead, ztail, dz;
    // double muPrevious, mindis, tempdis, SD;
    double SD = 0;
    _curve *Prov2, *Prov;
    if (checkfinalhead == false) {

        // check whether final->last can connect to other segments.
        z = complex(final->last->x1, final->last->x2); // last point in the last segment
        Prov2 = imageTracks->first;
        // see whether its image position overlaps with the head or tail of another segment, if yes, we connect it
        while (Prov2) {
            // fprintf(stderr, "can you here 5, imageTracks->length = %d\n", imageTracks->length);
            // see whether the image position overlaps with the head position
            zhead = complex(Prov2->first->x1, Prov2->first->x2);
            z = complex(final->last->x1, final->last->x2); // because final might be updated in the following, so z needs to be updated
            dz = zhead - z;
            if (abs(dz) < EPS ) {
                connectWithHead = 1;

                // // connect the tail of final with the head of Prov 2:

                imageTracks->drop(Prov2); // we are already connected, drop this curve from the linked list
                Prov = Prov2->next;
                Prov2->parity = 0; // remove the parity information of Prov2, since now it is part of "final"
                final->last->next = Prov2->first;
                Prov2->first->prev = final->last;
                final->last = Prov2->last;
                final->length += Prov2->length;
                Prov2 = Prov;
            }

            if (Prov2) {
                // see whether the image position overlaps with the tail image position
                ztail = complex(Prov2->last->x1, Prov2->last->x2);
                z = complex(final->last->x1, final->last->x2); // last point in the last segment

                dz = ztail - z;
                if (abs(dz) < EPS ) {
                    connectWithTail = 1;

                    // connect with the tail of Prov

                    Prov = Prov2->next;
                    imageTracks->drop(Prov2); // we are already connected, drop this curve from the linked list

                    Prov2->reverse();

                    final->last->next = Prov2->first;
                    Prov2->first->prev = final->last;
                    final->last = Prov2->last;
                    final->length += Prov2->length;

                    Prov2 = Prov;
                } else {
                    Prov2 = Prov2->next;
                }
            }

        }

    } else {
        // check head of final instead of check tail of final
        // z = complex(final->first->x1, final->first->x2); // first point in the last segment
        Prov2 = imageTracks->first;
        // see whether its image position overlaps with the head or tail of another segment, if yes, we connect it
        while (Prov2) {
            // see whether the image position overlaps with the head position
            zhead = complex(Prov2->first->x1, Prov2->first->x2);
            z = complex(final->first->x1, final->first->x2); // last point in the last segment
            dz = zhead - z;
            if (abs(dz) < EPS ) {
                connectWithHead = 1;

                imageTracks->drop(Prov2); // we are already connected, drop this curve from the linked list
                Prov = Prov2->next;
                Prov2->reverse();
                _final_parity
                Prov2->parity = final->parity;

                Prov2->last->next = final->first;
                final->first->prev = Prov2->last;
                Prov2->last = final->last;
                Prov2->length += final->length;

                final = Prov2;
                Prov2 = Prov;
            }

            if (Prov2) {
                // see whether the image position overlaps with the tail image position
                ztail = complex(Prov2->last->x1, Prov2->last->x2);
                z = complex(final->first->x1, final->first->x2); // last point in the last segment
                dz = ztail - z;
                if (abs(dz) < EPS ) {
                    connectWithTail = 1;
                    Prov = Prov2->next;

                    imageTracks->drop(Prov2); // we are already connected, drop this curve from the linked list

                    _final_parity
                    Prov2->parity = final->parity;

                    Prov2->last->next = final->first;
                    final->first->prev = Prov2->last;
                    Prov2->last = final->last;
                    Prov2->length += final->length;
                    final = Prov2;
                    Prov2 = Prov;

                } else {
                    Prov2 = Prov2->next;
                }
            }
        }

    }

#ifdef VERBOSE
    fprintf(stderr, "2444 connectWithHead = %d, connectWithTail = %d, imageTracks.length=%d,\n first and last image position difference:%.2e, final->length = %d\n", connectWithHead, connectWithTail , imageTracks->length, sqrt(*(final->first) - * (final->last)), final->length);
#endif

    // 2021.09.14, after previous connecting process, the imageTrack length will reduce, might be reduced to 0
    if (imageTracks->length == 0) {
        *ifkeep_connect = 0;
        *ifcontinue = 1; //2021.09.15
        return final;
    }

    if (connectWithHead || connectWithTail) {
        *ifkeep_connect = 1;
        return final;
    } else {
        *ifkeep_connect = 0;
        if (final->length > 2 ) {
            SD = *(final->first) - *(final->last);
            if (SD < EPS * EPS) {
                connected->append(final);
                if (imageTracks->length != 0) {
                    _newSegment
                    _final_parity
#ifdef VERBOSE
                    fprintf(stderr, "\n\nIn connectWithHead or connectWithTail, add a new closed curve, 'final' is updated\n\n");
#endif
                }
                *ifcontinue = 1; // this means final is updated (a new one), so should continue in the outer loop
            } else {
            }

        } else {
#ifdef VERBOSE
            fprintf(stderr, "\t 3340 in connect_head_or_tail, imageTracks->length = %d, no head and tail can connect \n", imageTracks->length);
#endif
        }
        *ifcontinue = 0;
        return final;
    }
}

_point *pnt_reach_head_or_tail(_point * p, _curve * Prov, _sols * imageTracks, int *reach)
{
    *reach = 0;
    _point *head_or_tail_pnt;
    double scale = 1e2, dzcomp, dz;
    complex z, znext, zhead, ztail;
    _curve *Prov2;
    z = complex(p->x1, p->x2);
    if (p->next) {
        znext = complex(p->next->x1, p->next->x2);
        dzcomp = abs(znext - z);
    }
    else {
        dzcomp = 0.0;
    }
    Prov2 = imageTracks->first;
    while (Prov2) {
        if (Prov2 != Prov) {
            zhead = complex(Prov2->first->x1, Prov2->first->x2);
            dz = abs(zhead - z) * scale;
            if (dz < dzcomp) {
                head_or_tail_pnt = Prov2->first;
                *reach = 1;
                return head_or_tail_pnt;
            }
            ztail = complex(Prov2->last->x1, Prov2->last->x2);
            dz = abs(ztail - z) * scale;
            if (dz < dzcomp) {
                head_or_tail_pnt = Prov2->last;
                *reach = 1;
                return head_or_tail_pnt;
            }
        }
        Prov2 = Prov2->next;
    }
    head_or_tail_pnt = NULL;
    return head_or_tail_pnt;
}

// // given the image position, find the source positon and Mu
// void lensEquation(double mlens[], complex zlens[], complex z, complex * zs, double * mu)
// {
//     complex dzsdz;

//     *zs = ( z - mlens[0] / conj(z - zlens[0]) - mlens[1] / conj(z - zlens[1]) - mlens[2] / conj(z - zlens[2]) );

//     dzsdz = mlens[0] / (z - zlens[0]) / (z - zlens[0]) + mlens[1] / (z - zlens[1]) / (z - zlens[1]) + mlens[2] / (z - zlens[2]) / (z - zlens[2]);
//     *mu = 1.0 / ( 1.0 - pow(abs(dzsdz), 2.0) );
// }


//
// check whether an image z is a true solution of the lens equation, and return the two eigenvalues of the Jacobian matrix and the angle of the ellipse
//
// int trueSolution(double mlens[], complex zlens[], double xs, double ys, complex z, double * mu, double * lambda1, double * lambda2, double * thetaJ, int nlens, complex *J1, complex *J2, complex * dJ, complex *J3)
int TripleLensing::trueSolution(double xs, double ys, complex z, double * mu)
{

    flag = 0; //

    Jxx = 1.0;
    Jyy = 0.0;
    Jxy = 0.0;
    x = z.re;
    y = z.im;

    switch (NLENS) {
    case 2:
        // binary lens
        rep = xs - x;// real part
        imp = ys - y;
        x_xj = x - zlens[0].re;
        y_yj = y - zlens[0].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        dzs = complex(rep, imp);

        break;

    case 3:
        // triple lens

        rep = xs - x;// real part
        imp = ys - y;
        x_xj = x - zlens[0].re;
        y_yj = y - zlens[0].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        x_xj = x - zlens[2].re;
        y_yj = y - zlens[2].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[2] * x_xj * xy_j2;
        imp += mlens[2] * y_yj * xy_j2;

        dzs = complex(rep, imp);
        break;

    default:

        zs = complex(xs, ys);
        dzs = zs - z;
        for (int i = 0; i < NLENS; i++) {
            dzs = dzs + mlens[i] / conj(z - zlens[i]);
        }

    }


    // check the difference between zs(real source position) and the position computed from lens equation using solved image position z
    tempzc = complex(0, 0);
    J1 = complex(0, 0);
    J2 = complex(0, 0);
    J3 = complex(0, 0);
    dJ = complex(0, 0);

    // 只有 在 满足 设定的 SOLEPS 时，才计算 mu，所以 之前 并没有算  "false solution" 的 magnification，那么算 PointSource 的 magnification 时，也依赖于你设定的 阈值 --> 我们想要的是不管 flag 是多少，都算出一个 mu 以备用
//     if (absdzs < SOLEPS)  {
//         flag = 1;
//         for (int i = 0; i < NLENS; i++) {
//             dx_db = x - zlens[i].re;
//             dy_db = y - zlens[i].im;
//             dx_db2 = dx_db * dx_db;
//             dy_db2 = dy_db * dy_db;

//             r2_1 = dx_db2 + dy_db2 + TINY; // avoid zero in the denominator
//             r2_2 = 1 / (r2_1 * r2_1);
//             Jxx += mlens[i] * (dx_db2 - dy_db2) * r2_2;
//             Jxy += 2.0 * mlens[i] * dx_db * dy_db * r2_2;
// #ifdef parabcorr
//             tempzc = complex(z - zlens[i]);
//             tempzc2 = tempzc * tempzc;
//             tempzc3 = tempzc2 * tempzc;
//             J1 = J1 + mlens[i] / tempzc2;
//             J2 = J2 - 2 * ( mlens[i] / tempzc3 );
// #endif

//         }
// #ifdef parabcorr
//         J1c = conj(J1);
//         dJ = 1 - (J1) * J1c;
// #endif



    for (int i = 0; i < NLENS; i++) {
        dx_db = x - zlens[i].re;
        dy_db = y - zlens[i].im;
        dx_db2 = dx_db * dx_db;
        dy_db2 = dy_db * dy_db;

        r2_1 = dx_db2 + dy_db2 + TINY; // avoid zero in the denominator
        r2_2 = 1.0 / (r2_1 * r2_1);
        Jxx += mlens[i] * (dx_db2 - dy_db2) * r2_2;
        Jxy += 2.0 * mlens[i] * dx_db * dy_db * r2_2;
#ifdef parabcorr
        tempzc = complex(z - zlens[i]);
        tempzc2 = tempzc * tempzc;
        tempzc3 = tempzc2 * tempzc;
        J1 = J1 + mlens[i] / tempzc2;
        J2 = J2 - 2 * ( mlens[i] / tempzc3 );
#endif
    }
#ifdef parabcorr
    J1c = conj(J1);
    dJ = 1 - (J1) * J1c;
#endif



    //
    //analytical results for the other components of the Jacobian
    //
    Jyy = 2.0 - Jxx;
    *mu = 1.0 / (Jxx * Jyy - Jxy * Jxy);

    absdzs = abs(dzs);
    if (absdzs < SOLEPS)  {
        flag = 1;
    }
    return (flag);
}



int TripleLensing::trueSolution_qtest(double xs, double ys, complex z)
{

    flag = 0; //
    Jxx = 1.0;
    Jyy = 0.0;
    Jxy = 0.0;

    // mu = -1.0e10;        // mark imaginary soltion;
    x = z.re;
    y = z.im;

    switch (NLENS) {
    case 2:
        // binary lens
        rep = xs - x;// real part
        imp = ys - y;
        x_xj = x - zlens[0].re;
        y_yj = y - zlens[0].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        dzs = complex(rep, imp);

        break;

    case 3:
        // triple lens
        rep = xs - x;// real part
        imp = ys - y;
        x_xj = x - zlens[0].re;
        y_yj = y - zlens[0].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        x_xj = x - zlens[2].re;
        y_yj = y - zlens[2].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[2] * x_xj * xy_j2;
        imp += mlens[2] * y_yj * xy_j2;

        dzs = complex(rep, imp);
        break;

    default:

        zs = complex(xs, ys);
        dzs = zs - z;
        for (int i = 0; i < NLENS; i++) {
            dzs = dzs + mlens[i] / conj(z - zlens[i]);
        }

    }


    // check the difference between zs(real source position) and the position computed from lens equation using solved image position z
    tempzc = complex(0, 0);
    J1 = complex(0, 0);
    J2 = complex(0, 0);
    J3 = complex(0, 0);
    dJ = complex(0, 0);

    if (abs(dzs) < SOLEPS)  {
        flag = 1;
        for (int i = 0; i < NLENS; i++) {
            dx_db = x - zlens[i].re;
            dy_db = y - zlens[i].im;
            dx_db2 = dx_db * dx_db;
            dy_db2 = dy_db * dy_db;

            r2_1 = dx_db2 + dy_db2 + TINY; // avoid zero in the denominator
            r2_2 = 1.0 / (r2_1 * r2_1);
            Jxx = Jxx + mlens[i] * (dx_db2 - dy_db2) * r2_2;
            Jxy = Jxy + 2.0 * mlens[i] * dx_db * dy_db * r2_2;

            tempzc = complex(z - zlens[i]);
            tempzc2 = tempzc * tempzc;
            tempzc3 = tempzc2 * tempzc;
            tempzc4 = tempzc2 * tempzc2; // add on 2021.06.08
            J1 = J1 + mlens[i] / tempzc2;
            J2 = J2 - 2.0 * ( mlens[i] / tempzc3 );
            J3 = J3 + 6.0 * ( mlens[i] / tempzc4 );
        }

        J1c = conj(J1);
        dJ = 1.0 - (J1) * J1c;

        //
        //analytical results for the other components of the Jacobian
        //
        Jyy = 2.0 - Jxx;
        mu = 1.0 / (Jxx * Jyy - Jxy * Jxy); // only when "abs(dzs) < SOLEPS"

    }

    return (flag);
}

int TripleLensing::trueSolution_nomu(double xs, double ys, complex z, double true_solution_threshold)
{

    flag = 0; //

    x = z.re;
    y = z.im;

    switch (NLENS) {
    case 2:
        // binary lens
        rep = xs - x;// real part
        imp = ys - y;
        x_xj = x - zlens[0].re;
        y_yj = y - zlens[0].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        dzs = complex(rep, imp);

        break;

    case 3:
        // triple lens
        rep = xs - x;// real part
        imp = ys - y;
        x_xj = x - zlens[0].re;
        y_yj = y - zlens[0].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        x_xj = x - zlens[2].re;
        y_yj = y - zlens[2].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[2] * x_xj * xy_j2;
        imp += mlens[2] * y_yj * xy_j2;

        dzs = complex(rep, imp);
        break;

    default:
        zs = complex(xs, ys);
        dzs = zs - z;
        for (int i = 0; i < NLENS; i++) {
            dzs = dzs + mlens[i] / conj(z - zlens[i]);
        }
    }

    if (abs(dzs) < true_solution_threshold) {
        flag = 1;
    }

    return (flag);
}


int TripleLensing::trueSolution_withmu(double xs, double ys, complex z, double true_solution_threshold, double * mu)
{
    flag = 0; //
    Jxx = 1.0;
    Jyy = 0.0;
    Jxy = 0.0;

    x = z.re;
    y = z.im;

    switch (NLENS) {
    case 2:
        // binary lens
        rep = xs - x;// real part
        imp = ys - y;
        x_xj = x - zlens[0].re;
        y_yj = y - zlens[0].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        dzs = complex(rep, imp);

        break;

    case 3:
        // triple lens
        rep = xs - x;// real part
        imp = ys - y;
        x_xj = x - zlens[0].re;
        y_yj = y - zlens[0].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        x_xj = x - zlens[2].re;
        y_yj = y - zlens[2].im;
        xy_j2 = 1.0 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[2] * x_xj * xy_j2;
        imp += mlens[2] * y_yj * xy_j2;

        dzs = complex(rep, imp);
        break;

    default:
        zs = complex(xs, ys);
        dzs = zs - z;
        for (int i = 0; i < NLENS; i++) {
            dzs = dzs + mlens[i] / conj(z - zlens[i]);
        }
    }

    // 只有 在 满足 设定的 SOLEPS 时，才计算 mu，所以 之前 并没有算  "false solution" 的 magnification，那么你算 PointSource 的 magnification 时，也依赖于你设定的 阈值
    if (abs(dzs) < true_solution_threshold)  {
        flag = 1;
        for (int i = 0; i < NLENS; i++) {
            dx_db = x - zlens[i].re;
            dy_db = y - zlens[i].im;
            dx_db2 = dx_db * dx_db;
            dy_db2 = dy_db * dy_db;

            r2_1 = dx_db2 + dy_db2 + TINY; // avoid zero in the denominator
            r2_2 = 1 / (r2_1 * r2_1);
            Jxx += mlens[i] * (dx_db2 - dy_db2) * r2_2;
            Jxy += 2.0 * mlens[i] * dx_db * dy_db * r2_2;
        }
    }

    Jyy = 2.0 - Jxx;
    *mu = 1.0 / (Jxx * Jyy - Jxy * Jxy);


    return (flag);
}


//print out one image track
void printOneTrack(_curve * curve)
{
    int i = 0;
    for (_point *p = curve->first; p; p = p->next) { //point
#ifdef VERBOSE
        printf("I am in printOneTrack, %f %f %f %f %d %d in track: x, y, phi, \n", p->x1, p->x2, p->phi, p->mu, p->flag, i);
#endif
        i++;
    }
}

// from  contour.cpp
//print out all image tracks
void printAllTracks(_sols * track)
{
    // int i = 0;
    // int itrack = 0;
    // FILE *fp;
    // fp = fopen("data/allImages.dat", "w");
    // printf("entered printAllTracks\n");

    // for (_curve *c = track->first; c; c = c->next) { //curves
    //     for (_point *p = c->first; p; p = p->next) { //point
    //         //      printf("%f %f %f %f #%d point in the %d-th image track\n", p->x1, p->x2, p->phi, p->mu, i, itrack);
    //         fprintf(fp, "%.17g %.17g %.17g %.17g %d #%d point in the %d-th image track\n", p->x1, p->x2, p->phi, p->mu, p->flag, i, itrack);
    //         i++;
    //     }
    //     itrack++;
    // }
    // fclose(fp);

    int i = 0;
    int itrack = 0;
    FILE *fp;

    fp = fopen("data/allImages.dat", "w");
    fprintf(fp, "%d\n", track->first->length);

    for (_curve *c = track->first; c; c = c->next) { //curves
        fprintf(fp, "%d ", c->length);
    }
    fprintf(fp, "\n");

    for (_curve *c = track->first; c; c = c->next) { //curves
        for (_point *p = c->first; p; p = p->next) { //point
            // DBL_DECIMAL_DIG =  17
            fprintf(fp, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,", p->x1, p->x2, p->phi, p->mu, p->zs.re, p->zs.im, p->flag);
            i++;
        }
        itrack++;
    }
    fclose(fp);


}

// void saveTracks_before(_sols * track)
// {
//     int i = 0;
//     int itrack = 0;
//     FILE *fp;

//     fp = fopen("data/allTracks_before.dat", "w");
//     // fp = fopen(filename, "w");
//     // printf("entered saveTracks\n\n\n");

//     for (_curve *c = track->first; c; c = c->next) { //curves
//         fprintf(fp, "%d ", c->length);
//     }
//     fprintf(fp, "\n");

//     for (_curve *c = track->first; c; c = c->next) { //curves
//         for (_point *p = c->first; p; p = p->next) { //point
//             fprintf(fp, "%f %f %f %f ", p->x1, p->x2, p->phi, p->mu);
//             i++;
//         }
//         itrack++;
//     }
//     fclose(fp);
// }

void saveTracks_before(_sols * track, int nphi)
{
    int i = 0;
    int itrack = 0;
    FILE *fp;

    fp = fopen("data/allTracks_before.dat", "w");
    // fp = fopen(filename, "w");
    // printf("entered saveTracks\n\n\n");

    fprintf(fp, "%d\n", nphi);

    for (_curve *c = track->first; c; c = c->next) { //curves
        // if (c->length >= 3) { // 下面也要改
        fprintf(fp, "%d ", c->length);
        // }
    }
    fprintf(fp, "\n");

    for (_curve *c = track->first; c; c = c->next) { //curves
        // if (c->length >= 3) {
        for (_point *p = c->first; p; p = p->next) { //point
            // DBL_DECIMAL_DIG =  17
            fprintf(fp, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,", p->x1, p->x2, p->phi, p->mu, p->zs.re, p->zs.im, p->flag);
            i++;
        }
        itrack++;
        // }
    }
    fclose(fp);
}

void saveTracks(_sols * track, int nphi)
{
    int i = 0;
    int itrack = 0;
    FILE *fp;

    fp = fopen("data/allTracks.dat", "w");
    // fp = fopen(filename, "w");
    // printf("entered saveTracks\n\n\n");

    fprintf(fp, "%d\n", nphi);

    for (_curve *c = track->first; c; c = c->next) { //curves
        if (c->length >= 3) { // 下面也要改
            fprintf(fp, "%d ", c->length);
        }
    }
    fprintf(fp, "\n");

    for (_curve *c = track->first; c; c = c->next) { //curves
        if (c->length >= 3) {
            for (_point *p = c->first; p; p = p->next) { //point
                // DBL_DECIMAL_DIG =  17
                fprintf(fp, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,", p->x1, p->x2, p->phi, p->mu, p->zs.re, p->zs.im, p->flag);
                i++;
            }
            itrack++;
        }
    }
    fclose(fp);
}


// void TripleLensing::outputImagesTriple(double xsCenter, double ysCenter, int nphi, double rs)
// {
//     _sols *imageTracks;
//     imageTracks = outputTracks(xsCenter, ysCenter, rs, nphi, 1.5);
//     if (!imageTracks)  return;
// #ifdef VERBOSE
//     printf("I am in outputImagesTriple, Number of image tracks: %d\n",  imageTracks->length);
// #endif
//     // Successful, print out image tracks
//     FILE *fileImage;
//     fileImage = fopen("data/images.dat", "w");
//     fprintf(fileImage, "%f %f %f\n", xsCenter, ysCenter, rs);
//     fprintf(fileImage, "%d\n", imageTracks->length);

//     for (_curve *c = imageTracks->first; c; c = c->next) { //curves
//         fprintf(fileImage, "%d\n", c->length);
//         for (_point *p = c->first; p; p = p->next) { //point
//             fprintf(fileImage, "%f %f %f %f \n", p->x1, p->x2, p->phi, p->mu);
//         }
//     }
//     fclose(fileImage);
// }

void get_crit_caus(double mlens[], complex zlens[], int nlens, int NPS, double * criticalx, double * criticaly, double * causticsx, double * causticsy, int *numcritical) {

    // double resxy[4*NPS+2]; //wrong
    double resxy[40 * NPS];
    double Zlens[nlens * 2];
    int i;
    for (i = 0; i < nlens; i++) {
        Zlens[i] = zlens[i].re;
        Zlens[i + nlens] = zlens[i].im;
    }
    outputCriticalTriple_list(resxy, mlens, Zlens, nlens, NPS);
    fprintf(stderr, "out of outputCriticalTriple_list\n");
    int numcrit = (int)resxy[0];
    *numcritical = numcrit;
    fprintf(stderr, "%d numcrit\n", numcrit);
    int offset = 2 * numcrit + 1;
    for (i = 0; i < numcrit; i++) {
        criticalx[i] = resxy[2 * i + 1];
        criticaly[i] = resxy[2 * i + 2];
        causticsx[i] = resxy[offset + 2 * i + 1];
        causticsy[i] = resxy[offset + 2 * i + 2];
    }
}


//output in two files the triple critical curves and caustics
void outputCriticalTriple_list(double allxys[], double mlens[], double zlens[], int nlens, int NPS)
{
    // VBBinaryLensing VBBL;
    // double allxys[4*NPS+2];
    complex Zlens[nlens];
    int i;
    for (i = 0; i < nlens; i++) {
        Zlens[i] = complex(zlens[i], zlens[i + nlens]);
    }

    int ncurves = 0;
    int ncritical = 0;
    _sols *criticalCurves;

    criticalCurves = VBBL.PlotCritTriple(mlens, Zlens, NPS, nlens);

    // check how many closed critical curves we have
    // first halfs are critical curves
    ncritical = criticalCurves->length / 2; // number of closed critical curves
    // allxys[0] = ncritical;
    // fprintf(stderr, "ncritical %d \n", ncritical);

    // int numcnt_at_each_single_caus[ncritical] = {0};
    int* numcnt_at_each_single_caus = new int[ncritical];

#ifdef VERBOSE
    printf("I am in outputCriticalTriple, Number of closed critical curves: %d\n", ncritical);
#endif

    // write out the critical curves and caustics separately

    ncurves = 0;
    int count_critical = 0;
    int count_caustic = 0;
    for (_curve *c = criticalCurves->first; c; c = c->next) {
        //int npoints = 0;
        ncurves++;

        // second half, caustics
        if (ncurves > ncritical) {      // second halfs are caustics
            for (_point *p = c->first; p; p = p->next) {
                //npoints++;

                count_caustic ++;
                allxys[2 * count_critical + 1 + 2 * count_caustic - 1] = p->x1;
                allxys[2 * count_critical + 1 + 2 * count_caustic] = p->x2;
            }
        } else {
            numcnt_at_each_single_caus[ncurves-1] = c->length;
            // first half, critical curves
            for (_point *p = c->first; p; p = p->next) { // first halfs are critical curves
                count_critical ++;
                allxys[count_critical * 2 - 1] = p->x1;
                allxys[count_critical * 2] = p->x2;
            }
        }

    }
    allxys[0] = count_critical;
    allxys[2 * count_critical + 1] = count_caustic;

    // save ncritical the and number of points at each critical points // 2021.10.05
    allxys[2 * count_critical + 1 + 2 * count_caustic + 1] = ncritical;
    for(int i = 0; i<ncritical; i++){
        allxys[2 * count_critical + 2 * count_caustic + 3 + i] = numcnt_at_each_single_caus[i];
    }

    delete[] numcnt_at_each_single_caus;

}


//output in two files the triple critical curves and caustics
void outputCriticalTriple(double mlens[], complex zlens[], int nlens, int npnt)
{
    fprintf(stderr, "Generating critical curves and caustics ...");
    // VBBinaryLensing VBBL;
    int ncurves = 0;
    int ncritical = 0;
    _sols *criticalCurves;
    FILE *fcritical, *fcaustics;
    criticalCurves = VBBL.PlotCritTriple(mlens, zlens, npnt, nlens);
    fcritical = fopen("./data/critical_curves.dat", "w");
    fcaustics = fopen("./data/caustics.dat", "w");

    // check how many closed critical curves we have
    // first halfs are critical curves
    ncritical = criticalCurves->length / 2; // number of closed critical curves
#ifdef VERBOSE
    printf("I am in outputCriticalTriple, Number of closed critical curves: %d\n", ncritical);
#endif

    // write out the critical curves and caustics separately

    ncurves = 0;
    for (_curve *c = criticalCurves->first; c; c = c->next) {
        int npoints = 0;

        ncurves++;

        if (ncurves > ncritical) {      // second halfs are caustics
            for (_point *p = c->first; p; p = p->next) {
                // fprintf(stderr,  "outputCriticalTriple, %.10lf %.10lf\n", p->x1, p->x2);
                fprintf(fcaustics, "%.10lf %.10lf\n", p->x1, p->x2);

                npoints++;
            }
        } else {
            for (_point *p = c->first; p; p = p->next) { // first halfs are critical curves
                fprintf(fcritical, "%.10lf %.10lf\n", p->x1, p->x2);
            }
        }
    }
    fclose(fcritical);  fclose(fcaustics);
    fprintf(stderr, " ... done.\n");
}



#define _addsubarea \
subArea +=  2*( (y2 + y1) * (x2 - x1) );

#define _addsubarea_parab \
subArea +=  (0.5*( (y2 + y1) * (x2 - x1) ) + one_24*(ds1+ds2)*dphi3);\

#define _addsubarea_parab2 \
subArea +=  (0.5*( (y2 + y1) * (x2 - x1) ) - one_24*(ds1+ds2)*dphi3);\


// double areaFunc(_sols * track, double rs,  int paracorr)
double TripleLensing::areaFunc(_sols * track, double rs, int finalnphi, double muPS, double muPrev, int *quality)
{
    area = 0.0;
    pos_area = 0.0;
    // int parity,
    trackcnt = 0;
    subArea = 0.0;
    // int special_flag = 0; // special_flag when dealing with the parity
    int parity_different_flag = 0;

    // double max_normal_area = -1; // normal means parity is normal
    // double max_abnormal_area = -1; //

#ifdef VERBOSE
    int cnti = 0;
    for (_curve *c = track->first; c; c = c->next) {

        if ( c->length > finalnphi * 0.25 ) { // just pring long tracks

            fprintf(stderr, "\t a new curve (length = %d, original parity = %d), p->phi, p->mu: \n", c->length, c->parity);

            cnti = 0;
            for (_point *p = c->last; p; p = p->prev) {
                fprintf(stderr, "%f %f, flag = %d, absdzs = %.5e\n", p->phi * 180 / M_PI, p->mu, p->flag, p->absdzs);
                cnti ++;
                if (cnti == 10) {
                    break;
                }
            }

            fprintf(stderr, "\t\t --- tail above, head below  -- \n");

            cnti = 0;
            for (_point *p = c->first; p; p = p->next) {
                fprintf(stderr, "%f %f, flag = %d, absdzs = %.5e\n", p->phi * 180 / M_PI, p->mu, p->flag, p->absdzs);
                cnti ++;
                if (cnti == 10) {
                    break;
                }
            }

        }


    }
#endif


    for (_curve *c = track->first; c; c = c->next) { //curves
        special_flag = 0;
        // before 2021.06.08
        // if ( c->length < 3) { // ignore tracks that are short than 3
        //     continue;
        // }

        // // after 2021.06.08
        if ( c->length < finalnphi * 0.25 ) { // ignore tracks that are short than the expected value
            continue;
        }


        subArea = 0.0;

        x1 = c->first->x1;
        y1 = c->first->x2;
        phi1 = c->first->phi;


        if (c->parity) {
            parity = c->parity; // original parity
        } else {
            fprintf(stderr, "wrong, a track with length %d has no parity assigned \n", c->length);
        }

        _c_parity;
        pscan = c->first;
        while ( pscan->next && pscan->next->next && ((pscan->mu * pscan->next->mu < 0) || ( abs(abs(pscan->phi - pscan->next->phi) - M_PI2) < EPS ) || ( abs(abs(pscan->next->phi - pscan->next->next->phi) - M_PI2) < EPS ) || abs(pscan->phi - pscan->next->phi) < EPS )) {
            pscan = pscan->next;
        }
        ph1 = pscan->phi; ph2 = pscan->next->phi; ph3 = pscan->next->next->phi;
        // if ( (special_flag==1) && (ph1 > ph2) && (ph2 > ph3)) {
        if ( (ph1 > ph2) && (ph2 > ph3)) {
            // fprintf(stderr, "phi1 > phi2 > phi3, c->parity before = %d \n", c->parity);
            c->parity *= -1;
        }

        for (_point *p = c->first->next; p; p = p->next) { //point

            x2 = p->x1; y2 = p->x2;
            _addsubarea

            x1 = x2; y1 = y2;
        }
        x1 = c->last->x1;
        y1 = c->last->x2;
        x2 = c->first->x1;
        y2 = c->first->x2;
        _addsubarea


        area += subArea * parity;
        pos_area += abs(subArea);
#ifdef VERBOSE
        fprintf(stderr, "track->length %d, c->length: %d, c->parity = %d, c->first->mu = %.5e, subArea = %f , subArea * parity: %.5e\n", track->length, c->length, parity, c->first->mu , subArea * 0.25 , subArea * parity * 0.25);
#endif


        if (parity != c->parity) {
            parity_different_flag = 1;
#ifdef VERBOSE
            fprintf(stderr, "be cautious, the parity of a track with length %d might be wrong\n", c->length);
#endif
        }

        // }

    }


    // somtimes parity may not work
    area = abs(area) * 0.25;
    pos_area *= 0.25;

    // the parity is mainly used for handling the situation when a ring like image is formed, in such case the area of two boundaries should be - instead of +
    // when a ring like structure is formed,

    // case 1, if the parity of one track is wrong, the magnification should be different, but not so different with
#ifdef VERBOSE
    fprintf(stderr, "pos_mag = %f, parity mag = %f\n", pos_area / areaSource, area / areaSource); // 5.239782
#endif


    // if (area < areaSource || ( (pos_area/area < muPS || area/areaSource < muPS) && parity_different_flag && muPS < 50) ) { // if muPS > 50, a ring probably formed
    if ( abs(area - pos_area) > EPS || (area < areaSource || ( parity_different_flag ) ))  { // if muPS > 50, a ring probably formed
        // fprintf(stderr, "use abs area, pos_area/areaSource = %f\n", pos_area/areaSource);
        // return pos_area * 0.25;
        *quality = 0;
        // return the area that are more close to the previous mu
        // abs(muPS - muPrev) < EPS means if muPrev == muPS, we trust the pos_area more
        if ( abs(area - areaSource * muPrev) > (pos_area - areaSource * muPrev) || (abs(muPS - muPrev) < EPS && parity_different_flag ) ) {
            return pos_area;
        } else {
            return area;
        }
    } else {
        *quality = 1;
        return area;
    }
}


#define _dphijump \
dphi = sqrt( ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) )/abs(dz1*dz2) );

double TripleLensing::areaFunc_parab(_sols * track, double rs)
{
    // double x1, x2, y1, y2, phi1, phi2, dphi, ds1, ds2;
    area = 0.0;
    pos_area = 0.0;
    // int parity,
    trackcnt = 0;
    subArea = 0.0;
#ifdef VERBOSE
    int cnti = 0;
    for (_curve *c = track->first; c; c = c->next) {
        fprintf(stderr, "\t a new curve, p->phi, p->mu: \n");
        cnti = 0;
        for (_point *p = c->first; p; p = p->next) {
            fprintf(stderr, "%f %f\n", p->phi * 180 / M_PI, p->mu);
            cnti ++;
            if (cnti == 10) {
                break;
            }
        }
    }
#endif


    for (_curve *c = track->first; c; c = c->next) { //curves
        if ( c->length < 3) {
            continue;
        }

        c->first->closepairparity = c->first->mu > 0 ? 1 : -1;
        c->last->closepairparity = c->last->mu > 0 ? 1 : -1;
#ifdef VERBOSE
        fprintf(stderr, "\t\t\t firtt mu, parity %f %d\n", c->first->mu, c->first->closepairparity);
        fprintf(stderr, "\t\t\t last mu, parity %f %d\n", c->last->mu, c->last->closepairparity);
#endif

        subArea = 0.0;

        x1 = c->first->x1;
        y1 = c->first->x2;
        phi1 = c->first->phi;
        c->first->phi = phi1 > 0 ? phi1 : M_PI2 + phi1;
        ds1 = c->first->ds;
        clpairparity1 = c->first->closepairparity;
        dz1 = c->first->dz;
        _c_parity;

        pscan = c->first;
        while ( pscan->next && pscan->next->next && ((pscan->mu * pscan->next->mu < 0) || ( abs(abs(pscan->phi - pscan->next->phi) - M_PI2) < 1e-10 ) || ( abs(abs(pscan->next->phi - pscan->next->next->phi) - M_PI2) < 1e-10 )  )) {
            pscan = pscan->next;
        }
        ph1 = pscan->phi; ph2 = pscan->next->phi; ph3 = pscan->next->next->phi;
        if ( ph1 > ph2 && ph2 > ph3) {
            c->parity *= -1;
        }
#ifdef VERBOSE
        fprintf(stderr, "ph1 ph2 ph3 %f %f %f, ph2 - ph3 = %f , c->parity %d\n", ph1 * 180 / M_PI, ph2 * 180 / M_PI, ph3 * 180 / M_PI, (abs(ph2 - ph3) - M_PI2) * 180 / M_PI, c->parity);
#endif
        parity = c->parity;
        for (_point *p = c->first->next; p; p = p->next) { //point

            x2 = p->x1; y2 = p->x2;
            phi2 = p->phi;
            p->phi = phi2 > 0 ? phi2 : M_PI2 + phi2;
            ds2 = p->ds;
            clpairparity2 = p->closepairparity;
            dz2 = p->dz;

            if (clpairparity1 == 1 && clpairparity2 == -1) {

                ds2 *= -1;
                _dphijump
                dphi3 = dphi;
                dphi3 *= dphi;
                dphi3 *= dphi;
                if (c->parity == -1) {
                    // _addsubarea_parab2
                    _addsubarea_parab
                } else {
                    _addsubarea_parab
                }

#ifdef VERBOSE
                fprintf(stderr, "\t\t\t jump over caustics, clpairparity1, clpairparity2 = %d %d\n", clpairparity1, clpairparity2);
#endif
                clpairparity1 = clpairparity2;
                phi1 = phi2;
                ds1 = p->ds;
                dz1 = dz2;
                x1 = x2;
                y1 = y2;


            } else if (clpairparity1 == -1 && clpairparity2 == 1) {
                ds1 *= -1;
                _dphijump
                dphi3 = dphi;
                dphi3 *= dphi;
                dphi3 *= dphi;

                if (c->parity == -1) {
                    // _addsubarea_parab2
                    _addsubarea_parab
                } else {
                    _addsubarea_parab
                }
#ifdef VERBOSE
                fprintf(stderr, "\t\t\t jump over caustics, clpairparity1, clpairparity2 = %d %d\n", clpairparity1, clpairparity2);
#endif
                clpairparity1 = clpairparity2;
                ds1 = p->ds;
                phi1 = phi2;
                dz1 = dz2;
                x1 = x2; y1 = y2;

            } else {
                dphi = (phi2 - phi1);
                if ( abs( abs(dphi) - M_PI2) < M_PI ) {
                    dphi = dphi + (dphi > 0 ? -1 : 1) * M_PI2;
                }
                dphi3 = abs(dphi);
                dphi3 *= dphi;
                dphi3 *= dphi; //

                _addsubarea_parab
                x1 = x2; y1 = y2;
                ds1 = ds2;
                phi1 = phi2;
                clpairparity1 = clpairparity2;
            }

        }
        x1 = c->last->x1;
        y1 = c->last->x2;
        x2 = c->first->x1;
        y2 = c->first->x2;

        phi1 = c->last->phi;
        phi2 = c->first->phi;
        ds1 = c->last->ds;
        ds2 = c->first->ds;
        clpairparity1 = c->last->closepairparity;
        clpairparity2 = c->first->closepairparity;
        dz1 = c->last->dz;
        dz2 = c->first->dz;


        if (clpairparity1 == 1 && clpairparity2 == -1) {
            // 1 - 2 --> 2 - 1
            ds1 *= -1;
            _dphijump
            dphi3 = dphi;
            dphi3 *= dphi;
            dphi3 *= dphi;
            if (c->parity == -1) {
                _addsubarea_parab2
            } else {
                _addsubarea_parab
            }
#ifdef VERBOSE
            fprintf(stderr, "\t\t\t jump over caustics, clpairparity1, clpairparity2 = %d %d\n", clpairparity1, clpairparity2);
#endif
        } else if (clpairparity1 == -1 && clpairparity2 == 1) {
            ds2 *= -1;
            _dphijump
            dphi3 = dphi;
            dphi3 *= dphi;
            dphi3 *= dphi;
            if (c->parity == -1) {
                _addsubarea_parab2
            } else {
                _addsubarea_parab
            }
#ifdef VERBOSE
            fprintf(stderr, "\t\t\t jump over caustics, clpairparity1, clpairparity2 = %d %d\n", clpairparity1, clpairparity2);
#endif
        } else {

            dphi = (phi2 - phi1);

            if ( abs( dphi + M_PI2 ) < M_PI ) {
                dphi = (dphi > 0 ? 1 : -1 ) * abs(dphi + M_PI2);
            }

            dphi3 = abs(dphi);
            dphi3 *= dphi;
            dphi3 *= dphi;

            _addsubarea_parab

        }

        area += subArea * parity;
        pos_area += abs(subArea);
#ifdef VERBOSE
        fprintf(stderr, "track->length %d, c->length: %d, c->parity = %d, c->first->mu = %f, subArea * parity: %f\n", track->length, c->length, c->parity, c->first->mu , subArea * parity);
#endif
    }
    // somtimes parity may not work
    area = abs(area);
    if (area < areaSource) {
        return pos_area;
    } else {
        return area;
    }


}

// _linkedarray *getphis_v3
_linkedarray *TripleLensing::getphis_v3(double xsCenter, double ysCenter, double rs) {
    _linkedarray *PHI;
    PHI = new _linkedarray;
#ifdef VERBOSE
    fprintf(stderr, "caustic crossing\n");
#endif
#ifdef verbose
    fprintf(stderr, "caustic crossing\n");
#endif
    distype = 2;
    double phi[secnum_priv];
    double mus[secnum_priv];
    double dphi, scanphi = 0;
    double minmu = 1e9;
    double maxmu = -1e9;
    dphi = 2 * M_PI / secnum_priv;
    for (int i = 0; i < secnum_priv; i++) {
        phi[i] = scanphi;
        xs = xsCenter + rs * cos(scanphi);
        ys = ysCenter + rs * sin(scanphi);

        mus[i] = TriplePS(xs, ys);

        if (mus[i] <= minmu) {
            minmu = mus[i];
            maxmuidx = i; // the name is not important
        }
        if (mus[i] >= maxmu) {
            maxmu = mus[i];
        }

        scanphi += dphi;
    }

#ifdef VERBOSE
    fprintf(stderr, "secnum_priv %d,basenum %d,  maxmuidx: %d, minmu = %f, maxmu = %f\n", secnum_priv, basenum_priv, maxmuidx, minmu, maxmu);
#endif
    for (int i = 0; i < secnum_priv; i++) {
        mus[i] = (mus[i] / minmu);
    }

    double psf[3] = {1, 2, 1};

    int npmus[secnum_priv];
    int secnum_privlist[secnum_priv];
    double offset[secnum_priv];
    npmus[0] = ceil(psf[0] * mus[secnum_priv - 1] + psf[1] * mus[0] + psf[2] * mus[1]) + 1;
    npmus[secnum_priv - 1] = ceil(psf[0] * mus[secnum_priv - 2] + psf[1] * mus[secnum_priv - 1] + psf[2] * mus[0]) + 1;


    if (0 < maxmuidx) {
        offset[0] = 2 * M_PI;
    } else {
        offset[0] = 0;
    }
    if (0 < secnum_priv - maxmuidx) {
        secnum_privlist[0] = maxmuidx;
    } else {secnum_privlist[0] = 0 - secnum_priv + maxmuidx;}

    for (int i = 1; i < secnum_priv - 1; i++) {
        npmus[i] = floor(psf[0] * mus[i - 1] + psf[1] * mus[i] + psf[2] * mus[i + 1]
                        ) + 1;

        if (i < maxmuidx) {
            offset[i] = 2 * M_PI;
        } else {
            offset[i] = 0;
        }
        if (i < secnum_priv - maxmuidx) {
            secnum_privlist[i] = i + maxmuidx;
        } else {secnum_privlist[i] = i - secnum_priv + maxmuidx;}


    }

    if (secnum_priv - 1 < maxmuidx) {
        offset[secnum_priv - 1] = 2 * M_PI;
    } else {
        offset[secnum_priv - 1] = 0;
    }
    if ( 1 >  maxmuidx) {
        secnum_privlist[secnum_priv - 1] = secnum_priv - 1 + maxmuidx;
    } else {secnum_privlist[secnum_priv - 1] = - 1 + maxmuidx;}
    dphi = M_PI / secnum_priv;
    int tempidx = 0;
    tempidx = secnum_privlist[0];
    PHI->linspace(offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, npmus[tempidx]*basenum_priv, 0);
    for (int i = 1; i < secnum_priv - 1; i++) {
        tempidx = secnum_privlist[i];
        PHI->linspace(offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, npmus[tempidx]*basenum_priv, 0);
    }
    tempidx = secnum_privlist[secnum_priv - 1];
    PHI->linspace(offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, npmus[tempidx]*basenum_priv, 1);

    return PHI;
}


double myatan(double x, double y) {
    // # let angle to be 0~2*M_PI
    if (x >= 0 && y == 0) {
        return 0.0;
    }
    if (x == 0 && y > 0) {
        return M_PI / 2;
    }
    if (x == 0 && y < 0) {
        return 3 * M_PI / 2;
    }
    if (y == 0 && x < 0) {
        return M_PI;
    }
    double ang = atan(y / x);
    if (ang > 0) {
        if (y > 0) {
            return ang;
        }
        else {
            return M_PI + ang;
        }
    }
    else {
        if (y < 0) {
            return 2 * M_PI + ang;
        }
        else {
            return M_PI + ang;
        }
    }
}

_linkedarray::_linkedarray(void) {
    length = 0;
    first = last = 0;
}


_linkedarray::~_linkedarray(void) {
    Node *scan1, *scan2;
    // scan1 = first;
    // for (int i = 0; i < length; i++) {
    //     scan2 = scan1->next;
    //     delete scan1;
    //     scan1 = scan2;
    // }
    // scan1 = NULL;
    // scan2 = NULL;
    // delete scan1;
    // delete scan2;

    scan1 = first;
    while (scan1) {
        scan2 = scan1->next;
        delete scan1;
        scan1 = scan2;
    }



}



void _linkedarray::append(double x) {
    Node *cc;
    cc = new Node(x);
    if (length == 0) {
        first = cc;
        last = cc;
        cc->prev = 0;
    }
    else {
        last->next = cc;
        cc->prev = last;
        last = cc;
    }
    cc->next = 0;
    length++;
}

void _linkedarray::append_behind(Node * place, double x) {
    Node *cc;
    cc = new Node(x);

    cc->next = place->next;
    place->next->prev = cc;
    cc->prev = place;
    place->next = cc;
    length++;
}

void _linkedarray::append_before(Node * place, double x) {
    Node *cc;
    cc = new Node(x);

    if (place->prev == 0 ) {
        cc->next = place;
        place->prev = cc;
        first = cc;
        cc->prev = 0;
        length ++;
    } else if (place->next == 0 ) {
        cc->prev = place->prev;
        place->prev->next = cc;

        cc->next = place;
        place->prev = cc;

    } else {
        Node *scan;
        scan = first;
        while (scan->next != place) {
            scan = scan->next;
        }

        cc->next =  scan->next;
        scan->next->prev = cc;
        scan->next = cc;
        cc->prev = scan;
        length++;
    }
}


void _linkedarray::linspace(double phi0, double phiend, int nphi,  int endpoint) {
    double dphi;
    double scanphi = phi0;
    if (endpoint) {
        dphi = (phiend - phi0) / (nphi - 1);
    } else {
        dphi = (phiend - phi0) / (nphi);
    }
    for (int i = 0; i < nphi; i++) {
        append(scanphi);
        scanphi += dphi;
    }
}

void _linkedarray::concatenate(_linkedarray * nc) {
    if (length > 0) {
        last->next = nc->first;
        nc->first->prev = last->next;
    }
    else {
        first = nc->first;
    };
    if (nc->length > 0) {
        // nc->first->prev = last;
        last = nc->last;
    }
    length += nc->length;
    nc->first = 0;
    nc->last = 0;
    nc->length = 0;
    delete nc;
}

void _linkedarray::extend() {

    if (length > 0) {
        Node *scan;
        Node *cc;
        // fprintf(stderr, "%d: ", length);
        double midvalues[length - 1];
        scan = first;
        for (int i = 0; i < length - 1; i++) {
            // fprintf(stderr, "%d, ", i);
            midvalues[i] = (scan->value + scan->next->value) / 2;
            scan = scan->next;
        }
        scan = first;
        int templen = length;
        for (int i = 0; i < templen - 1; i++) {
            cc = new Node( midvalues[i] );
            cc->next = scan->next;
            scan->next->prev = cc;

            scan->next = cc;
            cc->prev = scan;

            scan = cc->next;
            length++;
        }

    }
    else {
        fprintf(stderr, "length = 0\n");
    };
}

void _linkedarray::print() {
    if (length > 0) {
        Node *scan;
        scan = first;
        fprintf(stderr, "phis: \n");
        for (int i = 0; i < length; i++) {
            fprintf(stderr, "%f ,", scan->value / M_PI * 180);
            scan = scan->next;
        }
        fprintf(stderr, "\n");
    }
    else {
        fprintf(stderr, "length = 0\n");
    };
}

// #undef EPS_CLOSE
// #undef MAXIT_NEWTON

#define _toterr \
toterr = def_err;


// void bisec
void TripleLensing::bisec(_curve * lightkv, _point * scan2, _point * scan3, double errorTol_critical, double tE_inv, double u0, double salpha, double calpha, double t0, double rs) {
    int np1 = 1;
    double mag2, mag1, tm1, tm2, def_err, midtm, midtm_mag, midtm_mag_inter, timerr, toterr;
    double subdt, curr_t, curr_mag, tn, y1, y2;
    subdt = (scan3->x1 - scan2->x1) / (np1 + 1);
    _point *subscan, *pnt2beAdd, *scan;
    subscan = scan2;
    for (int i3 = 0; i3 < np1; i3++) {
        curr_t = scan2->x1 + subdt * (i3 + 1);
        tn = (curr_t - t0) * tE_inv;
        y1 = u0 * salpha + tn * calpha; //
        y2 = u0 * calpha - tn * salpha; //
        curr_mag = subscan->area;
        pnt2beAdd = new _point(curr_t, curr_mag, 0);
        pnt2beAdd->next = subscan->next;
        subscan->next->prev = pnt2beAdd;
        subscan->next = pnt2beAdd;
        pnt2beAdd->prev = subscan;
        lightkv->length += 1;
        subscan = subscan->next;
    }

    // update error
    scan = scan2;
    while (scan != scan3) {
        mag2 = scan->next->x2;
        mag1 = scan->x2;
        tm2 = scan->next->x1;
        tm1 = scan->x1;
        // compute the midtm_mag,
        midtm = (tm1 + tm2) * 0.5;
        tn = (midtm - t0) * tE_inv;
        y1 = u0 * salpha + tn * calpha; //
        y2 = u0 * calpha - tn * salpha; //
        midtm_mag = TripleMag( y1, y2, rs);
        midtm_mag = log10(midtm_mag - 1);
        scan->area = midtm_mag;
        midtm_mag_inter = (mag1 + mag2) * 0.5;
        midtm_mag_inter = pow(10, midtm_mag_inter) + 1;
        midtm_mag = pow(10, midtm_mag) + 1;
        def_err = fabs( (midtm_mag_inter - midtm_mag) / midtm_mag );
        _toterr
        scan->error = toterr;

        scan = scan->next;
    }

    scan = scan2;
    if (lightkv->length < MAXNPS) {
        while (scan != scan3) {
            if (scan->error >= adaerrTol) {
                timerr = 1e100;
                if (scan->prev && scan->next) {
                    timerr = fmax( fabs( (scan->x1 - scan->next->x1)),  fabs( (scan->x1 - scan->prev->x1)) );
                }

                if (timerr > timerrTol) {
                    bisec(lightkv, scan, scan->next, errorTol_critical, tE_inv, u0, salpha, calpha, t0, rs);
                }

            }
            scan = scan->next;
        }
    }
}

//tri_ada_sample_1D
void TripleLensing::TripleLkvAdap( double rs, _curve * lightkv, double t_start, double t_end, double alpha, double tE, double t0 , double u0) {
    double dt0 = (t_end - t_start) / (adanp0 - 1);
    double curr_t = t_start, curr_mag, y1, y2, mag1, mag2, def_err, tm1, tm2, midtm_mag, midtm_mag_inter , midtm, toterr;

    double errorTol_critical = adaerrTol;

    double salpha = sin(alpha), calpha = cos(alpha), tn, tE_inv = 1 / tE;
    _point *scan2;//error_s_pnt: smaller error pnt

    for (int i = 0; i < adanp0; i++) {
        curr_t = t_start + i * dt0;
        tn = (curr_t - t0) * tE_inv;
        y1 = u0 * salpha + tn * calpha; //
        y2 = u0 * calpha - tn * salpha; // Mine
        // curr_mag = tripleFS2python( y1, y2, rs);
        curr_mag = TripleMag( y1, y2, rs);
        curr_mag = log10(curr_mag - 1);

        lightkv->append(curr_t, curr_mag);
        if (i > 0) {
            mag2 = curr_mag;
            mag1 = lightkv->last->prev->x2;
            tm2 = curr_t;
            tm1 = lightkv->last->prev->x1;

            // compute the midtm_mag,
            midtm = (tm1 + tm2) * 0.5;
            tn = (midtm - t0) * tE_inv;
            y1 = u0 * salpha + tn * calpha; // Mine
            y2 = u0 * calpha - tn * salpha; // Mine
            // midtm_mag = tripleFS2python(y1, y2, rs);
            midtm_mag = TripleMag(y1, y2, rs);
            midtm_mag = log10(midtm_mag - 1);

            lightkv->last->prev->area = midtm_mag;

            // use real mu rather than log(mu - 1)
            midtm_mag = pow(10, midtm_mag) + 1;

            midtm_mag_inter = (mag1 + mag2) * 0.5;
            midtm_mag_inter = pow(10, midtm_mag_inter) + 1;
            def_err = fabs( (midtm_mag_inter - midtm_mag) / midtm_mag );
            // def_err = fabs( (midtm_mag_lin - midtm_mag_inter)/midtm_mag_lin );

            _toterr
            lightkv->last->prev->error = toterr;
            //using area to save midtm_mag
        }
    }

// add pnt where error is not satisfy the requirement
    // double subdt = 0;
    _point *scan3;

    scan2 = lightkv->first;
    while (scan2->next) {
        scan3 = scan2->next;
        if (lightkv->length < MAXNPS) {
            if (scan2->error > adaerrTol) {
                bisec(lightkv, scan2, scan2->next, errorTol_critical, tE_inv, u0, salpha, calpha, t0, rs);
            }
        }
        scan2 = scan3;
    }
}



// obtain the polynomical coefficients
// input: lens mass and positions, mlens[], zlens[], source positions: xs, ys
// output: c

void TripleLensing::polynomialCoefficients(double xs, double ys, complex c[])
{
    int i, j, k;
    // complex zs;
    zs = complex(xs, ys);
    zsc = conj(zs);

    for (i = 0; i < NLENS; i++) {
        /* denominator */
        for (j = 0; j <= NLENS; j++) { /* (zsc-conjugate(z[i])) * Product_{j=1}^N (z-z[j]) */
            q[j][i] = (zsc - zc[i]) * temp_const1[i][j];
        }

        /* Sum_{j=1}^n Product (z-z_k), k=1, n, but k !=j. */
        for (j = 0; j < NLENS; j++) {
            /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */

            /* doing the sum */
            for (k = 0; k < NLENS; k++) {
                q[k][i] = q[k][i] + mlens[j] * temp_const2[i][j][k];
                // q[k][i] = q[k][i] + mlens[j] * temp[k];
            }
        }
    }

    /* now clear the fractions, z-zs - Sum_i p_i/q_i */
    /* get the polynomial Product q_i, i=1, n */

    /* first term */
    qtemp[0] = 1.0;
    degrees = 0;
    for (i = 0; i < NLENS; i++) {
        for (j = 0; j <= NLENS; j++) {
            qtemp2[j] = q[j][i];
        }

        multiply(qtemp, degrees, qtemp2, NLENS, ctemp);

        degrees += NLENS;

        for (j = 0; j <= degrees; j++) {
            qtemp[j] = ctemp[j];
        }
    }

    /* get coefficients of (z-zs) Product_i=1^n q_i */
    multiply_z(ctemp, zs, degrees);

    /* copy the coefficients */
    for (i = 0; i < DEGREE + 1; i++) {
        c[i] = ctemp[i];
    }

    /* second term */
    for (i = 0; i < NLENS; i++) {
        degrees = 0;
        qtemp[0] = 1.0;
        for (j = 0; j < NLENS; j++) {
            if (j == i) continue;

            for (k = 0; k <= NLENS; k++) {
                qtemp2[k] = q[k][j];
            }

            multiply(qtemp, degrees, qtemp2, NLENS, ctemp);

            degrees += NLENS;

            for (k = 0; k <= degrees; k++) {
                qtemp[k] = ctemp[k];
            }
        }

        for (k = 0; k <= NLENS; k++) {
            ptemp[k] = p_const[k][i];
        }

        multiply(qtemp, degrees, ptemp, NLENS, ctemp);

        for (k = 0; k < DEGREE; k++) {
            c[k] = c[k] - ctemp[k];
        }
    }
}

/* multiply a polynomial of degree n by (z-a), returns a polynomial of degree n+1 */
void multiply_z(complex c[], complex a, int n)
{
    int j;

    c[n + 1] = c[n];
    for (j = n; j >= 1; j--) c[j] = c[j - 1] - c[j] * a;
    c[0] = c[0] * (-a);
}

void multiply_z_v2(complex c[][NLENS + 1], complex a, int n, int firstdim)
{
    int j;

    c[firstdim + 1][n + 1] = c[firstdim][n];
    if (firstdim == 0) {
        c[0][n + 1] = c[0][n];
        for (j = n; j >= 1; j--) {
            c[0][j] = c[0][j - 1] - c[0][j] * a;
        }
        c[0][0] = c[0][0] * (-a);
    } else
    {
        // for (int i = 0; i<=NLENS; i++) c[firstdim][i] = c[firstdim-1][i];
        // c[firstdim][0] = 1.0;

        c[firstdim][n + 1] = c[firstdim][n];
        for (j = n; j >= 1; j--) {
            c[firstdim][j] = c[firstdim][j - 1] - c[firstdim][j] * a;
        }
        c[firstdim][0] = c[firstdim][0] * (-a);

    }

}


/* multiply a polynomial of degree n by (z-a)^2, returns a polynomial of degree n+2 */
void multiply_zSquare(complex c[], complex a, int n)
{
    multiply_z(c, a, n);
    multiply_z(c, a, n + 1);
}

/* multiply two polynomials of degree of na and nb, input as a(na), b(nb) */
void multiply(complex a[], int na, complex b[], int nb, complex c[])
{
    int i, j;

    /* zero the array */
    for (i = 0; i <= na + nb; i++) {
        c[i] = 0.0;
    }

    /* now do the product */
    for (i = 0; i <= na; i++) {
        for (j = 0; j <= nb; j++) {
            c[i + j] = c[i + j] + a[i] * b[j];
        }
    }
}
