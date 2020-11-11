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

        // if(i == 0){
        //     temp_const1[i][0] = 1.0;
        //     temp_const2[i][0] = 1.0;
        // }

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

        // if(i == 0){
        //     temp_const1[i][0] = 1.0;
        //     temp_const2[i][0] = 1.0;
        // }

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
    // fprintf(stderr, "in reset2, mlens[0 1 2] %f %f %f\n", mlens[0], mlens[1], mlens[2] );

    // nphi = 32;
    // secnum = 45;
    // basenum = 1;
    // quad_err = 0;
    // quaderr_Tol = 1e-4;
    // CQ = 1;
    // relerr_mag = 1e-3;
    // TINY = 1.0e-20;


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

        // if(i == 0){
        //     temp_const1[i][0] = 1.0;
        //     temp_const2[i][0] = 1.0;
        // }

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
    // fprintf(stderr, "in reset2, mlens[0 1 2] %f %f %f\n", (this -> mlens)[0], (this -> mlens)[1], (this -> mlens)[2] );

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

        // if(i == 0){
        //     temp_const1[i][0] = 1.0;
        //     temp_const2[i][0] = 1.0;
        // }

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




// version 1, wrong :
    // for (int i = 0; i < NLENS; i++) {

    //     temp_const1[0] = 1.0;
    //     for (int j = 0; j < NLENS; j++) {
    //         multiply_z(temp_const1, zlens[j], j);
    //     }

    //     /* numerator, m[i] * (z-z1)*(z-z2) ... (z-zn) */
    //     for (int j = 0; j <= NLENS; j++) {
    //         p_const[j][i] = mlens[i] * temp_const1[j];
    //     }


    //     /* Sum_{j=1}^n Product (z-z_k), k=1, n, but k !=j. */
    //     for (int j = 0; j < NLENS; j++) {
    //         /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */
    //         temp_const2[0] = 1.0;
    //         degrees = 0;
    //         for (int k = 0; k < NLENS; k++) {
    //             if (k == j) continue;

    //             multiply_z(temp_const2, zlens[k], degrees);
    //             degrees++;
    //         }
    //     }
    // }

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

        // if(i == 0){
        //     temp_const1[i][0] = 1.0;
        //     temp_const2[i][0] = 1.0;
        // }

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



// double TripleLensing::TripleMag(double xsCenter, double ysCenter, double rs) {
//     // int secnum_adj = secnum;
//     // if (rs >= 1e-2) {
//     //     secnum_adj = fmin(360,  (int)(secnum_adj * rs / 1e-2) ); //if the source radius is large, we need more points on image track linking
//     // }
//     return tripleFS_v2_savehalf_quadtest(xsCenter, ysCenter, rs);
// }

void TripleLensing::tripleFSLimb2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, int secnum, int basenum, double quaderr_Tol, double relerr_mag, double mags[], int Np, double limba1, double RelTolLimb, double AbsTolLimb)
{
    // double Zlens[2 * NLENS] = { zlens[0].re, zlens[0].im, zlens[1].re, zlens[1].im, zlens[2].re, zlens[2].im };

    complex zlens[NLENS];
    for (int i = 0; i < NLENS ; i++) {
        zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
    }
    // TripleLensing TRIL(mlens, zlens);
    // TripleLensing(mlens, zlens);
    reset3(mlens, zlens);

    this->secnum = secnum;
    this->quaderr_Tol = quaderr_Tol;
    this->basenum = basenum;
    this->relerr_mag = relerr_mag;
    this->RelTolLimb = RelTolLimb;
    this->AbsTolLimb = AbsTolLimb;

    for (int i = 0; i < Np; i++) {
        // mags[i] = TripleMag(xsCenters[i], ysCenters[i], rs);

        mags[i] = TripleMagDark(xsCenters[i], ysCenters[i], rs, limba1);
        // fprintf(stderr, "i = %d, limb darkening mags[i] = %f\n", i, mags[i]);
    }

}

void TripleLensing::tripleFS2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, int secnum, int basenum, double quaderr_Tol, double relerr_mag, double mags[], int Np)
{
    // double Zlens[2 * NLENS] = { zlens[0].re, zlens[0].im, zlens[1].re, zlens[1].im, zlens[2].re, zlens[2].im };

    complex zlens[NLENS];
    for (int i = 0; i < NLENS ; i++) {
        zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
    }
    // TripleLensing TRIL(mlens, zlens);
    // TripleLensing(mlens, zlens);
    reset3(mlens, zlens);

    this->secnum = secnum;
    this->quaderr_Tol = quaderr_Tol;
    this->basenum = basenum;
    this->relerr_mag = relerr_mag;

    for (int i = 0; i < Np; i++) {
        mags[i] = TripleMag(xsCenters[i], ysCenters[i], rs);
    }

}



// double TripleLensing::tripleFS2python(double mlens[], double Zlens[], double xsCenter, double ysCenter, double rs) {
double TripleLensing::tripleFS2python(double xsCenter, double ysCenter, double rs) {
    // double Zlens[2 * NLENS] = { zlens[0].re, zlens[0].im, zlens[1].re, zlens[1].im, zlens[2].re, zlens[2].im };
    // complex zlens[NLENS];
    // for (int i = 0; i < NLENS ; i++) {
    //     zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
    // }
    // int finalnphi;
    // double quad_err, mu;
    mu = TripleMag(xsCenter, ysCenter, rs);
    // mu = tripleFS_v2_savehalf_quadtest(mlens, zlens, xsCenter, ysCenter, rs, nphi, &finalnphi, secnum, basenum, &quad_err, quaderr_Tol);
    return mu;
    // return log10(mu - 1);
}

double TripleLensing::tripleFS2python(double mlens[], double zlens[], double xsCenter, double ysCenter, double rs) {
    // double Zlens[2 * NLENS] = { zlens[0].re, zlens[0].im, zlens[1].re, zlens[1].im, zlens[2].re, zlens[2].im };
    // complex zlens[NLENS];
    // for (int i = 0; i < NLENS ; i++) {
    //     zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
    // }
    // int finalnphi;
    // double quad_err, mu;
    reset2(mlens, zlens);
    mu = TripleMag(xsCenter, ysCenter, rs);
    // mu = tripleFS_v2_savehalf_quadtest(mlens, zlens, xsCenter, ysCenter, rs, nphi, &finalnphi, secnum, basenum, &quad_err, quaderr_Tol);
    return mu;
    // return log10(mu - 1);
}



// double TripleLensing::TriplePS( double mlens[], complex zlens[], double xsCenter, double ysCenter) {
//     return triplePS(mlens, zlens, xsCenter, ysCenter, &nimages);
// }

// void TripleLensing::TripleLkvAdap(double mlens[], double zlens[], double rs, _curve *lightkv, int np0, int MAXNPS, double errTol, double t_start, double t_end, double alpha, double tE, double t0 , double u0, int np, double timerrTol) {


//     tri_ada_sample_1D(mlens, zlens, rs, nphi, secnum, basenum, quaderr_Tol, lightkv, np0, MAXNPS, errTol, t_start, t_end, alpha, tE, t0, u0, np, timerrTol);
// }

// add on 2020 Feb 04
void TripleLensing::solv_lens_equation(double zrxy[], double mlens[], double zlens[], double xs, double ys, int nlens) {
    // int DEGREE = NLENS*NLENS + 1;
    // complex zr[DEGREE];
    // // double zrxy[DEGREE*2];
    // complex coefficients[DEGREE + 1];

    // complex Zlens[NLENS];
    // int i;
    // for (i=0; i<NLENS; i++){
    //     Zlens[i] = complex(zlens[i], zlens[i+NLENS]);
    // }
    // polynomialCoefficients(mlens, Zlens, xs, ys, coefficients, NLENS, DEGREE);
    // cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

    // for (i=0; i<DEGREE;i++){
    //     zrxy[i] = zr[i].re;
    //     zrxy[i+DEGREE] = zr[i].im;
    // }
    // VBBinaryLensing VBBL;
    // VBBL.solv_lens_equation(zrxy, mlens, zlens,  xs,  ys, nlens);


    // int DEGREE = NLENS*NLENS + 1;
    complex zr[DEGREE];
    // double zrxy[DEGREE*2];
    complex coefficients[DEGREE + 1];

    complex Zlens[NLENS];
    int i;
    for (i = 0; i < NLENS; i++) {
        Zlens[i] = complex(zlens[i], zlens[i + NLENS]);
    }
    // polynomialCoefficients(mlens, Zlens, xs, ys, coefficients, NLENS, DEGREE);

    reset3(mlens, Zlens);
    polynomialCoefficients(xs, ys, coefficients);
    VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

    for (i = 0; i < DEGREE; i++) {
        zrxy[i] = zr[i].re;
        zrxy[i + DEGREE] = zr[i].im;
        // fprintf(stderr, "zr[%d] = %.15f, %.15f\n", i, zr[i].re, zr[i].im);
    }
    // return zrxy;
}


//pr:
// t0 = 7494.153;
// u0 = 0.021;
// tE = 74.62;
// s2 = 1.396;
// q2 = 0.029;
// alpha = 2.948; //rad
// s3 = 1.168;
// q3 = 3.27e-3;
// psi = 5.332; //rad
// rs = 0.22e-3;
//

void TripleLensing::TriLightCurve(double *pr, double *mags, double *y1s, double *y2s, int np) {
    // pr: t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs
    double t0 = pr[0], u0 = pr[1], tE = pr[2], s2 = pr[3], q2 = pr[4], alpha = pr[5], s3 = pr[6], q3 = pr[7], psi = pr[8], rs = pr[9];
    double salpha = sin(pr[3]), calpha = cos(pr[3]), tE_inv = 1 / tE;

    // double mlens[NLENS];
    double tn; // NLENS defined in .h file
    // complex zlens[NLENS];

    // parameter conversion, from s2, q2,, s3, q3, psi to mlens[i], zlens[i]
    double inv1andq2 = 1 / (1 + q2);
    double inv1andq2q3 = 1 / (1 + q2 + q3);
    mlens[0] = inv1andq2q3;
    mlens[1] = q2 * inv1andq2q3;
    mlens[2] = q3 * inv1andq2q3;
    zlens[0] = complex(- q2 * inv1andq2 * s2 , 0);
    zlens[1] = complex( inv1andq2 * s2 , 0);
    zlens[2] = complex(zlens[0].re + s3 * cos(psi), s3 * sin(psi));

    for (int i = 0; i < np; i++) {
        // tn = (ts[i] - t0) * tE_inv;
        // // Position of the source center
        // y1s[i] = u0 * salpha + tn * calpha;
        // y2s[i] = u0 * calpha - tn * salpha;
        mags[i] = TripleMag(y1s[i], y2s[i], rs);
    }
}


//output in two files the triple critical curves and caustics
void TripleLensing::outputCriticalTriple_list(double allxys[], double mlens[], double zlens[], int nlens, int NPS)
// void outputCriticalTriple_list(double allxys[], double mlens[], double zlens[], int NLENS, int NPS)
{
    VBBinaryLensing VBBL;
    VBBL.outputCriticalTriple_list(allxys, mlens, zlens, nlens, NPS);
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



// tripleFS_v2_savehalf_quadtest(mlens, zlens, xsCenter, ysCenter, rs, nphi, &finalnphi, secnum, basenum, &quad_err, quaderr_Tol);
double TripleLensing::TripleMag(double xsCenter, double ysCenter, double rs) {
    // VBBinaryLensing VBBL;

    // double quaderr = 0.0;//quaderr_Tol = 1e-2
    // muPS = tripleQuatrapoleTest( mlens, zlens, xsCenter, ysCenter, &nimages, rs, &quad_err );
    finalNPS = 0;
    rho2 = rs * rs;
    areaSource = M_PI * rho2;
    ftime = 1;
    relerr_priv = this->relerr_mag;
    secnum_priv = this->secnum;
    quad_err = this->quaderr_Tol;
    basenum = this->basenum;
    // fprintf(stderr, "in TripleMag, mlens[0 1 2] %f %f %f\n", (this->mlens)[0], (this->mlens)[1], (this->mlens)[2] );
    muPS = tripleQuatrapoleTest(xsCenter, ysCenter, rs);
    // fprintf(stderr, "in TripleMag, mlens[0 1 2] %f %f %f\n", (this->mlens)[0], (this->mlens)[1], (this->mlens)[2] );
    // fprintf(stderr, "quad_err %f, muPS = %f\n", quad_err, muPS);
    // *quad_err = quaderr;
    // if ( CQ * quad_err <= quaderr_Tol && muPS < 1.5) {
    if ( CQ * quad_err <= quaderr_Tol) {
#ifdef VERBOSE
        fprintf(stderr, "quad_err %f,using point source magnification, muPS = %f\n", quad_err, muPS);
#endif
        ifFinite = 0;
        return muPS;
    } else {
        ifFinite = 1;
#ifdef VERBOSE
        fprintf(stderr, "quad_err %f,using caustics crossing computation\n", quad_err);
#endif
        mu0 = muPS;
        // secnum_priv = fmin((int)(secnum_priv * fmax( fmax(1, quad_err), mu0 ) ), 90);
        secnum_priv = fmax(fmin((int)(2 * secnum_priv * ( abs(log10(quad_err * 1e4)) * mu0 * mu0) / abs(log10(rs))) , 45), secnum_priv);
        relerr_priv = fmax( relerr_priv , relerr_priv / fmax( quad_err, mu0 * mu0 ) * abs(log10(rs))  );


#ifdef VERBOSE
        printf("point magnification: %f, relerr_priv = %.3e\n", mu0, relerr_priv);
#endif

        _sols *imageTracks  = new _sols;
        _sols *prevstore  = new _sols;
        _linkedarray *phis = new _linkedarray;

        // phis = getphis_v3( mlens,  zlens,  xsCenter,  ysCenter,  rs,  nphi, secnum,  basenum,  &distype, &maxmuidx);
        phis = getphis_v3(  xsCenter,  ysCenter,  rs);
        // if (distype == 2) {relerr_priv /= 10.0;};
        // relerr_priv /= 10.0;

        // mu0 = 1;
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
                fprintf(stderr, "\n\nin tripleFS, i= %d, imageTracks = Null, nphi = %d, xsCenter=%f, ysCenter = %f\n\n", i, finalnphi, xsCenter, ysCenter);
#endif
                if (distype == 2) {
                    delete phis;
                    phis = new _linkedarray;
                    maxmuidx = (int)(absint(maxmuidx + 1)) % secnum;
                    basenum *= 2;
                    // phis = getphisNonfirst_caustic_crossing(mlens,  zlens,  xsCenter,  ysCenter,  rs,  nphi, secnum,  basenum, maxmuidx);
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
            area = areaFunc(imageTracks, rs);
#endif

            mu = area / areaSource;
#ifdef VERBOSE
            fprintf(stderr, "in tripleFS, i= %d, mu0= %f, mu= %f, nphi = %d, xsCenter=%f, ysCenter = %f, nimages = %d\n", i, mu0, mu, finalnphi, xsCenter, ysCenter, imageTracks->length);
#endif
#ifdef verbose
            fprintf(stderr, "in tripleFS, i= %d, mu0= %f, mu= %f, nphi = %d, xsCenter=%f, ysCenter = %f, nimages = %d\n", i, mu0, mu, *finalnphi, xsCenter, ysCenter, imageTracks->length);
#endif
            if (abs(mu - mu0) / mu < relerr_priv) break;
            mu0 = mu;
            phis->extend();
        }
#ifdef VERBOSE
        printf("mu = %f\n", mu);
        saveTracks(imageTracks);
#endif

        // show image track animation here:

        finalNPS = finalnphi;
        nimages = imageTracks->length; // add on 2020.10.23
        delete imageTracks;
        delete phis;
        delete prevstore;
        // fprintf(stderr, "\t\t\ttripleFS mu = %f, rs = %.25f\n", mu, rs);
        return (mu);
    }

}


//
//   implement Gould (2008) approximation, perhaps we can use static tables for phi to get rid of cos and sin
//
double TripleLensing::gould(double xsCenter, double ysCenter, double rs, double Gamma)
{
    // VBBinaryLensing VBBL;
    // TripleLensing TRIL;
    // TRIL.reset3(mlens, zlens);


    int i;
    double phi, dphi;
    double muRho[8];
    double muRhoHalf[4];
    double AhalfPlus, APlus, Across, A0;
    double A2rho2, A4rho4;// gould, 10%,
    double xs, ys;
    int nimages;
    double mag;

    // A0 = triplePS(mlens, zlens, xsCenter, ysCenter, &nimages);
    A0 = TriplePS(xsCenter, ysCenter);

    // 8 points at full radius
    dphi = M_PI * 0.25;
    for (i = 0; i < 8; i++) {
        phi = dphi * (i - 1);
        xs = xsCenter + rs * cos(phi);
        ys = ysCenter + rs * sin(phi);

        // muRho[i] = triplePS(mlens, zlens, xs, ys, &nimages);
        muRho[i] = TriplePS(xs, ys);
    }


    // four points at half radius
    dphi = M_PI * 0.5;
    for (i = 0; i < 4; i++) {
        phi = dphi * (i - 1);
        xs = xsCenter + 0.5 * rs * cos(phi);
        ys = ysCenter + 0.5 * rs * sin(phi);

        // muRhoHalf[i] = triplePS(mlens, zlens, xs, ys, &nimages);
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


double TripleLensing::TripleMagDark(double xsCenter, double ysCenter, double RSv, double a1) {
    double Mag = -1.0, Magold = 0., Tolv = this->AbsTolLimb;
    double tc, lb, rb, lc, rc, cb = 0, cc, r2, cr2, scr2;
    int c = 0, flag;
    double currerr, maxerr;
    annulus *first, *scan, *scan2;
    int nannuli, nannold, totNPS = 1, minannuli;
    // _sols *Images;

    double y_1 = xsCenter;
    double y_2 = ysCenter;
    while ((Mag < 0.9) && (c < 3)) {

        first = new annulus;
        first->bin = 0.;
        first->cum = 0.;// cumulative function F(r)
        first->Mag = TriplePS(y_1, y_2);
        // fprintf(stderr, "in TripleMagDark, TriplePS... = %f \n\n", first->Mag);
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
        // scan->Mag = BinaryMag(a, q, y_1, y_2, RSv, Tolv, &Images);
        scan->Mag = TripleMagFull(xsCenter, ysCenter, RSv);

        // fprintf(stderr, "TripleFS ??? = %f \n\n", scan->Mag);

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
        //          scan->err+=scan->Mag*Tolv*0.25; //Impose calculation of intermediate annulus at mag>4. Why?
        currerr = scan->err;
        flag = 0;
        // fprintf(stderr, "currerr %.5e Tolv %.5e\n", currerr, Tolv);
        nannuli = nannold = 1;
        while (((flag < nannold + 5) && (currerr > Tolv) && (currerr > this->RelTolLimb * Mag)) || (nannuli < minannuli)) {
            // update/find the annuli with the maximum error
            maxerr = 0;
            for (scan2 = first->next; scan2; scan2 = scan2->next) {
// #ifdef _PRINT_ERRORS_DARK
                // printf("\nnannuli=%d scan2->Mag=%le scan2->err=%le | Mag=%le currerr=%le, cb = %lf", nannuli, scan2->Mag, scan2->err, Mag, currerr, cb);
// #endif
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
                // cc = (3 * r2 * (1 - a1) - 2 * a1 * ((1-scr2) * cr2 - 1)) / (3 - a1);
                // cc = (3 * r2 * (- a1) - 2 * a1 * ((1-scr2) * cr2 - 1)) ;
                // -3 a1 r2 -2 a1 ( cr2 - 1)

                // cc = ( - a1 * (r2 - 2 * scr2*cr2)) ;
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

            scan->prev->Mag = TripleMagFull(xsCenter, ysCenter, RSv * cb);


            // fprintf(stderr, "nannuli = %d, Mag(cb) = %lf, total Mag = %lf, cb = %.20f, currerr = %.5e\n", nannuli, scan->prev->Mag, Mag, cb, currerr);

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
            // printf("\n%d", Images->length);
            printf("\n%d", nimages);
#endif
            // delete Images;

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

        // 这里 配合 BinaryMagMultiDark(), 当有多个 limb darkening，想计算不同 a 系数下的 magnification 时，第一次multidark = true，其他都是 false
        // if (multidark) {
        //     annlist = first;
        // } else {
        //     while (first) {
        //         scan = first->next;
        //         delete first;
        //         first = scan;
        //     }
        // }

        Tolv /= 10;
        c++;
    }
    NPS = totNPS;
    therr = currerr;
    // fprintf(stderr, "\t\t NPS = %d in VBBL\n", NPS);
    finalNPS = NPS;
    return Mag;
}



double TripleLensing::TripleMagFull(double xsCenter, double ysCenter, double rs) {
    // VBBinaryLensing VBBL;

    // double quaderr = 0.0;//quaderr_Tol = 1e-2
    // muPS = tripleQuatrapoleTest( mlens, zlens, xsCenter, ysCenter, &nimages, rs, &quad_err );
    finalNPS = 0;
    rho2 = rs * rs;
    areaSource = M_PI * rho2;
    ftime = 1;
    relerr_priv = relerr_mag;
    secnum_priv = secnum;
    // fprintf(stderr, "in TripleMag, mlens[0 1 2] %f %f %f\n", (this->mlens)[0], (this->mlens)[1], (this->mlens)[2] );
    // muPS = tripleQuatrapoleTest(xsCenter, ysCenter, rs);

    muPS = TriplePS(xsCenter, ysCenter);

    // fprintf(stderr, "in TripleMag, mlens[0 1 2] %f %f %f\n", (this->mlens)[0], (this->mlens)[1], (this->mlens)[2] );
    // fprintf(stderr, "quad_err %f, muPS = %f\n", quad_err, muPS);
    // *quad_err = quaderr;

    // if ( CQ * quad_err <= quaderr_Tol) {
    if ( 1 == 0 ) {
#ifdef VERBOSE
        fprintf(stderr, "quad_err %f,using point source magnification, muPS = %f\n", quad_err, muPS);
#endif
        ifFinite = 0;
        return muPS;
    } else {
        ifFinite = 1;
// #ifdef VERBOSE
//         fprintf(stderr, "quad_err %f,using caustics crossing computation\n", quad_err);
// #endif
        mu0 = muPS;
        // secnum_priv = fmin((int)(secnum_priv * fmax( fmax(1, quad_err), mu0 ) ), 90);
        // secnum_priv = fmax(fmin((int)(2 * secnum_priv * ( abs(log10(quad_err * 1e4)) * mu0 * mu0) / abs(log10(rs))) , 45), secnum_priv);
        // relerr_priv = fmax( relerr_priv , relerr_priv / fmax( quad_err, mu0 * mu0 ) * abs(log10(rs))  );


#ifdef VERBOSE
        printf("point magnification: %f, relerr_priv = %.3e\n", mu0, relerr_priv);
#endif

        _sols *imageTracks  = new _sols;
        _sols *prevstore  = new _sols;
        _linkedarray *phis = new _linkedarray;

        // phis = getphis_v3( mlens,  zlens,  xsCenter,  ysCenter,  rs,  nphi, secnum,  basenum,  &distype, &maxmuidx);
        phis = getphis_v3(xsCenter,  ysCenter,  rs);
        // if (distype == 2) {relerr_priv /= 10.0;};
        // relerr_priv /= 10.0;

        // mu0 = 1;
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
                    // phis = getphisNonfirst_caustic_crossing(mlens,  zlens,  xsCenter,  ysCenter,  rs,  nphi, secnum,  basenum, maxmuidx);
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
            area = areaFunc(imageTracks, rs);
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
        saveTracks(imageTracks);
#endif

        // show image track animation here:

        finalNPS = finalnphi;
        nimages = imageTracks->length; // add on 2020.10.23
        delete imageTracks;
        delete phis;
        delete prevstore;
        // fprintf(stderr, "before return in TripleMagFull mu = %f\n", mu);

        // fprintf(stderr, "\t\tin TripleMagFull, mu0= %f, mu= %f, nphi = %d, xsCenter=%f, ysCenter = %f, nimages = %d, rs = %f\n", mu0, mu, finalnphi, xsCenter, ysCenter, nimages, rs);

        return (mu);
    }

}
    

// tripleQuatrapoleTest Bozza 2018 equation 39 with c_Q = 1, add on
double TripleLensing::tripleQuatrapoleTest(double xs, double ys, double rs) {

    // polynomialCoefficients(mlens, zlens, xs, ys, coefficients, NLENS, DEGREE);

    // fprintf(stderr, "\t\t\t compare coefficients: BinaryCoefs and polynomialCoefficients \n");
//     // BinaryCoefs(xs, ys, coefficients);
//     clock_t tim0, tim1;
//     tim0 = clock();
//     for(int ti = 0; ti<100; ti++){
//     TripleCoefs(xs, ys, coefficients);
// }
// tim1 = clock();
// fprintf(stderr, "total time on TripleCoefs= %f s\n", (tim1 - tim0) / 1e6);

//     for(int i = 0; i<DEGREE; i++){
//         fprintf(stderr, "%f %f ; ", coefficients[i].re, coefficients[i].im);
//     }
//     fprintf(stderr, "\n");

//     tim0 = clock();
//     for(int ti = 0; ti<100; ti++){
//     polynomialCoefficients(xs, ys, coefficients);
// }
//     tim1 = clock();
// fprintf(stderr, "total time on polynomialCoefficients = %f s\n", (tim1 - tim0) / 1e6);
    // fprintf(stderr, "in tripleQuatrapoleTest, mlens[0 1 2] %f %f %f\n", mlens[0], mlens[1], mlens[2] );
    polynomialCoefficients(xs, ys, coefficients);

    // #ifdef VERBOSE
    // for(int i = 0; i<DEGREE; i++){
    //     fprintf(stderr, "%f %f ; ", coefficients[i].re, coefficients[i].im);
    // }
    // fprintf(stderr, "\n");
    // #endif

    VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

    muTotal = 0.0;
    eacherr = 0.0;
    // rho2 = rs * rs;

    nimages = 0;
    quad_err = 0.0;//quad_err
    for (int i = 0; i < DEGREE; i++) {
        // flag = trueSolution(mlens, zlens, xs, ys, zr[i], &mu, &lambda1, &lambda2, &thetaJ, NLENS, &J1, &J2, &dJ, &J3);
        flag = trueSolution_qtest(xs, ys, zr[i]);
        // fprintf(stderr, "flag = %d\n", flag);
        nimages += flag;
        if (flag) {
            muTotal += abs(mu);
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
            quad_err += eacherr / abs(dJ2 * dJ2 * dJ.re);

            // fprintf(stderr, "+= eacherr %f * rho2 %f / abs(dJ2 * dJ2 * dJ.re) dJ2 %f dJ.re %f; \n", eacherr , rho2 , dJ2, dJ.re );
        }
    }
    quad_err *= rho2;
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

    // add on 2020.04.10
    final->last->ds = p->ds;
    // add on 2020.04.14
    final->last->dz = p->dz;
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

        final->last->ds = p->ds;
        final->last->dz = p->dz;
        final->last->closepairparity = p->closepairparity;
    }

    // final = imageTracks->first;
    // imageTracks->first = imageTracks->first->next;
    // imageTracks->first->prev = 0;
    // imageTracks->length --;


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



//magnification for a triple lens point source
// double TripleLensing::triplePS(double mlens[], complex zlens[], double xs, double ys, int *nimages) {
double TripleLensing::TriplePS(double xs, double ys) {
    // VBBinaryLensing VBBL;
    // complex zr[DEGREE];
    // complex coefficients[DEGREE + 1];
    // double lambda1, lambda2, thetaJ;

    // polynomialCoefficients(mlens, zlens, xs, ys, coefficients, NLENS, DEGREE);
    polynomialCoefficients(xs, ys, coefficients);
    VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

    // double mu, muTotal = 0.0;
    muTotal = 0.0;
    // complex J1, J2, dJ, J3;
    nimages = 0;
    for (int i = 0; i < DEGREE; i++) {
        flag = trueSolution(xs, ys, zr[i]);
        nimages += flag;
        muTotal += ( flag ? abs(mu) : 0  );
        // if (flag) {
        //     muTotal += abs(mu);
        //     // printf("%f %f %f %f %f %f\n", real(zr[i]), imag(zr[i]), lambda1, lambda2, thetaJ, mu);
        // }
    }
    return (muTotal);
}

// double tripleFS2python(double mlens[], double Zlens[], double xsCenter, double ysCenter, double rs, int nphi, int secnum, int basenum, double quaderr_Tol) {
//     complex zlens[NLENS];
//     for (int i = 0; i < NLENS ; i++) {
//         zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
//     }
//     int finalnphi;
//     double quad_err, mu;
//     mu = tripleFS_v2_savehalf_quadtest(mlens, zlens, xsCenter, ysCenter, rs, nphi, &finalnphi, secnum, basenum, &quad_err, quaderr_Tol);
//     return log10(mu - 1);
// }






// out put pure image points around the edge of the source
// void TripleLensing::outImgPoints(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi)
void TripleLensing::outImgPoints(double xsCenter, double ysCenter, double rs, int nphi)
{
    // check if we need to use FS, if not, return 0; if yes, return 1;
    // (return the average mag on the boundary of edge of the source, this is not eligant)
    // VBBinaryLensing VBBL;
    fprintf(stderr, "Generating image positions corresponding to source boundary ...");
    /* now find the closed curves for all real and imaginary roots for a circular source */
    // double xs, ys,
    // double errormu = 0.2, relerrmu = 0.2;
    // double phi, dphi;
    // complex zr[DEGREE];
    // complex dzs;
    // complex coefficients[DEGREE + 1];

    // int nsolution = 0, FSflag = 0;
    nsolution = 0;
    FSflag = 0;
    int i, j, k;
    // double mu, lambda1, lambda2, thetaJ;//Sigma_mui = 0, absSigma_mui[nphi], muCenter
    // complex J1, J2, dJ, J3;

    // int nimages;          // number of true images
    // int flag;

    //double phi0 = 0.0;        // initial angle on the circle, occasionally it makes a difference

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
        flag = trueSolution(xs, ys, zr[i]);
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
        // Sigma_mui = 0;
        // absSigma_mui[j] = 0;
        for (int i = 0; i < DEGREE; i++) {
            // flag = trueSolution(mlens, zlens, xs, ys, zr[i], &mu, &lambda1, &lambda2, &thetaJ, NLENS, &J1, &J2, &dJ, &J3);
            flag = trueSolution( xs, ys, zr[i]);
            // printf("i,flag, mu : %d %d %f\n", i, flag, mu);
            nimages += flag;
            if (flag == 1) {        // true solution
                // Prov->last->flag = 1 ;
                // Sigma_mui += mu;
                // absSigma_mui[j] += absf(mu);
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



// void saveTracks(_sols * track)
// {
//     int i = 0;
//     int itrack = 0;
//     FILE *fp;

//     fp = fopen("data/allTracks.dat", "w");
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


//output in two files the triple critical curves and caustics 2020.09.08
void outputCaustics_for_crossSection(double mlens[], complex zlens[], int nlens, int npnt)
{
    // fprintf(stderr, "outputCaustics_for_crossSection, Generating caustics for cross section ...");
    // VBBinaryLensing VBBL;
    // int ncurves = 0;
    // int ncritical = 0;
    _sols *track;
    // FILE *fcritical, *fcaustics;
    // criticalCurves = VBBL.PlotCritTriple(mlens, zlens, npnt, nlens);
    track = VBBL.outCritTriple(mlens, zlens, npnt, nlens);
    // fcritical = fopen("./data/_critical_curves.dat", "w");
    // fcaustics = fopen("./data/_caustics.dat", "w");

    // check how many closed critical curves we have
    // first halfs are critical curves
    // ncritical = criticalCurves->length / 2; // number of closed critical curves
// #ifdef VERBOSE
//         printf("I am in outputCriticalTriple, Number of closed critical curves: %d\n", ncritical);
// #endif

    // write out the critical curves and caustics separately

    // ncurves = 0;
    // for (_curve *c = criticalCurves->first; c; c = c->next) {
    //         for (_point *p = c->first; p; p = p->next) {
    //             fprintf(fcritical, "%.10lf %.10lf\n", p->x1, p->x2);
    //             fprintf(fcaustics, "%.10lf %.10lf\n", p->ds, p->dJ);
    //             // npoints++;
    //         }
    // }
    // fclose(fcritical);  fclose(fcaustics);
    // fprintf(stderr, " ... done.\n");

    FILE *fp;

    fp = fopen("data/_allTracks.dat", "w");
    // fp = fopen(filename, "w");
    // printf("entered saveTracks\n\n\n");

    for (_curve *c = track->first; c; c = c->next) { //curves
        fprintf(fp, "%d ", c->length);
    }
    fprintf(fp, "\n");

    for (_curve *c = track->first; c; c = c->next) { //curves
        for (_point *p = c->first; p; p = p->next) { //point
            // fprintf(fp, "%f %f %f %f ", p->x1, p->x2, p->phi, p->mu);
            fprintf(fp, "%f %f %f %f ", p->ds, p->dJ, p->x1, p->x2);// caustics x, caustics y, critical x, critical y
            // i++;
        }
        // itrack++;
    }
    fclose(fp);


}




#define _c_parity \
if (c->length>1){\
            if (c->first->mu > 0 && c->first->next->mu > 0) {\
                c->parity = 1;\
            }\
            else if(c->first->mu < 0 && c->first->next->mu < 0) {\
                c->parity = -1;\
            }else if (c->first->mu > 0 && c->first->next->mu < 0 && c->first->next->next->mu < 0){\
                c->parity = -1;\
            }else if (c->first->mu < 0 && c->first->next->mu > 0 && c->first->next->next->mu > 0){\
                c->parity = 1;\
            }\
}else{\
    c->parity = (c->first->mu > 0) ? 1:-1;\
}

#define _final_parity \
// if (final->parity == 0) {\
//     if (final->length>1){\
//                 if (final->first->mu > 0 && final->first->next->mu > 0) {\
//                     final->parity = 1;\
//                 }\
//                 else if(final->first->mu < 0 && final->first->next->mu < 0) {\
//                     final->parity = -1;\
//                 }else if (final->first->mu > 0 && final->first->next->mu < 0 && final->first->next->next->mu < 0){\
//                     final->parity = -1;\
//                 }else if (final->first->mu < 0 && final->first->next->mu > 0 && final->first->next->next->mu > 0){\
//                     final->parity = 1;\
//                 }\
//     }else{\
//         final->parity = (final->first->mu > 0) ? 1:-1;\
//     }\
// }

// #define _final_parity \
// if (final->parity == 0) {\
// final->parity = (final->first->mu > 0) ? 1:-1;\
// }

#define _Prov_parity \
// if(Prov->length>1){\
// if (Prov->first->mu > 0 && Prov->first->next->mu > 0) {\
//     Prov->parity = 1;\
// }\
// else if (Prov->first->mu < 0 && Prov->first->next->mu < 0) {\
//     Prov->parity = -1;\
// }else if (Prov->first->mu > 0 && Prov->first->next->mu < 0 && Prov->first->next->next->mu < 0) {\
//     Prov->parity = -1;\
// }else if (Prov->first->mu < 0 && Prov->first->next->mu > 0 && Prov->first->next->next->mu > 0) {\
//     Prov->parity = 1;\
// }\
// }else{\
//     Prov->parity = (Prov->first->mu > 0) ? 1:-1;\
// }

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




// _sols *outputTracks_v2_savehalf(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, _linkedarray *PHI, _sols **prevstore, int ftime, unsigned short int *ifjump) {
_sols *TripleLensing::outputTracks_v2_savehalf(double xsCenter, double ysCenter, double rs, _linkedarray *PHI, _sols **prevstore) {
    // unsigned short int min_mass_idx = 0;
    // double tempmass = mlens[0];
    // for (int i = 1; i < NLENS; i++) {
    //     if ( mlens[i] < min_mass_idx ) {
    //         min_mass_idx = i;
    //     }
    // }

    // VBBinaryLensing VBBL;

    // now find the closed curves for all real and imaginary roots for a circular source
    // double xs, ys;
    // double phi, dphi;
    // complex zr[DEGREE];
    // complex dzs, J1, J2, dJ, dy, dz, J3;
    // double ds;
    // complex coefficients[DEGREE + 1];
    nsolution = 0;
    int i, j = 0, k;
    // double mu, lambda1, lambda2, thetaJ;

    _sols *allSolutions, *imageTracks;
    _curve *Prov = new _curve;
    _curve *Prov2 = new _curve;
    _curve *tempProv = new _curve;
    // Prov2 = NULL;
    // _curve *Prov;
    _point *pisso;
    _point *scanpisso;
    // double SD, MD, CD;

    // connect all images, including fake images, initialize ten _curve object to save solutions from polynomial solving
    allSolutions = new _sols;
    for (int i = 0; i < DEGREE; i++) {
        Prov = new _curve;
        allSolutions->append(Prov);
    }

    // int nimages;      // number of true images
    // int flag, nphi = PHI->length;
    pntnum = PHI->length;

    if (ftime) {
        *prevstore = NULL;
        // currstore = NULL;
        // delete *prevstore;
        // delete currstore;
        *prevstore = new _sols;
        Node *scan, *tempscan;
        scan = PHI->first;


        phi = scan->value;
        xs = xsCenter + rs * cos(phi);
        ys = ysCenter + rs * sin(phi);

        // polynomialCoefficients(mlens, zlens, xs, ys, coefficients, NLENS, DEGREE);
        polynomialCoefficients(xs, ys, coefficients);
        VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);
        // VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, false, true);

        // first point in all solution tracks
        tempProv = new _curve;
        (*prevstore)->append(tempProv);
        Prov = allSolutions->first;

        nimages = 0;
        for (int i = 0; i < DEGREE; i++) {
            // flag = trueSolution(mlens, zlens, xs, ys, zr[i], &mu, &lambda1, &lambda2, &thetaJ, NLENS, &J1, &J2, &dJ, &J3);
            flag = trueSolution(xs, ys, zr[i]);
            nimages += flag;

            Prov->append(zr[i].re, zr[i].im);
            Prov->last->phi = phi;
            Prov->last->mu = mu;
            Prov->last->zs = complex(xs, ys);
            Prov->last->thetaJ = thetaJ;

            (*prevstore)->last->append(zr[i].re, zr[i].im);
            (*prevstore)->last->last->phi = phi;
            (*prevstore)->last->last->mu = mu;
            (*prevstore)->last->last->zs = complex(xs, ys);
            (*prevstore)->last->last->thetaJ = thetaJ;

            if (flag == 1) {      // true solution
                Prov->last->flag = 1 ;
                (*prevstore)->last->last->flag = 1 ;

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

            Prov = Prov->next;
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
            Prov2 = new _curve;
            tempProv = new _curve;
            (*prevstore)->append(tempProv);
            nimages = 0;
            for (int i = 0; i < DEGREE; i++) {
                // flag = trueSolution(mlens, zlens, xs, ys, zr[i], &mu, &lambda1, &lambda2, &thetaJ, NLENS, &J1, &J2, &dJ, &J3);
                flag = trueSolution(xs, ys, zr[i]);
                nimages += flag;

                Prov2->append(zr[i].re, zr[i].im);
                Prov2->last->phi = phi;
                Prov2->last->mu = mu;
                Prov2->last->zs = complex(xs, ys);
                Prov2->last->thetaJ = thetaJ;

                (*prevstore)->last->append(zr[i].re, zr[i].im);
                (*prevstore)->last->last->phi = phi;
                (*prevstore)->last->last->mu = mu;
                (*prevstore)->last->last->zs = complex(xs, ys);
                (*prevstore)->last->last->thetaJ = thetaJ;


                if (flag == 1) {      // true solution
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

            scan->nimages = nimages;
            scan = scan->next;

            // attach the current solutions (both true and false) positions to last ones with the closest distance
            for (Prov = allSolutions->first; Prov;  Prov = Prov->next) {
                Prov2->closest(Prov->last, &pisso);
                Prov2->drop(pisso);
                Prov->append(pisso);
            }
        }
#ifdef VERBOSE
        fprintf(stderr, "*prevstore->length %d\n", (*prevstore)->length);
#endif

    } else { //firsttime = 0
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
        Prov = new _curve;
        Prov = allSolutions->first;
        nimages = 0;
        scanpisso = scanstorecurve->first;

        for (int i = 0; i < DEGREE; i++) {
            flag = (int)((scanpisso->flag + 1) / 2);
            nimages += flag;

            Prov->append(scanpisso->x1, scanpisso->x2);
            Prov->last->phi = scanpisso->phi;
            Prov->last->mu = scanpisso->mu;
            Prov->last->zs = scanpisso->zs; //complex(xs, ys);
            Prov->last->thetaJ = scanpisso->thetaJ;
            Prov->last->flag = scanpisso->flag;

            Prov->last->ds = scanpisso->ds;
            Prov->last->dz = scanpisso->dz;

            tempcurrstore->last->append(scanpisso->x1, scanpisso->x2);
            tempcurrstore->last->last->phi = scanpisso->phi;
            tempcurrstore->last->last->mu = scanpisso->mu;
            tempcurrstore->last->last->zs = scanpisso->zs; //complex(xs, ys);
            tempcurrstore->last->last->thetaJ = scanpisso->thetaJ;
            tempcurrstore->last->last->flag = scanpisso->flag;

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
                nimages = 0;
                for (int i = 0; i < DEGREE; i++) {
                    // flag = trueSolution(mlens, zlens, xs, ys, zr[i], &mu, &lambda1, &lambda2, &thetaJ, NLENS, &J1, &J2, &dJ, &J3);
                    flag = trueSolution(xs, ys, zr[i]);
                    nimages += flag;
                    Prov2->append(zr[i].re, zr[i].im);
                    Prov2->last->phi = phi;

                    Prov2->last->mu = mu;
                    Prov2->last->zs = complex(xs, ys);
                    Prov2->last->thetaJ = thetaJ;

                    tempcurrstore->last->append(zr[i].re, zr[i].im);
                    tempcurrstore->last->last->phi = phi;
                    tempcurrstore->last->last->mu = mu;
                    tempcurrstore->last->last->zs = complex(xs, ys);
                    tempcurrstore->last->last->thetaJ = thetaJ;

                    if (flag == 1) {      // true solution
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
                scan->nimages = nimages;
                scan = scan->next;

                // attach the current solutions positions to last ones with the closest distance
                for (Prov = allSolutions->first; Prov;  Prov = Prov->next) {
                    Prov2->closest(Prov->last, &pisso);
                    Prov2->drop(pisso);
                    Prov->append(pisso);
                }
            } else {
                delete Prov2;
                Prov2 = new _curve;
                tempProv = new _curve;
                tempcurrstore->append(tempProv);
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

                    Prov2->last->ds = scanpisso->ds;
                    Prov2->last->dz = scanpisso->dz;

                    tempcurrstore->last->append(scanpisso->x1, scanpisso->x2);
                    tempcurrstore->last->last->phi = scanpisso->phi;
                    tempcurrstore->last->last->mu = scanpisso->mu;
                    tempcurrstore->last->last->zs = scanpisso->zs; //complex(xs, ys);
                    tempcurrstore->last->last->thetaJ = scanpisso->thetaJ;
                    tempcurrstore->last->last->flag = scanpisso->flag;

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
                }
            }

        }
        delete (*prevstore);//200711
        (*prevstore) = new _sols;;//200711

        (*prevstore) = tempcurrstore;
        tempcurrstore = NULL;
        delete tempcurrstore;
    }

    //  printAllTracks(allSolutions);

    //allocate memory for potentially true image tracks, we make a factor four more segments just in case
    imageTracks = new _sols;
    for (int i = 0; i < (NLENS + 1) * DEGREE; i++) {
        Prov = new _curve;
        imageTracks->append(Prov);
    }


    //
    // Solution tracks may contain only true images, some only false solutions, while some may be mixed
    //
    complex z, zs;
    int mixedTracks = 0; // number of mixed tracks
    _curve *tempcurve;
    Prov2 = imageTracks->first;

    int trueImages = 0, falseImages = 0;
    int firstTime = 1;
    int disfirstTime = 1;
    int previousImage = 0;
    int Parity;  // within the same curve, we require the parity to be the same
    double ratio;
    double muPrevious, disPrevious;//mu
    for (Prov = allSolutions->first; Prov;  Prov = Prov->next) {
        trueImages = 0;
        falseImages = 0;
        firstTime = 1;
        disfirstTime = 1;
        previousImage = 0;
        for (_point *p = Prov->first; p; p = p->next) {
            z = complex(p->x1, p->x2);
            if (p->flag == 1) { // true solutions, should return identical source radius
                trueImages++;
                if (firstTime) {  // we enter the track the first time
                    firstTime = 0;
                    muPrevious = p->mu;
                } else {
                    // parity has changed, start a new curve
                    // actually there should be a furthur more test: test that whether p and p->next is close enough, sometimes they will be very far away from each other, in this case we should break the current Prov and start a new one.
                    // now this situation is done by adding a testing afterwards, see the "imageTracks_checked" below
                    if (muPrevious * p->mu < 0.0)  {
                        if (!Prov2->next) {
                            tempcurve = new _curve;
                            imageTracks->append(tempcurve);
                        }
                        Prov2 = Prov2->next;
                        muPrevious = p->mu;
                        disfirstTime = 1;
                    }

                    // begin ************************ add on 2020.07.14
                    if (disfirstTime && Prov2->length > 1) {
                        disPrevious = *(Prov2->last) - *(Prov2->last->prev); // sqrt(dis)
                        disfirstTime = 0;
                    }

                    if (Prov2->length > 1 && ((*(p) - * (Prov2->last)) >=  1e8 * disPrevious) ) {
                        if (!Prov2->next) {
                            tempcurve = new _curve;
                            imageTracks->append(tempcurve);
                        }
                        Prov2 = Prov2->next;
                        muPrevious = p->mu;
                        disfirstTime = 1;
                    }
                    // end ************************ add on 2020.07.14


                }

                addPoint(Prov2, p); // add point p to the track

                previousImage = TRUE_IMAGE;
            } else {

                falseImages++;

                //first time that we enter this
                if (firstTime) {
                    previousImage = FALSE_IMAGE;
                    firstTime = 0;
                }

                // already a true segment existing, we cut the segment, and start a new one
                if (previousImage == TRUE_IMAGE) { // previously we had true images, now we are encountering on false solutions
                    if (!Prov2->next) {
                        tempcurve = new _curve;
                        imageTracks->append(tempcurve);
                        disfirstTime = 1;
                    }

                    Prov2 = Prov2->next;  // we start a new image track

                    previousImage = FALSE_IMAGE;
                }
            }
        }

        if (trueImages) {
            if (!Prov2->next) {
                tempcurve = new _curve;
                imageTracks->append(tempcurve);
                disfirstTime = 1;
            }
            Prov2 = Prov2->next;
        }

        if (trueImages && falseImages) mixedTracks += 1; // there are tracks that are mixed
    }


    tempcurve = NULL;
    delete tempcurve;


    Prov = imageTracks->first;
    while (Prov) {
        Prov2 = Prov->next;
        if (!Prov->first || !Prov->last) { // add one more condition on 2020.03.19,
            // if (!Prov->first || !Prov->last || Prov->length < 2) { // need to remove on 2020.07.13
            // if (!Prov->first || !Prov->last || Prov->length <= fmax((int)(1e-2 * nphi), 10) ) {

            imageTracks->drop(Prov);
            delete Prov;
            Prov = Prov2;
        } else {Prov = Prov->next;}

    }
    // 到目前为止，得到了segment

#ifdef VERBOSE
    saveTracks_before(imageTracks);
#endif

#ifdef VERBOSE
    fprintf(stderr, "\nimageTracks->length %d before linking to closed continuous track\n", imageTracks->length);//for debug
#endif

    // show image track animation here:
    // saveTracks(imageTracks);

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
    // for (int i = 0; i < imageTracks->length; i++) {
    //     Prov = new _curve;
    //     connected->append(Prov);
    // }


    //
    //now attach self-closed, non-caustic crossing image tracks
    //
    int selfClosed = 0;
    Prov = imageTracks->first;
    // Prov2 = connected->first;
    // _curve *tempProv;
    tempProv = NULL;
    while (Prov) {
        SD = *(Prov->first) - *(Prov->last);
        if (SD < EPS * EPS && Prov->length > 3) { // closed tracks

            int pointInTrack = 0; // how many points in track

            // for (_point *p = Prov->first; p; p = p->next) {
            //     addPoint(Prov2, p); // add point p to the track
            //     pointInTrack++;
            // }

            tempProv = Prov->next;
            imageTracks->drop(Prov);
            // tempfinal = Prov;
            connected->append(Prov);



            selfClosed++;
#ifdef VERBOSE
            printf("entered, pointInTrack=%d imageTracks->length %d\n", Prov->length, imageTracks->length);
#endif
            //
            // such image tracks must have number of points == number of points on the source circle
            //
            // special case may be there is only one single point around a narrow cusp, / we probably
            // need to increase the number of points to solve the problem
            //
            if (Prov->length != pntnum) {
#ifdef VERBOSE
                fprintf(stderr, "the number of points is %d, should be equal to pntnum=%d!\n", Prov->length, pntnum);
#endif
            }

            // Prov2 = Prov2->next;  // go the next closed curve
            // drop the already connected curve from imageTracks
#ifdef VERBOSE
            fprintf(stderr, "imageTracks->length after drop of a self-closed track %d\n", imageTracks->length);
#endif
            // delete Prov;
            Prov = tempProv;

        }else{
            Prov = Prov->next;
        }

        // 2020.11.03 注释掉
//         else {
//             if (Prov->length == pntnum) {
//                 tempProv = Prov->next;
//                 dzs = Prov->first->zs - Prov->last->zs;
// #ifdef VERBOSE
//                 printf("first and last source position difference at place 00: %e\n", abs(dzs));
// #endif
//                 if (abs(dzs) < EPS) {
// #ifdef VERBOSE
//                     fprintf(stderr, "\t >>>connected->append(final) place 00, Prov->length = %d\n", Prov->length);
//                     // 2020.11.03 分析，虽然有可能头尾通过caustics crossing相连，但不能在第一次大循环时就抹除了所有的可能性，要设置一个flag，表示这不是第一次循环到这时，执行下面的。
// #endif
//                     imageTracks->drop(Prov);
//                     connected->append(Prov);
//                     Prov = tempProv;
//                 } else {
//                     Prov = Prov->next;
//                 }
//             } else {
//                 Prov = Prov->next;
//             }
//         }

 
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
            fprintf(stderr, "\t\t scan length %d\n", scan->length);
        }
#endif
    }

// begin **************************** remove on 2020.07.14
    // check if there is some segments actually needs to be separate to two segment
    // #添加对 imageTracks 中，是否存在某段连续的其实要分为两段的检查
//     _sols *imageTracks_checked;
//     imageTracks_checked = new _sols;
//     for (int i = 0; i < 2 * imageTracks->length; i++) {
//         Prov = new _curve;
//         imageTracks_checked->append(Prov);
//     }
//     int reach = 0;
//     _point *head_or_tail_pnt, *p;

//     Prov2 = imageTracks_checked->first;
//     Prov = imageTracks->first;
//     while (Prov) {
//         p = Prov->first;
//         while (p) {
//             // do not check p = Prov.first and Prov.last # 注意不要检查 p = Prov.first 和 Prov.last
//             if (p == Prov->first || p == Prov->last) {
//                 addPoint(Prov2, p);
//             }
//             else {
//                 head_or_tail_pnt = pnt_reach_head_or_tail(p, Prov, imageTracks, &reach);
//                 if (!reach) {addPoint(Prov2, p);}
//                 else {
// #ifdef VERBOSE
//                     fprintf(stderr, "special case!!!\n\n\n\n");
// #endif
//                     fprintf(stderr, "special case!!!\n\n\n\n");
//                     addPoint(Prov2, p);
//                     addPoint(Prov2, head_or_tail_pnt);
//                     Prov2 = Prov2->next;
//                 }
//             }
//             p = p->next;
//         }
//         Prov = Prov->next;
//         Prov2 = Prov2->next;
//     }

//     // prune empty curves in the connected image track array
//     Prov = imageTracks_checked->first;
//     while (Prov) {
//         if (!Prov->first || !Prov->last) { // drop any empty curve
//             Prov2 = Prov2->next;
//             imageTracks_checked->drop(Prov);
//             delete Prov;
//             Prov = Prov2;
//         } else {
//             Prov = Prov->next;
//         }

//     }
//     delete imageTracks;
//     imageTracks = new _sols;
//     imageTracks = imageTracks_checked;
// end **************************** remove on 2020.07.14


    // show image track animation here:
    // saveTracks(imageTracks);


    //
    // now we deal with the final image track that are non-closed, presumbly associated with caustic crossing
    //
    _curve *final;

    // final = new _curve;
    // final = newSegment(imageTracks);
    // imageTracks->drop(imageTracks->first, 1); // drop the first segment as it is already in final

    _newSegment


    //
    //we first connect identical image positions and then try identical source positions
    //
    complex zhead, ztail;//dz
    int ifcontinue = 0;

    // use another way, less reverse
    int head = 0, tail = 0;
    // double  muPrevious;
    _curve *markedCurve;
    while ( imageTracks->length >= 1) {
        _final_parity
        final = connect_head_or_tail(&ifcontinue, final, imageTracks, connected);
        if (ifcontinue) continue;
        else {

#ifdef VERBOSE
            // fprintf(stderr, "can not jump, reversing final and then check if connectWithHead or connectWithTail, or jump over caustics\n");
            fprintf(stderr, "can not connect, checking whether the head of final can connect or jump over caustics\n");
#endif
            // _final_parity
            // final->reverse();
            // final->parity *= -1;

            final = connect_head_or_tail(&ifcontinue, final, imageTracks, connected, true);
            if (ifcontinue) continue;
            else {
                markedCurve = jump_over_caustics(final, imageTracks, &head, &tail, markedCurve, true);
// here we jump over caustics
                if (head == 1) {               // jumped over caustics to the head of a curve
#ifdef VERBOSE
                    printf("jumped to head of a curve: %f %f %f %f, phi=%.15f\n", markedCurve->first->x1, markedCurve->first->x2, muPrevious, markedCurve->first->mu, markedCurve->first->phi * 180 / M_PI);
                    printf("jumped to head and going to tail at: %f %f\n", markedCurve->last->x1, markedCurve->last->x2);
#endif
                    // add on 2020.04.14

                    // markedCurve->first->closepairparity = (markedCurve->first->mu)>0?1:-1;
                    ifjump = 1;


                    // for (_point *p = markedCurve->first; p; p = p->next) { //  jumped over caustics to the head of a curve, attach the points in natural order
                    //     addPoint(final, p);
                    // }

#ifdef VERBOSE
                    fprintf(stderr, "imageTracks->length before drop %d\n", imageTracks->length);
#endif
                    imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list
// // delete markedCurve; // 2020.03.13
#ifdef VERBOSE
                    fprintf(stderr, "imageTracks->length after drop %d\n", imageTracks->length);
#endif

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

// if (imageTracks->length == 0){
//     connected->append(final);
//     break;
// }

                } else if (tail == 1) { // jumped over caustics to the tail of a curve
#ifdef VERBOSE
                    printf("jumped to tail of a curve: %f %f %f %f, phi=%.15f\n", markedCurve->last->x1, markedCurve->last->x2, muPrevious, markedCurve->last->mu, markedCurve->last->phi * 180 / M_PI);
                    printf("jumped to tail and going to head: %f %f\n", markedCurve->first->x1, markedCurve->first->x2);
#endif
                    // add on 2020.04.14
                    final->first->closepairparity = final->first->mu > 0 ? 1 : -1;
                    markedCurve->last->closepairparity = (markedCurve->last->mu) > 0 ? 1 : -1;
                    ifjump = 1;



                    // for (_point *p = markedCurve->last; p; p = p->prev) { //  jumped over caustics to the tail of a curve, attach the points in reverse order
                    //     addPoint(final, p);
                    // }

#ifdef VERBOSE
                    fprintf(stderr, "imageTracks->length before drop %d\n", imageTracks->length);
#endif
                    imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list
// // delete markedCurve; // 2020.03.13
#ifdef VERBOSE
                    fprintf(stderr, "imageTracks->length after drop %d\n", imageTracks->length);
#endif

                    markedCurve->last->next = final->first;
                    final->first->prev = markedCurve->last;
                    markedCurve->last = final->last;
                    markedCurve->length += final->length;
                    markedCurve->parity = final->parity;
                    final = markedCurve;
// if (imageTracks->length == 0){
//     connected->append(final);
//     break;
// }

                } else {      // we failed to jump at head


                    // we are not yet self closed, we need  to jump over caustics, we seek a point with the same source position but opposite parity, and smallest magnification ratio
                    muPrevious = final->last->mu;
                    markedCurve = jump_over_caustics(final, imageTracks, &head, &tail, markedCurve);
                    // here we jump over caustics
                    if (head == 1) {               // jumped over caustics to the head of a curve
#ifdef VERBOSE
                        printf("jumped to head of a curve: %f %f %f %f, final->last->mu %f, phi=%.15f\n", markedCurve->first->x1, markedCurve->first->x2, muPrevious, markedCurve->first->mu, final->last->mu, markedCurve->first->phi * 180 / M_PI);
                        printf("jumped to head and going to tail at: %f %f\n", markedCurve->last->x1, markedCurve->last->x2);
#endif
                        // add on 2020.04.14
                        final->last->closepairparity = (final->last->mu) > 0 ? 1 : -1;
                        markedCurve->first->closepairparity = (markedCurve->first->mu) > 0 ? 1 : -1;
                        ifjump = 1;

                        // for (_point *p = markedCurve->first; p; p = p->next) { //  jumped over caustics to the head of a curve, attach the points in natural order
                        //     addPoint(final, p);
                        // }

#ifdef VERBOSE
                        fprintf(stderr, "markedCurve->length %d\n", markedCurve->length);
                        fprintf(stderr, "imageTracks->length before drop %d, connected->length = %d\n", imageTracks->length, connected->length);
#endif
                        imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list
// // delete markedCurve; // 2020.03.13
#ifdef VERBOSE
                        fprintf(stderr, "imageTracks->length after drop %d, connected->length = %d\n", imageTracks->length, connected->length);
#endif

                        // 2020.07.17, attach two segments directly instead of for loop above, to speed up the code.
                        _final_parity
                        final->last->next = markedCurve->first;
                        markedCurve->first->prev = final->last;
                        final->last = markedCurve->last;
                        final->length += markedCurve->length;

// if (imageTracks->length == 0){
//     connected->append(final);
//     break;
// }

                    } else if (tail == 1) { // jumped over caustics to the tail of a curve
#ifdef VERBOSE
                        printf("jumped to tail of a curve: %f %f %f %f, final->last->mu %f, phi=%.15f\n", markedCurve->last->x1, markedCurve->last->x2, muPrevious, markedCurve->last->mu, final->last->mu, markedCurve->last->phi * 180 / M_PI);
                        printf("jumped to tail and going to head: %f %f\n", markedCurve->first->x1, markedCurve->first->x2);
#endif
                        // add on 2020.04.14
                        final->last->closepairparity = (final->last->mu) > 0 ? 1 : -1;

                        ifjump = 1;

                        // for (_point *p = markedCurve->last; p; p = p->prev) { //  jumped over caustics to the tail of a curve, attach the points in reverse order
                        //     addPoint(final, p);
                        // }

#ifdef VERBOSE
                        fprintf(stderr, "imageTracks->length before drop %d, connected->length = %d\n", imageTracks->length, connected->length);
#endif
                        imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list
// // delete markedCurve; // 2020.03.13
#ifdef VERBOSE
                        fprintf(stderr, "imageTracks->length after drop %d, connected->length = %d\n", imageTracks->length, connected->length);
#endif

                        // 2020.07.17
                        _final_parity
                        markedCurve->reverse();
                        markedCurve->first->closepairparity = (markedCurve->first->mu) > 0 ? 1 : -1;
                        final->last->next = markedCurve->first;
                        markedCurve->first->prev = final->last;
                        final->last = markedCurve->last;
                        final->length += markedCurve->length;

// if (imageTracks->length == 0){
//     connected->append(final);
//     break;
// }

                    } else {      // we failed to connect, something went wrong



#ifdef VERBOSE
                        fprintf(stderr, "connected with neither head nor tail, something went wrong, append final anyway\n");
#endif
                        _final_parity
#ifdef VERBOSE
                        fprintf(stderr, "\t >>>connected->append(final) place 5, final->length = %d\n", final->length);
#endif
                        connected->append(final);

#ifdef VERBOSE
                        fprintf(stderr, "due to may be small mass ratio (lack solution), some image tracks is not closed, just drop it, the area is small anyway, imageTracks.length = %d, connected.length = %d, final->length = %d\n", imageTracks->length, connected->length, final->length);
#endif

                        if (imageTracks->length > 0) {
                            connected->append(final);
                            final = NULL;
                            delete final;
#ifdef VERBOSE
                            fprintf(stderr, "deleting final\n");
#endif

                            // final = new _curve;
                            // final = newSegment(imageTracks);
                            // imageTracks->drop(imageTracks->first, 1);

                            _newSegment
                        }




                        continue;

                    }

                    if (head || tail) {

#ifdef VERBOSE
                        fprintf(stderr, "jumped, checking whether head and tail connects, ");
                        fprintf(stderr, "and continue the cycle to see weather we could attach to final a new subtrack\n");
#endif
                    } else {
#ifdef VERBOSE
                        fprintf(stderr, "not jumped, checking whether head and tail connects\n");
#endif
                    }

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
                        if (imageTracks->length == 0) {
                            delete allSolutions;
                            delete imageTracks;
                            Prov2 = NULL;
                            Prov = NULL;
                            delete Prov;
                            delete Prov2;
                            // delete markedCurve;
                            delete tempProv;
                            return (connected);
                        }

                        //we still have segments remaining, start a new curve
                        // final = new _curve;
                        // final = newSegment(imageTracks);
                        // imageTracks->drop(imageTracks->first, 1);


                        _newSegment
                        continue;
                    }
//                          else {
// #ifdef VERBOSE
//                             printf("check if current final is self-closed before continue the cycle place 1\n");
// #endif

//                             mindis = 1e9;
//                             for (_curve * Prov3 = imageTracks->first; Prov3; Prov3 = Prov3->next) {

//                                 tempdis = *(final->first) - *(Prov3->first);
//                                 if (tempdis <= mindis) {
//                                     mindis = tempdis;
//                                 }

//                                 tempdis = *(final->first) - *(Prov3->last);
//                                 if (tempdis <= mindis) {
//                                     mindis = tempdis;
//                                 }

//                                 tempdis = *(final->last) - *(Prov3->first);
//                                 if (tempdis <= mindis) {
//                                     mindis = tempdis;
//                                 }

//                                 tempdis = *(final->last) - *(Prov3->last);
//                                 if (tempdis <= mindis) {
//                                     mindis = tempdis;
//                                 }
//                             }
//                             if (mindis > SD) {
// #ifdef VERBOSE
//                                 fprintf(stderr, "No more track will be in final, start a new one and continue\n");
// #endif
//                                 // if (final->parity == 0) {
//                                 //     if (final->first->mu > 0) {
//                                 //         final->parity = 1;
//                                 //     }
//                                 //     else {
//                                 //         final->parity = -1;
//                                 //     }
//                                 // }
//                                 _final_parity
//                                 connected->append(final);
//                                 if (imageTracks->length != 0) {
//                                     // final = newSegment(imageTracks);
//                                     // imageTracks->drop(imageTracks->first, 1);
//                                     _newSegment
//                                 }
//                             }
//                         }
                    else {
                        // continue;   // not yet connected, we should continue the cycle


                        if ( ( abs(final->first->phi - final->last->phi) < EPS ||  abs( abs(final->first->phi - final->last->phi) - 2.0 * M_PI) < EPS) ) {

                            // we check whether the slope changes are smooth using four previous points
                            int i = 0;
                            double slope[2], slope_if_jumped;
                            double scatter;

#ifdef VERBOSE
                            fprintf(stderr, "number of points so far %d\n", final->length);
#endif

                            if (final->length < 3) {
#ifdef VERBOSE
                                fprintf(stderr, "too few points in determing slope ratios and scatters: needs to increase the number of points\n");
#endif

                                delete allSolutions;
                                delete imageTracks;
                                delete connected;
                                Prov2 = NULL;
                                Prov = NULL;

                                delete Prov;
                                delete tempProv;
                                delete Prov2;
                                delete final;
                                // delete markedCurve;
                                return (NULL);
                            }

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
#ifdef VERBOSE
                                printf("finally dropped %d\n", final->length);
#endif

                                if (final->length > 10 && imageTracks->length == 0) {
                                    addPoint(final, final->first); // we simply connect the first and last points
                                } else {
                                    continue;
                                }

                            } else {    // slopes changed too much, we should not jump; this is somewhat ambigious here
                                continue;
                            }
                        } else {
                            continue;
                        }


                    }


                    //ERROR correction begin. We forgot to check the case when the track is already closed when
                    //the first and last points share the same source position.
                    //
                    dzs = final->first->zs - final->last->zs;
#ifdef VERBOSE
                    printf("first and last source position difference %e\n", abs(dzs));
#endif
                    if (abs(dzs) < EPS) {

                        if (final->length >= PHI->length && imageTracks->length == 0) {
                            _final_parity

#ifdef VERBOSE
                            fprintf(stderr, "\t >>>connected->append(final) place 7, final->length = %d\n", final->length);
#endif
                            connected->append(final);

                            //we start a new curve
                            if (imageTracks->length != 0) {


                                // final = NULL;
                                // delete final;
                                // final = new _curve;
                                // final = newSegment(imageTracks);
                                // // temp = imageTracks->first;
                                // imageTracks->drop(imageTracks->first, 1);


                                _newSegment


                            } // we already incorporated the new curve, and so drop it
                            // delete temp;//add on 2020.03.13
                            continue;
                        } else {
                            continue;
                        }

                    } else {
                        if (imageTracks->length == 0) {
                            delete allSolutions;
                            delete imageTracks;
                            Prov2 = NULL;
                            Prov = NULL;
                            delete Prov;
                            delete tempProv;
                            delete Prov2;
                            // delete markedCurve;
                            return (connected); // we are done
                        } else {      // start a new curve

                            // final = new _curve;
                            // final = newSegment(imageTracks);
                            // // temp = imageTracks->first;
                            // imageTracks->drop(imageTracks->first, 1); // we already incorporated the new curve, and so drop it

                            _newSegment

                            continue;
                        }

                    }
                    //ERROR correction end


                }

                if (head || tail) {

#ifdef VERBOSE
                    fprintf(stderr, "jumped, checking whether head and tail connects, ");
                    fprintf(stderr, "and continue the cycle to see weather we could attach to final a new subtrack\n");
#endif
                } else {
#ifdef VERBOSE
                    fprintf(stderr, "not jumped, checking whether head and tail connects\n");
#endif
                }

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
                    if (imageTracks->length == 0) {
                        delete allSolutions;
                        delete imageTracks;
                        Prov2 = NULL;
                        Prov = NULL;
                        delete Prov;
                        delete Prov2;
                        // delete markedCurve;
                        delete tempProv;
                        return (connected);
                    }

                    //we still have segments remaining, start a new curve
                    // final = new _curve;
                    // final = newSegment(imageTracks);
                    // imageTracks->drop(imageTracks->first, 1);

                    _newSegment

                    continue;
                }
//                      else {
// #ifdef VERBOSE
//                         printf("check if current final is self-closed before continue the cycle place 2\n");
// #endif

//                         // continue;   // not yet connected, we should continue the cycle
//                         mindis = 1e9;
//                         for (_curve * Prov3 = imageTracks->first; Prov3; Prov3 = Prov3->next) {

//                             tempdis = *(final->first) - *(Prov3->first);
//                             if (tempdis <= mindis) {
//                                 mindis = tempdis;
//                             }

//                             tempdis = *(final->first) - *(Prov3->last);
//                             if (tempdis <= mindis) {
//                                 mindis = tempdis;
//                             }

//                             tempdis = *(final->last) - *(Prov3->first);
//                             if (tempdis <= mindis) {
//                                 mindis = tempdis;
//                             }

//                             tempdis = *(final->last) - *(Prov3->last);
//                             if (tempdis <= mindis) {
//                                 mindis = tempdis;
//                             }
//                         }
//                         if (mindis > SD) {
// #ifdef VERBOSE
//                             fprintf(stderr, "No more track will be in final, start a new one and continue 2, connected->length = %d\n", connected->length);
// #endif
//                             // if (final->parity == 0) {
//                             //     if (final->first->mu > 0) {
//                             //         final->parity = 1;
//                             //     }
//                             //     else {
//                             //         final->parity = -1;
//                             //     }
//                             // }
//                             _final_parity
//                             connected->append(final);
// #ifdef VERBOSE
//                             fprintf(stderr, "after append, connected->length = %d\n", connected->length);
// #endif

//                             if (imageTracks->length != 0) {
//                                 // final = newSegment(imageTracks);
//                                 // imageTracks->drop(imageTracks->first, 1);
//                                 _newSegment
//                                 continue;
//                             } else {
//                                 continue;
//                             }
//                         }
//                     }
                else {
                    // continue;   // not yet connected, we should continue the cycle


                    if ( ( abs(final->first->phi - final->last->phi) < EPS ||  abs( abs(final->first->phi - final->last->phi) - 2.0 * M_PI) < EPS) ) {

                        // we check whether the slope changes are smooth using four previous points
                        int i = 0;
                        double slope[2], slope_if_jumped;
                        double scatter;

#ifdef VERBOSE
                        fprintf(stderr, "number of points so far %d\n", final->length);
#endif

                        if (final->length < 3) {
#ifdef VERBOSE
                            fprintf(stderr, "too few points in determing slope ratios and scatters: needs to increase the number of points\n");
#endif

                            delete allSolutions;
                            delete imageTracks;
                            delete connected;
                            Prov2 = NULL;
                            Prov = NULL;

                            delete Prov;
                            delete tempProv;
                            delete Prov2;
                            delete final;
                            // delete markedCurve;
                            return (NULL);
                        }

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
#ifdef VERBOSE
                            printf("finally dropped %d\n", final->length);
#endif
                            if (final->length > 10 && imageTracks->length == 0) {
                                addPoint(final, final->first); // we simply connect the first and last points
                            }
                            else {
                                continue;
                            }

                        }
                        else {    // slopes changed too much, we should not jump; this is somewhat ambigious here
                            continue;
                        }
                    } else {
                        continue;
                    }


                }


                //ERROR correction begin. We forgot to check the case when the track is already closed when
                //the first and last points share the same source position.
                //
                dzs = final->first->zs - final->last->zs;
#ifdef VERBOSE
                printf("first and last source position difference %e\n", abs(dzs));
#endif
                if (abs(dzs) < EPS) {

                    if (final->length >= PHI->length && imageTracks->length == 0) {
                        _final_parity

#ifdef VERBOSE
                        fprintf(stderr, "\t >>>connected->append(final) place 7, final->length = %d\n", final->length);
#endif
                        connected->append(final);

                        //we start a new curve
                        if (imageTracks->length != 0) {

                            // final = NULL;
                            // delete final;
                            // final = new _curve;
                            // final = newSegment(imageTracks);
                            // // temp = imageTracks->first;
                            // imageTracks->drop(imageTracks->first, 1);

                            _newSegment

                        } // we already incorporated the new curve, and so drop it
                        // delete temp;//add on 2020.03.13
                        continue;
                    } else {
                        continue;
                    }

                } else {
                    if (imageTracks->length == 0) {
                        delete allSolutions;
                        delete imageTracks;
                        Prov2 = NULL;
                        Prov = NULL;
                        delete Prov;
                        delete tempProv;
                        delete Prov2;
                        // delete markedCurve;
                        return (connected); // we are done
                    } else {      // start a new curve
                        // final = new _curve;
                        // final = newSegment(imageTracks);
                        // // temp = imageTracks->first;
                        // imageTracks->drop(imageTracks->first, 1); // we already incorporated the new curve, and so drop it

                        _newSegment


                    }

                }
                //ERROR correction end





            }

        }






    }


//check whether the first and last points are connected, if not, join the first and last points
    SD = *(final->first) - *(final->last);

#ifdef VERBOSE
    printf("final->length = %d, joining the first and last point %f %f %f %f %f %f\n", final->length, final->first->x1, final->first->x2, final->first->phi,
           final->last->x1, final->last->x2, final->last->phi);
    fprintf(stderr, "SD = %.3e\n", SD);
#endif




    if (connected->last != final && final->length > 1) {
        _final_parity

        // this should comment out: on 2020.03.12 ?
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
    return (connected);
}


// _sols *TripleLensing::outputTracks(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, double phi0)
_sols *TripleLensing::outputTracks(double xsCenter, double ysCenter, double rs, int nphi, double phi0)
{
    // VBBinaryLensing VBBL;

    // now find the closed curves for all real and imaginary roots for a circular source
    // double xs, ys;
    // double phi, dphi;
    // complex zr[DEGREE];
    // complex dzs;
    // complex coefficients[DEGREE + 1];
    nsolution = 0;
    int i, j, k;
    // complex J1, J2, dJ, J3;
    // double mu, lambda1, lambda2, thetaJ;

    _sols *allSolutions, *imageTracks;
    _curve *Prov = new _curve;
    _curve *Prov2 = new _curve;
    // Prov2 = NULL;
    // _curve *Prov;
    _point *pisso;
    double SD, MD, CD;

    // connect all images, including fake images
    allSolutions = new _sols;
    for (int i = 0; i < DEGREE; i++) {
        Prov = new _curve;
        allSolutions->append(Prov);
    }

    // int nimages;      // number of true images
    // int flag;
    dphi = 2.0 * M_PI / (nphi - 1);
    phi = -dphi + phi0;

    for (int j = 0; j < nphi; j++) {
        phi += dphi;
        xs = xsCenter + rs * cos(phi);
        ys = ysCenter + rs * sin(phi);

        // polynomialCoefficients(mlens, zlens, xs, ys, coefficients, NLENS, DEGREE);
        polynomialCoefficients(xs, ys, coefficients);
        VBBL.cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

        if (j == 0) {   // first point in all solution tracks
            Prov = allSolutions->first;

            nimages = 0;
            for (int i = 0; i < DEGREE; i++) {
                // flag = trueSolution(mlens, zlens, xs, ys, zr[i], &mu, &lambda1, &lambda2, &thetaJ, NLENS, &J1, &J2, &dJ, &J3);
                flag = trueSolution(xs, ys, zr[i]);
                nimages += flag;

                Prov->append(zr[i].re, zr[i].im);
                Prov->last->phi = phi;
                Prov->last->mu = mu;
                Prov->last->zs = complex(xs, ys);
                Prov->last->thetaJ = thetaJ;

                if (flag == 1) {      // true solution
                    Prov->last->flag = 1 ;
                } else {    // false solution
                    Prov->last->flag = -1 ;
                }

                Prov = Prov->next;
            }
            // check image numbers
            if ((nimages - NLENS - 1) % 2 != 0) {
                fprintf(stderr, "warning: image number is wrong: %d\n", nimages);
            }
        } else {
            delete Prov2;
            Prov2 = new _curve;
            nimages = 0;
            for (int i = 0; i < DEGREE; i++) {
                // flag = trueSolution(mlens, zlens, xs, ys, zr[i], &mu, &lambda1, &lambda2, &thetaJ, NLENS, &J1, &J2, &dJ, &J3);
                flag = trueSolution(xs, ys, zr[i]);
                nimages += flag;

                Prov2->append(zr[i].re, zr[i].im);
                Prov2->last->phi = phi;
                Prov2->last->mu = mu;
                Prov2->last->zs = complex(xs, ys);
                Prov2->last->thetaJ = thetaJ;
                if (flag == 1) {      // true solution
                    Prov2->last->flag = +1 ;
                } else {
                    Prov2->last->flag = -1;
                }
            }

            //check image numbers
            if ((nimages - NLENS - 1) % 2 != 0) {
                fprintf(stderr, "warning: image number is wrong: %d\n", nimages);
            }

            // attach the current solutions positions to last ones with the closest distance
            for (Prov = allSolutions->first; Prov;  Prov = Prov->next) {
                Prov2->closest(Prov->last, &pisso);
                Prov2->drop(pisso);
                Prov->append(pisso);
                // delete pisso; //add on 2020.03.13
            }
        }
    }

    //allocate memory for potentially true image tracks, we make a factor four more segments just in case
    imageTracks = new _sols;
    for (int i = 0; i < 4 * DEGREE; i++) {
        Prov = new _curve;
        imageTracks->append(Prov);
    }
    //
    // Solution tracks may contain only true images, some only false solutions, while some may be mixed
    //
    complex z, zs;
    int mixedTracks = 0; // number of mixed tracks

    Prov2 = imageTracks->first;
    for (Prov = allSolutions->first; Prov;  Prov = Prov->next) {
        int trueImages = 0, falseImages = 0;

        int firstTime = 1;
        int previousImage = 0;

        int Parity;   // within the same curve, we require the parity to be the same

        double ratio;
        double mu, muPrevious;

        for (_point *p = Prov->first; p; p = p->next) {
            z = complex(p->x1, p->x2);

            if (p->flag == 1) { // true solutions, should return identical source radius
                trueImages++;

                if (firstTime) {  // we enter the track the first time
                    firstTime = 0;
                    muPrevious = p->mu;
                } else {
                    // parity has changed, start a new curve
                    if (muPrevious * p->mu < 0.0)  {
                        Prov2 = Prov2->next;
                        muPrevious = p->mu;
                    }
                }

                addPoint(Prov2, p); // add point p to the track

                previousImage = TRUE_IMAGE;
            } else {

                falseImages++;

                //first time that we enter this
                if (firstTime) {
                    previousImage = FALSE_IMAGE;

                    firstTime = 0;
                }

                // already a true segment existing, we cut the segment, and start a new one
                if (previousImage == TRUE_IMAGE) { // previously we had true images, now we are encountering on false solutions

                    Prov2 = Prov2->next;  // we start a new image track

                    previousImage = FALSE_IMAGE;
                }
            }
        }

        if (trueImages) Prov2 = Prov2->next;

        if (trueImages && falseImages) mixedTracks += 1; // there are tracks that are mixed
    }

    // prune empty image tracks first before proceeding
    for (Prov = imageTracks->first; Prov;  Prov = Prov->next) {
        // if (!Prov->first || !Prov->last) { // no points in curve and so drop the empty curve
        if (!Prov->first || !Prov->last || Prov->length < 2) { // no points in curve and so drop the empty curve
            imageTracks->drop(Prov);
            // delete Prov; //add on 2020.03.13 // do not use this
        }
    }


    // The easy case: every curve is self-closed without caustics-crossing, we are done
    if (mixedTracks == 0)  {
        if ( (imageTracks->length - (NLENS + 1)) % 2 == 0) { // double check the number of tracks should be 4 + 2*n in this case
            // delete allSolutions;
            fprintf(stderr, "return from plc1\n");

            for (Prov = imageTracks->first; Prov; Prov = Prov->next) {
                if (Prov->first->mu > 0) {
                    Prov->parity = 1;
                }
                else {
                    Prov->parity = -1;
                }

            }
            delete allSolutions;
            return (imageTracks);
        } else {
            fprintf(stderr, "image tracks should be %d, but is %d\n", NLENS + 1, imageTracks->length);
            delete allSolutions;
            delete imageTracks;
            return (NULL);
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
    connected = new _sols;
    for (int i = 0; i < imageTracks->length; i++) {
        Prov = new _curve;
        connected->append(Prov);
    }

    //
    //now attach self-closed, non-caustic crossing image tracks
    //
    int selfClosed = 0;

    for (Prov2 = connected->first, Prov = imageTracks->first; Prov;  Prov = Prov->next) {
        SD = *(Prov->first) - *(Prov->last);
        if (SD < EPS * EPS) { // closed tracks

            int pointInTrack = 0; // how many points in track
            for (_point *p = Prov->first; p; p = p->next) {
                addPoint(Prov2, p); // add point p to the track

                pointInTrack++;
            }

            selfClosed++;
#ifdef VERBOSE
            printf("entered: pointInTrack %d, imageTracks->length %d\n", pointInTrack, imageTracks->length);
#endif
            //
            // such image tracks must have number of points == number of points on the source circle
            //
            // special case may be there is only one single point around a narrow cusp, / we probably
            // need to increase the number of points to solve the problem
            //
            if (pointInTrack != nphi) {
#ifdef VERBOSE
                fprintf(stderr, "the number of points is %d, should be equal to nphi=%d!\n", pointInTrack, nphi);
#endif
            }

            Prov2 = Prov2->next;  // go the next closed curve
            imageTracks->drop(Prov);  // drop the already connected curve from imageTracks
        }
    }


    // prune empty curves in the connected image track array
    for (Prov = connected->first; Prov;  Prov = Prov->next) {
        if (!Prov->first || !Prov->last) { // drop any empty curve
            connected->drop(Prov);
            // delete Prov; //add on 2020.03.13
        }
        else {
            if (Prov->first->mu > 0) {
                Prov->parity = 1;
            }
            else {
                Prov->parity = -1;
            }

        }
    }


    // check is there some segments which actually need to separate into two subsegments
    //#添加对 imageTracks 中，是否存在某段连续的其实要分为两段的检查
    _sols *imageTracks_checked;
    imageTracks_checked = new _sols;
    for (int i = 0; i < 2 * imageTracks->length; i++) {
        Prov = new _curve;
        imageTracks_checked->append(Prov);
    }
    int reach = 0;
    _point *head_or_tail_pnt, *p;

    Prov2 = imageTracks_checked->first;
    Prov = imageTracks->first;
    while (Prov) {
        p = Prov->first;
        while (p) {
            // no need check p = Prov.first and Prov.last # 注意不要检查 p = Prov.first 和 Prov.last
            if (p == Prov->first || p == Prov->last) {
                addPoint(Prov2, p);
            }
            else {
                head_or_tail_pnt = pnt_reach_head_or_tail(p, Prov, imageTracks, &reach);
                if (!reach) {                                        addPoint(Prov2, p);}
                else {
                    fprintf(stderr, "special case!!!\n\n\n\n");
                    addPoint(Prov2, p);
                    addPoint(Prov2, head_or_tail_pnt);
                    Prov2 = Prov2->next;
                }
            }
            p = p->next;
        }
        Prov = Prov->next;
        Prov2 = Prov2->next;
    }

    // prune empty curves in the connected image track array
    for (Prov = imageTracks_checked->first; Prov;  Prov = Prov->next) {
        if (!Prov->first || !Prov->last) { // drop any empty curve
            imageTracks_checked->drop(Prov);
        }
    }

    // delete imageTracks;
    // _sols *imageTracks;
    imageTracks = imageTracks_checked;



    //
    // now we deal with the final image track that are non-closed, presumbly associated with caustic crossing
    //
    _curve *final;

    // final = new _curve;
    // final = newSegment(imageTracks);
    // imageTracks->drop(imageTracks->first); // drop the first segment as it is already in final

    _newSegment
    //
    //we first connect identical image positions and then try identical source positions
    //
    complex dz, zhead, ztail;
    int ifcontinue = 0;

    while ( imageTracks->length >= 1) {
        final = connect_head_or_tail(&ifcontinue, final, imageTracks, connected);

        if (ifcontinue) continue;
        else {
            if (!final->parity) {
                if (final->first->mu > 0) {
                    final->parity = 1;
                }
                else {
                    final->parity = -1;
                }
            }
#ifdef VERBOSE
            fprintf(stderr, "reversing final to see if connectWithHead or connectWithTail\n");
#endif

            final->reverse();
            final->parity *= -1;

            final = connect_head_or_tail(&ifcontinue, final, imageTracks, connected);
            if (ifcontinue) continue;
        }
        //
        //ERROR correction begin. We forgot to check the case when the track is already closed when
        //the first and last points share the same source position.
        //
        dzs = final->first->zs - final->last->zs;
#ifdef VERBOSE
        printf("first and last source position difference %e\n", abs(dzs));
#endif
        if (abs(dzs) < EPS) {
            if (!final->parity) {
                if (final->first->mu > 0) {
                    final->parity = 1;
                }
                else {
                    final->parity = -1;
                }
            }
            connected->append(final);

            printf("entered 99999e10\n");

            //we start a new curve
            if (imageTracks->length != 0) {

                // final = new _curve;
                // final = newSegment(imageTracks);
                // // temp = imageTracks->first;
                // imageTracks->drop(imageTracks->first);

                _newSegment

            } // we already incorporated the new curve, and so drop it
            // delete temp;//add on 2020.03.13
            continue;
        }
        //ERROR correction end


        // we are not yet self closed, we need  to jump over caustics, we seek a point with the same source position but opposite parity, and smallest magnification ratio
        int head = 0, tail = 0;
        double  muPrevious;
        muPrevious = final->last->mu;

        _curve *markedCurve;
        // markedCurve = new _curve;
        markedCurve = jump_over_caustics(final, imageTracks, &head, &tail, markedCurve);
        fprintf(stderr, "head tail outside func: %d %d\n", head, tail);

        // here we jump over caustics
        if (head == 1) {               // jumped over caustics to the head of a curve
#ifdef VERBOSE
            fprintf(stderr, "can you go here??\n");
            printf("jumped to head of a curve: %f %f %f %f\n", markedCurve->first->x1, markedCurve->first->x2, muPrevious, markedCurve->first->mu);
            printf("jumped to head and going to tail at: %f %f\n", markedCurve->last->x1, markedCurve->last->x2);
#endif

            for (_point *p = markedCurve->first; p; p = p->next) { //  jumped over caustics to the head of a curve, attach the points in natural order
                addPoint(final, p);
            }
        } else if (tail == 1) { // jumped over caustics to the tail of a curve
#ifdef VERBOSE
            printf("jumped to tail of a curve: %f %f %f %f\n", markedCurve->last->x1, markedCurve->last->x2, muPrevious, markedCurve->last->mu);
            printf("jumped to tail and going to head: %f %f\n", markedCurve->first->x1, markedCurve->first->x2);
#endif

            for (_point *p = markedCurve->last; p; p = p->prev) { //  jumped over caustics to the tail of a curve, attach the points in reverse order
                addPoint(final, p);
            }
        } else {      // we failed to connect, something went wrong
            fprintf(stderr, "connected with neither head nor tail, something went wrong\nreversing final and then jump over caustics\n");

            if (!final->parity) {
                if (final->first->mu > 0) {
                    final->parity = 1;
                }
                else {
                    final->parity = -1;
                }
            }
            final->reverse();
            final->parity *= -1;

            markedCurve = jump_over_caustics(final, imageTracks, &head, &tail, markedCurve);

// here we jump over caustics
            if (head == 1) {               // jumped over caustics to the head of a curve
#ifdef VERBOSE
                printf("jumped to head of a curve: %f %f %f %f\n", markedCurve->first->x1, markedCurve->first->x2, muPrevious, markedCurve->first->mu);
                printf("jumped to head and going to tail at: %f %f\n", markedCurve->last->x1, markedCurve->last->x2);
#endif

                for (_point *p = markedCurve->first; p; p = p->next) { //  jumped over caustics to the head of a curve, attach the points in natural order
                    addPoint(final, p);
                }
            } else if (tail == 1) { // jumped over caustics to the tail of a curve
#ifdef VERBOSE
                printf("jumped to tail of a curve: %f %f %f %f\n", markedCurve->last->x1, markedCurve->last->x2, muPrevious, markedCurve->last->mu);
                printf("jumped to tail and going to head: %f %f\n", markedCurve->first->x1, markedCurve->first->x2);
#endif

                for (_point *p = markedCurve->last; p; p = p->prev) { //  jumped over caustics to the tail of a curve, attach the points in reverse order
                    addPoint(final, p);
                }
            } else {      // we failed to connect, something went wrong
                fprintf(stderr, "connected with neither head nor tail, something went wrong\nreversing final and then jump over caustics\n");
                delete allSolutions;
                delete imageTracks;
                delete Prov;
                delete Prov2;
                delete final;
                return (NULL);
            }
        }

        imageTracks->drop(markedCurve); // drop the already connected, marked curve from the linked list
        fprintf(stderr, "dropping and deleting markedCurve\n");
        // delete markedCurve; // 2020.03.13
        if (head || tail) {
#ifdef VERBOSE
            fprintf(stderr, "jumped, checking whether head and tail connects\n");
            fprintf(stderr, "and continue the cycle to see weather we could attach to final a new subtrack\n");
#endif
            continue;
        }

        //    After jumping, we now check whether the tail and head of the current curve actually connect.
        SD = *(final->first) - *(final->last);

        if (SD < EPS * EPS) { // we are already connected

#ifdef VERBOSE
            fprintf(stderr, "jumped, head and tail connected!\n");
#endif

            if (!final->parity) {
                if (final->first->mu > 0) {
                    final->parity = 1;
                }
                else {
                    final->parity = -1;
                }
            }

            connected->append(final);

            //all segments connected, we are done
            if (imageTracks->length == 0) {
                fprintf(stderr, "I am deleting......\n");
                delete allSolutions;
                delete imageTracks;
                delete Prov;
                delete Prov2;
                delete markedCurve;
                return (connected);
            }

            //we still have segments remaining, start a new curve
            // final = new _curve;
            // final = newSegment(imageTracks);
            // // temp = imageTracks->first;
            // imageTracks->drop(imageTracks->first);

            _newSegment


            continue;
        }  else if ( ( abs(final->first->phi - final->last->phi) < EPS ||  abs( abs(final->first->phi - final->last->phi) - 2.0 * M_PI) < EPS) ) {

            // we check whether the slope changes are smooth using four previous points
            int i = 0;
            double slope[2], slope_if_jumped;
            double scatter;

#ifdef VERBOSE
            fprintf(stderr, "number of points so far %d\n", final->length);
#endif

            if (final->length < 3) {
                fprintf(stderr, "too few points in determing slope ratios and scatters: needs to increase the number of points\n");
                delete allSolutions;
                delete imageTracks;
                delete connected;
                delete Prov;
                delete Prov2;
                delete final;
                delete markedCurve;
                return (NULL);
            }

            for (_point *p = final->last;  i <= 1; i++, p = p->prev) {
                slope[i] = (p->prev->x2 - p->x2) / (p->prev->x1 -  p->x1);
#ifdef VERBOSE
                printf("slope: %f %f %d %f\n", p->x1, p->x2, i, slope[i]);
                printf("slope=%f\n", tan(p->thetaJ));
#endif
            }
            //      return(NULL);

            scatter = abs( log10(abs(slope[1] / slope[0])) );
            slope_if_jumped = (final->last->x2 - final->first->x2) / (final->last->x1 - final->first->x1);

            if (abs(log10(abs(slope[0] / slope_if_jumped))) < 2.0 * scatter) {
#ifdef VERBOSE
                printf("finally dropped %d\n", final->length);

                // for (_point *pp = final->last;  pp;  pp = pp->prev) {
                // printf("%f %f %e\n", pp->x1, pp->x2, pp->phi);
                // }
                //  return(NULL);
#endif

                if (final->length > 10) {
                    addPoint(final, final->first); // we simply connect the first and last points
                }
            } else {    // slopes changed too much, we should not jump; this is somewhat ambigious here
                continue;
            }
        }
        else {
#ifdef VERBOSE
            printf("continue the cycle\n");
#endif



            mindis = 1e9;
            for (_curve * Prov3 = imageTracks->first; Prov3; Prov3 = Prov3->next) {

                tempdis = *(final->first) - *(Prov3->first);
                if (tempdis <= mindis) {
                    mindis = tempdis;
                }

                tempdis = *(final->first) - *(Prov3->last);
                if (tempdis <= mindis) {
                    mindis = tempdis;
                }

                tempdis = *(final->last) - *(Prov3->first);
                if (tempdis <= mindis) {
                    mindis = tempdis;
                }

                tempdis = *(final->last) - *(Prov3->last);
                if (tempdis <= mindis) {
                    mindis = tempdis;
                }
            }
            if (mindis > SD) {
#ifdef VERBOSE
                fprintf(stderr, "No more track will be in final, start a new one and continue\n");
#endif
                // if (final->parity == 0) {
                //     if (final->first->mu > 0) {
                //         final->parity = 1;
                //     }
                //     else {
                //         final->parity = -1;
                //     }
                // }
                _final_parity
                connected->append(final);
                if (imageTracks->length != 0) {
                    // final = newSegment(imageTracks);
                    // imageTracks->drop(imageTracks->first, 1);
                    _newSegment
                }
            } else {
                continue;   // not yet connected, we should continue the cycle
            }

        }

        if (!final->parity) {
            if (final->first->mu > 0) {
                final->parity = 1;
            }
            else {
                final->parity = -1;
            }
        }

        connected->append(final);

        if (imageTracks->length == 0) {
            delete allSolutions;
            delete imageTracks;
            delete Prov;
            delete Prov2;
            delete markedCurve;
            return (connected); // we are done
        } else {      // start a new curve
            // final = new _curve;
            // final = newSegment(imageTracks);
            // // temp = imageTracks->first;
            // imageTracks->drop(imageTracks->first); // we already incorporated the new curve, and so drop it

            _newSegment

        }
    }

//check whether the first and last points are connected, if not, join the first and last points
#ifdef VERBOSE
    printf("joining the first and last point %f %f %f %f %f %f\n", final->first->x1, final->first->x2, final->first->phi,
           final->last->x1, final->last->x2, final->last->phi);
#endif

    SD = *(final->first) - *(final->last);
    final->append(final->first->x1, final->first->x2);


    if (connected->last != final) {
        if (!final->parity) {
            if (final->first->mu > 0) {
                final->parity = 1;
            }
            else {
                final->parity = -1;
            }
        }

        // this should comment out: on 2020.03.12 ?
        connected->append(final); // attach the final curve

    }
    printf("about to return\n");
    // delete allSolutions;
    // delete imageTracks;
    return (connected);
}

#undef TRUE_IMAGE
#undef FALSE_IMAGE



// head, tail, markedCurve = jump_over_caustics(final, imageTracks)
_curve *jump_over_caustics(_curve * final, _sols * imageTracks, int *head, int *tail, _curve * markedCurve, bool checkfinalhead)
{
    int jumpOverCaustics = 0;
    *tail = 0;//add on 2020.04.05
    *head = 0;//add on 2020.04.05
    double ratio, ratioMin = 1.0e10, muPrevious, disimg, mu;
    complex zs, dzs;
    complex zsPrevious, imgPrevious;
    // _curve *markedCurve;
    // int head = 0, tail = 0;


    if (checkfinalhead == false) {
        zsPrevious = final->last->zs;
        muPrevious = final->last->mu;
        imgPrevious = complex(final->last->x1, final->last->x2);
#ifdef VERBOSE
        printf("I am in jump_over_caustics(), need to jump over caustics at: %f %f %f, phi=%.15f, zs: %f, %f, imageTracks-length = %d\n", final->last->x1, final->last->x2,  muPrevious, final->last->phi * 180 / M_PI, zsPrevious.re, zsPrevious.im , imageTracks->length);
#endif


    } else {
        // check the head of final instead of tail of final. // 2020.07.17
        zsPrevious = final->first->zs;
        muPrevious = final->first->mu;
        imgPrevious = complex(final->first->x1, final->first->x2);
#ifdef VERBOSE
        printf("I am in jump_over_caustics(), need to jump over caustics at first point of final: %f %f %f, phi=%.15f, zs: %f, %f, imageTracks-length = %d\n", final->first->x1, final->first->x2,  muPrevious, final->first->phi * 180 / M_PI, zsPrevious.re, zsPrevious.im, imageTracks->length);
#endif
    }

    for (_curve *Prov2 = imageTracks->first; Prov2; Prov2 = Prov2->next) {
        // we first test whether the last end point connects with the head of the next curve
        zs = Prov2->first->zs;
        mu = Prov2->first->mu;
        dzs = zs - zsPrevious;
#ifdef VERBOSE
        printf("\t\t head x1 x2 mu: %f %f %f, zs: %f, %f\n", Prov2->first->x1, Prov2->first->x2, mu, zs.re, zs.im);
#endif
        disimg = abs(complex(Prov2->first->x1, Prov2->first->x2) - imgPrevious);
        if (abs(dzs) < EPS && muPrevious * mu < 0.0) { // connected with the head with opposite parity
            // if (abs(dzs) < EPS) { // connected with the head with opposite parity // remove on 2020.07.17
            ratio = abs( log10( abs(muPrevious / mu)) ) * disimg;
            //  printf("head - ratio, ratioMin=%f %f\n", ratio, ratioMin);
            if (ratio < ratioMin) {
                ratioMin = ratio;
                markedCurve = Prov2;
                // fprintf(stderr, "head???\n");
                *head = 1;   // head has the minimum
                *tail = 0;
            }
            // continue;
        }

        // we now test whether the last end point connects with the tail of the next curve
        zs = Prov2->last->zs;
        mu = Prov2->last->mu;
        dzs = zs - zsPrevious;

#ifdef VERBOSE
        printf("\t\t tail x1 x2 mu: %f %f %f, zs: %f, %f\n", Prov2->last->x1, Prov2->last->x2, mu, zs.re, zs.im);
#endif

        disimg = abs(complex(Prov2->last->x1, Prov2->last->x2) - imgPrevious);
        if (abs(dzs) < EPS && muPrevious * mu < 0.0) { // connected with the tail // remove on 2020.07.17
            // if (abs(dzs) < EPS) { // connected with the tail
            ratio = abs( log10( abs(muPrevious / mu)) ) * disimg;
            //  printf("tail - ratio, ratioMin=%f %f\n", ratio, ratioMin);
            if (ratio < ratioMin) {
                ratioMin = ratio;
                markedCurve = Prov2;
                // fprintf(stderr, "head???\n");
                *tail = 1;   // tail has the minimum
                *head = 0;
            }
            // continue;
        }
    }

    return markedCurve;
}

_curve *connect_head_or_tail(int *ifcontinue, _curve * final, _sols * imageTracks, _sols * connected, bool checkfinalhead)
{
    *ifcontinue = 0;
    int connectWithHead = 0, connectWithTail = 0;
    complex zsPrevious, z, zhead, ztail, dz;
    double muPrevious, mindis, tempdis, SD;
    _curve *Prov2, *Prov;
    if (checkfinalhead == false) {

        // zsPrevious = final->last->zs;
        // muPrevious = final->last->mu;
        z = complex(final->last->x1, final->last->x2); // last point in the last segment
        Prov2 = imageTracks->first;
        // see whether its image position overlaps with the head or tail of another segment, if yes, we connect it
        while (Prov2) {
            // fprintf(stderr, "can you here 5, imageTracks->length = %d\n", imageTracks->length);
            // see whether the image position overlaps with the head position
            zhead = complex(Prov2->first->x1, Prov2->first->x2);
            //add on 2020.04.05:
            z = complex(final->last->x1, final->last->x2); // last point in the last segment
            dz = zhead - z;
            if (abs(dz) < EPS) {
                connectWithHead = 1;

                // for (_point *p = Prov2->first; p; p = p->next)  addPoint(final, p);
                imageTracks->drop(Prov2); // we are already connected, drop this curve from the linked list

                Prov = Prov2->next;

                // 2020.07.17
                _final_parity
                final->last->next = Prov2->first;
                Prov2->first->prev = final->last;
                final->last = Prov2->last;
                final->length += Prov2->length;


                // delete Prov2; //add on 2020.03.13 // remove on 2020.07.17
                Prov2 = Prov;

                // continue;
            }
            // fprintf(stderr, "can you here 5.1, imageTracks->length = %d\n", imageTracks->length);

            if (Prov2) {
                // fprintf(stderr, "can you here 5.2, imageTracks->length = %d\n", imageTracks->length);
                // see whether the image position overlaps with the tail image position
                ztail = complex(Prov2->last->x1, Prov2->last->x2);

                //add on 2020.04.05:
                z = complex(final->last->x1, final->last->x2); // last point in the last segment

                dz = ztail - z;
                if (abs(dz) < EPS) {
                    connectWithTail = 1;

                    Prov = Prov2->next;
                    // for (_point *p = Prov2->last; p; p = p->prev)  addPoint(final, p); // attach the points in reverse order

                    imageTracks->drop(Prov2); // we are already connected, drop this curve from the linked list

                    Prov2->reverse();

                    _final_parity
                    // fprintf(stderr, "final->parity %d, final->length %d\n", final->parity, final->length);
                    // fprintf(stderr, "can you here 5.22, imageTracks->length = %d\n", imageTracks->length);
                    final->last->next = Prov2->first;
                    Prov2->first->prev = final->last;

                    final->last = Prov2->last;
                    final->length += Prov2->length;


                    // delete Prov2; //2020.03.13 // remove on 2020.07.17
                    Prov2 = Prov;
                    // continue;
                } else {
                    Prov2 = Prov2->next;
                }
                // fprintf(stderr, "can you here 5.3, imageTracks->length = %d\n", imageTracks->length);

            }

        }

    } else {
        // check head of final instead of check tail of final

        // zsPrevious = final->first->zs;
        // muPrevious = final->first->mu;
        z = complex(final->first->x1, final->first->x2); // first point in the last segment
        Prov2 = imageTracks->first;
        // see whether its image position overlaps with the head or tail of another segment, if yes, we connect it
        while (Prov2) {
            // see whether the image position overlaps with the head position
            zhead = complex(Prov2->first->x1, Prov2->first->x2);
            //add on 2020.04.05:
            z = complex(final->first->x1, final->first->x2); // last point in the last segment
            dz = zhead - z;
            if (abs(dz) < EPS) {
                connectWithHead = 1;

                // for (_point *p = Prov2->first; p; p = p->next)  addPoint(final, p);

                imageTracks->drop(Prov2); // we are already connected, drop this curve from the linked list
                Prov = Prov2->next;

                // 2020.07.17
                Prov2->reverse();
                _final_parity
                Prov2->parity = final->parity;

                // final->parity = 0;
                // _Prov2_parity
                Prov2->last->next = final->first;
                final->first->prev = Prov2->last;
                Prov2->last = final->last;
                Prov2->length += final->length;

                final = Prov2;

                // delete Prov2; //add on 2020.03.13 // remove on 2020.07.17
                Prov2 = Prov;

                // continue;
            }

            if (Prov2) {
                // see whether the image position overlaps with the tail image position
                ztail = complex(Prov2->last->x1, Prov2->last->x2);
                //add on 2020.04.05:
                z = complex(final->first->x1, final->first->x2); // last point in the last segment
                dz = ztail - z;
                if (abs(dz) < EPS) {
                    connectWithTail = 1;
                    Prov = Prov2->next;
                    // for (_point *p = Prov2->last; p; p = p->prev)  addPoint(final, p); // attach the points in reverse order

                    imageTracks->drop(Prov2); // we are already connected, drop this curve from the linked list

                    _final_parity
                    Prov2->parity = final->parity;

                    // final->parity = 0;
                    // _Prov2_parity
                    Prov2->last->next = final->first;
                    final->first->prev = Prov2->last;
                    Prov2->last = final->last;
                    Prov2->length += final->length;
                    final = Prov2;

                    // delete Prov2; //2020.03.13 // remove on 2020.07.17
                    Prov2 = Prov;

                    // continue;
                } else {
                    Prov2 = Prov2->next;
                }
            }

            // Prov2 = Prov2->next;
        }

    }


    // if ((connectWithHead || connectWithTail) ) {
    if ((connectWithHead || connectWithTail) && final->length > 2 ) {
        SD = *(final->first) - *(final->last);
#ifdef VERBOSE
        fprintf(stderr, "I am in connectWithHead = %d or connectWithTail = %d, imageTracks.length=%d,\n first and last image position difference:%f, final->length = %d\n", connectWithHead, connectWithTail , imageTracks->length, sqrt(SD), final->length);
#endif
        if (SD < EPS * EPS) {
            // if (final->parity == 0) {
            //     if (final->first->mu > 0) {
            //         final->parity = 1;
            //     }
            //     else {
            //         final->parity = -1;
            //     }
            // }
            _final_parity
#ifdef VERBOSE
            fprintf(stderr, "\n\nI am in connectWithHead or connectWithTail, add a new closed curve\n\n");
#endif
            connected->append(final);
            if (imageTracks->length != 0) {

                // final = newSegment(imageTracks);
                // imageTracks->drop(imageTracks->first, 1);

                _newSegment


            }
        } else {


//             mindis = 1e9;
//             for (_curve * Prov3 = imageTracks->first; Prov3; Prov3 = Prov3->next) {

//                 tempdis = *(final->first) - *(Prov3->first);
//                 if (tempdis <= mindis) {
//                     mindis = tempdis;
//                 }

//                 tempdis = *(final->first) - *(Prov3->last);
//                 if (tempdis <= mindis) {
//                     mindis = tempdis;
//                 }

//                 tempdis = *(final->last) - *(Prov3->first);
//                 if (tempdis <= mindis) {
//                     mindis = tempdis;
//                 }

//                 tempdis = *(final->last) - *(Prov3->last);
//                 if (tempdis <= mindis) {
//                     mindis = tempdis;
//                 }
//             }
//             if (mindis > SD) {
// #ifdef VERBOSE
//                 fprintf(stderr, "No more track will be in final, start a new one and continue\n");
// #endif
//                 // if (final->parity == 0) {
//                 //     if (final->first->mu > 0) {
//                 //         final->parity = 1;
//                 //     }
//                 //     else {
//                 //         final->parity = -1;
//                 //     }
//                 // }
//                 _final_parity
//                 connected->append(final);
//                 if (imageTracks->length != 0) {
//                     // final = newSegment(imageTracks);
//                     // imageTracks->drop(imageTracks->first, 1);
//                 _newSegment
//                 }
//             }



        }
        *ifcontinue = 1;
    } else {
#ifdef VERBOSE
        fprintf(stderr, "I am in connect_head_or_tail, imageTracks->length = %d, connectWithHead = %d, connectWithTail = %d \n", imageTracks->length, connectWithHead, connectWithTail);
#endif
    }
    return final;

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
int TripleLensing::trueSolution(double xs, double ys, complex z)
{
    // complex zs;
    // complex dzs;
    // complex dzsdz;

    flag = 0; //
    // double x, y;
    // double dx, dy;
    // complex tempzc, tempzc2, tempzc3, tempzc4, J1c;
    // double Jxx = 1.0, Jyy = 0.0, Jxy = 0.0, Jyx; // Jxx + Jyy = 2.0, Jxy = Jyx
    Jxx = 1.0;
    Jyy = 0.0;
    Jxy = 0.0;
    // Jyx; // Jxx + Jyy = 2.0, Jxy = Jyx
    // double sum2 = 0.0, sq = 0.0;
    // static double TINY = 1.0e-20;

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
        xy_j2 = 1 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        dzs = complex(rep, imp);

        break;

    case 3:
        // triple lens
        // binary lens

        rep = xs - x;// real part
        imp = ys - y;
        x_xj = x - zlens[0].re;
        y_yj = y - zlens[0].im;
        xy_j2 = 1 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        x_xj = x - zlens[2].re;
        y_yj = y - zlens[2].im;
        xy_j2 = 1 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[2] * x_xj * xy_j2;
        imp += mlens[2] * y_yj * xy_j2;

        dzs = complex(rep, imp);
        break;

    default:

        zs = complex(xs, ys);
        // dzs = zs - ( z - mlens[0] / conj(z - zlens[0]) - mlens[1] / conj(z - zlens[1]) - mlens[2] / conj(z - zlens[2]) );
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

    if (abs(dzs) < EPS)  {
        flag = 1;

        // x = real(z);
        // y = imag(z);

        for (int i = 0; i < NLENS; i++) {
            dx_db = x - zlens[i].re;
            dy_db = y - zlens[i].im;
            dx_db2 = dx_db * dx_db;
            dy_db2 = dy_db * dy_db;

            r2_1 = dx_db2 + dy_db2 + TINY; // avoid zero in the denominator
            r2_2 = 1 / (r2_1 * r2_1);
            Jxx += mlens[i] * (dx_db2 - dy_db2) * r2_2;
            Jxy += 2.0 * mlens[i] * dx_db * dy_db * r2_2;

            // add on 2020.04.09
#ifdef parabcorr
            // tempzc = (z - zlens[i]);
            // tempzc2 = tempzc * tempzc;
            // tempzc3 = tempzc2 * tempzc;
            // J1 = J1 + mlens[i] / tempzc2;
            // J2 = J2 - 2 * ( mlens[i] / tempzc3 );

            tempzc = complex(z - zlens[i]);
            tempzc2 = tempzc * tempzc;
            tempzc3 = tempzc2 * tempzc;
            J1 = J1 + mlens[i] / tempzc2;
            J2 = J2 - 2 * ( mlens[i] / tempzc3 );


#endif

        }
#ifdef parabcorr
        // if (qtest){
        J1c = conj(J1);
        dJ = 1 - (J1) * J1c;
        // J2c = conj(J2);
        // }
#endif

        //
        //analytical results for the other components of the Jacobian
        //
        Jyy = 2.0 - Jxx;
        mu = 1.0 / (Jxx * Jyy - Jxy * Jxy);


// #ifdef CALCULATE_ANGLE
//         sum2 = (Jxx + Jyy) / 2.0;
//         sq = sqrt(sum2 * sum2 - 1.0 / (*mu));

//         *lambda1 = sum2 + sq;
//         *lambda2 = sum2 - sq;

//         *thetaJ = 0.5 * atan(2.0 * Jxy / (Jyy - Jxx + TINY)); // avoid zero in the denomintor
//         if (*thetaJ < 0.0) *thetaJ += M_PI / 2.0;
// #endif

    }

    return (flag);
}



int TripleLensing::trueSolution_qtest(double xs, double ys, complex z)
{
    // complex zs;
    // complex dzs;
    // complex dzsdz;

    flag = 0; //
    // double x, y;
    // double dx, dy;
    // complex tempzc, tempzc2, tempzc3, tempzc4, J1c;
    // double Jxx = 1.0, Jyy = 0.0, Jxy = 0.0, Jyx; // Jxx + Jyy = 2.0, Jxy = Jyx
    Jxx = 1.0;
    Jyy = 0.0;
    Jxy = 0.0;
    // Jyx; // Jxx + Jyy = 2.0, Jxy = Jyx
    // double sum2 = 0.0, sq = 0.0;
    // static double TINY = 1.0e-20;

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
        xy_j2 = 1 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        dzs = complex(rep, imp);

        break;

    case 3:
        // triple lens
        // binary lens

        rep = xs - x;// real part
        imp = ys - y;
        x_xj = x - zlens[0].re;
        y_yj = y - zlens[0].im;
        xy_j2 = 1 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[0] * x_xj * xy_j2;
        imp += mlens[0] * y_yj * xy_j2;

        x_xj = x - zlens[1].re;
        y_yj = y - zlens[1].im;
        xy_j2 = 1 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[1] * x_xj * xy_j2;
        imp += mlens[1] * y_yj * xy_j2;

        x_xj = x - zlens[2].re;
        y_yj = y - zlens[2].im;
        xy_j2 = 1 / (x_xj * x_xj + y_yj * y_yj);
        rep += mlens[2] * x_xj * xy_j2;
        imp += mlens[2] * y_yj * xy_j2;

        dzs = complex(rep, imp);
        break;

    default:

        zs = complex(xs, ys);
        // dzs = zs - ( z - mlens[0] / conj(z - zlens[0]) - mlens[1] / conj(z - zlens[1]) - mlens[2] / conj(z - zlens[2]) );
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

    if (abs(dzs) < EPS)  {
        flag = 1;

        // x = real(z);
        // y = imag(z);

        for (int i = 0; i < NLENS; i++) {
            dx_db = x - zlens[i].re;
            dy_db = y - zlens[i].im;
            dx_db2 = dx_db * dx_db;
            dy_db2 = dy_db * dy_db;

            r2_1 = dx_db2 + dy_db2 + TINY; // avoid zero in the denominator
            r2_2 = 1 / (r2_1 * r2_1);
            Jxx += mlens[i] * (dx_db2 - dy_db2) * r2_2;
            Jxy += 2.0 * mlens[i] * dx_db * dy_db * r2_2;

            // add on 2020.04.09
// #ifdef parabcorr
            // if (qtest){
            // tempzc = (z - zlens[i]);
            // tempzc2 = tempzc * tempzc;
            // tempzc3 = tempzc2 * tempzc;
            // tempzc4 = tempzc2 * tempzc2;
            // J1 = J1 + mlens[i] / tempzc2;
            // J2 = J2 - 2 * ( mlens[i] / tempzc3 );
            // J3 = J3 + 6 * ( mlens[i] / tempzc4 );

            tempzc = complex(z - zlens[i]);
            tempzc2 = tempzc * tempzc;
            tempzc3 = tempzc2 * tempzc;
            J1 = J1 + mlens[i] / tempzc2;
            J2 = J2 - 2 * ( mlens[i] / tempzc3 );
            J3 = J3 + 6 * ( mlens[i] / tempzc4 ); 


            // }
// #endif

        }
// #ifdef parabcorr
        // if (qtest){
        J1c = conj(J1);
        dJ = 1 - (J1) * J1c;
        // J2c = conj(J2);
        // }
// #endif

        //
        //analytical results for the other components of the Jacobian
        //
        Jyy = 2.0 - Jxx;
        mu = 1.0 / (Jxx * Jyy - Jxy * Jxy);


// #ifdef CALCULATE_ANGLE
//         sum2 = (Jxx + Jyy) / 2.0;
//         sq = sqrt(sum2 * sum2 - 1.0 / (*mu));

//         *lambda1 = sum2 + sq;
//         *lambda2 = sum2 - sq;

//         *thetaJ = 0.5 * atan(2.0 * Jxy / (Jyy - Jxx + TINY)); // avoid zero in the denomintor
//         if (*thetaJ < 0.0) *thetaJ += M_PI / 2.0;
// #endif

    }

    return (flag);
}


// obtain the polynomical coefficients
// input: lens mass and positions, mlens[], zlens[], source positions: xs, ys
// output: c
// void polynomialCoefficients(double mlens[], complex zlens[], double xs, double ys, complex c[])
// {
//     int i, j, k;
//     int n;

//     complex zs, zsc;      /* source position and conjugate */
//     complex zc[NLENS];        /* conjugate of lens positions */

//     complex p[NLENS + 1][NLENS], q[NLENS + 1][NLENS];
//     complex temp[NLENS + 1];
//     complex ctemp[DEGREE + 1];
//     complex qtemp[DEGREE + 1], qtemp2[NLENS + 1];
//     complex ptemp[DEGREE + 1];
//     int degrees;

//     n = NLENS;
//     zs = complex(xs, ys);
//     zsc = conj(zs);

//     for (i = 0; i < n; i++) {
//         zc[i] = conj(zlens[i]);
//     }

//     for (i = 0; i < n; i++) {
//         temp[0] = 1.0;
//         for (j = 0; j < n; j++) {
//             multiply_z(temp, zlens[j], j);
//         }

//         /* numerator, m[i] * (z-z1)*(z-z2) ... (z-zn) */
//         for (j = 0; j <= n; j++) {
//             p[j][i] = mlens[i] * temp[j];
//         }

//         /* denominator */
//         for (j = 0; j <= n; j++) { /* (zsc-conjugate(z[i])) * Product_{j=1}^N (z-z[j]) */
//             q[j][i] = (zsc - zc[i]) * temp[j];
//         }

//         /* Sum_{j=1}^n Product (z-z_k), k=1, n, but k !=j. */
//         for (j = 0; j < n; j++) {
//             /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */
//             temp[0] = 1.0;
//             degrees = 0;
//             for (k = 0; k < n; k++) {
//                 if (k == j) continue;

//                 multiply_z(temp, zlens[k], degrees);
//                 degrees++;
//             }

//             /* doing the sum */
//             for (k = 0; k < n; k++) {
//                 q[k][i] = q[k][i] + mlens[j] * temp[k];
//             }
//         }
//     }

//     /* now clear the fractions, z-zs - Sum_i p_i/q_i */
//     /* get the polynomial Product q_i, i=1, n */

//     /* first term */
//     qtemp[0] = 1.0;
//     degrees = 0;
//     for (i = 0; i < n; i++) {
//         for (j = 0; j <= n; j++) {
//             qtemp2[j] = q[j][i];
//         }

//         multiply(qtemp, degrees, qtemp2, n, ctemp);

//         degrees += n;


//         for (j = 0; j <= degrees; j++) {
//             qtemp[j] = ctemp[j];
//         }
//     }

//     /* get coefficients of (z-zs) Product_i=1^n q_i */
//     multiply_z(ctemp, zs, degrees);

//     /* copy the coefficients */
//     for (i = 0; i < n * n + 2; i++) {
//         c[i] = ctemp[i];
//     }

//     /* second term */
//     for (i = 0; i < n; i++) {
//         degrees = 0;
//         qtemp[0] = 1.0;
//         for (j = 0; j < n; j++) {
//             if (j == i) continue;

//             for (k = 0; k <= n; k++) {
//                 qtemp2[k] = q[k][j];
//             }

//             multiply(qtemp, degrees, qtemp2, n, ctemp);

//             degrees += n;

//             for (k = 0; k <= degrees; k++) {
//                 qtemp[k] = ctemp[k];
//             }
//         }

//         for (k = 0; k <= n; k++) {
//             ptemp[k] = p[k][i];
//         }

//         multiply(qtemp, degrees, ptemp, n, ctemp);

//         for (k = 0; k < n * n + 1; k++) {
//             c[k] = c[k] - ctemp[k];
//         }
//     }
// }



// /* multiply a polynomial of degree n by (z-a), returns a polynomial of degree n+1 */
// void multiply_z(complex c[], complex a, int n)
// {
//     int j;

//     c[n + 1] = c[n];
//     for (j = n; j >= 1; j--) c[j] = c[j - 1] - c[j] * a;
//     c[0] = c[0] * (-a);
// }

// /* multiply a polynomial of degree n by (z-a)^2, returns a polynomial of degree n+2 */
// void multiply_zSquare(complex c[], complex a, int n)
// {
//     multiply_z(c, a, n);
//     multiply_z(c, a, n + 1);
// }

// /* multiply two polynomials of degree of na and nb, input as a(na), b(nb) */
// void multiply(complex a[], int na, complex b[], int nb, complex c[])
// {
//     int i, j;

//     /* zero the array */
//     for (i = 0; i <= na + nb; i++) {
//         c[i] = 0.0;
//     }

//     /* now do the product */
//     for (i = 0; i <= na; i++) {
//         for (j = 0; j <= nb; j++) {
//             c[i + j] = c[i + j] + a[i] * b[j];
//         }
//     }
// }

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
    int i = 0;
    int itrack = 0;
    FILE *fp;

    fp = fopen("data/allImages.dat", "w");
    printf("entered printAllTracks\n");


    for (_curve *c = track->first; c; c = c->next) { //curves
        for (_point *p = c->first; p; p = p->next) { //point
            //      printf("%f %f %f %f #%d point in the %d-th image track\n", p->x1, p->x2, p->phi, p->mu, i, itrack);
            fprintf(fp, "%f %f %f %f #%d point in the %d-th image track\n", p->x1, p->x2, p->phi, p->mu, i, itrack);
            i++;
        }
        itrack++;
    }

    fclose(fp);
}

void saveTracks_before(_sols * track)
{
    int i = 0;
    int itrack = 0;
    FILE *fp;

    fp = fopen("data/allTracks_before.dat", "w");
    // fp = fopen(filename, "w");
    // printf("entered saveTracks\n\n\n");

    for (_curve *c = track->first; c; c = c->next) { //curves
        fprintf(fp, "%d ", c->length);
    }
    fprintf(fp, "\n");

    for (_curve *c = track->first; c; c = c->next) { //curves
        for (_point *p = c->first; p; p = p->next) { //point
            fprintf(fp, "%f %f %f %f ", p->x1, p->x2, p->phi, p->mu);
            i++;
        }
        itrack++;
    }
    fclose(fp);
}

void saveTracks(_sols * track)
{
    int i = 0;
    int itrack = 0;
    FILE *fp;

    fp = fopen("data/allTracks.dat", "w");
    // fp = fopen(filename, "w");
    // printf("entered saveTracks\n\n\n");

    for (_curve *c = track->first; c; c = c->next) { //curves
        fprintf(fp, "%d ", c->length);
    }
    fprintf(fp, "\n");

    for (_curve *c = track->first; c; c = c->next) { //curves
        for (_point *p = c->first; p; p = p->next) { //point
            fprintf(fp, "%f %f %f %f ", p->x1, p->x2, p->phi, p->mu);
            i++;
        }
        itrack++;
    }
    fclose(fp);
}




void TripleLensing::outputImagesTriple(double xsCenter, double ysCenter, int nphi, double rs)
{
    _sols *imageTracks;
    // imageTracks = outputTracks(mlens, zlens, xsCenter, ysCenter, rs, nphi, 1.5);
    imageTracks = outputTracks(xsCenter, ysCenter, rs, nphi, 1.5);
    //we failed in returning image tracks
    if (!imageTracks)  return;
#ifdef VERBOSE
    printf("I am in outputImagesTriple, Number of image tracks: %d\n",  imageTracks->length);
#endif
    // Successful, print out image tracks
    FILE *fileImage;
    fileImage = fopen("data/images.dat", "w");
    fprintf(fileImage, "%f %f %f\n", xsCenter, ysCenter, rs);
    fprintf(fileImage, "%d\n", imageTracks->length);

    for (_curve *c = imageTracks->first; c; c = c->next) { //curves
        fprintf(fileImage, "%d\n", c->length);
        for (_point *p = c->first; p; p = p->next) { //point
            fprintf(fileImage, "%f %f %f %f \n", p->x1, p->x2, p->phi, p->mu);
        }
    }
    fclose(fileImage);
}

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
#ifdef VERBOSE
    printf("I am in outputCriticalTriple, Number of closed critical curves: %d\n", ncritical);
#endif

    // write out the critical curves and caustics separately
    // allxys:
    // num_critical_points
    // then critx_i, crity_i
    // numb_caustics_points
    //then causticsx_i, causticsy_i
    // [2, crx1,cry1,crx2,cry2, 1, cax1, cay1]
    //idx [0, 1,   2,   3,   4,    5, 6,    7]


    ncurves = 0;
    int count_critical = 0;
    int count_caustic = 0;
    for (_curve *c = criticalCurves->first; c; c = c->next) {
        int npoints = 0;
        ncurves++;

        // second half, caustics
        if (ncurves > ncritical) {      // second halfs are caustics
            for (_point *p = c->first; p; p = p->next) {
                npoints++;

                count_caustic ++;
                allxys[2 * count_critical + 1 + 2 * count_caustic - 1] = p->x1;
                allxys[2 * count_critical + 1 + 2 * count_caustic] = p->x2;
            }
        } else {
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

// #define _addsubarea \
// subArea +=  ( (y2 + y1) * (x2 - x1) - (y2 - y1) * (x2 + x1) );


#define _addsubarea_parab \
subArea +=  (0.5*( (y2 + y1) * (x2 - x1) ) + one_24*(ds1+ds2)*dphi3);\
// fprintf(stderr, "\tparabcorr %.15e, phi1, phi2 = %f %f, ds1, ds2 = %f, %f, dphi = %.3e\n", one_24*(ds1+ds2)*dphi3,phi1, phi2, ds1, ds2, dphi);

#define _addsubarea_parab2 \
subArea +=  (0.5*( (y2 + y1) * (x2 - x1) ) - one_24*(ds1+ds2)*dphi3);\
// fprintf(stderr, "\tparabcorr %.15e, phi1, phi2 = %f %f, ds1, ds2 = %f, %f, dphi = %.3e\n", -one_24*(ds1+ds2)*dphi3,phi1, phi2, ds1, ds2, dphi);


// double areaFunc(_sols * track, double rs,  int paracorr)
double TripleLensing::areaFunc(_sols * track, double rs)
{
    // double x1, x2, y1, y2, phi1, phi2, dphi, ds1, ds2;
    area = 0.0;
    pos_area = 0.0;
    // int parity,
    trackcnt = 0;
    subArea = 0.0;
    // double subareas[track->length];

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

        subArea = 0.0;

        x1 = c->first->x1;
        y1 = c->first->x2;
        phi1 = c->first->phi;

        // c->parity = c->first->mu > 0 ? 1 : -1;

        _c_parity;
        pscan = c->first;
        while ( pscan->next && pscan->next->next && ((pscan->mu * pscan->next->mu < 0) || ( abs(abs(pscan->phi - pscan->next->phi) - M_PI2) < EPS ) || ( abs(abs(pscan->next->phi - pscan->next->next->phi) - M_PI2) < EPS )  )) {
            pscan = pscan->next;
        }
        ph1 = pscan->phi; ph2 = pscan->next->phi; ph3 = pscan->next->next->phi;
        if ( ph1 > ph2 && ph2 > ph3) {
            c->parity *= -1;
        }
   

        parity = c->parity;
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
        fprintf(stderr, "track->length %d, c->length: %d, c->parity = %d, c->first->mu = %f, subArea * parity: %f\n", track->length, c->length, c->parity, c->first->mu , subArea * parity * 0.25);
#endif
    }
    // somtimes parity may not work
    area = abs(area) * 0.25;
    if (area < areaSource) {
        return pos_area * 0.25;
    } else {
        return area;
    }
}


// #define _dphijump \
dphi = 0;

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
    // double subareas[track->length];
    double tmpx1;

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
        // fprintf(stderr, "\t\t\t firtt mu %f\n", c->first->mu);
        // fprintf(stderr, "\t\t\t last mu %f\n", c->last->mu);
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

        // c->parity = c->first->mu>0? 1:-1;
        _c_parity;

        pscan = c->first;
        while ( pscan->next && pscan->next->next && ((pscan->mu * pscan->next->mu < 0) || ( abs(abs(pscan->phi - pscan->next->phi) - M_PI2) < 1e-10 ) || ( abs(abs(pscan->next->phi - pscan->next->next->phi) - M_PI2) < 1e-10 )  )) {
            pscan = pscan->next;
        }
        // fprintf(stderr, "\t\t\t phi1 phi2 phi3 = %f %f %f\n", ph1*180/M_PI, ph2*180/M_PI, ph3*180/M_PI );

        ph1 = pscan->phi; ph2 = pscan->next->phi; ph3 = pscan->next->next->phi;
        // ph1 = pscan->phi > 0 ? pscan->phi : M_PI2 + pscan->phi ;
        // ph2 = pscan->next->phi > 0 ? pscan->next->phi : M_PI2 + pscan->next->phi;
        // ph3 = pscan->next->next->phi > 0 ? pscan->next->next->phi : M_PI2 + pscan->next->next->phi;



        // if( (c->parity == 1) && (ph1<ph2) && (ph2<ph3)){
        //     c->parity = -1;
        // }else if( (c->parity == -1)  && ph1>ph2 && ph2>ph3){
        //     c->parity = 1;
        // }

        if ( ph1 > ph2 && ph2 > ph3) {
            c->parity *= -1;
        }
#ifdef VERBOSE
        fprintf(stderr, "ph1 ph2 ph3 %f %f %f, ph2 - ph3 = %f , c->parity %d\n", ph1 * 180 / M_PI, ph2 * 180 / M_PI, ph3 * 180 / M_PI, (abs(ph2 - ph3) - M_PI2)* 180 / M_PI, c->parity);
#endif
// if( (c->parity == 1) && (ph1>ph2) && (ph2>ph3)){
//             c->parity = -1;
//         }else if( (c->parity == -1)  && ph1<ph2 && ph2<ph3){
//             c->parity = 1;
//         }

// if( (ph1>ph2) && (ph2>ph3)){
//             c->parity = -1;
//         }else if( ph1<ph2 && ph2<ph3){
//             c->parity = 1;
//         }


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
                // ds1 *= -1;

                // tmpx1 = x1;
                // x1 = x2;
                // x2 = tmpx1;

                // dphi = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) )/sqrt(abs(dz1*dz2));
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
                // ds2 *= -1;

                // dphi = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) )/sqrt(abs(dz1*dz2));
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
                // blow is dphi^4 !!! not dphi3!!!
                // dphi3 = (phi2 - phi1);
                // dphi3 *= dphi3;
                // dphi3 *= dphi3;


                dphi = (phi2 - phi1);
                if ( abs( abs(dphi) - M_PI2) < M_PI ) {
                    // dphi = 0;
                    dphi = dphi + (dphi > 0 ? -1 : 1) * M_PI2;
                    // fprintf(stderr, "\t\t place 2 dphi = %f\n", dphi);
                    // dphi = (dphi>0?1:-1 )* abs(abs(dphi) - M_PI2);
                }
                // else{
                    // fprintf(stderr, "\t\t place 1 dphi = %f\n", dphi);
                // }

                dphi3 = abs(dphi);
                dphi3 *= dphi;
                // dphi3 *= -dphi; // why negative?
                dphi3 *= dphi; // 

                _addsubarea_parab
                x1 = x2; y1 = y2;
                ds1 = ds2;
                // ds1 = p->ds;
                phi1 = phi2;
                clpairparity1 = clpairparity2;
            }


            // fprintf(stderr, "\t\t\t phi1 = %f, phi2 = %f\n", phi1*180/M_PI, phi2*180/M_PI);


            // subArea += 0.25 * ( (y2 + y1) * (x2 - x1) - (y2 - y1) * (x2 + x1) );

            // subArea +=  ( (y2 + y1) * (x2 - x1) - (y2 - y1) * (x2 + x1) );


// #ifdef parabcorr
//             if (paracorr) {
//                 phi2 = p->phi;
//                 ds2 = p->ds;
//                 dphi = phi2 - phi1;
//                 subArea += one_24 * dphi * dphi * dphi * (ds1 + ds2);

//                 phi1 = phi2;
//                 ds1 = ds2;
//             }
// #endif

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
            // ds2 *= -1;
            ds1 *= -1;

            // dphi = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) )/sqrt(abs(dz1*dz2));
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
            // ds1 *= -1;
            ds2 *= -1;
            // dphi = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) )/sqrt(abs(dz1*dz2));
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


            // dphi = ( phi2 + M_PI2 - phi1);
            // fprintf(stderr, "\t\t phi1 = %f, phi2 = %f, dphi3 = %f\n", phi1*180/M_PI, phi2*180/M_PI, dphi3*180/M_PI);
            dphi = (phi2 - phi1);

            if ( abs( dphi + M_PI2 ) < M_PI ) {
                // dphi = 0;
                dphi = (dphi > 0 ? 1 : -1 ) * abs(dphi + M_PI2);
            }

            dphi3 = abs(dphi);
            dphi3 *= dphi;
            // dphi3 *= -dphi; // why negative??
            dphi3 *= dphi;

            _addsubarea_parab

        }






        // subArea += 0.25 * ( (y2 + y1) * (x2 - x1) - (y2 - y1) * (x2 + x1) );

        // subArea +=  ( (y2 + y1) * (x2 - x1) - (y2 - y1) * (x2 + x1) );
        // _addsubarea_parab

// #ifdef parabcorr
//         if (paracorr) {
//             phi1 = c->last->phi;
//             ds1 = c->last->ds;
//             phi2 = c->first->phi;
//             ds2 = c->first->ds;
//             if (phi2 < phi1) phi2 += 2 * M_PI;
//             dphi = phi2 - phi1;
//             subArea += one_24 * dphi * dphi * dphi * (ds1 + ds2);
//         }
// #endif
        // subareas[trackcnt] = subArea * parity;
        // trackcnt ++;
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



// _linkedarray *getphis_v3(double mlens[], complex zlens[], double xsCenter, double ysCenter, double rs, int nphi, int secnum, int basenum, int *distype, int maxmuidx) {
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
        // xs = xsCenter + rs*sin(scanphi); //bug
        ys = ysCenter + rs * sin(scanphi);
        // mus[i] = triplePS(mlens, zlens, xs, ys, &nimages);
        mus[i] = TriplePS(xs, ys);

        // mus[i] = tripleQuatrapoleTest(xs, ys, rs);
        // mus[i] *= quad_err;

        if (mus[i] <= minmu) {
            minmu = mus[i];
            maxmuidx = i;
        }
        if (mus[i] >= maxmu) {
            maxmu = mus[i];
        }

        scanphi += dphi;
    }

#ifdef VERBOSE
    fprintf(stderr, "secnum_priv %d,basenum %d,  maxmuidx: %d\n", secnum_priv, basenum, maxmuidx);
#endif
    for (int i = 0; i < secnum_priv; i++) {
        // mus[i] = log10(mus[i] / minmu);
        mus[i] = (mus[i] / minmu);
    }
    // double psf[3] = {0.5, 2, 0.5};
    // double psf[3] = {0.7, 1, 0.7};
    // double psf[3] = {0.8, 1.5, 0.8};
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
    PHI->linspace(offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, npmus[tempidx]*basenum, 0);
    for (int i = 1; i < secnum_priv - 1; i++) {
        tempidx = secnum_privlist[i];
        PHI->linspace(offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, npmus[tempidx]*basenum, 0);
    }
    tempidx = secnum_privlist[secnum_priv - 1];
    PHI->linspace(offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, npmus[tempidx]*basenum, 1);

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
    scan1 = first;
    for (int i = 0; i < length; i++) {
        scan2 = scan1->next;
        delete scan1;
        scan1 = scan2;
    }
    // scan1 = NULL;
    // scan2 = NULL;
    // delete scan1;
    // delete scan2;
}

// _thetas::~_thetas(void) {
//     _theta *scan, *scan2;
//     scan = first;
//     while (scan) {
//         scan2 = scan->next;
//         delete scan;
//         scan = scan2;
//     }
// }


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
    // 插值，相邻两点间插入该两点的平均值
    // PHI = PHI.repeat(2)[:-1]
    // PHI[1::2] += (PHI[2::2]-PHI[1::2])/2
    if (length > 0) {
        Node *scan;
        Node *cc;

        // wrong:
        // Node *temp;
        // for (int i =0; i<length-1;i++){
        //     temp = scan->next;
        //     cc = new Node( (scan->value+ temp->value)/2 );
        //     cc->next = temp;
        //     scan->next = cc;
        //     scan = scan->next;
        // }
        double midvalues[length - 1];
        scan = first;
        for (int i = 0; i < length - 1; i++) {
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

#undef EPS_CLOSE
#undef MAXIT_NEWTON

// // #define _LL (y-z)+coefs[21]/(zc-coefs[20])+coefs[22]/zc //Lens equation test
// complex _LL_Nlens(complex y, complex z, double mlens[], complex zlens[]) {
//     // y: source position
//     // z: image position
//     complex _LL_res;
//     _LL_res = y - z;
//     for (int i = 0; i < NLENS; i++) {
//         _LL_res = _LL_res + mlens[i] / conj( z - zlens[i]);
//     }
//     return _LL_res;
// }


// _linkedarray *thetas2phis(_thetas * thetas) {
//     _linkedarray *PHI;
//     PHI = new _linkedarray;
//     _theta *scan, *scan2;
//     scan = thetas->first;
//     _linkedarray *phis = new _linkedarray;
//     while (scan) {
//         PHI->append( scan->th );
//         scan2 = scan->next;
//         // delete scan;
//         scan = scan2;
//     }
//     return PHI;
// }


#define _toterr \
toterr = def_err;


// void TripleLensing::bisec(_curve *lightkv, _point *scan2, _point *scan3, double errTol, double errorTol_critical, double tE_inv, double u0, double salpha, double calpha, double t0, double mlens [], double zlens [], double rs, int nphi, int secnum, int basenum, double quaderr_Tol, int MAXNPS, double timerrTol) {
void TripleLensing::bisec(_curve * lightkv, _point * scan2, _point * scan3, double errorTol_critical, double tE_inv, double u0, double salpha, double calpha, double t0, double rs) {
    int np1 = 1, Xnp = 0;// np1 = 1, 3, 5
    double mag2, mag1, tm1, tm2, def_err, def_err3rd, def_err1d3d, midtm, midtm_mag, midtm_mag_inter, timerr, midtm_mag_lin, midtm_mag_inter3rd = 1, toterr;
    double subdt, curr_t, curr_mag, tn, y1, y2, dmag_dt_ratio = 1;
    subdt = (scan3->x1 - scan2->x1) / (np1 + 1);
    _point *subscan, *pnt2beAdd, *scan;
    subscan = scan2;
    for (int i3 = 0; i3 < np1; i3++) {
        curr_t = scan2->x1 + subdt * (i3 + 1);
        tn = (curr_t - t0) * tE_inv;
        y1 = u0 * salpha + tn * calpha; // Mine
        y2 = u0 * calpha - tn * salpha; // Mine
        // if (subscan->area){
        curr_mag = subscan->area;
        // }else{
        // curr_mag = tripleFS2python(mlens, zlens, y1, y2, rs, nphi, secnum, basenum, quaderr_Tol);
        // }
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
        midtm_mag = tripleFS2python( y1, y2, rs);
        scan->area = midtm_mag;
        midtm_mag_inter = (mag1 + mag2) * 0.5;
        midtm_mag_inter = pow(10, midtm_mag_inter) + 1;
        midtm_mag = pow(10, midtm_mag) + 1;
        def_err = fabs( (midtm_mag_inter - midtm_mag) / midtm_mag );
        _toterr
        scan->error = toterr;

        int testnp = 2;// 2, 4, 6
        for (int i4 = 0; i4 < 1; i4++) {
            midtm = tm1 + (tm2 - tm1) * (i4 + 1) / (testnp + 1);
            tn = (midtm - t0) * tE_inv;
            y1 = u0 * salpha + tn * calpha; //
            y2 = u0 * calpha - tn * salpha; //
            midtm_mag = tripleFS2python(y1, y2, rs);
            midtm_mag_inter = mag1 + (mag2 - mag1) * (i4 + 1) / (testnp + 1);
            midtm_mag_inter = pow(10, midtm_mag_inter) + 1;
            midtm_mag = pow(10, midtm_mag) + 1;
            def_err = fabs( (midtm_mag_inter - midtm_mag) / midtm_mag );
            _toterr
            // fprintf(stderr, "\t new scan->error %.3e, errTol = %.3e ,scan-length = %d\n", toterr,errTol, lightkv->length);
            if ( toterr  > scan->error) {
                scan->error =  toterr;
            }
        }


        scan = scan->next;
    }//while {
    scan = scan2;
    if (lightkv->length < MAXNPS) {
        while (scan != scan3) {
            if (scan->error >= adaerrTol) {

                timerr = 1e100;
                if (scan->prev && scan->next) {
                    timerr = fmax( fabs( (scan->x1 - scan->next->x1)),  fabs( (scan->x1 - scan->prev->x1)) );
                }

                if (timerr > timerrTol) {
                    // bisec(lightkv, scan, scan->next, errTol, errorTol_critical, tE_inv, u0, salpha, calpha, t0, mlens, zlens, rs, nphi, secnum, basenum, quaderr_Tol, MAXNPS, timerrTol);
                    bisec(lightkv, scan2, scan2->next, errorTol_critical, tE_inv, u0, salpha, calpha, t0, rs);
                }

            } else if ( scan->error > errorTol_critical ) {

            }

            scan = scan->next;
        }
    }

}

// void tri_ada_sample_1D( double mlens[], double zlens[], double rs, int nphi, int secnum, int basenum, double quaderr_Tol, _curve *lightkv, int np0, int MAXNPS, double errTol, double t_start, double t_end, double alpha, double tE, double t0 , double u0, int np, double timerrTol);
void TripleLensing::TripleLkvAdap( double rs, _curve * lightkv, double t_start, double t_end, double alpha, double tE, double t0 , double u0) {
    double dt0 = (t_end - t_start) / (adanp0 - 1);
    double curr_t = t_start, curr_mag, y1, y2, mag1, mag2, def_err, tm1, tm2, mag0, tm0, k1, k2, def_curerr, midtm_mag, midtm_mag_lin, midtm_mag_inter , midtm, toterr;

    double errorTol_critical = adaerrTol * 0.1;

    double salpha = sin(alpha), calpha = cos(alpha), tn, tE_inv = 1 / tE;
    _point *pnt2beAdd, *scan, *scan2, *err_l_pnt, *err_s_pnt;//error_s_pnt: smaller error pnt
    for (int i = 0; i < adanp0; i++) {
        curr_t = t_start + i * dt0;
        tn = (curr_t - t0) * tE_inv;
        y1 = u0 * salpha + tn * calpha; // Mine
        y2 = u0 * calpha - tn * salpha; // Mine
        curr_mag = tripleFS2python( y1, y2, rs);
        lightkv->append(curr_t, curr_mag);
        if (i > 0) {
            // error estimator, can be changed

            mag2 = curr_mag;
            mag1 = lightkv->last->prev->x2;
            tm2 = curr_t;
            tm1 = lightkv->last->prev->x1;

            // compute the midtm_mag,
            midtm = (tm1 + tm2) * 0.5;
            tn = (midtm - t0) * tE_inv;
            y1 = u0 * salpha + tn * calpha; // Mine
            y2 = u0 * calpha - tn * salpha; // Mine
            midtm_mag = tripleFS2python(y1, y2, rs);
            lightkv->last->prev->area = midtm_mag;

            // use real mu rather than log(mu - 1)
            midtm_mag = pow(10, midtm_mag) + 1;

            midtm_mag_inter = (mag1 + mag2) * 0.5;
            midtm_mag_inter = pow(10, midtm_mag_inter);
            def_err = fabs( (midtm_mag_inter - midtm_mag) / midtm_mag );
            // def_err = fabs( (midtm_mag_lin - midtm_mag_inter)/midtm_mag_lin );

            _toterr
            lightkv->last->prev->error = toterr;
            //using area to save midtm_mag


        }
    }

// add pnt where error is not satisfy the requirement
    double subdt = 0;
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


// void TripleLensing::polynomialCoefficients(double xs, double ys, complex c[])
// {
//     int i, j, k;
//     int n;

//     complex zs, zsc;      /* source position and conjugate */
//     complex zc[NLENS];        /* conjugate of lens positions */

//     complex p[NLENS + 1][NLENS], q[NLENS + 1][NLENS];
//     complex temp[NLENS + 1];
//     complex ctemp[DEGREE + 1];
//     complex qtemp[DEGREE + 1], qtemp2[NLENS + 1];
//     complex ptemp[DEGREE + 1];
//     int degrees;

//     n = NLENS;
//     zs = complex(xs, ys);
//     zsc = conj(zs);

//     for (i = 0; i < n; i++) {
//         zc[i] = conj(zlens[i]);
//     }

//     for (i = 0; i < n; i++) {
//         temp[0] = 1.0;
//         for (j = 0; j < n; j++) {
//             multiply_z(temp, zlens[j], j);
//         }

//         /* numerator, m[i] * (z-z1)*(z-z2) ... (z-zn) */
//         for (j = 0; j <= n; j++) {
//             p[j][i] = mlens[i] * temp[j];
//         }

//         /* denominator */
//         for (j = 0; j <= n; j++) { /* (zsc-conjugate(z[i])) * Product_{j=1}^N (z-z[j]) */
//             q[j][i] = (zsc - zc[i]) * temp[j];
//         }

//         /* Sum_{j=1}^n Product (z-z_k), k=1, n, but k !=j. */
//         for (j = 0; j < n; j++) {
//             /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */
//             temp[0] = 1.0;
//             degrees = 0;
//             for (k = 0; k < n; k++) {
//                 if (k == j) continue;

//                 multiply_z(temp, zlens[k], degrees);
//                 degrees++;
//             }

//             /* doing the sum */
//             for (k = 0; k < n; k++) {
//                 q[k][i] = q[k][i] + mlens[j] * temp[k];
//             }
//         }
//     }

//     /* now clear the fractions, z-zs - Sum_i p_i/q_i */
//     /* get the polynomial Product q_i, i=1, n */

//     /* first term */
//     qtemp[0] = 1.0;
//     degrees = 0;
//     for (i = 0; i < n; i++) {
//         for (j = 0; j <= n; j++) {
//             qtemp2[j] = q[j][i];
//         }

//         multiply(qtemp, degrees, qtemp2, n, ctemp);

//         degrees += n;


//         for (j = 0; j <= degrees; j++) {
//             qtemp[j] = ctemp[j];
//         }
//     }

//     /* get coefficients of (z-zs) Product_i=1^n q_i */
//     multiply_z(ctemp, zs, degrees);

//     /* copy the coefficients */
//     for (i = 0; i < n * n + 2; i++) {
//         c[i] = ctemp[i];
//     }

//     /* second term */
//     for (i = 0; i < n; i++) {
//         degrees = 0;
//         qtemp[0] = 1.0;
//         for (j = 0; j < n; j++) {
//             if (j == i) continue;

//             for (k = 0; k <= n; k++) {
//                 qtemp2[k] = q[k][j];
//             }

//             multiply(qtemp, degrees, qtemp2, n, ctemp);

//             degrees += n;

//             for (k = 0; k <= degrees; k++) {
//                 qtemp[k] = ctemp[k];
//             }
//         }

//         for (k = 0; k <= n; k++) {
//             ptemp[k] = p[k][i];
//         }

//         multiply(qtemp, degrees, ptemp, n, ctemp);

//         for (k = 0; k < n * n + 1; k++) {
//             c[k] = c[k] - ctemp[k];
//         }
//     }
// }


// obtain the polynomical coefficients
// input: lens mass and positions, mlens[], zlens[], source positions: xs, ys
// output: c

void TripleLensing::polynomialCoefficients(double xs, double ys, complex c[])
{
    int i, j, k;
    // int n;
    // complex zsc;      /* source position and conjugate */
    // complex zc[NLENS];        /* conjugate of lens positions */

    // complex p[NLENS + 1][NLENS], q[NLENS + 1][NLENS];
    // complex temp[NLENS + 1];
    // complex ctemp[DEGREE + 1];
    // complex qtemp[DEGREE + 1], qtemp2[NLENS + 1];
    // complex ptemp[DEGREE + 1];
    // int degrees;

    // n = NLENS;

    // complex zs;
    zs = complex(xs, ys);
    zsc = conj(zs);

    // fprintf(stderr, "in polynomialCoefficients, mlens[0 1 2] %f %f %f\n", mlens[0], mlens[1], mlens[2] );



    for (i = 0; i < NLENS; i++) {

        // temp[0] = 1.0;
        // for (j = 0; j < NLENS; j++) {
        //     multiply_z(temp, zlens[j], j);
        // }

        /* numerator, m[i] * (z-z1)*(z-z2) ... (z-zn) */
        // for (j = 0; j <= NLENS; j++) {
        //     p_const[j][i] = mlens[i] * temp[j];
        // }

        /* denominator */
        for (j = 0; j <= NLENS; j++) { /* (zsc-conjugate(z[i])) * Product_{j=1}^N (z-z[j]) */
            q[j][i] = (zsc - zc[i]) * temp_const1[i][j];
        }

        /* Sum_{j=1}^n Product (z-z_k), k=1, n, but k !=j. */
        for (j = 0; j < NLENS; j++) {
            /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */

            // temp[0] = 1.0;
            // degrees = 0;
            // for (k = 0; k < NLENS; k++) {
            //     if (k == j) continue;

            //     multiply_z(temp, zlens[k], degrees);
            //     degrees++;
            // }

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

// void TripleLensing::polynomialCoefficients(double xs, double ys, complex c[])
// {
//     int i, j, k;
//     // int n;
//     // complex zsc;      /* source position and conjugate */
//     // complex zc[NLENS];        /* conjugate of lens positions */

//     // complex p[NLENS + 1][NLENS], q[NLENS + 1][NLENS];
//     // complex temp[NLENS + 1];
//     // complex ctemp[DEGREE + 1];
//     // complex qtemp[DEGREE + 1], qtemp2[NLENS + 1];
//     // complex ptemp[DEGREE + 1];
//     // int degrees;

//     // n = NLENS;

//     // complex zs;
//     zs = complex(xs, ys);
//     zsc = conj(zs);



//     for (i = 0; i < NLENS; i++) {

//         temp[0] = 1.0;
//         for (j = 0; j < NLENS; j++) {
//             multiply_z(temp, zlens[j], j);
//         }

//         /* numerator, m[i] * (z-z1)*(z-z2) ... (z-zn) */
//         for (j = 0; j <= NLENS; j++) {
//             p[j][i] = mlens[i] * temp[j];
//         }

//         /* denominator */
//         for (j = 0; j <= NLENS; j++) { /* (zsc-conjugate(z[i])) * Product_{j=1}^N (z-z[j]) */
//             q[j][i] = (zsc - zc[i]) * temp[j];
//         }

//         /* Sum_{j=1}^n Product (z-z_k), k=1, n, but k !=j. */
//         for (j = 0; j < NLENS; j++) {
//             /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */
//             temp[0] = 1.0;
//             degrees = 0;
//             for (k = 0; k < NLENS; k++) {
//                 if (k == j) continue;

//                 multiply_z(temp, zlens[k], degrees);
//                 degrees++;
//             }

//              doing the sum
//             for (k = 0; k < NLENS; k++) {
//                 q[k][i] = q[k][i] + mlens[j] * temp[k];
//             }
//         }
//     }

//     /* now clear the fractions, z-zs - Sum_i p_i/q_i */
//     /* get the polynomial Product q_i, i=1, n */

//     /* first term */
//     qtemp[0] = 1.0;
//     degrees = 0;
//     for (i = 0; i < NLENS; i++) {
//         for (j = 0; j <= NLENS; j++) {
//             qtemp2[j] = q[j][i];
//         }

//         multiply(qtemp, degrees, qtemp2, NLENS, ctemp);

//         degrees += NLENS;


//         for (j = 0; j <= degrees; j++) {
//             qtemp[j] = ctemp[j];
//         }
//     }

//     /* get coefficients of (z-zs) Product_i=1^n q_i */
//     multiply_z(ctemp, zs, degrees);

//     /* copy the coefficients */
//     for (i = 0; i < DEGREE + 1; i++) {
//         c[i] = ctemp[i];
//     }

//     /* second term */
//     for (i = 0; i < NLENS; i++) {
//         degrees = 0;
//         qtemp[0] = 1.0;
//         for (j = 0; j < NLENS; j++) {
//             if (j == i) continue;

//             for (k = 0; k <= NLENS; k++) {
//                 qtemp2[k] = q[k][j];
//             }

//             multiply(qtemp, degrees, qtemp2, NLENS, ctemp);

//             degrees += NLENS;

//             for (k = 0; k <= degrees; k++) {
//                 qtemp[k] = ctemp[k];
//             }
//         }

//         for (k = 0; k <= NLENS; k++) {
//             ptemp[k] = p[k][i];
//         }

//         multiply(qtemp, degrees, ptemp, NLENS, ctemp);

//         for (k = 0; k < DEGREE; k++) {
//             c[k] = c[k] - ctemp[k];
//         }
//     }
// }



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


// From Tianshu, should be further simplified, otherwise it is slow
// void TripleLensing::TripleCoefs(double m1, double m2, double m3, complex z1, complex z2, complex z3, complex zs, complex *coefs)
// void TripleLensing::TripleCoefs(double xs, double ys, complex * coefs)
// {}