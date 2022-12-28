// #include "VBBinaryLensingLibrary.h"
#include "TripleLensingLibrary.h"
#include <stdio.h>
#include <math.h>
#include <time.h>

// VBBinaryLensing VBBL;

TripleLensing::TripleLensing() {
    nphi = 32;
    secnum = 45;
    basenum = 1;
    quad_err = 0;
    quaderr_Tol = 1e-4;
    CQ = 1;
    relerr_mag = 1e-3;
    TINY = 1.0e-20;
    NLENS = 3;
    DEGREE = NLENS * NLENS + 1;
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


void TripleLensing::setnlens(int nlens) {
    // later you need a destructor for the arrays created here
    NLENS = nlens;
    DEGREE = NLENS * NLENS + 1;

    // https://blog.csdn.net/r_e_d_h_a_t/article/details/78178196

    // complex p[NLENS + 1][NLENS], q[NLENS + 1][NLENS], p_const[NLENS + 1][NLENS], temp[NLENS + 1], temp_const1[NLENS + 1][NLENS + 1], temp_const2[NLENS + 1][NLENS + 1][NLENS + 1], temp_const22[NLENS + 1], ctemp[DEGREE + 1], qtemp[DEGREE + 1], qtemp2[NLENS + 1], ptemp[DEGREE + 1], zr[DEGREE], coefficients[DEGREE + 1];

    temp  = new   complex [NLENS + 1];
    temp_const22  = new   complex [NLENS + 1];
    ctemp  = new   complex [DEGREE + 1];
    qtemp  = new   complex [DEGREE + 1];
    qtemp2  = new   complex [NLENS + 1];
    ptemp  = new   complex [DEGREE + 1];
    zr  = new   complex [DEGREE];
    zc  = new   complex [NLENS];
    coefficients  = new   complex [DEGREE + 1];

    p   =   new   complex*[NLENS + 1];
    q   =   new   complex*[NLENS + 1];
    p_const   =   new   complex*[NLENS + 1];
    temp_const1   =   new   complex*[NLENS + 1];

    for (int  i = 0;   i < NLENS + 1;   ++i) {
        p[i]   =   new   complex[NLENS];
        q[i]   =   new   complex[NLENS];
        p_const[i]   =   new   complex[NLENS];
        temp_const1[i]   =   new   complex[NLENS + 1];
    }

    temp_const2 = new   complex**[NLENS + 1];
    for (int  i = 0;   i < NLENS + 1;   ++i) {
        temp_const2[i] = new   complex*[NLENS + 1];
        for (int  j = 0;   j < NLENS + 1;   ++j) {
            temp_const2[i][j] = new   complex[NLENS + 1];
        }
    }
    // complex p[NLENS + 1][NLENS], q[NLENS + 1][NLENS], p_const[NLENS + 1][NLENS], temp[NLENS + 1], temp_const1[NLENS + 1][NLENS + 1], temp_const2[NLENS + 1][NLENS + 1][NLENS + 1], temp_const22[NLENS + 1], ctemp[DEGREE + 1], qtemp[DEGREE + 1], qtemp2[NLENS + 1], ptemp[DEGREE + 1], zr[DEGREE], coefficients[DEGREE + 1];
}

TripleLensing::~TripleLensing() {
    // destructor
    delete [] temp;
    delete [] temp_const22;
    delete [] ctemp;
    delete [] qtemp;
    delete [] qtemp2;
    delete [] ptemp;
    delete [] zr;
    delete [] zc;
    delete [] coefficients;

    for (int  i = 0;   i < NLENS + 1;   ++i) {
        delete [] p[i];
        delete [] q[i] ;
        delete [] p_const[i];
        delete [] temp_const1[i];
    }
    delete [] p;
    delete [] q;
    delete [] p_const;
    delete [] temp_const1;

    for (int  i = 0;   i < NLENS + 1;   ++i) {
        for (int  j = 0;   j < NLENS + 1;   ++j) {
            delete [] temp_const2[i][j];
        }
        delete [] temp_const2[i];
    }
    delete [] temp_const2;

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



void TripleLensing::reset3(double mlens[], complex zlens[]) {
    // just update mlens and zlens
    this -> mlens = mlens;
    this -> zlens = zlens;


    // for(int i = 0; i<NLENS; i++){
    //     printf("mlens, x, y = %f, %f, %f\n", mlens[i], zlens[i].re, zlens[i].im);
    // }

    //used in polynomialCoefficients
    for (int i = 0; i < NLENS; i++) {
        zc[i] = conj(zlens[i]);
        // printf("313 %d, zc[i] = %f %f\n",i, zc[i].re, zc[i].im);
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
    // for (int i = 0; i < NLENS+1; i++) {
    //     for (int j = 0; j < NLENS+1; j++) {
    //         for (int k = 0; k < NLENS + 1; k++) {
    //             printf("i, j, k, temp_const2 = %d, %d, %d, %f %f\n", i, j, k, temp_const2[i][j][k].re, temp_const2[i][j][k].im);
    //         }
    //     }
    // }


    // for (int i = 0; i < NLENS+1; i++) {
    //     for (int j = 0; j < NLENS+1; j++) {
    //             printf("i, j, temp_const1 = %d, %d, %f %f\n", i, j, temp_const1[i][j].re, temp_const1[i][j].im);
    //     }
    // }


}


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

    if (rs > 1e-10) {
        for (int i = 0; i < Np; i++) {
            mags[i] = arrTripleMag(xsCenters[i], ysCenters[i], rs);
        }
    } else {
        for (int i = 0; i < Np; i++) {
            mags[i] = TriplePS(xsCenters[i], ysCenters[i]);
        }
    }

}

// Gould approximation calculate interface to Python
void TripleLensing::tripleGould2py(double mlens[], double Zlens[], double xsCenters[], double ysCenters[], double rs, double Gamma, double mags[], int Np)
{


    double A2rho2, A4rho4;
    int Np2 = 2 * Np;
    complex zlens[NLENS];
    for (int i = 0; i < NLENS ; i++) {
        zlens[i] = complex( Zlens[i * 2], Zlens[i * 2 + 1] );
    }
    reset3(mlens, zlens);
    for (int i = 0; i < Np; i++) {
        mags[i] = gould(xsCenters[i], ysCenters[i], rs, Gamma, &A2rho2, &A4rho4);
        mags[Np + i] = A2rho2;
        mags[Np2 + i] = A4rho4;
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
        cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

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
    mu = arrTripleMag(xsCenter, ysCenter, rs);
    return mu;
}

double TripleLensing::tripleFS2python(double mlens[], double zlens[], double xsCenter, double ysCenter, double rs) {
    reset2(mlens, zlens);
    mu = arrTripleMag(xsCenter, ysCenter, rs);
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
    cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

    for (i = 0; i < DEGREE; i++) {
        zrxy[i] = zr[i].re;
        zrxy[i + DEGREE] = zr[i].im;
    }
}


void TripleLensing::outsys(double mlens[], complex zlens[], double t0, double u0, double tE, double s2, double q2, double alpha, double s3, double q3, double psi, double rs, double xsCenter, double ysCenter) {
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


void TripleLensing::arrareaFunc(double *area){
        // nfinal_closed_image = closed_image_info.shape[0]
        // area_array = np.zeros(nclosed_image)
        *area = 0;
        double tmparea, x1, y1, x2, y2, tmpdb;
        unsigned short int jj , kk , offset, nsegs;
        int parity, NPT, j, npts, isgn, scanpnt, hj, tj;

        for (jj=0; jj<nfinal_closed_image; jj++){
            parity = closed_image_info[jj][0]; nsegs = closed_image_info[jj][1];
            #ifdef VERBOSE
                printf("%d, nsegs = %d, \n", jj, nsegs);
                for (kk=0; kk< nsegs; kk++){
                    offset = kk*3;
                    j = closed_image_info[jj][2+offset]; hid = closed_image_info[jj][3+offset]; tid = closed_image_info[jj][4+offset];
                    printf("\t, %d, %d, %d\n", j, hid, tid);
                }
            #endif

            NPT = 0;
            tmparea = 0;
            for (kk=0; kk< nsegs; kk++){
                offset = kk*3;
                // j, hid, tid = closed_image_info[jj, 2+offset: 5+offset]
                j = closed_image_info[jj][2+offset]; hid = closed_image_info[jj][3+offset]; tid = closed_image_info[jj][4+offset];
                npts = fabs(tid - hid) + 1;
                if (npts < 2){
                    // print(jj, kk, "this segment has only %d point"%npts)
                    continue; // # just in case 
                }
                NPT += npts;
                if (hid < tid){
                    isgn = 1;
                }else{
                    isgn = -1;                    
                }
                scanpnt = hid;
                x1 = allSolutions_x[j][scanpnt]; y1 = allSolutions_y[j][scanpnt];
                scanpnt += isgn;
                while (scanpnt != tid+isgn){
                    x2 = allSolutions_x[j][scanpnt]; y2 = allSolutions_y[j][scanpnt];
                    tmparea += 0.5 * (y1 + y2) * (x2 - x1);
                    x1 = x2; y1 = y2; // x1, y1 = x2, y2
                    scanpnt += isgn;
                }
            }
            if (nsegs > 1){ //# if nsegs > 1, you need to add the area between two segments
                for (kk = 0; kk < nsegs - 1; kk++){
                    offset = kk*3;
                    // j, hid, tid = closed_image_info[jj, 2+offset: 5+offset]
                    j = closed_image_info[jj][2+offset]; hid = closed_image_info[jj][3+offset]; tid = closed_image_info[jj][4+offset];
                    x1 = allSolutions_x[j][tid]; y1 = allSolutions_y[j][tid];

                    offset = kk*3 + 3;
                    // j, hid, tid = closed_image_info[jj, 2+offset: 5+offset]
                    j = closed_image_info[jj][2+offset]; hid = closed_image_info[jj][3+offset]; tid = closed_image_info[jj][4+offset];
                    x2 = allSolutions_x[j][hid]; y2 = allSolutions_y[j][hid];
                    tmparea += 0.5 * (y1 + y2) * (x2 - x1);
                }
            }
            //# the tail and head
            hj = closed_image_info[jj][2]; hid = closed_image_info[jj][3];
            tj = closed_image_info[jj][2+3*nsegs-3]; tid = closed_image_info[jj][2+3*nsegs-1];            
            x1 = allSolutions_x[tj][tid]; y1 = allSolutions_y[tj][tid];
            x2 = allSolutions_x[hj][hid]; y2 = allSolutions_y[hj][hid];
            tmparea += 0.5 * (y1 + y2) * (x2 - x1);
#ifdef VERBOSE
            printf("%d, track total length = %d, parity = %d, nsegs = %d, area = %f\n", jj, NPT, parity, nsegs, tmparea);
#endif
            tmparea *= parity;
            if (NPT > finalnphi*0.25){
                *area += tmparea;
            }
        }
        *area = fabs(*area);
        // #area_array[jj] = area
        // return abs(total_area/ (np.pi * rs * rs) )
}


double TripleLensing::arrTripleMag(double xsCenter, double ysCenter, double rs) {
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

    for (int jj = 0; jj < DEGREE; jj++) {
        zr[jj] = complex(0, 0);
    }

    muPS = tripleQuatrapoleTest(xsCenter, ysCenter, rs); // this is not robust, didn't check the number of true solution is 4, 6, 8, or 10
    // scanf("%d", &ftime);

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

        bool prevstore;
        double mindphi = M_PI, tmpdouble;
        get_arrphi(xsCenter, ysCenter, rs, &finalnphi); // ok
        for (int i = 1; i < finalnphi; i++) {
            tmpdouble = ARRPHI[i] - ARRPHI[i - 1];
            if (mindphi > tmpdouble) mindphi = tmpdouble;
        }

        if (0) {
            printf("nphi = %d, mindphi = %f\n", finalnphi, mindphi * 180 / M_PI); // nphi = 704, mindphi = 0.166667
            for (int i = 0; i < finalnphi; i++) {
                printf("i = %d, phi = %f\n", i, ARRPHI[i] * 180 / M_PI);
            }
            extend_arrphi(&finalnphi);
            printf("after extend nphi = %d, mindphi = %f\n", finalnphi, mindphi * 180 / M_PI); // nphi = 704, mindphi = 0.166667
            for (int i = 0; i < finalnphi; i++) {
                printf("i = %d, phi = %f\n", i, ARRPHI[i] * 180 / M_PI);
            }
            // nphi = 16, mindphi = 22.500000
            // i = 0, phi = 270.000000
            // i = 1, phi = 292.500000
            // i = 2, phi = 315.000000
            // after extend nphi = 31,
            // i = 0, phi = 270.000000
            // i = 1, phi = 281.250000
            // i = 2, phi = 292.500000
            // i = 3, phi = 303.750000
        }

        // fprintf(stderr, "initial nphi = %d\n", finalnphi);
        // for (int i = 0; i < finalnphi + 5; i++){
        //     if (i<5 || ((i<finalnphi+2) && (i>finalnphi - 5))){
        //         printf("i = %d, phi = %f\n", i, ARRPHI[i]*180/M_PI);
        //     }
        // }
        // initial nphi = 1322
        // i = 0, phi = 710.000000
        // i = 1, phi = 710.400000
        // i = 2, phi = 710.800000
        // i = 1319, phi = 1069.111111
        // i = 1320, phi = 1069.555556
        // i = 1321, phi = 1070.000000
        // i = 1322, phi = 0.000000

//         for (int i = 0; i < 10; i++) {
        for (int i = 0; i < 7; i++) { // 2022.06.05
            if (i > 0) { // 2021.06.08
#ifdef VERBOSE
                fprintf(stderr, "i = %d, finalnphi = %d\n", i, finalnphi);
#endif
            }
            prevstore = true;
            if (i == 0) {prevstore = false;};
            // imageTracks = outputTracks_v2_savehalf(xsCenter,  ysCenter,  rs, phis, &prevstore);
            arroutputTracks(xsCenter, ysCenter, rs, prevstore, finalnphi, mindphi);
            arrareaFunc(&area);

// #ifdef VERBOSE
//         char arr[1024];
//         fprintf(stderr, "please input, area = %f >>> ", area);
//         scanf("%c%*c", &arr[0]);
// #endif

            //finalnphi = phis->length;
#ifdef VERBOSE
            //saveTracks(imageTracks, finalnphi);// temp
            fprintf(stderr, "\t\t\t imageTracks saved (temp). \n");
#endif


            // if (!imageTracks) {
            if (0) {
#ifdef VERBOSE
                fprintf(stderr, "\n\nin tripleFS, i= %d, imageTracks = Null (be careful), nphi = %d, xsCenter=%f, ysCenter = %f\n\n", i, finalnphi, xsCenter, ysCenter);
#endif
                if (distype == 2) {
                    maxmuidx = (int)(absint(maxmuidx + 1)) % secnum;
                    basenum_priv *= 2;
                    //phis = getphis_v3( xsCenter,  ysCenter, ys);
                    get_arrphi(xsCenter, ysCenter, rs, &finalnphi);
                } else {
                    //phis->extend();
                    extend_arrphi(&finalnphi); mindphi *= 0.5;
                }
                i -= 1;
                ftime = 1;
                continue;
            }
            ftime = 0;

            mu = area / areaSource;

            area_quality = area_quality_local;
            if (area_quality_local == 0) {
                _curr_relerr_priv *= 2;
            }

#ifdef VERBOSE
            fprintf(stderr, "in tripleFS, i= %d, mu0= %f, mu= %f, nphi = %d, xsCenter=%f, ysCenter = %f, nimages = %d, abs(mu - mu0) / mu = %.3e, errTol = %.3e, \n", i, mu0, mu, finalnphi, xsCenter, ysCenter, nfinal_closed_image, abs(mu - mu0) / mu, _curr_relerr_priv);
#endif
#ifdef verbose
            fprintf(stderr, "in tripleFS, i= %d, mu0= %f, mu= %f, nphi = %d, xsCenter=%f, ysCenter = %f, nimages = %d\n", i, mu0, mu, finalnphi, xsCenter, ysCenter, imageTracks->length);
#endif
            if (abs(mu - mu0) / mu < _curr_relerr_priv) {
                break;
            } else if (abs(mu - mu0) / mu < _curr_relerr_priv * i  ||  finalnphi > 5e5 ) { // add "|  imageTracks->length > 1e6" on 2022.06.04, when logq is too small <~ -6, the memory usage increases due to he amount of sampling points on the source boundary
#ifdef VERBOSE
                fprintf(stderr, "return mu = %f, because mu0 = %f, mu = %f, relerr_priv*i = %.3e this might be a hard case, area_quality_local = %d\n", mu, mu0, mu, _curr_relerr_priv * i, area_quality_local);

#endif

                area_quality = 2; // return because of looser threshold, quality might be low
                return (mu);
//             } else if ( (i > 2 && (0.5 * ( abs(mu / mu0) + abs(mu0 / mu) > 2.1)) ) ) {
            } else if ( (i > 2 && ( 0.5 * (abs(mu / mu0) + abs(mu0 / mu)) > 2.1 ) ) ) { // modified on 2022.06.04
#ifdef VERBOSE
                fprintf(stderr, "return muPS = %f, because mu0 = %f, mu = %f are too arbitrary, this might be a hard case, area_quality_local=%d\n", muPS, mu0, mu, area_quality_local);
#endif

                area_quality = 3; // return point source magnification
                return (muPS);
            }

            mu0 = mu;
            //phis->extend();
            extend_arrphi(&finalnphi);  mindphi *= 0.5;
        }
#ifdef VERBOSE
        printf("mu = %f\n", mu);
        // saveTracks(imageTracks);
#endif
        finalNPS = finalnphi;
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
    printf("I am in gould, second and fourth order corrections: %f %f\n", *A2rho2 / 2.0, *A4rho4 / 3.0);
#endif

    mag = A0 + *A2rho2 / 2.0 * (1.0 - 0.2 * Gamma)  + *A4rho4 / 3.0 * (1.0 - 11.0 / 35.0 * Gamma);
    return (mag);
}


// tripleQuatrapoleTest Bozza 2018 equation 39 with c_Q = 1, add on
double TripleLensing::tripleQuatrapoleTest(double xs, double ys, double rs) {

    // for(int i = 0; i<NLENS; i++) printf("1199 %d, zc[i] = %f %f\n",i, zc[i].re, zc[i].im);


    // printf("in tripleQuatrapoleTest, xs, ys = %f, %f\n", xs, ys);
    polynomialCoefficients(xs, ys, coefficients);
    // for (int i = 0; i<DEGREE; i++) printf("coe %d = %f %f\n", i, coefficients[i].re, coefficients[i].im);
    cmplx_roots_gen(zr, coefficients, DEGREE, true, true);
    // for (int i = 0; i<DEGREE; i++) printf("zr %d = %f %f\n", i, zr[i].re, zr[i].im);



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


double TripleLensing::TriplePS(double xs, double ys) {
    double mu;
    double mulist[DEGREE];
    int total_parity = 0, flaglist[DEGREE], absdzslist[DEGREE];

    polynomialCoefficients(xs, ys, coefficients);
    cmplx_roots_gen(zr, coefficients, DEGREE, true, true);
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
    cmplx_roots_gen(zr, coefficients, DEGREE, true, true);
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
        cmplx_roots_gen(zr, coefficients, DEGREE, true, true);
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


void TripleLensing::get_arrphi(double xsCenter, double ysCenter, double rs, int *retnphi) {
    // calculate phis and save to ARRPHI
    secnum_priv = this->secnum;
    basenum_priv = this->basenum;

    distype = 2;

    double phi[secnum_priv]; // secnum_priv is specified by user, can also be fixed at compile
    double mus[secnum_priv];
    double offset[secnum_priv];
    int npmus[secnum_priv];
    int secnum_privlist[secnum_priv];
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
    fprintf(stderr, "at get_arrphi secnum_priv %d,basenum %d,  maxmuidx: %d, minmu = %f, maxmu = %f\n", secnum_priv, basenum_priv, maxmuidx, minmu, maxmu);
#endif
    for (int i = 0; i < secnum_priv; i++) {
        mus[i] = (mus[i] / minmu);
    }

    double psf[3] = {1, 2, 1};
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

    //PHI->linspace(offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, npmus[tempidx]*basenum_priv, 0);
    *retnphi = 0;
    arrlinspace(ARRPHI, offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, *retnphi, npmus[tempidx]*basenum_priv, 0);
    *retnphi += npmus[tempidx] * basenum_priv;
    for (int i = 1; i < secnum_priv - 1; i++) {
        tempidx = secnum_privlist[i];
        // PHI->linspace(offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, npmus[tempidx]*basenum_priv, 0);
        arrlinspace(ARRPHI, offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, *retnphi, npmus[tempidx]*basenum_priv, 0);
        *retnphi += npmus[tempidx] * basenum_priv;
    }
    tempidx = secnum_privlist[secnum_priv - 1];
    // PHI->linspace(offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, npmus[tempidx]*basenum_priv, 1);
    arrlinspace(ARRPHI, offset[tempidx] + phi[tempidx] - dphi + M_PI2, offset[tempidx] + phi[tempidx] + dphi + M_PI2, *retnphi, npmus[tempidx]*basenum_priv, 1);
    *retnphi += npmus[tempidx] * basenum_priv;

    // return PHI;
}


void TripleLensing::arrlinspace(double *arr, double phi0, double phiend, int insertidx, int nphi, int endpoint) {
    // insertidx: the index to be insert point
    double dphi;
    double scanphi = phi0;
    if (endpoint) {
        dphi = (phiend - phi0) / (nphi - 1);
    } else {
        dphi = (phiend - phi0) / (nphi);
    }
    for (int i = insertidx; i < nphi + insertidx; i++) {
        arr[i] = scanphi;
        scanphi += dphi;
    }
}

void TripleLensing::extend_arrphi(int *retnphi) {

    for (int i = *retnphi - 1; i > -1; i--) {
        ARRPHI[i * 2] = ARRPHI[i];
    }

    *retnphi = (*retnphi + *retnphi - 1);

    for (int i = 1; i < *retnphi - 1; i = i+2) {
        ARRPHI[i] = (ARRPHI[i - 1] + ARRPHI[i + 1]) * 0.5;
    }

}

void TripleLensing::arroutputTracks(double xsCenter, double ysCenter, double rs, bool prevstore, int nphi, double mindphi) {
        double mindsource;
        mindsource = rs*2*sin((0.5*mindphi));// # minimum distance between two sampled points at the source boundary
        mindsource = mindsource*mindsource;
        int curr_total_parity = 0, i, j, offset, saveindex;
        unsigned short int curr_nimages = 0; //# how many true solutions for a certain source position
        bool ifsave;

        for(i = 0; i<DEGREE; i++){posflagnums[i] = 0;}
        for(i=0; i<DEGREE; i++){
            for (j=0; j< 2+3*(DEGREE*DEGREE+1); j++){
                closed_image_info[i][j] = 0;
            }
        }
        for(i=0; i<DEGREE*DEGREE; i++){
            for (j=0; j< 3; j++){
                true_segments_info[i][j] = 0;
            }
        }
            

        if (!prevstore){ //# equavalent to cpp code: if (ftime)
            // # solve lens equations for the first time 
            for (j = 0; j<nphi; j++){
                phi = ARRPHI[j];
                xs = xsCenter + rs * cos(phi);
                ys = ysCenter + rs * sin(phi);

                //res = sol_len_equ_cpp(self.mlens, self.zlens, xs, ys, self.nlens, DEGREE)
                polynomialCoefficients(xs, ys, coefficients);
                cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

                if (j == 0){
                    curr_total_parity = 0; curr_nimages = 0;
                    for(i=0; i<DEGREE;i++){
                        //flag = trueSolution(self.mlens, self.zlens, xs, ys, res[i], cal_ang = False, NLENS = self.nlens);// # [flag, mu, absdzs, xxx]
                        flag = trueSolution(xs, ys, zr[i], &mu);
                        allSolutions_x[i][j] = zr[i].re;
                        allSolutions_y[i][j] = zr[i].im;
                        allSolutions_flag[i][j] = flag==1?true:false;
                        allSolutions_mu[i][j] = mu;
                        allSolutions_srcx[i][j] = xs;
                        allSolutions_srcy[i][j] = ys;
                        allSolutions_absdzs[i][j] = absdzs;
                        posflagnums[i] += flag;

                        if (flag){curr_total_parity += one_parity(mu);}
                        curr_nimages += flag;
                    }
                    if (((curr_nimages - NLENS)%2 != 1) || (curr_total_parity != CORRECT_PARITY) ) {
                        // # there might be some true images being classified as false images
                        saveimage(curr_nimages, curr_total_parity, j, &ifsave, &saveindex);
                        if (ifsave){
                            allSolutions_flag[saveindex][j] = 1;
                            posflagnums[saveindex] += 1;
                        }
                    }

                }else{
                    // # for j > 0, you need to attached the closest solution to the existing allSolutions
                    curr_total_parity = 0; curr_nimages = 0;
                    for(i=0; i<DEGREE;i++){
                        // flag = trueSolution(self.mlens, self.zlens, xs, ys, res[i], cal_ang = False, NLENS = self.nlens) #
                        flag = trueSolution(xs, ys, zr[i], &mu);
                        Prov_x[i] = zr[i].re;
                        Prov_y[i] = zr[i].im;
                        Prov_flag[i] = flag==1?true:false;
                        Prov_mu[i] = mu;
                        Prov_srcx[i] = xs;
                        Prov_srcy[i] = ys;
                        Prov_absdzs[i] = absdzs;
                        preProv_x[i] = allSolutions_x[i][j-1];// # the last point in the i-th image track
                        preProv_y[i] = allSolutions_y[i][j-1];
                        if (flag){ curr_total_parity += one_parity(mu); }
                        curr_nimages += flag;
                    }

                    //attach_idx = self.get_closest_idx(preProv_x, preProv_y, Prov_x, Prov_y)
                    get_closest_idx(preProv_x, preProv_y, Prov_x, Prov_y, attach_idx);

                    for(i=0; i<DEGREE;i++){
                        allSolutions_x[i][j] = Prov_x[attach_idx[i]];
                        allSolutions_y[i][j] = Prov_y[attach_idx[i]];
                        allSolutions_flag[i][j] = Prov_flag[attach_idx[i]];
                        allSolutions_mu[i][j] = Prov_mu[attach_idx[i]];
                        allSolutions_srcx[i][j] = Prov_srcx[attach_idx[i]];
                        allSolutions_srcy[i][j] = Prov_srcy[attach_idx[i]];
                        allSolutions_absdzs[i][j] = Prov_absdzs[attach_idx[i]];
                        posflagnums[i] += Prov_flag[attach_idx[i]]?1:0;
                    }
                    if ( ((curr_nimages - NLENS)%2 != 1) || (curr_total_parity != CORRECT_PARITY) ){
                        // # there might be some true images being classified as false images
                        saveimage(curr_nimages, curr_total_parity, j, &ifsave, &saveindex);
                        if (ifsave){
                            allSolutions_flag[saveindex][j] = 1;
                            posflagnums[saveindex] += 1;
                        }
                    }
                }
            }
        }else{ //# re-use previously obtained solutions
            //## allSolutions_#[i, k*2] <-- prevstore_#[i, k], where k = 0, 1, 2, ...
            //for j in range( (nphi-1)//2, -1, -1 ):
            for (j = (int)((nphi-1)*0.5); j>-1; j--){
                for (i=0; i<DEGREE; i++){
                    allSolutions_x[i][j*2] = allSolutions_x[i][j];
                    allSolutions_y[i][j*2] = allSolutions_y[i][j];
                    allSolutions_flag[i][j*2] = allSolutions_flag[i][j];
                    allSolutions_mu[i][j*2] = allSolutions_mu[i][j];
                    allSolutions_srcx[i][j*2] = allSolutions_srcx[i][j];
                    allSolutions_srcy[i][j*2] = allSolutions_srcy[i][j];
                    allSolutions_absdzs[i][j*2] = allSolutions_absdzs[i][j];
                    posflagnums[i] += allSolutions_flag[i][j];
                }
            }

            // for j in range(1, nphi, 2):
            for (j=1; j<nphi; j=j+2){
                phi = ARRPHI[j];
                xs = xsCenter + rs * cos(phi);
                ys = ysCenter + rs * sin(phi);             
                // res = sol_len_equ_cpp(self.mlens, self.zlens, xs, ys, self.nlens, DEGREE)
                polynomialCoefficients(xs, ys, coefficients);
                cmplx_roots_gen(zr, coefficients, DEGREE, true, true);

                //# for j > 0, you need to attached the closest solution to the existing allSolutions
                curr_total_parity = 0; curr_nimages = 0;
                for(i=0; i<DEGREE;i++){
                    // flag = trueSolution(self.mlens, self.zlens, xs, ys, res[i], cal_ang = False, NLENS = self.nlens)
                    flag = trueSolution(xs, ys, zr[i], &mu);
                    Prov_x[i] = zr[i].re;
                    Prov_y[i] = zr[i].im;
                    Prov_flag[i] = flag==1?true:false;
                    Prov_mu[i] = mu;
                    Prov_srcx[i] = xs;
                    Prov_srcy[i] = ys;
                    Prov_absdzs[i] = absdzs;
                    preProv_x[i] = allSolutions_x[i][j-1];// # the last point in the i-th image track
                    preProv_y[i] = allSolutions_y[i][j-1];
                    if (flag){ curr_total_parity += one_parity(mu); }
                    curr_nimages += flag;
                }
                // attach_idx = self.get_closest_idx(preProv_x, preProv_y, Prov_x, Prov_y)
                get_closest_idx(preProv_x, preProv_y, Prov_x, Prov_y, attach_idx);

                for(i=0; i<DEGREE;i++){
                    allSolutions_x[i][j] = Prov_x[attach_idx[i]];
                    allSolutions_y[i][j] = Prov_y[attach_idx[i]];
                    allSolutions_flag[i][j] = Prov_flag[attach_idx[i]];
                    allSolutions_mu[i][j] = Prov_mu[attach_idx[i]];
                    allSolutions_srcx[i][j] = Prov_srcx[attach_idx[i]];
                    allSolutions_srcy[i][j] = Prov_srcy[attach_idx[i]];
                    allSolutions_absdzs[i][j] = Prov_absdzs[attach_idx[i]];
                    posflagnums[i] += Prov_flag[attach_idx[i]]?1:0;
                }
                if (((curr_nimages - NLENS)%2 != 1) || (curr_total_parity != CORRECT_PARITY)){
                    // # there might be some true images being classified as false images
                    saveimage(curr_nimages, curr_total_parity, j, &ifsave, &saveindex);
                    if (ifsave){
                        allSolutions_flag[saveindex][j] = 1;
                        posflagnums[saveindex] += 1;
                    }
                }
            }
        }
        // # return allSolutions_x, allSolutions_y, allSolutions_srcx, allSolutions_srcy, allSolutions_mu, allSolutions_flag, allSolutions_absdzs, None, None

        // # now, select and connect true segments
        // # initialize array to store the information of true image segments:
        // # i.e., the i-th track, the index of head, the index of tail, (the length of the segments, this is not necessary)
        ntrue_segments = 0;
        // # true_segments_info = np.zeros((10 * DEGREE, 3 )).astype(int)
        // true_segments_info = np.zeros((DEGREE * DEGREE, 3 )).astype(int) // ## 这个也可以事先定义好一定长度的 array
        // # each row contains where this segments is, a, b, c,
        npure_close_segments = 0;

        //imgTrack_type = np.zeros(DEGREE) // # 0: pure false, 1: pure true, 2: mixed
        for (i=0; i<DEGREE; i++){
            imgTrack_type[i] = 0; 
            #ifdef VERBOSE
            printf("i = %d, posflagnums[i] = %d\n", i, posflagnums[i]);
            #endif
        }

        for (i=0; i<DEGREE; i++){
            // # iter over each image Track, which may contains only true images, false images, or mixed
            if (posflagnums[i] <= 1 || ( (1.0 * posflagnums[i] / nphi < 0.5) && fabs(allSolutions_mu[i][0]) < EPS  && fabs(allSolutions_mu[i][nphi-1]) < EPS ) ){
                    // if VERBOSE: print(i, "pure false");
                            }
                // pass # pure false image Track
            else if (posflagnums[i] == nphi) { //: #(1.0 * posflagnums[i] / nphi) > 0.99){
                imgTrack_type[i] = 1;
                //if VERBOSE: print(i, "pure true")
                true_segments_info[ntrue_segments][0] = i;
                true_segments_info[ntrue_segments][1] = 0;
                true_segments_info[ntrue_segments][2] = nphi - 1;// # index 0, 1, ..., nphi - 1 are true image points
                ntrue_segments += 1;
                npure_close_segments += 1;
            }else{
                imgTrack_type[i] = 2;
            }
        }

        // # now for each of the type == 2 image Tracks, select out all true image segments
        for (i=0; i<DEGREE; i++){
            if (imgTrack_type[i] == 2){
                // # select true image tracks
                //if VERBOSE: print(i, "now you only need to proceed this posflagnum %d, posflagnum/length %f, first->flag %d, mu %f, x1 x2 = %f, %f"%(posflagnums[i], 1.0 * posflagnums[i] / nphi, allSolutions_flag[i][0], allSolutions_mu[i][0], allSolutions_x[i][0], allSolutions_y[i][0]))
                // # iter over this mixed track, record the indexes
                scan_true_idx = 0;
                while (scan_true_idx < nphi){
                    if (allSolutions_flag[i][scan_true_idx]){
                        true_segments_info[ntrue_segments][0] = i;
                        true_segments_info[ntrue_segments][1] = scan_true_idx;
                        while (scan_true_idx < nphi && allSolutions_flag[i][scan_true_idx] == 1){
                            scan_true_idx += 1;   
                        }
                        true_segments_info[ntrue_segments][2] = scan_true_idx - 1;
                        ntrue_segments += 1;                        
                    }else{
                        while (scan_true_idx < nphi && allSolutions_flag[i][scan_true_idx] == 0){
                            scan_true_idx += 1;
                        }
                    }
                }
            }
        }
#ifdef VERBOSE
        fprintf(stderr, "can you 6200, ntrue_segments = %d\n", ntrue_segments);
#endif
        // # now, connect segments into closed track
        //# if a segment has length == npoints, then pass, do not need to handle this, we just need to add this segment to the closed_image_info
        // # otherwise, you need to connect the current segment
        // # by the way, you need to record how many segments left, to be connected.

        // # initialize array to store the information of true-close-image boundaries
        // # there are at most DEGREE true-close-image boundaries
        nclosed_image = 0;
        // # closed_image_info = np.zeros(( min(DEGREE, ntrue_segments), 2 + 3 * (ntrue_segments - npure_close_segments + 1) )) # at most DEGREE closed-images, #### npure_close_segments is wrong number, not all segments with length nphi is closed
        //closed_image_info = np.zeros(( min(DEGREE, ntrue_segments), 2 + 3 * (ntrue_segments + 1) ));// # at most DEGREE closed-images, 
        // # each row contains how this closed-image is built from segments, i.e.,
        // # (parity, from_n_segments, a, b, c, ...)
        // # a is the allSolutions index where the segments belongs to
        // # b, c is the head_index, tail_index (if b>c, then means this segment is actually reversed)
        nfinal_closed_image = 0;
        //already_done_segments = np.zeros(ntrue_segments)
        for (i=0;i<DEGREE*DEGREE;i++){already_done_segments[i] = false;}

        open_seg_leftover = ntrue_segments;

        //if VERBOSE: print("true_segments_info before sort by length = ", true_segments_info[:ntrue_segments,:])
        // # do we need to order by the length of the segments in ntrue_segments?
        // # bubble sort? // pros: when create new image track, you always start with the longest one among the remaining segments
        //tmp_segments_length = np.zeros(ntrue_segments);
        for (i=0;i<DEGREE*DEGREE;i++){tmp_segments_length[i] = 0;}


        for (i=0;i<ntrue_segments;i++){
            tmp_segments_length[i] = (true_segments_info[i][2] - true_segments_info[i][1]) + 1;
        }


        //################################### rewrite below two lines
        //sorted_lenidx = np.argsort(tmp_segments_length)[::-1]
        // true_segments_info = true_segments_info[sorted_lenidx]
        
        // below is not neccessary
        // for index from 0 to ntrue_segments-1, argsort by tmp_segments_length
        if(false){
        myargsort(tmp_segments_length, sorted_lenidx, ntrue_segments, 1);
        for (i=0; i<ntrue_segments; i++){ temp_true_segments_info[i] = true_segments_info[sorted_lenidx[i]]; }
        // for (i=0; i<ntrue_segments; i++){ true_segments_info[i] = temp_true_segments_info[i]; } // array type 'unsigned short [3]' is not assignable
        }//###################################


        // if VERBOSE: print("ntrue_segments = ", ntrue_segments)
        // if VERBOSE: print(true_segments_info[:ntrue_segments,:])
        // if VERBOSE: # print info of head tail at each segments
        //     for iv in range(ntrue_segments):
        //         print(">>>", iv, "true_segments_info: ", true_segments_info[iv, :], "x, y, mu, flag, absdzs")
        //         jv, hidv, tidv = true_segments_info[iv, :]
        //         print("head ", allSolutions_x[jv, hidv], allSolutions_y[jv, hidv], allSolutions_mu[jv, hidv], allSolutions_flag[jv, hidv], allSolutions_absdzs[jv, hidv])
        //         print("tail ", allSolutions_x[jv, tidv], allSolutions_y[jv, tidv], allSolutions_mu[jv, tidv], allSolutions_flag[jv, tidv], allSolutions_absdzs[jv, tidv])


        // ### first find all closed tracks
        if_creat_new = true;
        for (i=0; i<ntrue_segments; i++){
            if (already_done_segments[i] || open_seg_leftover <= 0) {continue;}
            j = true_segments_info[i][0]; hid = true_segments_info[i][1]; tid = true_segments_info[i][2];//j, hid, tid = true_segments_info[i, :]
            head1[0] = j; head1[1] = hid; tail1[0] = j; tail1[1] = tid; //head1, tail1 = [j, hid], [j, tid]

            if (hid == 0 && tid == nphi - 1 && head_tail_close(head1, tail1)<EPS2){
                closed_image_info[nfinal_closed_image][0] = allSolutions_mu[j][hid] > 0?1:-1; //1 if allSolutions_mu[j, hid] > 0 else -1 // # final->first->mu > 0 ? 1 : -1;
                closed_image_info[nfinal_closed_image][1] = 1; //# this closed image is originate from 1 segment
                closed_image_info[nfinal_closed_image][2] = j; //# which allSolutions_x index this image belongs to
                closed_image_info[nfinal_closed_image][3] = hid;
                closed_image_info[nfinal_closed_image][4] = tid;

                already_done_segments[i] = 1;
                open_seg_leftover -= 1;
                nfinal_closed_image += 1;

                // if VERBOSE: print(i, ">>>>>> initialize new segments, start from seg ", i, "open_seg_leftover = ", open_seg_leftover, "hid, tid = ", hid, tid)
                // // # to check whether a track is closed or not
                // if VERBOSE: print(i, "close track", "head", allSolutions_x[head1[0], head1[1]], allSolutions_y[head1[0], head1[1]], allSolutions_x[tail1[0], tail1[1]], allSolutions_y[tail1[0], tail1[1]],)
            }
        }


        if_creat_new = true;
        continue_left = 1;
        while (open_seg_leftover > 0){
            connectEPS = EPS2;
            for (i=0;i<ntrue_segments;i++){
                if (already_done_segments[i] || open_seg_leftover <= 0){continue;}
                j = true_segments_info[i][0]; hid = true_segments_info[i][1]; tid = true_segments_info[i][2]; //j, hid, tid = true_segments_info[i, :]

                if (if_creat_new){
                    if (hid == tid){ //# and open_seg_leftover > 1:
                        //if VERBOSE: print("<<<<<< segment %d has only 1 data point, we do not start from this, continue_left = %d"%(i, continue_left))
                        if (continue_left > 0){
                            continue_left -= 1;
                            break; //# because you have sorted segments by length, so length == 1 means there are only
                        }else{ //# there are only length = 1 segments, you only need to attach them to the head/tail of the existing head1, tail1
                            connectEPS = pow(allSolutions_x[head1[0]][head1[1]] - allSolutions_x[tail1[0]][tail1[1]], 2)+pow(allSolutions_y[head1[0]][head1[1]] - allSolutions_y[tail1[0]][tail1[1]], 2); }
                    }else{
                        head1[0] = j; head1[1] = hid; tail1[0] = j; tail1[1] = tid; //head1, tail1 = [j, hid], [j, tid]

                        closed_image_info[nfinal_closed_image][0] = allSolutions_mu[j][hid] > 0?1:-1;// 1 if allSolutions_mu[j, hid] > 0 else -1 //# final->first->mu > 0 ? 1 : -1;
                        closed_image_info[nfinal_closed_image][1] = 1;// # this closed image is originate from 1 segment
                        closed_image_info[nfinal_closed_image][2] = j;// # which allSolutions_x index this image belongs to
                        closed_image_info[nfinal_closed_image][3] = hid;
                        closed_image_info[nfinal_closed_image][4] = tid;

                        already_done_segments[i] = 1;
                        open_seg_leftover -= 1;

                        if (if_creat_new){
                            #ifdef VERBOSE
                            printf(">>>>>> initialize new segments, start from seg %d, open_seg_leftover = %d, hid, tid = %d %d", i, open_seg_leftover, hid, tid);
                            #endif
                        }
                        if_creat_new = false;
                    }
                }
                // # if hid == 0 and tid == nphi - 1 and self.head_tail_close(head1, tail1, allSolutions_x, allSolutions_y)<EPS_2:
                if ((head_tail_close(head1, tail1)<EPS2) && (connectEPS == EPS2)){
                    // # to check whether a track is closed or not
                    // printf("%d, close track head %f %f %f %f", i, allSolutions_x[head1[0]][head1[1]], allSolutions_y[head1[0]][head1[1]], allSolutions_x[tail1[0]][tail1[1]], allSolutions_y[tail1[0]][tail1[1]]);
                    nfinal_closed_image += 1;
                    if_creat_new = true;
                    continue;
                }else{
                    ifcontinue = false;
                    // # now, test whether we can connect segment i with "other segments"
                    bestconnect_type = 0;
                    bestconnect_dis = 1e10;
                    bestconnect_i2 = 0;
                    canweconnect = false;
                    for (i2 =0; i2<ntrue_segments; i2++){
                        if (open_seg_leftover > 0 && (!already_done_segments[i2]) ){
                            // # # judge whether two segments (i, i2) can be connected together
                            // j2, hid2, tid2 = true_segments_info[i2, :]
                            j2 = true_segments_info[i2][0]; hid2 = true_segments_info[i2][1]; tid2 = true_segments_info[i2][2];
                            head2[0] = j2; head2[1] = hid2; tail2[0] = j2; tail2[1] = tid2;// head2, tail2 = [j2, hid2], [j2, tid2]
                            // # call a function to test whether these two segments can be connected
                            // itype, connectdis = self.if_two_segments_connect(head1, tail1, head2, tail2, allSolutions_x, allSolutions_y, connectEPS)
                            // printf("in 6341 before check connect\n");
                            if_two_segments_connect(head1, tail1, head2, tail2, connectEPS, &itype, &connectdis);
                            // printf("in 6343 after check connect, %d, %f \n", itype, connectdis);
                            if (itype > 0){
                                bestconnect_type = itype;
                                bestconnect_dis = connectdis;
                                bestconnect_i2 = i2;
                                canweconnect = true;
                                break;
                            }else if (itype < 0){
                                canweconnect = true;
                                if (bestconnect_dis > connectdis){
                                    bestconnect_dis = connectdis;
                                    bestconnect_type = itype;
                                    bestconnect_i2 = i2;
                                }
                            }
                        }
                    }
                    
                    // if canweconnect and ((bestconnect_type > 0) or (bestconnect_type<0) and open_seg_leftover == 1 ):
                    if (canweconnect && ((bestconnect_type > 0) || ((bestconnect_type<0) && open_seg_leftover == 1 ))){
                        //j2, hid2, tid2 = true_segments_info[bestconnect_i2, :]
                        j2 = true_segments_info[bestconnect_i2][0]; hid2 = true_segments_info[bestconnect_i2][1]; tid2 = true_segments_info[bestconnect_i2][2];
                        head2[0] = j2; head2[1] = hid2; tail2[0] = j2; tail2[1] = tid2;

#ifdef VERBOSE
                    printf("\t canweconnect %d, bestconnect_type %d, bestconnect_dis = %e, bestconnect_i2 = %d\n", canweconnect, bestconnect_type, bestconnect_dis, bestconnect_i2);
#endif
                        // if VERBOSE: print("\t head2 = ", allSolutions_x[head2[0], head2[1]], allSolutions_y[head2[0], head2[1]], "tail2 = ", allSolutions_x[tail2[0], tail2[1]], allSolutions_y[tail2[0], tail2[1]])

                        // # we can connect seg_i with seg_i2
                        existing_seg_n = closed_image_info[nfinal_closed_image][1];
                        offset = 3*existing_seg_n;

                        if (fabs(bestconnect_type) == 3){ //# tail connect with head
                            closed_image_info[nfinal_closed_image][2 + offset] = j2;
                            closed_image_info[nfinal_closed_image][3 + offset] = hid2;
                            closed_image_info[nfinal_closed_image][4 + offset] = tid2;
                            tail1[0] = tail2[0]; tail1[1] = tail2[1]; //tail1 = tail2;
                        }else if (fabs(bestconnect_type) == 4){ //# tail connect with tail
                            closed_image_info[nfinal_closed_image][2 + offset] = j2;
                            closed_image_info[nfinal_closed_image][3 + offset] = tid2;
                            closed_image_info[nfinal_closed_image][4 + offset] = hid2;
                            tail1[0] = head2[0]; tail1[1] = head2[1]; //tail1 = head2;
                        }else if (fabs(bestconnect_type) == 2){ //# head connect with tail
                            // # move prev segments behind
                            // closed_image_info[nfinal_closed_image][5: 5+offset] = closed_image_info[nfinal_closed_image][2: 2+offset]
                            for(i4 = 1+offset; i4>1; i4--){closed_image_info[nfinal_closed_image][i4+3] = closed_image_info[nfinal_closed_image][i4];}
                            closed_image_info[nfinal_closed_image][2] = j2;
                            closed_image_info[nfinal_closed_image][3] = hid2;
                            closed_image_info[nfinal_closed_image][4] = tid2;
                            head1[0] = head2[0]; head1[1] = head2[1]; //head1 = head2;
                        }else if (fabs(bestconnect_type) == 1){ //# head connect with head
                            // # move prev segments behind
                            // closed_image_info[nfinal_closed_image][5: 5+offset] = closed_image_info[nfinal_closed_image][2: 2+offset]
                            for(i4 = 1+offset; i4>1; i4--){closed_image_info[nfinal_closed_image][i4+3] = closed_image_info[nfinal_closed_image][i4];}
                            closed_image_info[nfinal_closed_image][2] = j2;
                            closed_image_info[nfinal_closed_image][3] = tid2;
                            closed_image_info[nfinal_closed_image][4] = hid2;
                            head1[0] = tail2[0]; head1[1] = tail2[1]; //head1 = tail2;
                        }

                        closed_image_info[nfinal_closed_image][1] += 1;
                        already_done_segments[i2] = 1;
                        open_seg_leftover -= 1;
                        // #continue # you have connect seg2 with seg1, now proceed to the next segment
                        // # check whether current is already closed
                        if (head_tail_close(head1, tail1)<EPS2 || (open_seg_leftover<=0 && connectEPS == EPS2)){
                            nfinal_closed_image += 1;
                            if_creat_new = true;
                            break;
                        }else{
                            //if VERBOSE: print("break 322, open_seg_leftover, nfinal_closed_image, if_creat_new = ", open_seg_leftover, nfinal_closed_image, if_creat_new)
                            ifcontinue = true;// # connect once, then try whether we can connect again
                            // # break
                        }
                    }

                    if (ifcontinue){continue;} //# if connected above, then continue, otherwise try jump

                    // if VERBOSE: print("330 continue")
                    // # for i2 in range(ntrue_segments):
                    // #     if open_seg_leftover > 0 and (not already_done_segments[i2]):
                    // # now, test whether we can jump from segment i to "other segments"
                    // # different from connect, jump is more trivial, need to find the best place to jump
                    canwejump = false;
                    bestjumptype = 0;
                    bestjump_i2 = 0;
                    bestjump_fac_mu = 1e10;// # find the minimum how_close, shoule be very close to 2.0
                    bestjump_dis = 1e10;

                    // if VERBOSE:
                    //     print("\t before test jump, already_done_segments = ", already_done_segments)
                    //     print("head1, tail1", head1, tail1)
                    //     print("\tcurrent head: x, y, mu, flag, srcx, srcy = ", allSolutions_x[head1[0],head1[1]], allSolutions_y[head1[0],head1[1]],allSolutions_mu[head1[0],head1[1]],allSolutions_flag[head1[0],head1[1]],allSolutions_srcx[head1[0],head1[1]], allSolutions_srcy[head1[0],head1[1]])
                    //     print("\tcurrent tail: x, y, mu, flag, srcx, srcy = ", allSolutions_x[tail1[0],tail1[1]], allSolutions_y[tail1[0],tail1[1]],allSolutions_mu[tail1[0],tail1[1]],allSolutions_flag[tail1[0],tail1[1]],allSolutions_srcx[tail1[0],tail1[1]], allSolutions_srcy[tail1[0],tail1[1]])  
                    
                    for (i3=0; i3<ntrue_segments;i3++){ // in range(ntrue_segments):
                        if (open_seg_leftover > 0 && (!already_done_segments[i3])){
                            // # # judge whether two segments (i, i2) can be connected together
                            // j2, hid2, tid2 = true_segments_info[i3, :]
                            // head2, tail2 = [j2, hid2], [j2, tid2]
                            j2 = true_segments_info[i3][0]; hid2 = true_segments_info[i3][1]; tid2 = true_segments_info[i3][2];
                            head2[0] = j2; head2[1] = hid2; tail2[0] = j2; tail2[1] = tid2;

                            // # call a function to test whether these two segments can be connected
                            // if VERBOSE: print(">>>>>> i3 = %d"%(i3))
                            // itype, how_close, jumpdis = self.if_two_segments_jump(head1, tail1, head2, tail2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu)
                            if_two_segments_jump(head1, tail1, head2, tail2, &itype, &how_close, &jumpdis);
                            if (bestjump_fac_mu > how_close && bestjump_dis >= jumpdis){
                                bestjump_fac_mu = how_close;
                                bestjump_dis = jumpdis;
                                bestjump_i2 = i3;
                                bestjumptype = itype;
                            }
                            if (itype > 0 && jumpdis == 0){// # jumpdis has an upper limit: the mimimum distances between two sample points in the source limb, dphi
                                bestjump_fac_mu = how_close;
                                bestjump_dis = jumpdis;
                                bestjump_i2 = i3;
                                bestjumptype = itype;

                                canwejump = true;
                                break;
                                // # if bestjump_fac_mu > how_close and bestjump_dis >= jumpdis:
                                // #     bestjump_fac_mu = how_close
                                // #     bestjump_dis = jumpdis
                                // #     bestjump_i2 = i3
                                // #     bestjumptype = itype
                            }else{
                                // if VERBOSE:
                                //     print("\t440, we can not jump, itype = %d the bestjump_fac_mu = %f, bestjump_dis = %e"%(itype, bestjump_fac_mu, bestjump_dis))
                                //     print("\thead2, tail2", head2, tail2)                               
                                //     print("\ttmp seg %d head: x, y, mu, flag, srcx, srcy = "%i3, allSolutions_x[head2[0],head2[1]], allSolutions_y[head2[0],head2[1]],allSolutions_mu[head2[0],head2[1]],allSolutions_flag[head2[0],head2[1]],allSolutions_srcx[head2[0],head2[1]], allSolutions_srcy[head2[0],head2[1]])
                                //     print("\ttmp seg %d tail: x, y, mu, flag, srcx, srcy = "%i3, allSolutions_x[tail2[0],tail2[1]], allSolutions_y[tail2[0],tail2[1]],allSolutions_mu[tail2[0],tail2[1]],allSolutions_flag[tail2[0],tail2[1]], allSolutions_srcx[tail2[0],tail2[1]], allSolutions_srcy[tail2[0],tail2[1]])
                            }
                        }
                    }
                    if (canwejump){
                        // #print(i, "best jump to ", bestjump_i2, "bestjump_fac_mu = ", bestjump_fac_mu)
                        // # can find some place to jump, 
                        j2 = true_segments_info[bestjump_i2][0]; hid2 = true_segments_info[bestjump_i2][1]; tid2 = true_segments_info[bestjump_i2][2];
                        head2[0] = j2; head2[1] = hid2; tail2[0] = j2; tail2[1] = tid2;


                        // if VERBOSE:
                        //     print(">>> best jump to ", bestjump_i2, "bestjump_fac_mu = ", bestjump_fac_mu, 'type = ', bestjumptype, allSolutions_x[head2[0], head2[1]], allSolutions_y[head2[0], head2[1]], allSolutions_x[tail2[0], tail2[1]], allSolutions_y[tail2[0], tail2[1]], "bestjump_dis = ", bestjump_dis)           

                        //     print("\thead2, tail2", head2, tail2)                               
                        //     print("\ttmp seg %d head: x, y, mu, flag, srcx, srcy = "%i3, allSolutions_x[head2[0],head2[1]], allSolutions_y[head2[0],head2[1]],allSolutions_mu[head2[0],head2[1]],allSolutions_flag[head2[0],head2[1]],allSolutions_srcx[head2[0],head2[1]], allSolutions_srcy[head2[0],head2[1]])
                        //     print("\ttmp seg %d tail: x, y, mu, flag, srcx, srcy = "%i3, allSolutions_x[tail2[0],tail2[1]], allSolutions_y[tail2[0],tail2[1]],allSolutions_mu[tail2[0],tail2[1]],allSolutions_flag[tail2[0],tail2[1]], allSolutions_srcx[tail2[0],tail2[1]], allSolutions_srcy[tail2[0],tail2[1]])


                        // # we can connect seg_i with seg_i2
                        existing_seg_n = closed_image_info[nfinal_closed_image][1];
                        offset = 3*existing_seg_n;

                        if (fabs(bestjumptype) == 3){ //# tail connect with head
                            closed_image_info[nfinal_closed_image][2 + offset] = j2;
                            closed_image_info[nfinal_closed_image][3 + offset] = hid2;
                            closed_image_info[nfinal_closed_image][4 + offset] = tid2;
                            tail1[0] = tail2[0]; tail1[1] = tail2[1]; //tail1 = tail2;
                        }else if (fabs(bestjumptype) == 4){ //# tail connect with tail
                            closed_image_info[nfinal_closed_image][2 + offset] = j2;
                            closed_image_info[nfinal_closed_image][3 + offset] = tid2;
                            closed_image_info[nfinal_closed_image][4 + offset] = hid2;
                            tail1[0] = head2[0]; tail1[1] = head2[1]; //tail1 = head2;
                        }else if (fabs(bestjumptype) == 2){ //# head connect with tail
                            // # move prev segments behind
                            // closed_image_info[nfinal_closed_image][5: 5+offset] = closed_image_info[nfinal_closed_image][2: 2+offset]
                            for(i4 = 1+offset; i4>1; i4--){closed_image_info[nfinal_closed_image][i4+3] = closed_image_info[nfinal_closed_image][i4];}
                            closed_image_info[nfinal_closed_image][2] = j2;
                            closed_image_info[nfinal_closed_image][3] = hid2;
                            closed_image_info[nfinal_closed_image][4] = tid2;
                            head1[0] = head2[0]; head1[1] = head2[1]; //head1 = head2;
                        }else if (fabs(bestjumptype) == 1){ //# head connect with head
                            // # move prev segments behind
                            // closed_image_info[nfinal_closed_image][5: 5+offset] = closed_image_info[nfinal_closed_image][2: 2+offset]
                            for(i4 = 1+offset; i4>1; i4--){closed_image_info[nfinal_closed_image][i4+3] = closed_image_info[nfinal_closed_image][i4];}
                            closed_image_info[nfinal_closed_image][2] = j2;
                            closed_image_info[nfinal_closed_image][3] = tid2;
                            closed_image_info[nfinal_closed_image][4] = hid2;
                            head1[0] = tail2[0]; head1[1] = tail2[1]; //head1 = tail2;
                        }
                        closed_image_info[nfinal_closed_image][1] += 1;
                        already_done_segments[bestjump_i2] = 1;
                        open_seg_leftover -= 1;


                        // # check whether current is already closed
                        if (head_tail_close(head1, tail1)<EPS2 || open_seg_leftover<=0){
                            nfinal_closed_image += 1;
                            if_creat_new = true;
                        }else{
                            // if VERBOSE: print("break 398, open_seg_leftover, nfinal_closed_image", open_seg_leftover, nfinal_closed_image)
                            // if VERBOSE: print(already_done_segments)
                            continue;// # after jump, not close, and open_seg_leftover > 0
                        }
                    }else{
                        // # we can not connect, and we can not jump, we need to test whether current
                        // # we need to test whether it is closed, other wise, we regard this as a close image anyway
                        // # new test whether head1, tail1 is close, and whether there is no open segments left
                        // #if self.head_tail_close(head1, tail1, allSolutions_x, allSolutions_y, EPS_2):
                        // if VERBOSE: print("410, open_seg_leftover", open_seg_leftover)
                        if (open_seg_leftover > 0){
                            nfinal_closed_image += 1;
                            if_creat_new = true;
                        }
                    }
                }
            }
        }
        //if VERBOSE: print("nfinal_closed_image = ", nfinal_closed_image)
        //if VERBOSE: print(closed_image_info[:nfinal_closed_image,:].astype(int))

        //return true_segments_info[:ntrue_segments,:].astype(int), closed_image_info[:nfinal_closed_image,:].astype(int)

}


void TripleLensing::if_head_tail_jump(int head1[], int tail1[], bool *ifjumpbool, double *how_close, double *dis){
        srcxprev = allSolutions_srcx[head1[0]][head1[1]];
        srcyprev = allSolutions_srcy[head1[0]][head1[1]];
        srcxpost = allSolutions_srcx[tail1[0]][tail1[1]];
        srcypost = allSolutions_srcy[tail1[0]][tail1[1]];
        muprev = allSolutions_mu[head1[0]][head1[1]];
        mupost = allSolutions_mu[tail1[0]][tail1[1]];

        *how_close = ( fabs(muprev / mupost) + fabs(mupost/muprev) );
        *dis = pow(srcxprev - srcxpost, 2) + pow(srcyprev - srcypost, 2); //self.if_dis_close(srcxprev, srcyprev, srcxpost, srcypost);
        *ifjumpbool = false;

        if (muprev * mupost < 0.0){
            if (head1[1] == tail1[1]){
                *dis = 0;
                *ifjumpbool = true;
            }else{
                *ifjumpbool = (*how_close < 2.5);
            }
        }
        // return ifjump, how_close, dis
}

void TripleLensing::if_two_segments_jump(int head1[], int tail1[],int head2[], int tail2[], int *itype, double *bestjump_fac_mu, double *bestjump_dis){
    //printf("in 6581, if_two_segments_jump \n");

        *itype = 0;
        // # we prefer connect, rather than jump, so, first test image connectivity; if we cannot connect, then try whether we can "jump"
        // # besides, you need to have a seperate function for jump
        // # we prefer to find the easiest scenario
        *itype = 0;
        *bestjump_fac_mu = 1e10;
        *bestjump_dis= 1e10;
        double dis;

        // ifjump, how_close, dis = self.if_head_tail_jump( tail1, head2, allSolutions_srcx, allSolutions_srcy, allSolutions_mu)
        if_head_tail_jump(tail1, head2, &ifjumpbool, &how_close, &dis);
        if (dis == 0 || (ifjumpbool && (*bestjump_fac_mu > how_close) && (*bestjump_dis >= dis))){
            *bestjump_fac_mu = how_close;
            *bestjump_dis = dis;
            *itype = 3;// # tail connect with head
        }

        if_head_tail_jump(tail1, tail2, &ifjumpbool, &how_close, &dis);
        if (dis == 0 || (ifjumpbool && (*bestjump_fac_mu > how_close) && (*bestjump_dis >= dis))){
            *bestjump_fac_mu = how_close;
            *bestjump_dis = dis;
            *itype = 4;// # tail connect with tail
        }

        if_head_tail_jump(head1, tail2, &ifjumpbool, &how_close, &dis);
        if (dis == 0 || (ifjumpbool && (*bestjump_fac_mu > how_close) && (*bestjump_dis >= dis))){
            *bestjump_fac_mu = how_close;
            *bestjump_dis = dis;
            *itype = 2;// # head connect with tail
        }

        if_head_tail_jump(head1, head2, &ifjumpbool, &how_close, &dis);
        if (dis == 0 || (ifjumpbool && (*bestjump_fac_mu > how_close) && (*bestjump_dis >= dis))){
            *bestjump_fac_mu = how_close;
            *bestjump_dis = dis;
            *itype = 1;// # head connect with head
        }

        // if besttype == 3 and VERBOSE: print("tail jump to head", besttype, tail1, head2, bestjump_fac_mu, bestjump_dis)
        // if besttype == 4 and VERBOSE: print("tail jump to tail", besttype, tail1, tail2, bestjump_fac_mu,bestjump_dis)
        // if besttype == 2 and VERBOSE: print("head jump to tail", besttype, head1, tail2, bestjump_fac_mu,bestjump_dis)
        // if besttype == 1 and VERBOSE: print("head jump to head", besttype, head1, head2, bestjump_fac_mu,bestjump_dis)
        // return besttype, bestjump_fac_mu, bestjump_dis
}

double TripleLensing::head_tail_close(int head1[], int tail1[]){
    double x1, y1, x2, y2;
    x1 = allSolutions_x[head1[0]][head1[1]];
    y1 = allSolutions_y[head1[0]][head1[1]];
    x2 = allSolutions_x[tail1[0]][tail1[1]];
    y2 = allSolutions_y[tail1[0]][tail1[1]];
    //#print("x1, y1, x2, y2", x1, y1, x2, y2)
    return pow(x1 - x2, 2) + pow(y1 - y2, 2);
}

void TripleLensing::if_two_segments_connect(int head1[], int tail1[],int head2[], int tail2[], double EPS_2, int *itype, double *dis){
    // printf("in 6637, if_two_segments_connect \n");
        // # head1 = [head1_j, head1_idx]
        // #print(head1, tail1, head2, tail2)
        // # connect seg1 and seg2
        // # itype 0 = can not connect
        // # itype 1 = head head --> how to handle this? --> you need to insert seg2 into the front of seg1, and reverse seg2
        // # itype 2 = head tail --> how to handle this? --> you need to insert seg2 into the front of seg1, do not need to reverse seg2
        // # itype 3 = tail head --> how to handle this? --> the easiest scenerio, just add seg2 behind seg1
        // # itype 4 = tail tail --> how to handle this? --> add seg2 behind seg1, and reverse seg2
        *itype = 0;
        // # we prefer connect, rather than jump, so, first test image connectivity; if we cannot connect, then try whether we can "jump"
        // # besides, you need to have a seperate function for jump
        int besttype = 0;
        double bestdis = 1e10;

        // # we prefer to find the easiest scenario
        *dis = head_tail_close(tail1, head2);
        if (*dis < EPS_2){ //# either they are overlapping or they are continuous in terms of phi
            *itype = 3;// # tail connect with head
            // if VERBOSE: print("tail connect with head", itype, tail1, head2, allSolutions_x[tail1[0], tail1[1]], allSolutions_y[tail1[0], tail1[1]] )
            return;// itype, dis
        }else if (fabs(tail1[1] - head2[1]) == 1){ //: ## if you can not connect by distance, you need to check whether the phi is continuous
            if (bestdis > *dis){
                bestdis = *dis;
                besttype = -3;
            }
        }

        *dis = head_tail_close( tail1, tail2);
        if (*dis < EPS_2){ //# either they are overlapping or they are continuous in terms of phi
            *itype = 4;// # tail connect with head
            // if VERBOSE: print("tail connect with head", itype, tail1, head2, allSolutions_x[tail1[0], tail1[1]], allSolutions_y[tail1[0], tail1[1]] )
            return;// itype, dis
        }else if (fabs(tail1[1] - tail2[1]) == 1){ //: ## if you can not connect by distance, you need to check whether the phi is continuous
            if (bestdis > *dis){
                bestdis = *dis;
                besttype = -4;
            }
        }

        *dis = head_tail_close( head1, tail2);
        if (*dis < EPS_2){ //# either they are overlapping or they are continuous in terms of phi
            *itype = 2;// # tail connect with head
            // if VERBOSE: print("tail connect with head", itype, tail1, head2, allSolutions_x[tail1[0], tail1[1]], allSolutions_y[tail1[0], tail1[1]] )
            return;// itype, dis
        }else if (fabs(head1[1] - tail2[1]) == 1){ //: ## if you can not connect by distance, you need to check whether the phi is continuous
            if (bestdis > *dis){
                bestdis = *dis;
                besttype = -2;
            }
        }

        *dis = head_tail_close( head1, head2);
        if (*dis < EPS_2){ //# either they are overlapping or they are continuous in terms of phi
            *itype = 1;// # tail connect with head
            // if VERBOSE: print("tail connect with head", itype, tail1, head2, allSolutions_x[tail1[0], tail1[1]], allSolutions_y[tail1[0], tail1[1]] )
            return;// itype, dis
        }else if (fabs(head1[1] - head2[1]) == 1){ //: ## if you can not connect by distance, you need to check whether the phi is continuous
            if (bestdis > *dis){
                bestdis = *dis;
                besttype = -1;
            }
        }

        if (besttype == 0){
            //# if we can not connect by overlapping, we still have chance to connect two points with abs(phi_index1 - phi_index2) == 1
            *itype = 0; *dis = 1e10;
            return;// itype, 1e10
        }else{
            //if VERBOSE: print("possible connect by continuous phi, besttype = %d, bestdis = %e"%(besttype, bestdis))
            *itype = besttype; *dis = bestdis;
            return;// besttype, bestdis
        }
}


void myargsort(int arr[], int index[], int n, int mode){
    //# mode 1 降序
    //# n = len(arr)
    // index = list(range(n))
    int i, j;
    for (i=0; i<n; i++){
        index[i] = i;
    }
    for (i=0; i<n-1; i++){
        for (j=i+1; j<n; j++){
            if (mode == 1){
                if (arr[i] < arr[j]){
                    myswap(&arr[i], &arr[j]); //arr[i], arr[j] = arr[j], arr[i]
                    myswap(&index[i], &index[j]); //index[i], index[j] = index[j], index[i]
                }
            } else if (mode == 0){
                if (arr[i] > arr[j]){
                    myswap(&arr[i], &arr[j]);             
                    myswap(&index[i], &index[j]);             
                }
            }
        }
    }
}

void mydbswap(double *a, double *b){
    double t;
    t = *a;
    *a = *b;
    *b = t;
}

void myswap(int *a, int *b){
    int t;
    t = *a;
    *a = *b;
    *b = t;
}

void TripleLensing::saveimage(int curr_nimages, int curr_total_parity, int j, bool *ifsave, int *saveindex){
        // # find one true solutions that has been classified as false solution (due to the simple criterion absdzs > EPS)
        *ifsave = false;
        *saveindex = 0;
        double minabsdzs = 1e10; // this variable can be defined as a class member
        for (int i=0; i<DEGREE; i++){
            if (allSolutions_flag[i][j] == 0){ // and allSolutions_mu[i, j] != -1e10:
                if (allSolutions_absdzs[i][j] < minabsdzs){
                    minabsdzs = allSolutions_absdzs[i][j];
                    *saveindex = i;
                }
            }
        }
        // # after save the solution, whether the nimages and parity are correct
        if (((curr_nimages + 1 - NLENS)%2 == 1) && ((curr_total_parity+one_parity(allSolutions_mu[*saveindex][j])) == CORRECT_PARITY)){
            *ifsave = true;
        }
}


void TripleLensing::get_closest_idx(double preProv_x[], double preProv_y[], double Prov_x[], double Prov_y[], unsigned short int res_idx[]){
    // # for each point i (Prov_x[i], Prov_y[i]), find the index in preProv_x, preProv_y so that the point j is closest to i
    // has_removed = np.zeros(self.DEGREE)
    // res_idx = np.zeros(self.DEGREE).astype(int)
    // #return np.arange(0, 10).astype(int)
    for(int i=0; i<DEGREE; i++){
        has_removed[i] = false;
    }
    double mindis, dis;
    unsigned short int minidx;

    for (int i=0; i<DEGREE; i++){
        currx = preProv_x[i];
        curry = preProv_y[i];
        mindis = 1e100;
        minidx = 0;
        for (int j=0; j<DEGREE; j++){
            if (!has_removed[j]){
                dis = pow(currx - Prov_x[j], 2) + pow(curry - Prov_y[j], 2);
                if (dis<mindis){
                    mindis = dis;
                    minidx = j;
                }
            }
        }
        res_idx[i] = minidx;
        has_removed[minidx] = 1;
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
    // printf("zsc = %f %f\n", zsc.re, zsc.im);
    for (i = 0; i < NLENS; i++) {
        // printf("%d, zc[i] = %f %f\n",i, zc[i].re, zc[i].im);

        /* denominator */
        for (j = 0; j <= NLENS; j++) { /* (zsc-conjugate(z[i])) * Product_{j=1}^N (z-z[j]) */
            q[j][i] = (zsc - zc[i]) * temp_const1[i][j];
            // printf("j, i = %d %d, q[j][i] = %f %f\n",j, i, q[j][i].re, q[j][i].im);
        }
        //printf("q[j][i] = %f %f\n", q[j][i].re, q[j][i].im);

        /* Sum_{j=1}^n Product (z-z_k), k=1, n, but k !=j. */
        for (j = 0; j < NLENS; j++) {
            /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */

            /* doing the sum */
            for (k = 0; k < NLENS; k++) {
                q[k][i] = q[k][i] + mlens[j] * temp_const2[i][j][k];
                // q[k][i] = q[k][i] + mlens[j] * temp[k];
            }
        }
        // printf("k, i = %d, %d, q[k][i] = %f %f, temp_const2[i][j][k] %f %f\n", k, i, q[k][i].re, q[k][i].im, temp_const2[i][j][k].re, temp_const2[i][j][k].re);

    }
    // printf("can you 5961\n"); // ok


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
                qtemp2[k] = q[k][j]; // mark
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

// void multiply_z_v2(complex c[][NLENS + 1], complex a, int n, int firstdim)
void multiply_z_v2(complex **c, complex a, int n, int firstdim)
{
    int j;
    // printf("firstdim, n, a = %d, %d, %f, %f\n", firstdim, n, a.re, a.im);

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


int one_parity(double mu){
    if (mu == -1e10){
        return 0;}
    return mu>0?1:-1;
}


void TripleLensing::cmplx_roots_gen(complex *roots, complex *poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points) {
    //roots - array which will hold all roots that had been found.
    //If the flag 'use_roots_as_starting_points' is set to
    //.true., then instead of point(0, 0) we use value from
    //this array as starting point for cmplx_laguerre

    //poly - is an array of polynomial cooefs, length = degree + 1,
    //poly[0] x ^ 0 + poly[1] x ^ 1 + poly[2] x ^ 2 + ...

    //degree - degree of the polynomial and size of 'roots' array

    //polish_roots_after - after all roots have been found by dividing
    //original polynomial by each root found,
    //you can opt in to polish all roots using full
    //polynomial

    //use_roots_as_starting_points - usually we start Laguerre's
    //method from point(0, 0), but you can decide to use the
    //values of 'roots' array as starting point for each new
    //root that is searched for.This is useful if you have
    //very rough idea where some of the roots can be.
    //

    complex poly2[MAXM];
    static int i, j, n, iter;
    bool success;
    complex coef, prev;

    if (!use_roots_as_starting_points) {
        for (int jj = 0; jj < degree; jj++) {
            roots[jj] = complex(0, 0);
        }
    }

    for (j = 0; j <= degree; j++) poly2[j] = poly[j];

    // Don't do Laguerre's for small degree polynomials
    if (degree <= 1) {
        if (degree == 1) roots[0] = -poly[0] / poly[1];
        return;
    }

    for (n = degree; n >= 3; n--) {
        cmplx_laguerre2newton(poly2, n, &roots[n - 1], iter, success, 2);
        if (!success) {
            roots[n - 1] = complex(0, 0);
            cmplx_laguerre(poly2, n, &roots[n - 1], iter, success);
        }

        // Divide by root
        coef = poly2[n];
        for (i = n - 1; i >= 0; i--) {
            prev = poly2[i];
            poly2[i] = coef;
            coef = prev + roots[n - 1] * coef;
        }
    }


    //Find the to last 2 roots
    solve_quadratic_eq(roots[1], roots[0], poly2);
    //cmplx_laguerre2newton(poly2, 2, &roots[1], iter, success, 2);
    //if (!success) {
    //  solve_quadratic_eq(roots[1], roots[0], poly2);
    //}
    //else {
    //  roots[0] = -(roots[1] + poly2[1] / poly2[2]); // Viete's Formula for the last root
    //}



    if (polish_roots_after) {
        for (n = 0; n < degree; n++) {
            cmplx_newton_spec(poly, degree, &roots[n], iter, success); // Polish roots with full polynomial
        }
    }

    return;
}

void TripleLensing::solve_quadratic_eq(complex &x0, complex &x1, complex *poly) {
    complex a, b, c, b2, delta;
    a = poly[2];
    b = poly[1];
    c = poly[0];
    b2 = b * b;
    delta = sqrt(b2 - 4 * a * c);
    if (real(conj(b)*delta) >= 0) {
        x0 = -0.5 * (b + delta);
    }
    else {
        x0 = -0.5 * (b - delta);
    }
    if (x0 == complex(0., 0.)) {
        x1 = complex(0., 0.);
    }
    else { //Viete's formula
        x1 = c / x0;
        x0 = x0 / a;
    }
    return;

}

void TripleLensing::solve_cubic_eq(complex &x0, complex &x1, complex &x2, complex *poly) {
    //Cubic equation solver for comples polynomial (degree=3)
    //http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
    // poly is an array of polynomial cooefs, length = degree+1, poly[0] is constant
    //  0               1               2           3
    //poly[0] x^0 + poly[1] x^1 + poly[2] x^2 + poly[3] x^3
    complex zeta = complex(-0.5, 0.8660254037844386);
    complex zeta2 = complex(-0.5, -0.8660254037844386);
    double third = 0.3333333333333333;
    complex s0, s1, s2;
    complex E1; //x0+x1+x2
    complex E2; //x0*x1+x1*x2+x2*x0
    complex E3; //x0*x1*x2
    complex A, B, a_1, E12, delta, A2;

    complex val, x;
    a_1 = 1 / poly[3];
    E1 = -poly[2] * a_1;
    E2 = poly[1] * a_1;
    E3 = -poly[0] * a_1;

    s0 = E1;
    E12 = E1 * E1;
    A = 2.0 * E1 * E12 - 9.0 * E1 * E2 + 27.0 * E3;
    B = E12 - 3.0 * E2;
    //quadratic equation z^2 - A * z + B^3 where roots are equal to s1^3 and s2^3
    A2 = A * A;
    delta = sqrt(A2 - 4.0 * (B * B * B));
    if (real(conj(A) * delta) >= 0.0) { // scalar product to decide the sign yielding bigger magnitude
        s1 = cbrt(0.5 * (A + delta));
    }
    else
    {
        s1 = cbrt(0.5 * (A - delta));
    }
    if (s1.re == 0.0 && s1.im == 0.0) {
        s2 = complex(0, 0);
    }
    else {
        s2 = B / s1;
    }

    x0 = third * (s0 + s1 + s2);
    x1 = third * (s0 + s1 * zeta2 + s2 * zeta);
    x2 = third * (s0 + s1 * zeta + s2 * zeta2);

    return;

}

void TripleLensing::cmplx_laguerre(complex *poly, int degree, complex *root, int &iter, bool &success) {
    //Subroutine finds one root of a complex polynomial using
    //Laguerre's method. In every loop it calculates simplified
    //Adams' stopping criterion for the value of the polynomial.
    //
    //Uses 'root' value as a starting point(!!!!!)
    //Remember to initialize 'root' to some initial guess or to
    //point(0, 0) if you have no prior knowledge.
    //
    //poly - is an array of polynomial cooefs
    //
    //length = degree + 1, poly(1) is constant
    //  1              2                3
    //poly(1) x ^ 0 + poly(2) x ^ 1 + poly(3) x ^ 2 + ...
    //
    //degree - a degree of the polynomial
    //
    //root - input: guess for the value of a root
    //output : a root of the polynomial
    //iter - number of iterations performed(the number of polynomial
    //evaluations and stopping criterion evaluation)
    //
    //success - is false if routine reaches maximum number of iterations
    //
    //For a summary of the method go to :
    //http://en.wikipedia.org/wiki/Laguerre's_method
    //
    static int FRAC_JUMP_EVERY = 10;
    const int FRAC_JUMP_LEN = 10;
    double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297,
                                         0.91577881, 0.25921289, 0.50487203,
                                         0.08177045, 0.13653241, 0.306162,
                                         0.37794326, 0.04618805, 0.75132137
                                       }; // some random numbers

    double faq; //jump length
    // double FRAC_ERR = 2.0e-15; //Fractional Error for double precision
    complex p, dp, d2p_half; //value of polynomial, 1st derivative, and 2nd derivative
    // static int i, j, k;
    static int i, k;
    bool good_to_go;
    complex denom, denom_sqrt, dx, newroot;
    double ek, absroot, abs2p;
    complex fac_newton, fac_extra, F_half, c_one_nth;
    double one_nth, n_1_nth, two_n_div_n_1;
    complex c_one = complex(1, 0);
    complex zero = complex(0, 0);
    double stopping_crit2;

    //--------------------------------------------------------------------------------------------

    // //EXTREME FAILSAFE! not usually needed but kept here just to be on the safe side. Takes care of first coefficient being 0
    // if (false) {
    //  if (degree < 0) {
    //      printf("Error: cmplx_laguerre: degree<0");
    //      return;
    //  }
    //  if (poly[degree] == complex(0, 0)) {
    //      if (degree == 0) return;
    //      cmplx_laguerre(poly, degree - 1, root, iter, success);
    //  }
    //  if (degree <= 1) {
    //      if (degree == 0) {
    //          success = false; // we just checked if poly[0] is zero and it isnt
    //          printf("Warning: cmplx_laguerre: degree = 0 and poly[0] does not equal zero, no roots");
    //          return;
    //      }
    //      else {
    //          *root = -poly[0] / poly[1];
    //          return;
    //      }
    //  }
    // } // End of EXTREME failsafe

    good_to_go = false;
    one_nth = 1.0 / degree;
    n_1_nth = (degree - 1.0) * one_nth;
    two_n_div_n_1 = 2.0 / n_1_nth;
    c_one_nth = complex(one_nth, 0.0);
    for (i = 1; i <= MAXIT; i++) {
        ek = abs(poly[degree]); // Preparing stopping criterion
        absroot = abs(*root);
        // Calculate the values of polynomial and its first and second derivatives
        p = poly[degree];
        dp = zero;
        d2p_half = zero;
        for (k = degree - 1; k >= 0; k--) {
            d2p_half = dp + d2p_half * (*root);
            dp = p + dp * *root;
            p = poly[k] + p * (*root); // b_k
            //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
            //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
            //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
            //Eq 8.
            ek = absroot * ek + abs(p);
        }
        iter += 1;

        abs2p = real(conj(p) * p);
        if (abs2p == 0) return;
        stopping_crit2 = pow(FRAC_ERR * ek, 2.0);
        if (abs2p < stopping_crit2) {
            //(simplified a little Eq. 10 of Adams 1967)
            //do additional iteration if we are less than 10x from stopping criterion
            if (abs2p < 0.01 * stopping_crit2) {
                return; // we are at a good place!
            }
            else {
                good_to_go = true;
            }
        }
        else {
            good_to_go = false;
        }

        faq = 1.0;
        denom = zero;
        if (dp != zero) {
            fac_newton = p / dp;
            fac_extra = d2p_half / dp;
            F_half = fac_newton * fac_extra;
            denom_sqrt = sqrt(c_one - two_n_div_n_1 * F_half);

            //NEXT LINE PROBABLY CAN BE COMMENTED OUT. Check if compiler outputs positive real
            if (real(denom_sqrt) >= 0.0) {
                denom = c_one_nth + n_1_nth * denom_sqrt;
            }
            else {
                denom = c_one_nth - n_1_nth * denom_sqrt;
            }
        }

        if (denom == 0) {
            dx = (absroot + 1.0) * expcmplx(complex(0.0, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * M_PI));
        }
        else {
            dx = fac_newton / denom;
        }


        newroot = *root - dx;
        if (newroot == *root) return; //nothing changes so return
        if (good_to_go) {
            *root = newroot;
            return;
        }
        if (i % FRAC_JUMP_EVERY == 0) { //decide whether to do a jump of modified length (to break cycles)
            faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
            newroot = *root - faq * dx; // do jump of semi-random length
        }
        *root = newroot;
    }
    success = false; // too many iterations here
    return;
}

void TripleLensing::cmplx_newton_spec(complex *poly, int degree, complex *root, int &iter, bool &success) {
    //Subroutine finds one root of a complex polynomial
    //Newton's method. It calculates simplified Adams' stopping
    //criterion for the value of the polynomial once per 10 iterations (!),
    //after initial iteration. This is done to speed up calculations
    //when polishing roots that are known preety well, and stopping
    // criterion does significantly change in their neighborhood.

    //Uses 'root' value as a starting point (!!!!!)
    //Remember to initialize 'root' to some initial guess.
    //Do not initilize 'root' to point (0,0) if the polynomial
    //coefficients are strictly real, because it will make going
    //to imaginary roots impossible.

    // poly - is an array of polynomial cooefs
    //  length = degree+1, poly(1) is constant
    //0                 1               2
    //poly[0] x^0 + poly[1] x^1 + poly[2] x^2 + ...
    //degree - a degree of the polynomial
    // root - input: guess for the value of a root
    //        output: a root of the polynomial
    //iter - number of iterations performed (the number of polynomial evaluations)
    //success - is false if routine reaches maximum number of iterations

    //For a summary of the method go to:
    //http://en.wikipedia.org/wiki/Newton's_method

    int FRAC_JUMP_EVERY = 10;
    const int FRAC_JUMP_LEN = 10;
    double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045, 0.13653241, 0.306162, 0.37794326, 0.04618805, 0.75132137 }; //some random numbers
    double faq; //jump length
    // double FRAC_ERR = 2e-15;
    complex p; //value of polynomial
    complex dp; //value of 1st derivative
    int i, k;
    bool good_to_go;
    complex dx, newroot;
    double ek, absroot, abs2p;
    complex zero = complex(0, 0);
    double stopping_crit2;

    iter = 0;
    success = true;

    // //the next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version
    // // if (false) { //change false to true if you would like to use caustion about haveing first coefficient == 0
    //  if (degree < 0) {
    //      printf("Error: cmplx_newton_spec: degree<0");
    //      return;
    //  }
    //  if (poly[degree] == zero) {
    //      if (degree == 0) return;
    //      cmplx_newton_spec(poly, degree, root, iter, success);
    //      return;
    //  }
    //  if (degree <= 1) {
    //      if (degree == 0) {
    //          success = false;
    //          printf("Warning: cmplx_newton_spec: degree=0 and poly[0]!=0, no roots");
    //          return;
    //      }
    //      else {
    //          *root = -poly[0] / poly[1];
    //          return;
    //      }
    //  }
    // }
    //end EXTREME Failsafe

    good_to_go = false;

    stopping_crit2 = 0.0; //value not important, will be initialized anyway on the first loop
    for (i = 1; i <= MAXIT; i++) {
        faq = 1.0;
        //prepare stoping criterion
        //calculate value of polynomial and its first two derivatives
        p = poly[degree];
        dp = zero;
        if (i % 10 == 1) { //calculate stopping criterion every tenth iteration
            ek = abs(poly[degree]);
            absroot = abs(*root);
            for (k = degree - 1; k >= 0; k--) {
                dp = p + dp * (*root);
                p = poly[k] + p * (*root); //b_k
                //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                //Communications of ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                //Eq. 8
                ek = absroot * ek + abs(p);
            }
            stopping_crit2 = pow(FRAC_ERR * ek, 2);
        }
        else { // calculate just the value and derivative
            for (k = degree - 1; k >= 0; k--) { //Horner Scheme, see for eg. Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                dp = p + dp * (*root);
                p = poly[k] + p * (*root);
            }
        }

        iter = iter + 1;

        abs2p = real(conj(p) * p);
        if (abs2p == 0.0) return;
        if (abs2p < stopping_crit2) { //simplified a little Eq. 10 of Adams 1967
            if (dp == zero) return; //if we have problem with zero, but we are close to the root, just accept
            //do additional iteration if we are less than 10x from stopping criterion
            if (abs2p < 0.01 * stopping_crit2) return; //return immediatley because we are at very good place
            else {
                good_to_go = true; //do one iteration more
            }
        }

        else {
            good_to_go = false; //reset if we are outside the zone of the root
        }
        if (dp == zero) {
            //problem with zero
            dx = (abs(*root) + 1.0) * expcmplx(complex(0.0, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * M_PI));
        }
        else {
            dx = p / dp; // Newton method, see http://en.wikipedia.org/wiki/Newton's_method
        }
        newroot = *root - dx;
        if (newroot == *root) return; //nothing changes -> return
        if (good_to_go) {//this was jump already after stopping criterion was met
            *root = newroot;
            return;
        }
        if (i % FRAC_JUMP_EVERY == 0) { // decide whether to do a jump of modified length (to break cycles)
            faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
            newroot = *root - faq * dx;
        }
        *root = newroot;
    }
    success = false;
    return;
    //too many iterations here
}

void TripleLensing::cmplx_laguerre2newton(complex *poly, int degree, complex *root, int &iter, bool &success, int starting_mode) {
    //Subroutine finds one root of a complex polynomial using
    //Laguerre's method, Second-order General method and Newton's
    //method - depending on the value of function F, which is a
    //combination of second derivative, first derivative and
    //value of polynomial [F=-(p"*p)/(p'p')].

    //Subroutine has 3 modes of operation. It starts with mode=2
    //which is the Laguerre's method, and continues until F
    //becames F<0.50, at which point, it switches to mode=1,
    //i.e., SG method (see paper). While in the first two
    //modes, routine calculates stopping criterion once per every
    //iteration. Switch to the last mode, Newton's method, (mode=0)
    //happens when becomes F<0.05. In this mode, routine calculates
    //stopping criterion only once, at the beginning, under an
    //assumption that we are already very close to the root.
    //If there are more than 10 iterations in Newton's mode,
    //it means that in fact we were far from the root, and
    //routine goes back to Laguerre's method (mode=2).

    //Uses 'root' value as a starting point (!!!!!)
    //Remember to initialize 'root' to some initial guess or to
    //point (0,0) if you have no prior knowledge.

    //poly - is an array of polynomial cooefs
    //  0                   1               2
    //  poly[0] x^0 + poly[1] x^1 + poly[2] x^2
    //degree - a degree of the polynomial
    //root - input: guess for the value of a root
    //      output: a root of the polynomial
    //iter - number of iterations performed (the number of polynomial
    //       evaluations and stopping criterion evaluation)
    //success - is false if routine reaches maximum number of iterations
    //starting_mode - this should be by default = 2. However if you
    //                choose to start with SG method put 1 instead.
    //                Zero will cause the routine to
    //                start with Newton for first 10 iterations, and
    //                then go back to mode 2.

    //For a summary of the method see the paper: Skowron & Gould (2012)

    int FRAC_JUMP_EVERY = 10;
    const int FRAC_JUMP_LEN = 10;
    double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045, 0.13653241, 0.306162, 0.37794326, 0.04618805, 0.75132137 }; //some random numbers

    double faq; //jump length
    // double FRAC_ERR = 2.0e-15;

    complex p; //value of polynomial
    complex dp; //value of 1st derivative
    complex d2p_half; //value of 2nd derivative
    int i, j, k;
    bool good_to_go;
    //complex G, H, G2;
    complex denom, denom_sqrt, dx, newroot;
    double ek, absroot, abs2p, abs2_F_half;
    complex fac_netwon, fac_extra, F_half, c_one_nth;
    double one_nth, n_1_nth, two_n_div_n_1;
    int mode;
    complex c_one = complex(1, 0);
    complex zero = complex(0, 0);
    double stopping_crit2;

    iter = 0;
    success = true;
    stopping_crit2 = 0; //value not important, will be initialized anyway on the first loop

    // //next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
    // if (false) {//change false to true if you would like to use caution about having first coefficent == 0
    //  if (degree < 0) {
    //      printf("Error: cmplx_laguerre2newton: degree < 0");
    //      return;
    //  }
    //  if (poly[degree] == zero) {
    //      if (degree == 0) return;
    //      cmplx_laguerre2newton(poly, degree, root, iter, success, starting_mode);
    //      return;
    //  }
    //  if (degree <= 1) {
    //      if (degree == 0) {//// we know from previous check that poly[0] not equal zero
    //          success = false;
    //          printf("Warning: cmplx_laguerre2newton: degree = 0 and poly[0] = 0, no roots");
    //          return;
    //      }
    //      else {
    //          *root = -poly[0] / poly[1];
    //          return;
    //      }
    //  }
    // }
    // //end EXTREME failsafe

    j = 1;
    good_to_go = false;

    mode = starting_mode; // mode = 2 full laguerre, mode = 1 SG, mode = 0 newton

    for (;;) { //infinite loop, just to be able to come back from newton, if more than 10 iteration there

        ////////////
        ///mode 2///
        ////////////

        if (mode >= 2) {//Laguerre's method
            one_nth = 1.0 / (degree); ///
            n_1_nth = (degree - 1) * one_nth; ////
            two_n_div_n_1 = 2.0 / n_1_nth;
            c_one_nth = complex(one_nth, 0.0);

            for (i = 1; i <= MAXIT; i++) {
                faq = 1.0;

                //prepare stoping criterion
                ek = abs(poly[degree]);
                absroot = abs(*root);
                //calculate value of polynomial and its first two derivative
                p = poly[degree];
                dp = zero;
                d2p_half = zero;
                for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                    d2p_half = dp + d2p_half * (*root);
                    dp = p + dp * (*root);
                    p = poly[k - 1] + p * (*root); // b_k
                    //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                    //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                    //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                    //Eq 8.
                    ek = absroot * ek + abs(p);
                }
                abs2p = real(conj(p) * p); // abs(p)
                iter = iter + 1;
                if (abs2p == 0) return;

                stopping_crit2 = pow(FRAC_ERR * ek, 2);
                if (abs2p < stopping_crit2) {//(simplified a little Eq. 10 of Adams 1967)
                    //do additional iteration if we are less than 10x from stopping criterion
                    if (abs2p < 0.01 * stopping_crit2) return; // ten times better than stopping criterion
                    //return immediately, because we are at very good place
                    else {
                        good_to_go = true; //do one iteration more
                    }
                }
                else {
                    good_to_go = false; //reset if we are outside the zone of the root
                }

                denom = zero;
                if (dp != zero) {
                    fac_netwon = p / dp;
                    fac_extra = d2p_half / dp;
                    F_half = fac_netwon * fac_extra;

                    abs2_F_half = real(conj(F_half) * F_half);
                    if (abs2_F_half <= 0.0625) {//F<0.50, F/2<0.25
                        //go to SG method
                        if (abs2_F_half <= 0.000625) {//F<0.05, F/2<0.02
                            mode = 0; //go to Newton's
                        }
                        else {
                            mode = 1; //go to SG
                        }
                    }

                    denom_sqrt = sqrt(c_one - two_n_div_n_1 * F_half);

                    //NEXT LINE PROBABLY CAN BE COMMENTED OUT
                    if (real(denom_sqrt) > 0.0) {
                        //real part of a square root is positive for probably all compilers. You can ù
                        //test this on your compiler and if so, you can omit this check
                        denom = c_one_nth + n_1_nth * denom_sqrt;
                    }
                    else {
                        denom = c_one_nth - n_1_nth * denom_sqrt;
                    }
                }
                if (denom == zero) {//test if demoninators are > 0.0 not to divide by zero
                    dx = (abs(*root) + 1.0) + expcmplx(complex(0.0, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * M_PI)); //make some random jump
                }
                else {
                    dx = fac_netwon / denom;
                }
                newroot = *root - dx;
                if (newroot == *root) return; // nothing changes -> return
                if (good_to_go) {//this was jump already after stopping criterion was met
                    *root = newroot;
                    return;
                }
                if (mode != 2) {
                    *root = newroot;
                    j = i + 1; //remember iteration index
                    break; //go to Newton's or SG
                }
                if ((i % FRAC_JUMP_EVERY) == 0) { //decide whether to do a jump of modified length (to break cycles)
                    faq = FRAC_JUMPS[((i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN)];
                    newroot = *root - faq * dx; // do jump of some semi-random length (0 < faq < 1)
                }
                *root = newroot;
            } //do mode 2

            if (i >= MAXIT) {
                success = false;
                return;
            }
        }

        ////////////
        ///mode 1///
        ////////////

        if (mode == 1) {//SECOND-ORDER GENERAL METHOD (SG)

            for (i = j; i <= MAXIT; i++) {
                faq = 1.0;
                //calculate value of polynomial and its first two derivatives
                p = poly[degree];
                dp = zero;
                d2p_half = zero;
                if ((i - j) % 10 == 0) {
                    //prepare stopping criterion
                    ek = abs(poly[degree]);
                    absroot = abs(*root);
                    for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                        d2p_half = dp + d2p_half * (*root);
                        dp = p + dp * (*root);
                        p = poly[k - 1] + p * (*root); //b_k
                        //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                        //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                        //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                        //Eq 8.
                        ek = absroot * ek + abs(p);
                    }
                    stopping_crit2 = pow(FRAC_ERR * ek, 2);
                }
                else {
                    for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                        d2p_half = dp + d2p_half * (*root);
                        dp = p + dp * (*root);
                        p = poly[k - 1] + p * (*root); //b_k
                    }
                }
                abs2p = real(conj(p) * p); //abs(p)**2
                iter = iter + 1;
                if (abs2p == 0.0) return;

                if (abs2p < stopping_crit2) {//(simplified a little Eq. 10 of Adams 1967)
                    if (dp == zero) return;
                    //do additional iteration if we are less than 10x from stopping criterion
                    if (abs2p < 0.01 * stopping_crit2) return; //ten times better than stopping criterion
                    //ten times better than stopping criterion
                    else {
                        good_to_go = true; //do one iteration more
                    }
                }
                else {
                    good_to_go = false; //reset if we are outside the zone of the root
                }
                if (dp == zero) {//test if denominators are > 0.0 not to divide by zero
                    dx = (abs(*root) + 1.0) * expcmplx(complex(0.0, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * M_PI)); //make some random jump
                }
                else {
                    fac_netwon = p / dp;
                    fac_extra = d2p_half / dp;
                    F_half = fac_netwon * fac_extra;

                    abs2_F_half = real(conj(F_half) * F_half);
                    if (abs2_F_half <= 0.000625) {//F<0.05, F/2<0.025
                        mode = 0; //set Newton's, go there after jump
                    }
                    dx = fac_netwon * (c_one + F_half); //SG
                }
                newroot = *root - dx;
                if (newroot == *root) return; //nothing changes -> return
                if (good_to_go) {
                    *root = newroot; //this was jump already after stopping criterion was met
                    return;
                }
                if (mode != 1) {
                    *root = newroot;
                    j = i + 1; //remember iteration number
                    break; //go to Newton's
                }
                if ((i % FRAC_JUMP_EVERY) == 0) { // decide whether to do a jump of modified length (to break cycles)
                    faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
                    newroot = *root - faq * dx; //do jump of some semi random lenth (0 < faq < 1)
                }
                *root = newroot;
            }
            if (i >= MAXIT) {
                success = false;
                return;
            }

        }
        //------------------------------------------------------------------------------- mode 0
        if (mode == 0) { // Newton's Method

            for (i = j; i <= j + 10; i++) { // Do only 10 iterations the most then go back to Laguerre
                faq = 1.0;

                //calc polynomial and first two derivatives
                p = poly[degree];
                dp = zero;
                if (i == j) { // Calculating stopping criterion only at the beginning
                    ek = abs(poly[degree]);
                    absroot = abs(*root);
                    for (k = degree; k >= 1; k--) {
                        dp = p + dp * (*root);
                        p = poly[k - 1] + p * (*root);
                        ek = absroot * ek + abs(p);
                    }
                    stopping_crit2 = pow(FRAC_ERR * ek, 2.0);
                }
                else {
                    for (k = degree; k >= 1; k--) {
                        dp = p + dp * (*root);
                        p = poly[k - 1] + p * (*root);
                    }
                }
                abs2p = real(conj(p) * p);
                iter = iter + 1;
                if (abs2p == 0.0) return;

                if (abs2p < stopping_crit2) {
                    if (dp == zero) return;
                    // do additional iteration if we are less than 10x from stopping criterion
                    if (abs2p < 0.01 * stopping_crit2) {
                        return; // return immediately since we are at a good place
                    }
                    else {
                        good_to_go = true; // do one more iteration
                    }
                }
                else {
                    good_to_go = false;
                }

                if (dp == zero) {
                    dx = (abs(*root) + 1.0) * expcmplx(complex(0.0, 2 * M_PI * FRAC_JUMPS[i % FRAC_JUMP_LEN])); // make a random jump
                }
                else {
                    dx = p / dp;
                }

                newroot = *root - dx;
                if (newroot == *root) return;
                if (good_to_go) {
                    *root = newroot;
                    return;
                }
                *root = newroot;
            }
            if (iter >= MAXIT) {
                //too many iterations
                success = false;
                return;
            }
            mode = 2; //go back to Laguerre's. Happens when could not converge with 10 steps of Newton
        }

    }/// end of infinite loop
}


//////////////////////////////
//////////////////////////////
////////complex methods and operators
//////////////////////////////
//////////////////////////////


complex::complex(double a, double b) {
    re = a;
    im = b;
}

complex::complex(double a) {
    re = a;
    im = 0;
}

complex::complex(void) {
    re = 0;
    im = 0;
}

double absf(double z) {
    if (z < 0) {
        return -z;
    }
    else {
        return z;
    }
    // return sqrt(z.re * z.re + z.im * z.im);
}

double absint(int z) {
    if (z < 0) {
        return -z;
    }
    else {
        return z;
    }
}

double abs(complex z) {
    return sqrt(z.re * z.re + z.im * z.im);
}

complex conj(complex z) {
    return complex(z.re, -z.im);
}

complex sqrt(complex z) {
    double md = sqrt(z.re * z.re + z.im * z.im);
    return (md > 0) ? complex((sqrt((md + z.re) / 2) * ((z.im > 0) ? 1 : -1)) , sqrt((md - z.re) / 2)) : 0.0;
}

double real(complex z) {
    return z.re;
}

double imag(complex z) {
    return z.im;
}

complex operator+(complex p1, complex p2) {
    return complex(p1.re + p2.re, p1.im + p2.im);
}

complex operator-(complex p1, complex p2) {
    return complex(p1.re - p2.re, p1.im - p2.im);
}

complex operator*(complex p1, complex p2) {
    return complex(p1.re * p2.re - p1.im * p2.im, p1.re * p2.im + p1.im * p2.re);
}

complex operator/(complex p1, complex p2) {
    double md = p2.re * p2.re + p2.im * p2.im;
    return complex((p1.re * p2.re + p1.im * p2.im) / md, (p1.im * p2.re - p1.re * p2.im) / md);
}

complex operator+(complex z, double a) {
    return complex(z.re + a, z.im);
}

complex operator-(complex z, double a) {
    return complex(z.re - a, z.im);
}

complex operator*(complex z, double a) {
    return complex(z.re * a, z.im * a);
}

complex operator/(complex z, double a) {
    return complex(z.re / a, z.im / a);
}

complex operator+(double a, complex z) {
    return complex(z.re + a, z.im);
}

complex operator-(double a, complex z) {
    return complex(a - z.re, -z.im);
}

complex operator*(double a, complex z) {
    return complex(a * z.re, a * z.im);
}

complex operator/(double a, complex z) {
    double md = z.re * z.re + z.im * z.im;
    return complex(a * z.re / md, -a * z.im / md);
}


complex operator+(complex z, int a) {
    return complex(z.re + a, z.im);
}

complex operator-(complex z, int a) {
    return complex(z.re - a, z.im);
}

complex operator*(complex z, int a) {
    return complex(z.re * a, z.im * a);
}

complex operator/(complex z, int a) {
    return complex(z.re / a, z.im / a);
}

complex operator+(int a, complex z) {
    return complex(z.re + a, z.im);
}

complex operator-(int a, complex z) {
    return complex(a - z.re, -z.im);
}

complex operator*(int a, complex z) {
    return complex(a * z.re, a * z.im);
}

complex operator/(int a, complex z) {
    double md = z.re * z.re + z.im * z.im;
    return complex(a * z.re / md, -a * z.im / md);
}

complex operator-(complex z) {
    return complex(-z.re, -z.im);
}

bool operator==(complex p1, complex p2) {
    if (p1.re == p2.re && p1.im == p2.im) return true;
    return false;
}

bool operator!=(complex p1, complex p2) {
    if (p1.re == p2.re && p1.im == p2.im) return false;
    return true;
}

complex expcmplx(complex p1) {
    double r = exp(p1.re);
    double theta = atan2(p1.im, p1.re);
    return complex(r * cos(theta), r * sin(theta));
}

complex cbrt(complex z) {
    complex zout;
    double r, r_cube, theta, theta_cube;
    r = abs(z);
    r_cube = pow(r, 0.333333333333);
    theta = atan2(z.im, z.re);
    theta_cube = theta / 3.;
    return  complex(r_cube * cos(theta_cube), r_cube * sin(theta_cube));
}