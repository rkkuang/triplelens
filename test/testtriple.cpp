/*
Testing the usage of codes for finite source star lensed by triple lens system.
including:
  basic usage like single magnification computation, output critical curves and caustics, and image tracks
  light curve generating
  magnification map generating
  adaptive sampling of light curve

run ./bin/testtriple after make

The output file can further be processed with provided python scripts (inside test folder) for visualization

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "../src/VBBinaryLensingLibrary.h"
#include "../src/TripleLensingLibrary.h"
#include <time.h>
#include <fstream>
#include <string.h>

// #define TESTLKV // flag controls whether test light curve generating
// #define TESTMAP // flag controls whether test magnification map generating
// #define TESTADA // flag controls whether test adaptive sampling of light curve
using namespace std;

int main()
{
  // Declaration of an instance of TripleLensing class.
  TripleLensing TRIL;
  // TRIL.secnum = 90;
  TRIL.basenum = 2;
  double mlens[NLENS]; // NLENS defined in .h file
  complex zlens[NLENS];
  double u0, t0, tE, alpha, s2, q2, s3, q3, psi, rs;

  // we adopt the triple lens geometry of Han et al., 2017 https://doi.org/10.3847/1538-3881/aa9179
  // Table 3, sol C wide
  // event parameters
  t0 = 7494.153;
  u0 = 0.021;
  tE = 74.62;
  s2 = 1.396;
  q2 = 0.029;
  alpha = 2.948; //rad
  s3 = 1.168;
  q3 = 3.27e-3;
  psi = 5.332; //rad
  rs = 0.22e-3;


  // parameter conversion, from s2, q2,, s3, q3, psi to mlens[i], zlens[i]
  double inv1andq2 = 1 / (1 + q2);
  double inv1andq2q3 = 1 / (1 + q2 + q3);
  mlens[0] = inv1andq2q3;
  mlens[1] = q2 * inv1andq2q3;
  mlens[2] = q3 * inv1andq2q3;
  zlens[0] = complex(- q2 * inv1andq2 * s2 , 0);
  zlens[1] = complex( inv1andq2 * s2 , 0);
  zlens[2] = complex(zlens[0].re + s3 * cos(psi), s3 * sin(psi));

  //*******************************************************************************************
  // part 1, basic usage, //after runing this cpp test code, run python3 test/plotcaus.py to visualize
  double muTri;
  // double xsCenter = -0.022, ysCenter = 0.031; // source position
  // double xsCenter = -0.039744764784384, ysCenter = -0.027607534315906; // source position
  double xsCenter = -0.034747426672208, ysCenter = -0.026627816352184; // source position

  // for better visulization of image tracks, we set a larger source radius:
  rs = 0.005;

  // simply call TripleMag function to get the magnification of the source at (xsCenter, ysCenter)
  muTri = TRIL.TripleMag(mlens, zlens, xsCenter, ysCenter, rs);
  fprintf(stderr, "muTri: %f\n", muTri);

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


  // output critical curves and caustics, and image positions corresponding to the source boundary
  outputCriticalTriple(mlens, zlens, NLENS);
  outImgPoints(mlens, zlens, xsCenter, ysCenter, rs, 3001);


  //*******************************************************************************************
  // part 2, light curve generating, //after runing this cpp test code, run python3 test/plotlkv.py to visualize
// full lightcurve
  double salpha = sin(alpha), calpha = cos(alpha), tn, tE_inv = 1 / tE;
#ifdef TESTLKV

  rs = 2.2e-4; // we set back the source radius to generate the lightcurve
  rs = 0.005;

  double range = 1;
  double Zlens[2 * NLENS] = { zlens[0].re, zlens[0].im, zlens[1].re, zlens[1].im, zlens[2].re, zlens[2].im };

// Now let us calculate the light curve on np points equally spaced between t0-3tE and t0+3tE:
  int np = 10000; // ~ 10m cadence
  double *t_array = NULL, *mag_array = NULL, *y1_array = NULL, *y2_array = NULL;
  t_array = new double[np];
  mag_array = new double[np];
  y1_array = new double[np];
  y2_array = new double[np];
  // double t_array[np], mag_array[np], y1_array[np], y2_array[np];
  double t_start = 7470, t_end = 7510;
  double dt = (t_end - t_start) / np;
  clock_t tim0, tim1;
  for (int i = 0; i < np; i++) {
    t_array[i] = t_start + i * dt;
  }
  tim0 = clock();
  std::cout.width(3);//i的输出为3位宽
  for (int i = 0; i < np; i++) {
    std::cout << 100.0 * i / np << "%";
    cout << flush << '\r' << string(i * 100.0 / np, '#');
    // std::cout << "\b\b\b\b";
    tn = (t_array[i] - t0) * tE_inv;
    // Position of the source center
    y1_array[i] = u0 * salpha + tn * calpha;
    y2_array[i] = u0 * calpha - tn * salpha;
    mag_array[i] = TRIL.TripleMag(mlens, zlens, y1_array[i], y2_array[i], rs);
    // mag_array[i] = TRIL.TriplePS(mlens, zlens, y1_array[i], y2_array[i]);
  }
  std::cout << std::endl;
  tim1 = clock();
  fprintf(stderr, "total time on lightkv generating = %f s\n", (tim1 - tim0) / 1e6);

  FILE *ftrilkv; //save light curve to .dat file
  ftrilkv = fopen("./data/ftrilkv_rs5e-3_adap.dat", "w");

  for (int j = 0; j < np; j++) {
    fprintf(ftrilkv, "%.15f %.15f %.15f %.15f ", t_array[j], y1_array[j], y2_array[j], mag_array[j]);
  }
  fclose(ftrilkv);
#endif


  // // test at certain time
  double test_t = 7495.58, mu;
  test_t = 7496.488;//7496.044 54.714842 & muFS = 54.734133; //7496.868;
  tn = (test_t - t0) * tE_inv;
  xsCenter = -u0 * salpha + tn * calpha;
  ysCenter = u0 * calpha - tn * salpha;
  // mu = tripleFS_v2_savehalf_gould(mlens, zlens, xsCenter, ysCenter, rs, nphi, &finalnphi, secnum, basenum, &quad_err);
  mu = TRIL.TripleMag(mlens, zlens, xsCenter, ysCenter, rs);
  fprintf(stderr, "\t\t muFS = %f, xsCenter = %.15f, ysCenter = %.15f\n", mu, xsCenter, ysCenter);

  //*******************************************************************************************
  // part 3, magnification map generating, //after runing this cpp test code, run python3 test/plotmap.py to visualize
#ifdef TESTMAP
  rs = 2.2e-4; // we set back the source radius to generate the lightcurve
  rs = 5e-2;
  //set the magnification map boundary
  float xlim0, xlim1, ylim0, ylim1;
  // xlim0 = -0.5; xlim1 = 1; ylim0 = -0.75; ylim1 = 0.75;
  xlim0 = -0.2; xlim1 = 1; ylim0 = -0.6; ylim1 = 0.6;
  int ImgSize = 128;
  int ImgSize2 = ImgSize * ImgSize;
  double dx = (xlim1 - xlim0) / (ImgSize - 1);
  double dy = (ylim1 - ylim0) / (ImgSize - 1);
  double px = xlim0, py = ylim0;
  double *mus = NULL, *y1_array = NULL, *y2_array = NULL, muFS;
  mus = new double[ImgSize2];
  y1_array = new double[ImgSize2];
  y2_array = new double[ImgSize2];
  // double mus = new [ImgSize2], y1_array[ImgSize2], y2_array[ImgSize2], muFS;
  clock_t tim0, tim1;
  std::cout.width(3);//i的输出为3位宽
  tim0 = clock();
  for (int i = 0; i < ImgSize; i++) {
    py = ylim0;
    for (int j = 0; j < ImgSize; j++) {
      std::cout << 100.0 * (i * ImgSize + j) / ImgSize2 << "%";
      cout << flush << '\r' << string((i * ImgSize + j) * 100.0 / ImgSize2, '#');
      muFS = TRIL.TripleMag(mlens, zlens, px, py, rs);
      mus[i * ImgSize + j] = muFS;
      y1_array[i * ImgSize + j] = px;
      y2_array[i * ImgSize + j] = py;
      py += dy;
    }
    px += dx;
  }
  std::cout << std::endl;
  tim1 = clock();
  fprintf(stderr, "total time on map generating = %f s\n", (tim1 - tim0) / 1e6);

  FILE *fmumap;
  fmumap = fopen("./data/magmap.dat", "w");
  for (int i = 0; i < ImgSize; i++) {
    for (int j = 0; j < ImgSize; j++) {
      fprintf(fmumap, "%.15f %.15f %.15f ", y1_array[i * ImgSize + j], y2_array[i * ImgSize + j], mus[i * ImgSize + j]);
    }
  }
  fclose(fmumap);
#endif

//*******************************************************************************************
// part 4, adaptive light curve sampling,
// after runing this cpp test code, run python3 test/plotlkvadap.py to visualize

#ifdef TESTADA
  _curve *lightkv = new _curve;
  rs = 2.2e-4;
  rs = 0.005;
  double t_start = 7470, t_end = 7510;
  clock_t tim0, tim1;
  int np0 = 2, MAXNPS = 2000;
  double quad_err, quaderr_Tol = 1e-4, errTol = 0.001, timerrTol = 1e-5;
  double Zlens[2 * NLENS] = { zlens[0].re, zlens[0].im, zlens[1].re, zlens[1].im, zlens[2].re, zlens[2].im };
  fprintf(stderr, "sampling ... \n");
  tim0 = clock();
  TRIL.TripleLkvAdap(mlens, Zlens, rs, lightkv, np0, MAXNPS, errTol, t_start, t_end, alpha, tE, t0, u0, 1, timerrTol);
  tim1 = clock();
  fprintf(stderr, "total time on lightkv adaptive sampling = %f s\n", (tim1 - tim0) / 1e6);
  //rs 2.2e-4, adaptive sampling = 32.138063 s,  pntcnt = 477, total time on lightkv generating = 12.054293 s
  //rs 1e-2, total time on lightkv adaptive sampling =50.124823  s, pntcnt = 226 ,total time on lightkv generating = 254.249186 s
  //rs 5e-2, 1215.239477 s,
  // save data
  FILE *ftrilkv_adaptive;
  ftrilkv_adaptive = fopen("./data/ftrilkv_adaptive_rs5e-3.dat", "w");
  fprintf(stderr, "open ftrilkv_adaptive.dat, pntcnt = %d\n", lightkv->length);
  for (_point *scan = lightkv->first; scan; scan = scan->next) {
    fprintf(ftrilkv_adaptive, "%.15f ", scan->x1);
  }
  fprintf(ftrilkv_adaptive, "\n");
  for (_point *scan = lightkv->first; scan; scan = scan->next) {
    fprintf(ftrilkv_adaptive, "%.15f ", scan->x2);
  }
  fclose(ftrilkv_adaptive);
#endif

  return (0);
}