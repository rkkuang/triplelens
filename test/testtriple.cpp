/*
Testing the usage of codes for finite source star lensed by triple lens system.
including:
  basic usage like single magnification computation, output critical curves and caustics, and image tracks
  light curve generating
  magnification map generating
  adaptive sampling of light curve

run ./bin/testtriple after make

For applying to Binary lens, set NLENS to 2 in src/TripleLensingLibrary.h and then redo compile

The output file can further be processed with provided python scripts (inside the test folder) for visualization

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


// main tests
// #define TESTLKV // flag controls whether test light curve generating
// #define TESTADA // flag controls whether test adaptive sampling of light curve
// #define TESTMAP // flag controls whether test magnification map generating
// #define LMBLKV // flag controls test limb darkening light curve calculation


using namespace std;
VBBinaryLensing vbbl;

int main()
{
  int NLENS = 4; // how many number of point lenses
  // Declaration of an instance of TripleLensing class.
  fprintf(stderr, "NLENS = %d\n", NLENS);
  double mlens[NLENS]; // NLENS defined in .h file
  complex zlens[NLENS];
  double u0, t0, tE, alpha, s2, q2, s3, q3, psi, rs;
  double muTri, xsCenter, ysCenter;

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
  zlens[0] = complex(- q2 * inv1andq2 * s2 , 0);
  zlens[1] = complex( inv1andq2 * s2 , 0);

  TripleLensing TRIL;
  TRIL.setnlens(NLENS);

  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  // part 1, basic usage, //after runing this cpp test code, run python3 test/plotcaus.py to visualize
  // for better visulization of image tracks, we set a larger source radius:
  rs = 0.01;
  xsCenter = 0; // source position
  ysCenter = 0;

  if (NLENS == 3) {
    mlens[2] = q3 * inv1andq2q3;
    zlens[2] = complex(zlens[0].re + s3 * cos(psi), s3 * sin(psi));
    // // TripleLensing TRIL(mlens, zlens); // add the third paramter nlens on 2022.11.01, for the flexibility to change the number of lenses
    // TRIL.reset3(mlens, zlens);
    // TRIL.secnum = 90; // divide the source bondary into how many parts
    // TRIL.quaderr_Tol = 1e-3; // the Quadrupole test tolerance
    // TRIL.basenum = 2; // the number density of sampled dots among each part
  } else if (NLENS == 4) { //end if NLENS == 3
    // TripleLensing TRIL(mlens, zlens);
    mlens[2] = q3 * inv1andq2q3 * 0.5;
    zlens[2] = complex(zlens[0].re + s3 * cos(psi), s3 * sin(psi));
    mlens[3] = mlens[2];
    zlens[3] = complex(-zlens[2].re, -zlens[2].im);
  }

  TRIL.reset3(mlens, zlens);
  TRIL.secnum = 90;
  TRIL.quaderr_Tol = 1e-3;
  TRIL.basenum = 2;
  // simply call TripleMag function to get the magnification of the source at (xsCenter, ysCenter)
  muTri = TRIL.TripleMag(xsCenter, ysCenter, rs);
  fprintf(stderr, "muTri: %f\n", muTri);
  //// output triple lensing event parameters
  TRIL.outsys(mlens, zlens, t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs, xsCenter, ysCenter);
  // outputCriticalTriple(mlens, zlens,NLENS, 10000);

  // output critical curves and caustics, and image positions corresponding to the source boundary
  outputCriticalTriple(mlens, zlens, NLENS, 10000);
  TRIL.outImgPoints(xsCenter, ysCenter, rs, 10001);
  // outputCaustics_for_crossSection(mlens, zlens, NLENS,  3000);


  double salpha = sin(alpha), calpha = cos(alpha), tn, tE_inv = 1 / tE;

  // // test at certain time on a light curve
  double test_t = 7490, mu;
  rs = 0.005;
  test_t = 7497.100830078125;//7496.488;//7496.044 54.714842 & muFS = 54.734133; //7496.868;
  tn = (test_t - t0) * tE_inv;
  xsCenter = u0 * salpha + tn * calpha;
  ysCenter = u0 * calpha - tn * salpha;
  mu = TRIL.TripleMag(xsCenter, ysCenter, rs);
  fprintf(stderr, "\t\t muFS = %f, xsCenter = %.15f, ysCenter = %.15f\n", mu, xsCenter, ysCenter);

//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
// part 2, light curve generating, //after runing this cpp test code, run python3 test/plotlkv.py to visualize
// full lightcurve
#ifdef TESTLKV
  rs = 0.005;
// Now let us calculate the light curve on np points equally spaced between t0-3tE and t0+3tE:
  int np = 1000;
  double *t_array = NULL, *mag_array = NULL, *y1_array1 = NULL, *y2_array1 = NULL;
  t_array = new double[np];
  mag_array = new double[np];
  y1_array1 = new double[np];
  y2_array1 = new double[np];
  double t_start = 7470, t_end = 7510;
  double dt = (t_end - t_start) / np;
  clock_t tim0, tim1;
  for (int i = 0; i < np; i++) {
    t_array[i] = t_start + i * dt;
  }
  fprintf(stderr, "generating light curve (%d points) ... \n", np);
  tim0 = clock();
  std::cout.width(3);
  for (int i = 0; i < np; i++) {
    std::cout << 100.0 * i / np << "%";
    cout << flush << '\r' << string(i * 100.0 / np, '#');

    tn = (t_array[i] - t0) * tE_inv;
    // Position of the source center
    y1_array1[i] = u0 * salpha + tn * calpha;
    y2_array1[i] = u0 * calpha - tn * salpha;
    mag_array[i] = TRIL.TripleMag(y1_array1[i], y2_array1[i], rs);
    // fprintf(stderr, "\t\t tn = %f, muFS = %f \n", t_array[i], mag_array[i] );
    // mag_array[i] = TRIL.TriplePS(mlens, zlens, y1_array[i], y2_array[i]);
  }
  std::cout << std::endl;
  tim1 = clock();
  fprintf(stderr, "total time on lightkv generating = %f s\n", (tim1 - tim0) / 1e6);

  FILE *ftrilkv; //save light curve to .dat file
  ftrilkv = fopen("./data/fulltrilkv_rs5e-3.dat", "w");

  for (int j = 0; j < np; j++) {
    fprintf(ftrilkv, "%.15f %.15f %.15f %.15f ", t_array[j], y1_array1[j], y2_array1[j], mag_array[j]);
  }
  fclose(ftrilkv);
#endif


//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
// part 3, magnification map generating, //after runing this cpp test code, run python3 test/plotmap.py to visualize
#ifdef TESTMAP
  rs = 5e-2; // we set back the source radius to generate the magnification map
  //set the magnification map boundary
  float xlim0, xlim1, ylim0, ylim1;
  xlim0 = -0.2; xlim1 = 1; ylim0 = -0.6; ylim1 = 0.6;
  int ImgSize = 128;
  int ImgSize2 = ImgSize * ImgSize;
  double dx = (xlim1 - xlim0) / (ImgSize - 1);
  double dy = (ylim1 - ylim0) / (ImgSize - 1);
  double px = xlim0, py = ylim0;
  double *mus = NULL, *y1_array2 = NULL, *y2_array2 = NULL, muFS;
  mus = new double[ImgSize2];
  y1_array2 = new double[ImgSize2];
  y2_array2 = new double[ImgSize2];
  // double mus = new [ImgSize2], y1_array[ImgSize2], y2_array[ImgSize2], muFS;
  clock_t tim02, tim12;
  std::cout.width(3);//i的输出为3位宽
  fprintf(stderr, "generating magnification map (%d points) ... \n", ImgSize2);
  tim02 = clock();
  for (int i = 0; i < ImgSize; i++) {
    py = ylim0;
    for (int j = 0; j < ImgSize; j++) {
      std::cout << 100.0 * (i * ImgSize + j) / ImgSize2 << "%";
      cout << flush << '\r' << string((i * ImgSize + j) * 100.0 / ImgSize2, '#');
      muFS = TRIL.TripleMag(px, py, rs);
      mus[i * ImgSize + j] = muFS;
      y1_array2[i * ImgSize + j] = px;
      y2_array2[i * ImgSize + j] = py;
      py += dy;
    }
    px += dx;
  }
  std::cout << std::endl;
  tim12 = clock();
  fprintf(stderr, "total time on map generating = %f s\n", (tim12 - tim02) / 1e6);

  FILE *fmumap;
  fmumap = fopen("./data/magmap0.05.dat", "w");
  for (int i = 0; i < ImgSize; i++) {
    for (int j = 0; j < ImgSize; j++) {
      fprintf(fmumap, "%.15f %.15f %.15f ", y1_array2[i * ImgSize + j], y2_array2[i * ImgSize + j], mus[i * ImgSize + j]);
    }
  }
  fclose(fmumap);
#endif

//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
// part 4, adaptive light curve sampling,
// after runing this cpp test code, run python3 test/plotlkvadap.py to visualize

#ifdef TESTADA
  _curve *lightkv = new _curve;
  rs = 0.005;
  double t_start_adap = 7470, t_end_adap = 7510;
  clock_t tim0_adp, tim1_adap;

  fprintf(stderr, "sampling ... \n");
  // TRIL.reset(mlens, zlens);
  tim0_adp = clock();
  TRIL.TripleLkvAdap(rs, lightkv, t_start_adap, t_end_adap, alpha, tE, t0, u0);
  tim1_adap = clock();
  fprintf(stderr, "total time on lightkv adaptive sampling = %f s\n", (tim1_adap - tim0_adp) / 1e6);
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


//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
// test limb darkening light curve
#ifdef LMBLKV
  rs = 0.1; // use a larger source radius to see the limb-darkening effect
  double limba1 = 0.51; // limb darkening coefficient
  int np3 = 100; //
  double *t_array3 = NULL, *mag_array3 = NULL, *mag_array_limb = NULL, *y1_array3 = NULL, *y2_array3 = NULL;
  t_array3 = new double[np3];
  mag_array3 = new double[np3];
  mag_array_limb = new double[np3];
  y1_array3 = new double[np3];
  y2_array3 = new double[np3];
  double t_start3 = 7470, t_end3 = 7510;
  double dt3 = (t_end3 - t_start3) / np3;
  clock_t tim03, tim13;
  for (int i = 0; i < np3; i++) {
    t_array3[i] = t_start3 + i * dt3;
  }

  TRIL.secnum = 90; // divide the source bondary into how many parts
  TRIL.quaderr_Tol = 1e-3; // the Quadrupole test tolerance
  TRIL.basenum = 2;
  TRIL.relerr_mag = 1e-3;

  // parameters controls the accuracy of limb-darkening calculation
  TRIL.RelTolLimb = 1e-2;
  TRIL.AbsTolLimb = 1e-2;

  fprintf(stderr, "generating uniform brightness light curve (%d points) ... \n", np3);
  tim03 = clock();
  std::cout.width(3);
  for (int i = 0; i < np3; i++) {
    std::cout << 100.0 * i / np3 << "%";
    cout << flush << '\r' << string(i * 100.0 / np3, '#');

    tn = (t_array3[i] - t0) * tE_inv;
    // Position of the source center
    y1_array3[i] = u0 * salpha + tn * calpha;
    y2_array3[i] = u0 * calpha - tn * salpha;

    mag_array3[i] = TRIL.TripleMag(y1_array3[i], y2_array3[i], rs);

  }
  std::cout << std::endl;
  tim13 = clock();
  fprintf(stderr, "total time on lightkv generating = %f s\n", (tim13 - tim03) / 1e6);

  fprintf(stderr, "generating limb darkening light curve (%d points) ... \n", np3);
  tim03 = clock();
  std::cout.width(3);
  for (int i = 0; i < np3; i++) {
    std::cout << 100.0 * i / np3 << "%";
    cout << flush << '\r' << string(i * 100.0 / np3, '#');
    mag_array_limb[i] = TRIL.TripleMagDark(y1_array3[i], y2_array3[i], rs, limba1); // not exactly, due to y1_array, y2_array are already obtained
  }
  std::cout << std::endl;
  tim13 = clock();
  fprintf(stderr, "total time on limb darkening lightkv generating = %f s\n", (tim13 - tim03) / 1e6);

  FILE *ftrilkv3, *ftrilkv_limb; //save light curve to .dat file
  ftrilkv3 = fopen("./data/unilkv.dat", "w");
  ftrilkv_limb = fopen("./data/limblkv.dat", "w");

  for (int j = 0; j < np3; j++) {
    fprintf(ftrilkv3, "%.15f %.15f %.15f %.15f ", t_array3[j], y1_array3[j], y2_array3[j], mag_array3[j]);
    fprintf(ftrilkv_limb, "%.15f %.15f %.15f %.15f ", t_array3[j], y1_array3[j], y2_array3[j], mag_array_limb[j]);
  }
  fclose(ftrilkv3);
  fclose(ftrilkv_limb);
#endif



  return (0);
}
