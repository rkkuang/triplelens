/*
Testing the usage of codes for finite source star lensed by triple lens system.
including:
  basic usage like single magnification computation, output critical curves and caustics, and image tracks
  light curve generating
  magnification map generating
  adaptive sampling of light curve

run ./bin/testtriple after make

For applying to Binary lens, set NLENS to 2 in src/TripleLensingLibrary.h

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
// #define MAPCOMPARE // flag controls whether compare with VBBL using binary system - magnification map

// #define LKVCOMPARE // flag controls whether compare with VBBL using binary system - light curve
// #define LMBCOMPARE // flag controls whether compare with VBBL on limb darkening - light curve


using namespace std;

VBBinaryLensing vbbl;

int main()
{
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

  // s2 = 1;
  // q2 = 1;

  // parameter conversion, from s2, q2,, s3, q3, psi to mlens[i], zlens[i]
  double inv1andq2 = 1 / (1 + q2);
  double inv1andq2q3 = 1 / (1 + q2 + q3);
  mlens[0] = inv1andq2q3;
  mlens[1] = q2 * inv1andq2q3;
  zlens[0] = complex(- q2 * inv1andq2 * s2 , 0);
  zlens[1] = complex( inv1andq2 * s2 , 0);
  if (NLENS == 3) {
    mlens[2] = q3 * inv1andq2q3;
    zlens[2] = complex(zlens[0].re + s3 * cos(psi), s3 * sin(psi));

    //*******************************************************************************************
    //*******************************************************************************************
    //*******************************************************************************************
    //*******************************************************************************************
    //*******************************************************************************************
    // part 1, basic usage, //after runing this cpp test code, run python3 test/plotcaus.py to visualize

    // double xsCenter = -0.022, ysCenter = 0.031; // source position
    // double xsCenter = -0.039744764784384, ysCenter = -0.027607534315906; // source position
    xsCenter = -0.034747426672208;
    ysCenter = -0.026627816352184; // source position

    // for better visulization of image tracks, we set a larger source radius:
    rs = 0.005;

    // test by directly assign mlens and zlens
    // mlens[0] = 0.779062;
    // mlens[1] = 0.198325;
    // mlens[2] = 0.0226121;
    // zlens[0] = complex(-0.187325, -0.0208104);
    // zlens[1] = complex(0.802407, -0.0208104);
    // zlens[2] = complex(-0.583751, 0.89951);


    xsCenter = 0.25;
    ysCenter = 0;
    rs = 0.3;
    TripleLensing TRIL(mlens, zlens);
    TRIL.secnum = 45;
    TRIL.quaderr_Tol = 1e-2;
    TRIL.basenum = 1;


    // simply call TripleMag function to get the magnification of the source at (xsCenter, ysCenter)
    muTri = TRIL.TripleMag(xsCenter, ysCenter, rs);
    fprintf(stderr, "muTri: %f\n", muTri);

    //// output triple lensing event parameters
    outsys(mlens, zlens, t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs, xsCenter, ysCenter);
    // outputCriticalTriple(mlens, zlens,NLENS, 10000);
  }

  // //end if NLENS == 3
  TripleLensing TRIL(mlens, zlens);
  TRIL.secnum = 45;
  TRIL.quaderr_Tol = 1e-2;
  TRIL.basenum = 1;
  TRIL.reset(mlens, zlens);


  // output critical curves and caustics, and image positions corresponding to the source boundary
  outputCriticalTriple(mlens, zlens, NLENS, 10000);
  TRIL.outImgPoints(xsCenter, ysCenter, rs, 10001);
  // outputCaustics_for_crossSection(mlens, zlens, NLENS,  3000);


  double salpha = sin(alpha), calpha = cos(alpha), tn, tE_inv = 1 / tE;

  // // test at certain time on a light curve
  double test_t = 7495.58, mu;
  rs = 0.005;
  test_t = 7496.776;//7496.488;//7496.044 54.714842 & muFS = 54.734133; //7496.868;
  tn = (test_t - t0) * tE_inv;
  xsCenter = -u0 * salpha + tn * calpha;
  ysCenter = u0 * calpha - tn * salpha;
  mu = TRIL.TripleMag(xsCenter, ysCenter, rs);
  fprintf(stderr, "\t\t muFS = %f, xsCenter = %.15f, ysCenter = %.15f\n", mu, xsCenter, ysCenter);

  // to test whether tripleFS2python works, it works
  // it should work on python
  double mu2;
  mu2 = TRIL.tripleFS2python(xsCenter, ysCenter, rs);
  double Zlens[2 * NLENS] = { zlens[0].re, zlens[0].im, zlens[1].re, zlens[1].im};
  double mu3[2], xs[2], ys[2];
  xs[0] = xsCenter;
  xs[1] = 2*xsCenter;
  ys[0] = ysCenter;
  ys[1] = 2*ysCenter;

  TRIL.tripleFS2py(mlens, Zlens, xs, ys, rs, 45, 1, 1e-2,1e-3, mu3, 2);
  fprintf(stderr, "\t\t muFS2 = %f, err = %le \n", mu2, (mu2-mu)/mu);
  fprintf(stderr, "\t\t muFS3 = %f, err = %le \n", mu3[0], (mu3[0]-mu)/mu);
  fprintf(stderr, "\t\t muFS3 = %f \n", mu3[1]);




//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
// part 2, light curve generating, //after runing this cpp test code, run python3 test/plotlkv.py to visualize
// full lightcurve
#ifdef TESTLKV

  rs = 2.2e-4; // we set back the source radius to generate the lightcurve
  rs = 0.005;
  rs = 0.001;

  double range = 1;
  double Zlens[2 * NLENS] = { zlens[0].re, zlens[0].im, zlens[1].re, zlens[1].im, zlens[2].re, zlens[2].im };

// Now let us calculate the light curve on np points equally spaced between t0-3tE and t0+3tE:
  int np = 10000; // ~ 10m cadence
  double *t_array = NULL, *mag_array = NULL, *y1_array = NULL, *y2_array = NULL;
  t_array = new double[np];
  mag_array = new double[np];
  y1_array = new double[np];
  y2_array = new double[np];
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
    mag_array[i] = TRIL.TripleMag(y1_array[i], y2_array[i], rs);
    // mag_array[i] = TRIL.TriplePS(mlens, zlens, y1_array[i], y2_array[i]);
  }
  std::cout << std::endl;
  tim1 = clock();
  fprintf(stderr, "total time on lightkv generating = %f s\n", (tim1 - tim0) / 1e6);

  FILE *ftrilkv; //save light curve to .dat file
  ftrilkv = fopen("./data/ftrilkv_rs1e-3_adap_v2.dat", "w");

  for (int j = 0; j < np; j++) {
    fprintf(ftrilkv, "%.15f %.15f %.15f %.15f ", t_array[j], y1_array[j], y2_array[j], mag_array[j]);
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
      muFS = TRIL.TripleMag(px, py, rs);
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
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
// part 4, adaptive light curve sampling,
// after runing this cpp test code, run python3 test/plotlkvadap.py to visualize



#ifdef TESTADA
  _curve *lightkv = new _curve;
  rs = 2.2e-4;
  rs = 0.005;
  rs = 0.001;
  double t_start = 7470, t_end = 7510;
  clock_t tim0, tim1;


  // double Zlens[2 * NLENS] = { zlens[0].re, zlens[0].im, zlens[1].re, zlens[1].im, zlens[2].re, zlens[2].im };
  fprintf(stderr, "sampling ... \n");
  TRIL.reset(mlens, zlens);
  tim0 = clock();
  TRIL.TripleLkvAdap(rs, lightkv, t_start, t_end, alpha, tE, t0, u0);
  tim1 = clock();
  fprintf(stderr, "total time on lightkv adaptive sampling = %f s\n", (tim1 - tim0) / 1e6);
  //rs 2.2e-4, adaptive sampling = 32.138063 s,  pntcnt = 477, total time on lightkv generating = 12.054293 s
  //rs 1e-2, total time on lightkv adaptive sampling =50.124823  s, pntcnt = 226 ,total time on lightkv generating = 254.249186 s
  //rs 5e-2, 1215.239477 s,
  // save data
  FILE *ftrilkv_adaptive;
  ftrilkv_adaptive = fopen("./data/ftrilkv_adaptive_rs1e-3_v2.dat", "w");
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
// comparing limb darkening results from VBBL and Triple lensing
#ifdef LMBCOMPARE
  FILE *ftrilkv, *fbinlkv, *ftimratio; //save light curve to .dat file
  vbbl.Tol = 1e-3;
  double s, q;
  s = 1;
  q = 0.5;
  alpha = 0.95;
  tE = 100.3;
  t0 = 7154;

  rs = 1e-3;
  string strrs = "1e-3";

  double limba1 = 0.51; // limb darkening coefficient

  string stru0;
  string fnametri = "./data/limblkvcmp_tri_rs";
  string fnameb = "./data/limblkvcmp_bin_rs";
  string fnametim = "./data/limblkvcmp_tim_rs";

  int which = 1;
  switch (which) {
  case 1:
    u0 = -0.24;// -0.24, 0.48, 0.8
    stru0 = "-0.24";
    break;

  case 2:
    u0 = 0.48;// -0.24, 0.48, 0.8
    stru0 = "0.48";
    break;

  case 3:
    u0 = 0.8;// -0.24, 0.48, 0.8
    stru0 = "0.8";
    break;
  }


  string final_name1 = fnametri + strrs + "u0" + stru0 + ".dat";
  string final_name2 = fnameb + strrs + "u0" + stru0 + ".dat";
  string final_name3 = fnametim + strrs + "u0" + stru0 + ".dat";

  ftrilkv = fopen(final_name1.c_str(), "w");
  fbinlkv = fopen(final_name2.c_str(), "w");
  ftimratio = fopen(final_name3.c_str(), "w");


  mlens[0] = 1.0 / (1.0 + q);
  mlens[1] = q / (1.0 + q);//small
  zlens[0] = complex(-q * s / (1 + q), 0);
  zlens[1] = complex(s / (1 + q), 0);
  salpha = sin(alpha);
  calpha = cos(alpha);

  TRIL.reset(mlens, zlens);
  TRIL.quaderr_Tol = 1e-5;
  TRIL.secnum = 15;
  TRIL.relerr_mag = 5e-4;

  int np = 5000; // ~ 10m cadence
  double *t_array = NULL, *mag_array = NULL, *y1_array = NULL, *y2_array = NULL, *mag_array_VBBL = NULL;
  int *ifFinate_array_tri = NULL, *ifFinate_array_vbb = NULL;
  t_array = new double[np];
  mag_array = new double[np];
  mag_array_VBBL = new double[np];
  y1_array = new double[np];
  y2_array = new double[np];
  ifFinate_array_tri = new int[np];
  ifFinate_array_vbb = new int[np];
  double t_start = 7000, t_end = 7300;
  double dt = (t_end - t_start) / np;
  double dtvbb, dttri;
  clock_t tim0, tim1;

  int Nexp = 10, Navg = 5, i3 = 0, i2 = 0;
  double dlogrho = 2.0 / (Nexp - 1);

  if (which == 3) {
    outsys(mlens, zlens, t0, u0, tE, s, q, alpha, s3, q3, psi, rs, xsCenter, ysCenter);
  }

  for (int i = 0; i < np; i++) {
    t_array[i] = t_start + i * dt;
    tn = (t_array[i] - t0) * tE_inv;
    // Position of the source center
    y1_array[i] = u0 * salpha + tn * calpha;
    y2_array[i] = u0 * calpha - tn * salpha;
  }

  tim0 = clock();
  for (int i = 0; i < np; i++) {
    // mag_array_VBBL[i] = vbbl.BinaryMag2(s, q, y1_array[i], y2_array[i], rs);
    mag_array_VBBL[i] = vbbl.BinaryMagDark(s, q, y1_array[i], y2_array[i], rs, limba1, vbbl.Tol);

    // ifFinate_array_vbb[i] = vbbl.ifFinite;
    ifFinate_array_vbb[i] = 1;
  }
  tim1 = clock();
  dtvbb = (tim1 - tim0) / 1.0e6;

  fprintf(stderr, "i3 = %d, i2 = %d, total time on map generating by vbbl = %f s\n", i3, i2, (tim1 - tim0) / 1e6);


  tim0 = clock();
  // for (int i = 0; i < np; i++) {
  for (int i = 3291; i < 3292; i++) {
    mag_array[i] = TRIL.TripleMagDark(y1_array[i], y2_array[i], rs, limba1);

    fprintf(stderr, "i = %d, trimag = %lf, binmag = %lf, relerr = %le\n\n", i, mag_array[i], mag_array_VBBL[i], (mag_array[i]-mag_array_VBBL[i])/mag_array_VBBL[i]);

    // ifFinate_array_tri[i] = TRIL.ifFinite;
    ifFinate_array_tri[i] = 1;
  }
  tim1 = clock();
  dttri = (tim1 - tim0) / 1.0e6;

  fprintf(stderr, "i3 = %d, i2 = %d, total time on map generating by TRIL = %f s\n", i3, i2, (tim1 - tim0) / 1e6);



  fprintf(ftimratio, "%.15f %.15f %.15f \n", dttri, dtvbb, dttri / dtvbb);

  for (int j = 0; j < np; j++) {
    fprintf(ftrilkv, "%.15f %.15f %.15f %.15f %d ", t_array[j], y1_array[j], y2_array[j], mag_array[j], ifFinate_array_tri[j]);
    fprintf(fbinlkv, "%.15f %.15f %.15f %.15f %d ", t_array[j], y1_array[j], y2_array[j], mag_array_VBBL[j], ifFinate_array_vbb[j]);
  }
  fclose(ftrilkv);
  fclose(fbinlkv);
  fclose(ftimratio);

  outputCriticalTriple(mlens, zlens, NLENS, 10000);
  // TRIL.outImgPoints(xsCenter, ysCenter, rs, 10001);
#endif




//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
#ifdef LKVCOMPARE
  FILE *ftrilkv, *fbinlkv, *ftimratio; //save light curve to .dat file
  vbbl.Tol = 1e-3;
  double s, q;
  s = 1;
  q = 0.5;
  alpha = 0.95;
  tE = 100.3;
  t0 = 7154;

  rs = 1e-3;
  string strrs = "1e-3";

  string stru0;
  string fnametri = "./data/flkvcmp_tri_rs";
  string fnameb = "./data/flkvcmp_bin_rs";
  string fnametim = "./data/flkvcmp_tim_rs";

  int which = 1;
  switch (which) {
  case 1:
    u0 = -0.24;// -0.24, 0.48, 0.8
    stru0 = "-0.24";
    break;

  case 2:
    u0 = 0.48;// -0.24, 0.48, 0.8
    stru0 = "0.48";
    break;

  case 3:
    u0 = 0.8;// -0.24, 0.48, 0.8
    stru0 = "0.8";
    break;
  }


  string final_name1 = fnametri + strrs + "u0" + stru0 + ".dat";
  string final_name2 = fnameb + strrs + "u0" + stru0 + ".dat";
  string final_name3 = fnametim + strrs + "u0" + stru0 + ".dat";

  ftrilkv = fopen(final_name1.c_str(), "w");
  fbinlkv = fopen(final_name2.c_str(), "w");
  ftimratio = fopen(final_name3.c_str(), "w");


  mlens[0] = 1.0 / (1.0 + q);
  mlens[1] = q / (1.0 + q);//small
  zlens[0] = complex(-q * s / (1 + q), 0);
  zlens[1] = complex(s / (1 + q), 0);
  salpha = sin(alpha);
  calpha = cos(alpha);

  TRIL.reset(mlens, zlens);
  TRIL.quaderr_Tol = 1e-5;
  TRIL.secnum = 10;
  TRIL.relerr_mag = 5e-4;

  // y1s[i] = pr[2] * salpha - tn * calpha;
  // y2s[i] = -pr[2] * calpha - tn * salpha;
  // mags[i] = BinaryMag2(s, q, y1s[i], y2s[i], rho);

  int np = 5000; // ~ 10m cadence
  double *t_array = NULL, *mag_array = NULL, *y1_array = NULL, *y2_array = NULL, *mag_array_VBBL = NULL;
  int *ifFinate_array_tri = NULL, *ifFinate_array_vbb = NULL;
  t_array = new double[np];
  mag_array = new double[np];
  mag_array_VBBL = new double[np];
  y1_array = new double[np];
  y2_array = new double[np];
  ifFinate_array_tri = new int[np];
  ifFinate_array_vbb = new int[np];
  double t_start = 7000, t_end = 7300;
  double dt = (t_end - t_start) / np;
  double dtvbb, dttri;
  clock_t tim0, tim1;

  int Nexp = 10, Navg = 5, i3 = 0, i2 = 0;
  double dlogrho = 2.0 / (Nexp - 1);

  if (which == 3) {
    outsys(mlens, zlens, t0, u0, tE, s, q, alpha, s3, q3, psi, rs, xsCenter, ysCenter);
  }

  for (int i = 0; i < np; i++) {
    t_array[i] = t_start + i * dt;
    tn = (t_array[i] - t0) * tE_inv;
    // Position of the source center
    y1_array[i] = u0 * salpha + tn * calpha;
    y2_array[i] = u0 * calpha - tn * salpha;
  }

  tim0 = clock();
  for (int i = 0; i < np; i++) {
    mag_array_VBBL[i] = vbbl.BinaryMag2(s, q, y1_array[i], y2_array[i], rs);
    ifFinate_array_vbb[i] = vbbl.ifFinite;
  }
  tim1 = clock();
  dtvbb = (tim1 - tim0) / 1.0e6;

  fprintf(stderr, "i3 = %d, i2 = %d, total time on map generating by vbbl = %f s\n", i3, i2, (tim1 - tim0) / 1e6);

  tim0 = clock();
  for (int i = 0; i < np; i++) {
  // for (int i = 3291; i < 3292; i++) {
    mag_array[i] = TRIL.TripleMag(y1_array[i], y2_array[i], rs);
    ifFinate_array_tri[i] = TRIL.ifFinite;

    // fprintf(stderr, "i = %d, trimag = %lf, binmag = %lf, relerr = %le\n\n", i, mag_array[i], mag_array_VBBL[i], (mag_array[i]-mag_array_VBBL[i])/mag_array_VBBL[i]);

  }
  tim1 = clock();
  dttri = (tim1 - tim0) / 1.0e6;

  fprintf(stderr, "i3 = %d, i2 = %d, total time on map generating by TRIL = %f s\n", i3, i2, (tim1 - tim0) / 1e6);



  fprintf(ftimratio, "%.15f %.15f %.15f \n", dttri, dtvbb, dttri / dtvbb);

  for (int j = 0; j < np; j++) {
    fprintf(ftrilkv, "%.15f %.15f %.15f %.15f %d ", t_array[j], y1_array[j], y2_array[j], mag_array[j], ifFinate_array_tri[j]);
    fprintf(fbinlkv, "%.15f %.15f %.15f %.15f %d ", t_array[j], y1_array[j], y2_array[j], mag_array_VBBL[j], ifFinate_array_vbb[j]);
  }
  fclose(ftrilkv);
  fclose(fbinlkv);
  fclose(ftimratio);

  outputCriticalTriple(mlens, zlens, NLENS, 10000);
  // TRIL.outImgPoints(xsCenter, ysCenter, rs, 10001);
#endif



//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************


// #undef MAPCOMPARE
//   fprintf(stderr, "\n\n\n\n\n\n\n\n\n\n\n");
//   vbbl.Tol = 1e-3;

//   double s, q;
//   s = 0.8;
//   q = 0.1;
//   mlens[0] = 1.0 / (1.0 + q);
//   mlens[1] = q / (1.0 + q);//small
//   zlens[0] = complex(-q * s / (1 + q), 0);
//   zlens[1] = complex(s / (1 + q), 0);
//   rs = 1e-2;

//   TRIL.reset(mlens, zlens);
//   TRIL.quaderr_Tol = 1e-3;
//   TRIL.secnum = 10;
//   TRIL.relerr_mag = 1e-3;
//   double muFSvbbl, muFSTRIL;

//   xsCenter = 0.460629917681217; ysCenter = -0.5;
//   xsCenter = 0; ysCenter = 0;
//   xsCenter = 0.0275590531528; ysCenter = 0.003937005996704;
//   xsCenter = 0.059055116027594; ysCenter = 0.011811021715403;
//   xsCenter = -0.405511811375618; ysCenter = -0.437007874250412;
//   xsCenter = -0.035433072596788; ysCenter = 0.003937005996704;
//   xsCenter = -0.043307088315487; ysCenter = 0.003937005996704;
//   // xsCenter = -0.043307088315487; ysCenter = 0.011811021715403;
//   xsCenter = -0.460629921406507; ysCenter = 0.492125980556011;
//   xsCenter = -0.161417324095964; ysCenter = 0.263779524713755;
//   xsCenter = -0.145669292658567; ysCenter = 0.248031493276358;
//   xsCenter = -0.02755905687809; ysCenter = -0.011811025440693;
//   xsCenter = -0.5; ysCenter = -0.492125984281301;
//   xsCenter = -0.043307088315487; ysCenter = 0.003937005996704;
//   xsCenter = -0.043307088315487; ysCenter = -0.003937009721994;
//   // xsCenter = -0.5; ysCenter = -0.484251968562603;
//   // xsCenter = -0.035433072596788; ysCenter = -0.003937009721994;
//   xsCenter = -0.035433072596788; ysCenter = 0.003937005996704;

//   // muFSvbbl = vbbl.BinaryMag2(s, q, xsCenter, ysCenter, rs);
//   muFSvbbl = vbbl.BinaryMag(s, q, xsCenter, ysCenter, rs, vbbl.Tol);
//   // muFSTRIL = TRIL.TripleMag(xsCenter, ysCenter, rs);
//   muFSTRIL = TRIL.TripleMag_vbbl(xsCenter, ysCenter, rs, vbbl.Tol);
//   fprintf(stderr, "muFSvbbl %f (%d points), muFSTRIL %.15f (%d points), relerr = %.3e\n", muFSvbbl, vbbl.finalNPS , muFSTRIL, TRIL.finalNPS ,(muFSTRIL - muFSvbbl) / muFSvbbl);
//   // output critical curves and caustics, and image positions corresponding to the source boundary
//   outputCriticalTriple(mlens, zlens, NLENS);
//   TRIL.outImgPoints(xsCenter, ysCenter, rs, 3001);
//   outsys(mlens, zlens, t0, u0, tE, s2, q2, alpha, s3, q3, psi, rs, xsCenter, ysCenter);

#ifdef MAPCOMPARE
  // compare speed with Bozza's code using binary lens system
  // VBBinaryLensing vbbl;
  vbbl.Tol = 1e-3;
  double s, q;
  s = 0.8;
  q = 0.1;
  mlens[0] = 1.0 / (1.0 + q);
  mlens[1] = q / (1.0 + q);//small
  zlens[0] = complex(-q * s / (1 + q), 0);
  zlens[1] = complex(s / (1 + q), 0);
  rs = 1e-2;

  TRIL.reset(mlens, zlens);
  TRIL.quaderr_Tol = 1e-3;
  TRIL.secnum = 5;
  TRIL.relerr_mag = 1e-3;

// #define MULTIEXP

  FILE *fmumap, *fmumapvbbl, *ftimecomp;

#ifdef MULTIEXP
  ftimecomp = fopen("./data/magmapTimeRatio6_10_2.dat", "w");
#else
  fmumap = fopen("./data/magmapTRIL_s0.8q0.1rs1e-2_new5.dat", "w");
  fmumapvbbl = fopen("./data/magmapvbbl_s0.8q0.1rs1e-2_new5.dat", "w");
#endif

  //set the magnification map boundary
  float xlim0, xlim1, ylim0, ylim1;
  // xlim0 = -0.5; xlim1 = 1; ylim0 = -0.75; ylim1 = 0.75;
  xlim0 = -0.5; xlim1 = 0.5; ylim0 = -0.5; ylim1 = 0.5;
  int ImgSize = 128;
  int ImgSize2 = ImgSize * ImgSize;
  double dx = (xlim1 - xlim0) / (ImgSize - 1);
  double dy = (ylim1 - ylim0) / (ImgSize - 1);
  double px = xlim0, py = ylim0;
  double *mus = NULL, *musvbbl = NULL , *y1_array = NULL, *y2_array = NULL, muFS;
  mus = new double[ImgSize2];
  musvbbl = new double[ImgSize2];
  y1_array = new double[ImgSize2];
  y2_array = new double[ImgSize2];

  double dtvbb, dttri;
  clock_t tim0, tim1;
  std::cout.width(3);//i的输出为3位宽

  int Nexp = 10, Navg = 5, i3 = 0, i2 = 0;
  double dlogrho = 2.0 / (Nexp - 1); // -4 to -2

#ifdef MULTIEXP
  fprintf(ftimecomp, "%d %d \n", Nexp, Navg);
  for (i3 = 0; i3 < Nexp; i3++) {
    // TRIL.secnum = (int)((i3 * dlogrho)*30 + 90);
    rs = pow( 10, -4 + i3 * dlogrho );
    for (i2 = 0; i2 < Navg; i2++) { // ten time average
#endif

      tim0 = clock();
      px = xlim0;
      for (int i = 0; i < ImgSize; i++) {
        py = ylim0;
        for (int j = 0; j < ImgSize; j++) {
          muFS = vbbl.BinaryMag2(s, q, px, py, rs);
          // muFS = vbbl.BinaryMag(s, q, px, py, rs, vbbl.Tol);
          musvbbl[i * ImgSize + j] = muFS;
          y1_array[i * ImgSize + j] = px;
          y2_array[i * ImgSize + j] = py;
          py += dy;
        }
        px += dx;
      }
      tim1 = clock();
      dtvbb = (tim1 - tim0) / 1.0e6;
      fprintf(stderr, "i3 = %d, i2 = %d, total time on map generating by vbbl = %f s\n", i3, i2, (tim1 - tim0) / 1e6);

      tim0 = clock();
      px = xlim0;
      for (int i = 0; i < ImgSize; i++) {
        py = ylim0;
        for (int j = 0; j < ImgSize; j++) {
          // fprintf(stderr, "xsCenter = %.15f; ysCenter = %.15f;\n", px, py);
          muFS = TRIL.TripleMag(px, py, rs);
          mus[i * ImgSize + j] = muFS;
          y1_array[i * ImgSize + j] = px;
          y2_array[i * ImgSize + j] = py;
          py += dy;
        }
        px += dx;
      }
      tim1 = clock();
      dttri = (tim1 - tim0) / 1.0e6;

      fprintf(stderr, "i3 = %d, i2 = %d, total time on map generating by TRIL = %f s\n", i3, i2, (tim1 - tim0) / 1e6);

#ifdef MULTIEXP
      fprintf(ftimecomp, "%.15f %.15f %.15f ", rs, dtvbb, dttri);
    }
  }
  fclose(ftimecomp);
#else
      for (int i = 0; i < ImgSize; i++) {
        for (int j = 0; j < ImgSize; j++) {
          fprintf(fmumap, "%.15f %.15f %.15f ", y1_array[i * ImgSize + j], y2_array[i * ImgSize + j], mus[i * ImgSize + j]);
          fprintf(fmumapvbbl, "%.15f %.15f %.15f ", y1_array[i * ImgSize + j], y2_array[i * ImgSize + j], musvbbl[i * ImgSize + j]);
        }
      }
      fclose(fmumap);
      fclose(fmumapvbbl);
#endif
#endif





  return (0);
}