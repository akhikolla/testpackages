#include <Rcpp.h>
#include "htdp.h"
using namespace Rcpp;

// Reimplementation of interactive subroutine GETMDY in htdp.for
void c_getmdy(const int month, const int iday, const int iyear, double *date, int *mins)
{
  // Calculate modified Julian date
  int mjd = 0, mjd0 = 0, mjd1 = 0, i = 1;
  iymdmj_(&iyear, &month, &iday, &mjd);

  // Calculate time in decimal years (date)
  iymdmj_(&iyear, &i, &i, &mjd0);
  int iyear1 = iyear + 1;
  iymdmj_(&iyear1, &i, &i, &mjd1);
  int day = mjd - mjd0;
  int denom = mjd1 - mjd0;
  *date = iyear + (day / denom);

  // calculate Julian time in minutes (mins)
  *mins = mjd * 24 * 60;

  return;
}

// [[Rcpp::export(name = ".htdpinit")]]
void htdpinit()
{
  // Obtain parameters defining crustal motion model
  model_();

  // Initialize transformation parameters from ITRF94 to other reference frames
  settp_();

  // Initialize conversion table between reference frame identifiers
  setrf_();

  return;
}

// [[Rcpp::export]]
DataFrame displace(NumericMatrix xy, Date t0, Date t1, int iopt)
{
  int nrows = xy.nrow();

  if (iopt < 1 || iopt > 23 || iopt == 4) {
    stop("Invalid reference frame");
  }

  // Extract MDY for GETMDY subroutine
  int d0 = t0.getDay();
  int m0 = t0.getMonth();
  int y0 = t0.getYear();
  int d1 = t1.getDay();
  int m1 = t1.getMonth();
  int y1 = t1.getYear();

  if (y0 <= 1906 || y1 <= 1906) {
    stop("The model is not valid for dates prior to 1906");
  }

  // Calculate decimal years and Julian time in minutes
  int min0 = 0, min1 = 0;
  double date0 = 0, date1 = 0;
  c_getmdy(m0, d0, y0, &date0, &min0);
  c_getmdy(m1, d1, y1, &date1, &min1);

  // Ensure lat/lon is in valid range
  for (int i = 0; i < nrows; i++) {
    if (xy(i,1) < -90 || xy(i,1) > 90) {
      stop("Invalid latitude");
    }
    if (xy(i,0) < -180 || xy(i,0) > 180) {
      stop("Invalid longitude");
    }
  }

  // Main HTDP routines
  double lat0 = 0, lon0 = 0, eht0 = 0;
  double lat1 = 0, lon1 = 0, eht1 = 0;
  int jregn = 0;
  double vn = 0, ve = 0, vu = 0;
  double dn = 0, de = 0, du = 0;

  // Displacement vectors
  std::vector<double> dx;
  std::vector<double> dy;
  std::vector<double> dz;

  // Velocity vectors
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> vz;

  for (int i = 0; i < nrows; i++) {
    // lat/lon in positive N/W, radians
    lat0 = xy(i,1) * (M_PI / 180.0);
    lon0 = -xy(i,0) * (M_PI / 180.0);

    // Predict velocity in iopt reference frame
    predv_(&lat0, &lon0, &eht0, &date0, &iopt, &jregn, &vn, &ve, &vu);

    vx.push_back(ve);
    vy.push_back(vn);
    vz.push_back(vu);

    // Predict coordinates and displacements from time MIN1 to time MIN2
    newcor_(&lat0, &lon0, &eht0, &min0, &min1, &lat1, &lon1, &eht1, &dn, &de, &du, &vn, &ve, &vu);

    dx.push_back(de);
    dy.push_back(dn);
    dz.push_back(du);
  }

  Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::Named("de")=dx,
                                               Rcpp::Named("dn")=dy,
                                               Rcpp::Named("du")=dz,
                                               Rcpp::Named("ve")=vx,
                                               Rcpp::Named("vn")=vy,
                                               Rcpp::Named("vu")=vz);

  return df;
}
