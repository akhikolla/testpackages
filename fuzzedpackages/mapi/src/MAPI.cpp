#include <iostream>
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <math.h>

/*
NOTE: "using namespace" are commented out and replaced by the syntax namespace::object in order to avoid ambiguities between std and Rcpp.
WARNING: DO NOT UNCOMMENT!
using namespace std;
using namespace Rcpp;
*/


// //' Function countMatches_cpp
// //' 
// //' From the list of integer vectors which reprresents the intersection of grid cells and ellipses, this function counts the total number of matches (items) in this list of vectors.
// //'
// //' @param inter A list of integer vector containing intersections between cells and ellipses
// //' @return the total number of items (sum of lengths of all vectors in the list).
// [[Rcpp::export]]
long countMatches_cpp(Rcpp::List inter) {
  int n = inter.size();
  long resu = 0L;
  //Rcout << "n=" << n << "\n";
  for (int i=0; i<n; i++) {
    Rcpp::List ells = inter(i);
    long ne = ells.size();
	resu += ne;
  }
  return(resu);
}

// [[Rcpp::export]]
Rcpp::List getValues_cpp(Rcpp::NumericVector cells, Rcpp::List inter, Rcpp::DoubleVector weights, Rcpp::DoubleVector values) {
  int n = inter.size();
  Rcpp::List resu(0);
  for (int i=0; i<n; i++) {
    Rcpp::List ells = inter(i);
    int ne = ells.size();
    if (ne == 0) {
        Rcpp::NumericMatrix elem(0,0);
        resu.push_back(elem);
    } else {
      Rcpp::NumericMatrix elem(ne,2);
      for (int j=0; j<ne; j++) {
        int ie = int(ells[j]) - 1;
        if (ie < weights.size()) {
          double w = weights(ie);
          double v = values(ie);
	  elem(j, 0) = v;
	  elem(j, 1) = w;
        } else {
          Rcpp::Rcout << "overflow: ie="<<ie<<"\n";
          break;
        }
      }
      resu.push_back(elem);
    }
  }
  return(resu);
}

// //' Function parseInter_cpp
// //' 
// //' From the list of integer vectors which represents the intersection of grid cells and ellipses, this function returns a numeric matrix with one row per grid cell and five columns : cell gid, the number of intersecting ellipses, the weighted-averaged value of intesecting ellipses values, the sum of intersecting ellipses weights and the weighted standard deviation of intesecting ellipses values.
// //'
// //' @param cells An integer vector containing cells ids
// //' @param inter A list of integer vector containing intersections between cells and ellipses
// //' @param weights A double-precision vector containing weights of the ellipses
// //' @param values A double-precision vector containing values of the ellipses
// //' @return a numeric matrix with one row per grid cell and five columns : cell gid, the number of intersecting ellipses, the weighted-averaged value of intesecting ellipses values, the sum of intersecting ellipses weights and the weighted standard deviation of intesecting ellipses values.
// [[Rcpp::export]]
Rcpp::NumericMatrix parseInter_cpp(Rcpp::NumericVector cells, Rcpp::List inter, Rcpp::DoubleVector weights, Rcpp::DoubleVector values) {
  int n = inter.size();
  Rcpp::NumericMatrix resu(n,5);
  //Rcout << "n=" << n << "\n";
  for (int i=0; i<n; i++) {
    int gid = cells(i);
    Rcpp::List ells = inter(i);
    int ne = ells.size();
    //Rcpp::Rcout << "i="<< i << "\t gid=" << gid << "\t ne=" << ne << "\n";
    if (ne == 0) {
      resu(i,0) = gid;
      resu(i,1) = NA_REAL;
      resu(i,2) = NA_REAL;
      resu(i,3) = NA_REAL;
      resu(i,4) = NA_REAL;
      //cout << i << "<-NA\n";
    } else {
      // Weighted mean
      double valuesSum = 0.0;
      double weightsSum = 0.0;
      double squareSum = 0.0;
      for (int j=0; j<ne; j++) {
        int ie = int(ells[j]) - 1;
        if (ie < weights.size()) {
          double w = weights(ie);
          double v = values(ie);
          if (! (std::isnan(w) || std::isnan(v))) {
            valuesSum  += w * v;
            squareSum  += w * pow(v, 2);
            weightsSum += w;
          }
        } else {
          Rcpp::Rcout << "overflow: ie="<<ie<<"\n";
          break;
        }
      }
      double avg = valuesSum / weightsSum;
      double var = (squareSum / weightsSum) - pow(avg,2) ;
      double stdv = sqrt(var);
      //Rcpp::Rcout << gid<<"\t"<<avg<<"\t"<<weightsSum<<"\t"<<ne<<"\n";
      resu(i,0) = gid;
      resu(i,1) = ne;
      resu(i,2) = avg;
      resu(i,3) = weightsSum;
      resu(i,4) = stdv;
    }
  }
  //Rcpp::Rcout << "\n Fonction optimisée avec Konig-Huygens : stdv = " << resu(1,4)<<"\n";
  return(resu);
}



// //' Function parseInterPerm_cpp
// //' 
// //' From the list of integer vectors which reprresents the intersection of grid cells and ellipses, this function returns a numeric vector containing the weighted-averaged value of intesecting ellipses values.
// //'
// //' @param cells An integer vector containing cells ids
// //' @param inter A list of integer vector containing intersections between cells and ellipses
// //' @param weights A numeric vector containing weights of the ellipses
// //' @param values A numeric vector containing values of the ellipses
// //' @return a numeric vector with the weighted-averaged value of intesecting ellipses values.
// [[Rcpp::export]]
Rcpp::DoubleVector parseInterPerm_cpp(Rcpp::NumericVector cells, Rcpp::List inter, Rcpp::DoubleVector weights, Rcpp::DoubleVector values) {
  int n = inter.size();
  Rcpp::DoubleVector resu(n);
  //Rcout << "n=" << n << "\n";
  for (int i=0; i<n; i++) {
    Rcpp::List ells = inter(i);
    int ne = ells.size();
    //Rcout << "i="<< i << "\t gid=" << gid << "\t ne=" << ne << "\n";
    if (ne == 0) {
      resu(i) = NA_REAL;
      //cout << i << "<-NA\n";
    } else {
      // Weighted mean
      double valuesSum = 0.0;
      double weightsSum = 0.0;
      for (int j=0; j<ne; j++) {
        int ie = int(ells[j]) - 1;
        if (ie < weights.size()) {
          double w = weights(ie);
          double v = values(ie);
          if (! (std::isnan(w) || std::isnan(v))) {
            valuesSum  += w * v;
            weightsSum += w;
          }
        } else {
          Rcpp::Rcout << "overflow: ie="<<ie<<"\n";
          break;
        }
      }
      resu(i) = valuesSum / weightsSum;
    }
  }
  return(resu);
}

// Function for sorting Rcomplex as if CPoint in original source
bool RcomplexSorter (Rcomplex i, Rcomplex j) { 
	//return (i.r < j.r || i.r == j.r && i.i < j.i);
	return ( (i.r < j.r) || ((i.r == j.r) && (i.i < j.i)) );
}
// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
double c_cross(const Rcomplex &O, const Rcomplex &A, const Rcomplex &B) { return (A.r - O.r) * (B.i - O.i) - (A.i - O.i) * (B.r - O.r); }

// //' @title Function convex_hull
// //' 
// //' @description Computes convex hull of cloud of points.
// //'   Implementation of Andrew's monotone chain 2D convex hull algorithm, 
// //'     modified by S. Piry for Rcomplex points, march 2018.
// //'   Asymptotic complexity: O(n log n).
// //' 
// //'   Note: the last point in the returned list is not the same as the first one (unclosed polygon).
// //' 
// //' @references https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#C++
// //' 
// //' @param P complex vector, the point coordiantes
// //'
// //' @return a complex vector as coordinates of summits of the convex hull in counter-clockwise order.
// //'
std::vector<Rcomplex> convex_hull(Rcpp::ComplexVector P) {
  int n = P.size(), k = 0;
  if (n == 1) {
    std::vector<Rcomplex> H(1);
    H[0] = P[0];
    return H;
  } else {
    std::vector<Rcomplex> H(2*n);
    // Sort points lexicographically
    std::sort(P.begin(), P.end(), RcomplexSorter);
    // Build lower hull
    for (int i = 0; i < n; ++i) {
      while (k >= 2 && c_cross(H[k-2], H[k-1], P[i]) <= 0) k--;
      H[k++] = P[i];
    }
    // Build upper hull
    for (int i = n-2, t = k+1; i >= 0; i--) {
      while (k >= t && c_cross(H[k-2], H[k-1], P[i]) <= 0) k--;
      H[k++] = P[i];
    }
    H.resize(k-1);
    return H;
  }
}

// //' @title Function mkCc_cpp
// //' 
// //' @description Builds points for an circle.
// //'
// //' @param e double, the radius length
// //' @param c0 complex, coordinates of the center
// //' @param fic double vector, the N+1 angles for a quarter circle
// //'
// //' @return a complex vector as coordinates of the points
Rcpp::ComplexVector mkCc_cpp(double e, Rcomplex c0, Rcpp::DoubleVector fic) {
  int n = fic.size();
  Rcpp::ComplexVector eco(n);
  Rcpp::ComplexVector r((n-1) * 4 + 1);
  Rcomplex p, cr;
  int ip = 0;
  for (int i=0; i < n; i++) {
    p.r = e*cos(fic[i]);
    p.i = e*sin(fic[i]);
    eco[i] = p;
    cr = c0 + p;
    r[ip++] = cr;
  }
  for (int i=n-2; i > 0; i--) {
    p.r = - eco[i].r;
    p.i =   eco[i].i;
    cr = c0 + p;
    r[ip++] = cr;
  }
  for (int i=0; i < n; i++) {
    p.r = - eco[i].r;
    p.i = - eco[i].i;
    cr = c0 + p;
    r[ip++] = cr;
  }
  for (int i=n-2; i > 0; i--) {
    p.r =  eco[i].r;
    p.i = - eco[i].i;
    cr = c0 + p;
    r[ip++] = cr;
  }
  r[ip++] = r[0];
  return(r);
}

// //' @title Function mkEc_cpp
// //' 
// //' @description Builds points for an ellipse.
// //'
// //' @param a double, the half-minor axis length
// //' @param b double, the half-major axis length
// //' @param c0 complex, coordinates of the center
// //' @param rot complex, the rotation
// //' @param fic double vector, the N+1 angles for a quarter circle
// //'
// //' @return a complex vector as coordinates of the points
Rcpp::ComplexVector mkEc_cpp(double a, double b, Rcomplex c0, Rcomplex rot, Rcpp::DoubleVector fic) {
  int n = fic.size();
  Rcpp::ComplexVector eco(n);
  Rcpp::ComplexVector r((n-1) * 4 + 1);
  Rcomplex p, cr;
  int ip = 0;
  for (int i=0; i < n; ++i) {
    p.r = a*cos(fic[i]);
    p.i = b*sin(fic[i]);
    eco[i] = p;
    cr = c0 + p * rot;
    r[ip++] = cr;
  }
  for (int i=n-2; i > 0; i--) {
    p.r = - eco[i].r;
    p.i =   eco[i].i;
    cr = c0 + p * rot;
    r[ip++] = cr;
  }
  for (int i=0; i < n; ++i) {
    p.r = - eco[i].r;
    p.i = - eco[i].i;
    cr = c0 + p * rot;
    r[ip++] = cr;
  }
  for (int i=n-2; i > 0; i--) {
    p.r =   eco[i].r;
    p.i = - eco[i].i;
    cr = c0 + p * rot;
    r[ip++] = cr;
  }
  r[ip++] = r[0];
  return(r);
}

// //' Function mkP4st_cpp
// //' 
// //' Builds the convex hull of the ellipse joining two spatial points and their error circles
// //'
// //' @param r An numeric vector containing x1, y1, x2, y2, e1, e2 respectively the coordinates x,y and the error circle radius of point 1 (2).
// //' @param N The number of segments per quarter-circle
// //' @param ecc The eccentricity of the ellipse
// //' @return a numeric matrix of coordinates (x,y) in colums with one row per point of the convex hull
// [[Rcpp::export]]
Rcpp::NumericMatrix mkP4st_cpp(Rcpp::DoubleVector r, Rcpp::IntegerVector N, Rcpp::DoubleVector ecc) {
  int ic = 0;
  // 	int id = r[ic++];
  double x1 = r[ic++];
  double y1 = r[ic++];
  double x2 = r[ic++];
  double y2 = r[ic++];
  double e1 = r[ic++];
  double e2 = r[ic++];
  std::complex<double> d = std::complex<double>(x1-x2, y1-y2);
  Rcomplex c0; c0.r = (x1+x2)/2.0; c0.i = (y1+y2)/2.0;
  double f = abs(d) / 2.0;
  int N0 = N[0];
  double Nf = static_cast<double>(N0);
  Rcpp::DoubleVector fic(N0+1);
  for (int i=0; i<=N0; i++) {
    fic[i] = (M_PI / 2.0L) * static_cast<double>(i) / Nf;
  }
  if (f > 0.0) { // different points
    double rot = arg(d);
    std::complex<double> crot = exp(std::complex<double>(0,1) * rot); // missing sugar...
    Rcomplex rcrot; rcrot.r=crot.real(); rcrot.i=crot.imag();
    double a = f / ecc[0];
    double b = sqrt(pow(a,2) - pow(f,2));
    Rcpp::ComplexVector ellipse = mkEc_cpp(a, b, c0, rcrot, fic);
    Rcomplex c1; c1.r=x1; c1.i=y1;
    Rcpp::ComplexVector ec1 = mkCc_cpp(e1, c1, fic);
    Rcomplex c2; c2.r=x2; c2.i=y2;
    Rcpp::ComplexVector ec2 = mkCc_cpp(e2, c2, fic);
    // Now let's build the convex hull of the three set of points
    int npE = ellipse.size();
    int npC1 = ec1.size();
    int npC2 = ec2.size();
    int ipt = 0;
    Rcpp::ComplexVector pts(npE+npC1+npC2);
    for(int i=0; i<npE;  i++) { pts[ipt++]=ellipse[i]; }
    for(int i=0; i<npC1; i++) { pts[ipt++]=ec1[i]; }
    for(int i=0; i<npC2; i++) { pts[ipt++]=ec2[i]; }
    std::vector<Rcomplex> ch1 = convex_hull(pts);
    int np = ch1.size();
    Rcpp::NumericMatrix resu(np+1, 2);
    for(int i=0; i<np; i++) {
      Rcomplex pp = ch1[i];
      resu(i,0) = pp.r;
      resu(i,1) = pp.i;
    }
    // already closed?? Not that sure...
    // close convex hull polygon
    resu(np,0) = resu(0,0);
    resu(np,1) = resu(0,1);
    return(resu);
  } else { // same point
    double e3 = (e1 > e2 ? e1 : e2);
    if (e3 > 0.0) {
      // Return error circle
      Rcomplex c1; c1.r=x1; c1.i=y1;
      Rcpp::ComplexVector ec = mkCc_cpp(e3, c1, fic);
      int np = ec.size();
      Rcpp::NumericMatrix resu(np, 2);
      for (int i=0; i < np; i++) {
        resu(i,0) = ec[i].r;
        resu(i,1) = ec[i].i;
      }
      return(resu);
    } else {
      // Return Polygon = point
      Rcpp::NumericMatrix resu(2, 2);
      resu(0,0) = r[1];
      resu(0,1) = r[2];
      return(resu);
    }
  }
}
