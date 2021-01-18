#ifndef STRUCTS_H
#define STRUCTS_H

#include <vector>
#include <RcppArmadillo.h>


/* DYNAMIC CREATION OF LOCAL VECTORS AND MATRICES */
typedef std::vector<double> DoubleVec;
typedef std::vector<int> IntVec;
typedef std::vector<IntVec> IntMatrix;
typedef std::vector<DoubleVec> DoubleMatrix;

using namespace std;
using namespace Rcpp;
using namespace arma;




class HelperVariables {
public:
  DoubleVec cumy;
  DoubleMatrix cumymat;
  DoubleMatrix cumx;
  DoubleMatrix cumxsq;
  DoubleVec cumysq;
  IntVec cumksize;

  double ybar;
  double ysqall;
  colvec Y;
  mat X;
  uvec pred_cols;


  // Constructors
  HelperVariables(SEXP y, SEXP x) {
    Y = as<colvec>(y);
    X = as<mat>(x);
    pred_cols = zeros<uvec>(X.n_cols-1); // this feels like something that ought not to be needed
    for (int i = 0; i < X.n_cols-1; i++) {
      pred_cols[i] = i+1;
    }
  }
  HelperVariables(SEXP y, SEXP x, SEXP pid) {
    Y = as<colvec>(y);
    X = as<mat>(x);
    IntegerVector id(pid);

    int N = id[id.size()-1]+1; // the number of locations
    int curr_id;

    cumy.push_back(Y[0]);
    cumysq.push_back(pow(Y[0], 2));
    cumksize.push_back(1); // the number of obs per location (I feel this is buggy... check)

    for (int i = 1; i < X.n_cols; i++) {
      DoubleVec cumxvec(N);
      DoubleVec cumxsqvec(N);
      cumxvec[0] = X(0,i);
      cumxsqvec[0] = pow(X(0,i), 2);
      curr_id = 0;
      for (int j = 1; j < Y.n_rows; j++) {
        if (id[j] > curr_id) {
          if (i == 1){
            cumy.push_back(cumy[curr_id] + Y[j]);
            cumysq.push_back(cumysq[curr_id] + pow(Y[j], 2));
            cumksize.push_back(cumksize[curr_id] + 1);
          }
          cumxvec[curr_id+1] = cumxvec[curr_id] + X(j,i);
          cumxsqvec[curr_id+1] = cumxsqvec[curr_id] + pow(X(j,i), 2);
          curr_id++;
        } else {
          if (i == 1) {
            cumy[curr_id] += Y[j];
            cumysq[curr_id] += pow(Y[j], 2);
            cumksize[curr_id]++;
          }
          cumxvec[curr_id] += X(j,i);
          cumxsqvec[curr_id] += pow(X(j,i), 2);
        }
      }
      cumx.push_back(cumxvec);
      cumxsq.push_back(cumxsqvec);
    }
    ybar = cumy[N-1]/Y.n_rows;
  }
  HelperVariables(NumericMatrix data, SEXP pid)
  {
    IntegerVector id(pid);
    int mm = data.ncol();
    int nn2 = data.nrow();
    int N = id[id.size()-1]+1; // number of locations
    // Rprintf("N:%d, idsize:%d, idval1:%d, idval22:%d\n", N, id.size(), id[1], id[22]);

    int curr_id = 0;
    cumksize.push_back(1);

    DoubleVec cumyc(N);
    cumymat.assign(mm, cumyc);
    ysqall = 0.0;
    ybar = 0.0;

    for (int i = 0; i < mm; i++) {
      cumymat[i][0] = data(0, i);
      ysqall += pow(data(0, i), 2);
    }
    for (int j = 1; j < nn2; j++) {
      if (id[j] > curr_id) {
        for (int i = 0; i < mm; i++) {
          cumymat[i][id[j]] = cumymat[i][curr_id] + data(j,i);
          ysqall += pow(data(j, i), 2);
        }
        cumksize.push_back(cumksize[curr_id]+1);
        curr_id++;
      } else {
        for (int i = 0; i < mm; i++) {
          cumymat[i][curr_id] += data(j, i);
          ysqall += pow(data(j, i), 2);
        }
        cumksize[curr_id]++;
      }
    }
    for (int i = 0; i < mm; i++)
      ybar += cumymat[i][N-1];
    ybar /= (nn2 * mm);
    // Rprintf("mm: %d, nn2: %d, N:%d, cumksize.size: %d\n", mm, nn2, N, cumksize.size());
  }
  // other methods
  void print() {
    Rprintf("Helper Variables Print ----\n");
    Rprintf("ybar:%0.2f, cumy[last]:%0.2f", ybar, cumy[Y.n_rows-1]);
    for (int i = 0; i < cumy.size(); i++) {
      Rprintf("i:%d, k:%0.2d, Y:%0.2f, Ysq:%0.2f, X:%0.2f, Xsq:%0.2f\n",
              i, cumksize[i], cumy[i], cumysq[i], cumx[0][i], cumxsq[0][i]);
    }
  }
};

#endif
