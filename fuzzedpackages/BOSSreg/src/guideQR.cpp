// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

// givens rotation
// algorithm 1.1 in the reference: Hammarling (2008) 'Updating the QR Factorization and the Least Squares Problem'
void givens(double a, double b, double &c, double &s){
  if(b == 0){
    c = 1;
    s = 0;
  }else{
    if(fabs(b) >= fabs(a)){
      double t = -a / b;
      s = 1 / sqrt(1 + pow(t,2));
      c = s * t;
    }else{
      double t = -b / a;
      c = 1 / sqrt(1 + pow(t,2));
      s = c * t;
    }
  }
}

// update Qr and v
void updateQv(arma::vec &v, arma::mat &Qr){
  int r = Qr.n_cols;
  double c, s;
  for (int i=r-1; i>=1; i--){
    givens(v(i-1), v(i), c, s);
    // update v
    v(i-1) = c * v(i-1) - s * v(i);
    // update Q
    arma::mat rotate = { {c, s},
                         {-s, c} };
    Qr.cols(i-1,i) =  Qr.cols(i-1,i) * rotate;
  }
}

// update QR decomposition
void updateQR(arma::mat &Ql, arma::mat &Qr, arma::mat &R, arma::vec newcol){
  int n = Ql.n_rows;
  int l = Ql.n_cols;
  arma::mat Q = join_rows(Ql,Qr);
  arma::vec v = Q.t() * newcol;
  arma::vec v_toupdate = v.subvec(l,n-1);
  updateQv(v_toupdate, Qr);
  // update v
  v.subvec(l,n-1) = v_toupdate;
  // update Q
  Ql = join_rows(Ql, Qr.col(0));
  Qr.shed_col(0);
  // update R
  R = join_cols(R, arma::zeros(l).t());
  R = join_rows(R, v.subvec(0,l));
}

// guided QR decomposition, based on semi-partial correlation
// [[Rcpp::export]]
List guideQR(arma::mat x, arma::vec y, int maxstep){
  // int n = x.n_rows;
  int p = x.n_cols;
  // which_remain stores the indices of remaining preditors
  arma::vec which_remain = arma::vec(maxstep);
  std::iota (std::begin(which_remain), std::end(which_remain), 0);

  arma::mat Q, R, x_remain = x;
  arma::vec steps = arma::vec(maxstep);

  // first variable to step in
  arma::vec cor_tmp = arma::abs(x.t() * y);
  arma::uword i = cor_tmp.index_max();
  steps(0) = i;
  which_remain.shed_row(i);
  x_remain.shed_col(i);
  arma::qr(Q, R, x.col(i));
  arma::mat Ql = Q, Qr = Q;
  Ql = Ql.col(0);
  Qr.shed_col(0);
  R = R.row(0);

  arma::mat resid_tmp = x_remain;
  int j_remain, j_tmp;
  // step 2 and beyond
  //int i_step = 1;
  for (int i_step=1; i_step<=maxstep-1; i_step++){
    // calculate the semi-partial correlations for each remaining predictor
    if(i_step < p-1){
      j_tmp = Ql.n_cols - 1; // -1 to match index in C
      resid_tmp = resid_tmp - Ql.col(j_tmp) * Ql.col(j_tmp).t() * x_remain;
      j_remain = arma::abs(y.t() * arma::normalise(resid_tmp)).index_max();
    }else{
      j_remain = 0;
    }

    steps(i_step) = which_remain(j_remain);
    // update QR
    updateQR(Ql, Qr, R, x_remain.col(j_remain));

    // update others
    if(i_step < p-1){
      resid_tmp.shed_col(j_remain);
      x_remain.shed_col(j_remain);
      which_remain.shed_row(j_remain);
    }
  }

  // add 1 to match the index in R
  steps = steps + 1;
  return List::create(Named("Q") = Ql,
                      Named("R") = R,
                      Named("steps") = steps);

}
