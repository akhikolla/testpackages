// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using std::endl;

// [[Rcpp::interfaces(r, cpp)]]

// Augmented Lagrangian method (ALM) nested with Groupwise majorization descent algorithm (GMD)
// for linearly constrained sparse group regression problem. GMD was proposed by Yang and Zou (2014) <doi:10.1007/s11222-014-9498-5>
// The algorithm was designed for sparse log-contrast regression with functional compositional
// predictors proposed by Sun et al. (2020) <<arXiv:1808.02403>.
// [[Rcpp::export]]

Rcpp::List ALM_GMD(arma::vec y,
                   arma::mat Z,
                   arma::mat Zc,
                   arma::mat Zc_proj,
                   arma::vec beta,
                   arma::vec lambda,
                   arma::vec pf,
                   int dfmax,
                   int pfmax,
                   arma::mat A,
                   arma::vec b,
                   arma::umat group_index,
                   int inner_maxiter,
                   int outer_maxiter,
                   double inner_eps,
                   double outer_eps,
                   double mu_ratio,
                   double u_ini,
                   double tol
) {

  int inner_inter, outer_inter, i, j; /* index for loop */
  double inner_err, outer_err; /* convergence for loop */

  int n = y.size();
  int p = group_index.n_cols;
  int k = group_index.n_rows;
  int p1 = p * k;
  int m = Zc.n_cols;
  int pp = p1 + m;
  int nlam = lambda.size();

  int npass = 0;
  double S_norm, tt;
  double u;
  //double upper = 1 ;//+ 1e-10;

  arma::vec dd(p+1);
  arma::vec gamma(p), ita(p);
  arma::vec r(n), C_linear(k), alpha(k);
  arma::vec S(k), xi(p+1), laml(p);
  arma::vec beta_old_inner(pp), beta_old_outer(pp);
  arma::vec diff_group(k), diff_control(m), Ldiff(pp);
  arma::vec alpha_mu(p1);
  arma::uvec ids;

  arma::vec err_flag = arma::zeros(nlam);
  arma::mat path(pp, nlam);
  arma::vec df(nlam), Ni(p);

  group_index -= 1;

  if(mu_ratio < 1) {
    u_ini = 0;
    outer_maxiter = 1;
    for(i = 0; i < p; i++) {
      gamma(i) = max(eig_sym(Z.cols(group_index.col(i)).t() * Z.cols(group_index.col(i))));
    }
    ita.zeros();
  } else {
    for(i = 0; i < p; i++) {
      gamma(i) = max(eig_sym(Z.cols(group_index.col(i)).t() * Z.cols(group_index.col(i))));
      ita(i) = max(eig_sym(A.cols(group_index.col(i)).t() * A.cols(group_index.col(i))));
    }
  }

  gamma = gamma / n; //gamma = gamma * upper / n;
  //ita  *= upper;

  Ni.zeros();
  df.zeros();
  xi(p) = 1;
  /*----------- outer loop for path starts-------------*/

  for(j = 0; j < nlam; j++) {
    /*----------- middle loop ALM starts-------------*/

    laml = pf * lambda(j);
    alpha.zeros();
    r = y - join_rows(Z, Zc) * beta;

    u = u_ini;
    C_linear = A * beta.subvec(0, p1 - 1) - b;
    outer_inter = 1;
    outer_err = 1;
    alpha_mu = A.t() * alpha;
    xi.subvec(0,p-1) = gamma + u * ita;

    while((outer_inter <= outer_maxiter) && (outer_err > outer_eps)) {
      /*----------- inner loop GMD group-lasso stats-------------*/

      beta_old_outer = beta;
      inner_inter = 1;
      inner_err = 1;

      while( (inner_inter <= inner_maxiter) && (inner_err > inner_eps)) {

        npass++;
        beta_old_inner = beta;
        dd.zeros(); //dd = 0;

        //update coefficients for composition variables
        for(i = 0; i < p; i++) {

          diff_group = beta.elem(group_index.col(i));

          S = Z.cols(group_index.col(i)).t() * r / n - alpha_mu.rows(group_index.col(i)) - u * A.cols(group_index.col(i)).t() * C_linear;
          S += xi(i) * diff_group;
          S_norm = sqrt(sum(S % S));

          tt = S_norm - laml(i);

          if( tt > tol ) {
            beta.elem(group_index.col(i)) = S * tt / (xi(i) * S_norm);
          } else {
            beta.elem(group_index.col(i)).fill(0.0);
          }

          diff_group -= beta.elem(group_index.col(i));
          if(any(diff_group != 0.0)) {
            dd(i) = sum(diff_group % diff_group);
            r += Z.cols(group_index.col(i)) * diff_group;
            C_linear -= A.cols(group_index.col(i)) * diff_group;
          }
        }

        //update coefficients for control variables
        diff_control = beta.subvec(p1, p1 + m - 1);
        beta.subvec(p1, p1 + m - 1) = Zc_proj * r + diff_control;
        diff_control -= beta.subvec(p1, p1 + m - 1);
        if(any(diff_control != 0.0)){
          r += Zc * diff_control;
          dd(p) = sum(diff_control % diff_control); //std::max(dd, sum(diff_control % diff_control));
        }

        inner_inter++;
        //check convergence for inner loop
        inner_err = max(dd % (xi % xi));
        if(inner_err <= inner_eps) {
          inner_err =  sum(dd) / std::max(sum(beta_old_inner % beta_old_inner), 1e-30);
        }
      }

      /*----------- inner loop GMD group-lasso ends-------------*/
      outer_inter++;
      //update lagrange multipliers and penalty term
      alpha_mu = alpha_mu + u * A.t() * C_linear;
      u = u * mu_ratio;
      xi.subvec(0,p-1) = gamma + u * ita;

      //check convergence for middle loop
      outer_err = max(abs(C_linear));
      if(outer_err <= outer_eps) {
        Ldiff = beta - beta_old_outer;
        outer_err = sum(Ldiff % Ldiff) / std::max(sum(beta_old_outer % beta_old_outer), 1e-30);
      }

    }

    /*----------- middle loop ALM ends-------------*/

    /*
     QUESTION
     Current: after convergence,
     absolute value of each element of long-coefficient vector beta < tol then set to 0
     Feature: in inner loop
     maximun absolute value of each group coefficient vector < tol then set the whole group to 0
     for(i = 0; i < p; i++) {
        if(max(abs(beta.elem(group_index.col(i)))) < tol ) {
            beta.elem(group_index.col(i)).fill(0)
        }
     }
     */
    ids = find(abs(beta.subvec(0, p1 - 1)) < tol);
    beta.elem(ids).fill(0);

    for(i = 0; i < p; i++) {
      if( any(abs(beta(group_index.col(i))) > 0)) {
        Ni(i) = 1;
        df(j) = df(j) + 1;
      }
    }

    path.col(j) = beta;

    if(outer_err > outer_eps) {
      err_flag(j) = 1;
    } else if (path.col(j).has_nan() || path.col(j).has_inf()) {
      err_flag(j) = 2;
      beta.zeros();
    }

    if(df(j) > dfmax || sum(Ni) > pfmax) break;

  }

  /*----------- outer loop for path ends-------------*/

  //formulate output
  arma::vec output_lambda;
  arma::vec output_df;
  arma::mat output_path(pp, j);
  arma::vec output_error(j);
  int output_ll;
  output_lambda = lambda.subvec(0, j-1);
  output_df = df.subvec(0, j-1);
  output_path = path.cols(0, j-1);
  output_ll = npass;
  output_error = err_flag.subvec(0,j-1);

  return Rcpp::List::create(Rcpp::Named("beta") = output_path,
                            Rcpp::Named("lam") = output_lambda,
                            Rcpp::Named("df") = output_df,
                            Rcpp::Named("npass") = output_ll,
                            //Rcpp::Named("u") = u,
                            Rcpp::Named("error") = output_error);
}



// [[Rcpp::interfaces(r, cpp)]]


// Augmented Lagrangian method (ALM) nested with coordinate descent (CD) for linearly constrained
// lasso regression problem. This algorithm follows the paper Lin et al. (2014) Variable Selection in
// Regression with Compostional <doi:10.1093/biomet/asu031>.
// [[Rcpp::export]]
Rcpp::List ALM_CD(arma::vec y,
                  arma::mat Z,
                  arma::mat Zc,
                  arma::mat Zc_proj,
                  arma::vec beta,
                  arma::vec lambda,
                  arma::vec pf,
                  double b,
                  arma::vec A,
                  int dfmax,
                  int pfmax,
                  int inner_maxiter,
                  int outer_maxiter,
                  double inner_eps,
                  double outer_eps,
                  double mu_ratio,
                  double u_ini,
                  double tol
) {
  int inner_inter, outer_inter, i, j; /* index for loop */
  double inner_err, outer_err; /* convergence for loop */

  int n = y.size();
  int p = Z.n_cols;
  int m = Zc.n_cols;
  int nlam = lambda.size();

  double C_linear;
  double u;
  double S;
  double tt;
  int npass = 0;

  arma::vec gamma(p);
  arma::vec xi(p);
  arma::vec ita(p);
  arma::vec r(n);
  arma::vec Aalpha_mu(p);
  //double alpha_mu = 0;
  arma::vec laml(p);
  arma::vec diff_control(m);
  arma::vec beta_old_inner(p+m), beta_old_outer(p+m), Ldiff(p+m);
  arma::vec err_flag(nlam);
  arma::mat path(p+m, nlam);
  arma::vec df(nlam), Ni(p);
  arma::uvec ids;

  for(i = 0; i < p; i++) {
    gamma(i) = sum(Z.col(i) % Z.col(i));
  }
  gamma = gamma / n;

  if(mu_ratio < 1) {
    u_ini = 0;
    outer_maxiter = 1;
    ita.zeros();
  } else {
    ita = A % A;
  }

  Ni.zeros();
  df.zeros();
  Aalpha_mu.zeros();
  err_flag.zeros();

  /*----------- outer loop for path starts -------------*/

  for(j = 0; j < nlam; j++) {

    /*----------- middle loop ALM starts-------------*/

    laml = pf * lambda(j);
    u = u_ini;

    r = y - join_rows(Z, Zc) * beta;
    //C_linear = sum(beta.subvec(0, p - 1));
    C_linear = sum(A % beta.subvec(0, p - 1))- b;
    xi = gamma + u * ita;

    outer_inter = 1;
    outer_err = 1;

    while((outer_inter <= outer_maxiter) && (outer_err > outer_eps)) {

      /*----------- inner loop Block-wise-lasso stats-------------*/

      beta_old_outer = beta;
      inner_inter = 1;
      inner_err = 1;

      while( (inner_inter <= inner_maxiter) && (inner_err > inner_eps)) {
        npass++;
        beta_old_inner = beta;

        //update coefficients for composition variables
        for(i = 0; i < p; i++) {

          r = r + Z.col(i) * beta(i);
          C_linear = C_linear - A(i) * beta(i);
          S = sum(Z.col(i) % r) / n - u * A(i) * C_linear - Aalpha_mu(i);
          tt = std::abs(S) - laml(i);
          if( tt > 0.0 ) {
            if(S > 0.0) {
              beta(i) =  tt / xi(i);
            } else {
              beta(i) =  -tt / xi(i);
            }
            r = r - Z.col(i) * beta(i);
            C_linear = C_linear + A(i) * beta(i);
          } else {
            beta(i) = 0.0;
          }
        }

        //update coeffients for control variables
        diff_control = beta.subvec(p, p + m - 1);
        beta.subvec(p, p + m - 1) = Zc_proj * r + diff_control;
        diff_control -= beta.subvec(p, p + m - 1);
        if(any(diff_control != 0.0)){
          r = r + Zc * diff_control;
        }

        inner_inter++;
        //check convergence for inner loop
        Ldiff = beta - beta_old_inner;
        inner_err = sum((beta_old_inner - beta) % (beta_old_inner - beta));
        if( inner_err < (inner_eps * (p + 1)) ) {
          inner_err = inner_err / std::max(sum(beta_old_inner % beta_old_inner), 1e-30);
        }

        // internal check point output
        //Rcout << "inner_int: " << inner_inter - 1 << "; inner_err: " << inner_err << endl;
      }

      /*----------- inner loop Block-wise-lasso ends-------------*/

      outer_inter++;
      //update lagrange multiplier and penalty term
      Aalpha_mu = Aalpha_mu +  u * C_linear * A;
      u = u * mu_ratio;
      xi = gamma + u * ita;

      //check convergence for middle loop
      outer_err = std::abs(C_linear);
      if(outer_err <= outer_eps) {
        Ldiff = beta - beta_old_outer;
        outer_err = sum(Ldiff % Ldiff) / std::max(sum(beta_old_outer % beta_old_outer), 1e-30);
      }

      // internal check point output
      //Rcout << "outer_int: " << outer_inter - 1 << "; outer_err: " << outer_err << endl;

    }

    /*----------- middle loop ALM ends-------------*/

    ids = find(abs(beta.subvec(0, p - 1)) < tol);
    beta.elem(ids).fill(0.0);

    for(i = 0; i < p; i++) {
      if(std::abs(beta(i)) > 0.0) {
        Ni(i) = 1;
        df(j) = df(j) + 1;
      }
    }

    path.col(j) = beta;

    if(outer_err > outer_eps) {
      err_flag(j) = 1;
    } else if (path.col(j).has_nan() || path.col(j).has_inf()) {
      err_flag(j) = 2;
      beta.zeros();
    }

    if(df(j) > dfmax || sum(Ni) > pfmax) break;

  }

  /*----------- outer loop for path ends-------------*/

  //formulate output
  arma::vec output_lambda;
  arma::vec output_df;
  arma::mat output_path(p+m, j);
  arma::vec output_error(j);
  int output_ll;
  output_lambda = lambda.subvec(0, j-1);
  output_df = df.subvec(0, j-1);
  output_path = path.cols(0, j-1);
  output_ll = npass;
  output_error = err_flag.subvec(0,j-1);

  return Rcpp::List::create(Rcpp::Named("beta") = output_path,
                            Rcpp::Named("lam") = output_lambda,
                            Rcpp::Named("df") = output_df,
                            Rcpp::Named("npass") = output_ll,
                            //Rcpp::Named("u") = u,
                            Rcpp::Named("error") = output_error);
}


// Augmented Lagrangian method (ALM) nested with coordinate descent (CD) for zero-sum constraints
// with lasso regression problem. This algorithm follows the paper Lin et al. (2014) Variable
// Selection in Regression with Compostional <doi:10.1093/biomet/asu031>.
// [[Rcpp::export]]
Rcpp::List ALM_CD_comp(arma::vec y,
                       arma::mat Z,
                       arma::mat Zc,
                       arma::mat Zc_proj,
                       arma::vec beta,
                       arma::vec lambda,
                       arma::vec pf,
                       int dfmax,
                       int pfmax,
                       int inner_maxiter,
                       int outer_maxiter,
                       double inner_eps,
                       double outer_eps,
                       double mu_ratio,
                       double u_ini,
                       double tol
) {
  int inner_inter, outer_inter, i, j; /* index for loop */
  double inner_err, outer_err; /* convergence for loop */

  int n = y.size();
  int p = Z.n_cols;
  int m = Zc.n_cols;
  int nlam = lambda.size();
  int npass = 0;

  double u;
  double S;
  double tt;
  double C_linear;
  double alpha_mu = 0;

  arma::vec laml(p);
  arma::vec gamma(p);
  arma::vec xi(p);
  arma::vec r(n);

  arma::vec diff_control(m);
  arma::vec beta_old_inner(p+m), beta_old_outer(p+m), Ldiff(p+m);
  arma::vec err_flag = arma::zeros(nlam);
  arma::mat path(p+m, nlam);
  arma::vec df(nlam), Ni(p);
  arma::uvec ids;

  for(i = 0; i < p; i++) {
    gamma(i) = sum(Z.col(i) % Z.col(i));
  }
  gamma = gamma / n;

  if(mu_ratio < 1) {
    u_ini = 0;
    outer_maxiter = 1;
  }

  Ni.zeros();
  df.zeros();

  /*----------- outer loop for path starts-------------*/

  for(j = 0; j < nlam; j++) {

    /*----------- middle loop ALM starts-------------*/
    laml = pf * lambda(j);
    u = u_ini;

    r = y - join_rows(Z, Zc) * beta;
    C_linear = sum(beta.subvec(0, p - 1));
    xi = gamma + u;

    outer_inter = 1;
    outer_err = 1;

    while((outer_inter <= outer_maxiter) && (outer_err > outer_eps)) {

      /*----------- inner loop Block-wise-lasso starts-------------*/

      beta_old_outer = beta;
      inner_inter = 1;
      inner_err = 1;

      while( (inner_inter <= inner_maxiter) && (inner_err > inner_eps)) {

        npass++;
        beta_old_inner = beta;

        //update coeffients for composition variables
        for(i = 0; i < p; i++) {

          r = r + Z.col(i) * beta(i);
          C_linear = C_linear - beta(i);
          S = sum(Z.col(i) % r) / n - u * C_linear - alpha_mu;
          tt = std::abs(S) - laml(i);
          if( tt > 0.0 ) {
            if(S > 0.0) {
              beta(i) =  tt / xi(i);
            } else {
              beta(i) =  -tt / xi(i);
            }

            r = r - Z.col(i) * beta(i);
            C_linear = C_linear + beta(i);
          } else {
            beta(i) = 0.0;
          }
        }

        //update coeffients for control variables
        diff_control = beta.subvec(p, p + m - 1);
        beta.subvec(p, p + m - 1) = Zc_proj * r + diff_control;
        diff_control -= beta.subvec(p, p + m - 1);
        if(any(diff_control != 0.0)){
          r = r + Zc * diff_control;
        }

        inner_inter++;

        //check convergence for inner loop
        Ldiff = beta - beta_old_inner;
        inner_err = sum((beta_old_inner - beta) % (beta_old_inner - beta));
        if( inner_err < (inner_eps * (p + 1)) ) {
          inner_err = inner_err / std::max(sum(beta_old_inner % beta_old_inner), 1e-30);
        }

        //internal check point output
        //Rcout << "inner_int: " << inner_inter << "; inner_err: " << inner_err << endl;

      }

      /*----------- inner loop Block-wise-lasso ends-------------*/

      outer_inter++;

      //update lagrange multiplier and penalty term
      alpha_mu = alpha_mu + u * C_linear;
      u = u * mu_ratio;
      xi = gamma + u;

      //check convergence for middle loop
      outer_err = std::abs(C_linear);
      if(outer_err <= outer_eps) {
        Ldiff = beta - beta_old_outer;
        outer_err = sum(Ldiff % Ldiff) / std::max(sum(beta_old_outer % beta_old_outer), 1e-30);
      }

    }

    /*----------- middle loop ALM ends-------------*/

    ids = find(abs(beta.subvec(0, p - 1)) < tol);
    beta.elem(ids).fill(0.0);

    for(i = 0; i < p; i++) {
      if(std::abs(beta(i)) > 0.0) {
        Ni(i) = 1;
        df(j) = df(j) + 1;
      }
    }

    path.col(j) = beta;

    if(outer_err > outer_eps) {
      err_flag(j) = 1;
    } else if (path.col(j).has_nan() || path.col(j).has_inf()) {
      err_flag(j) = 2;
      beta.zeros();
    }

    if(df(j) > dfmax || sum(Ni) > pfmax) break;

  }

  /*----------- outer loop for path ends-------------*/

  //formulate output
  arma::vec output_lambda;
  arma::vec output_df;
  arma::mat output_path(p+m, j);
  arma::vec output_error(j);
  int output_ll;
  output_lambda = lambda.subvec(0, j-1);
  output_df = df.subvec(0, j-1);
  output_path = path.cols(0, j-1);
  output_ll = npass;
  output_error = err_flag.subvec(0,j-1);

  return Rcpp::List::create(Rcpp::Named("beta") = output_path,
                            Rcpp::Named("lam") = output_lambda,
                            Rcpp::Named("df") = output_df,
                            Rcpp::Named("npass") = output_ll,
                            //Rcpp::Named("u") = u,
                            Rcpp::Named("error") = output_error);
}
