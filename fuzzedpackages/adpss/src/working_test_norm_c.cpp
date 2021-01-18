#include <Rcpp.h>
#include <R.h>
// using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]] 
//' @useDynLib adpss, .registration = TRUE
//' @importFrom Rcpp sourceCpp

////static int glb_sgn = 0;

template <class c1>
static double bisection_inverse(double (*fx)(double, c1 *), double y, c1 *info,
  double sol_l = 0, double sol_u = 10, 
  bool larger = false, bool smaller = false, bool exact = false,
  double prec = 1e-8) {
    //Rcpp::Rcout << "# bisection_inverse # START" << std::endl;

  // Search target range //
  double itv = sol_u - sol_l;

    //Rcpp::Rcout << "# bisection_inverse # fx(lower_lim) = fx(" << sol_l << ")" << std::endl;
  double fx_l = (*fx)(sol_l, info);
    //Rcpp::Rcout << "# bisection_inverse # fx(lower_lim) = fx(" << sol_l << ") = " << fx_l << std::endl;
    //Rcpp::Rcout << "# bisection_inverse # fx(upper_lim) = fx(" << sol_u << ")" << std::endl;
  double fx_u = (*fx)(sol_u, info);
    //Rcpp::Rcout << "# bisection_inverse # fx(upper_lim) = fx(" << sol_u << ") = " << fx_u << std::endl;
  int sgn = (fx_l < fx_u) - (fx_l > fx_u);
    //glb_sgn = sgn;
    //Rcpp::Rcout << "# bisection_inverse # ..... sol_l = " << sol_l << "; fx_l = " << fx_l 
    //  << "; sol_u = " << sol_u << "; fx_u = " << fx_u << "; sign = " << sgn << "." << std::endl;

  fx_l = sgn * fx_l;
  fx_u = sgn * fx_u;
  y = sgn * y;

  // Sub Settings //
  if ( exact ) { prec = 0; smaller = false; larger = false; }
  if ( prec == 0 ) { exact = TRUE; smaller = FALSE; larger = FALSE; }
  if ( larger ) { smaller = FALSE; exact = FALSE; }
  if ( smaller ) { larger = FALSE; exact = FALSE; }
  if ( sgn == -1 && !exact ) { larger = !larger; smaller = !smaller; }

    //Rcpp::Rcout << "# bisection_inverse # ..... Passed initial check" << std::endl;

  // Search //
  int ii = 0;
  while( 1 ){
    R_CheckUserInterrupt();
    //Rcpp::Rcout << "# bisection_inverse # ..... Search solution range [" << sol_l << ", " << sol_u << "]" << std::endl;
    //Rcpp::Rcout << "# bisection_inverse # ..... sol_l = " << sol_l << "; fx_l = " << fx_l 
    //  << "; sol_u = " << sol_u << "; fx_u = " << fx_u << "; sign = " << sgn << "." << std::endl;
    if ( fx_l <= y && fx_u >= y ) {
      // Found //
      break;
    } else if ( fx_l < y && fx_u < y ) {
      // x > sol_u //
      sol_u = sol_u + itv;
      sol_l = sol_l + itv;
      fx_l = fx_u;
      fx_u = sgn * (*fx)(sol_u, info);
    } else if ( fx_l > y && fx_u > y ) {
      // x < sol_l //
      sol_u = sol_u - itv;
      sol_l = sol_l - itv;
      fx_u = fx_l;
      fx_l = sgn * (*fx)(sol_l, info);
    }
    ii += 1;
    //Rcpp::Rcout << "ii: " << ii << std::endl;
  }
  //Rcpp::Rcout << "# bisection_inverse # Solution lies between: [" << sol_l << ", " << sol_u << "]" << std::endl;

  // Search solution //
  //Rcpp::Rcout << "# bisection_inverse # Is exact solution on the bound?" << std::endl;
  double sol;
  // Exact solution is on range limit //
  if ( fx_l == y ){
    sol = sol_l;
    //Rcpp::Rcout << "# bisection_inverse # Lower bound is exact." << std::endl;
    //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
    return sol;
  } else if ( fx_u == y ){
    sol = sol_u;
    //Rcpp::Rcout << "# bisection_inverse # Upper bound is exact." << std::endl;
    //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
    return sol;
  }

  // Search //
  //Rcpp::Rcout << "# bisection_inverse # Search solution." << std::endl;
  double step = itv / 2., fx_m, fx_mw;
  sol = sol_l + step;
  while(1){
    R_CheckUserInterrupt();
    fx_m = (*fx)(sol, info);
    fx_mw = sgn * fx_m;
    //Rcpp::Rcout << "# bisection_inverse # ..... x: " << sol << "; fx: " << fx_m << "." << std::endl;  //%%%%%%%%
    if ( fx_mw == y ) {
      // Exact solution //
        //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
      return sol;
    } else if ( fx_mw < y ) {
      // fx_l < x < fx_m //
      if ( smaller && step <= prec ) {
        //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
        return sol;
      } else if ( larger && step  <= prec ) {
        //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
        return sol_u;
      } else {
        fx_l = fx_mw;
        sol_l = sol;
        step = step / 2.;
        sol = sol + step;
      }
    } else if ( fx_mw > y ) {
      // fx_m < x < fx_u //
      if ( larger && step <= prec ) {
        //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
        return sol;
      } else if ( smaller && step  <= prec ) {
        //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
        return sol_l;
      } else {
        fx_u = fx_mw;
        sol_u = sol;
        step = step / 2.;
        sol = sol - step;
      }
    }
  }
}

/*
template <class c1>
static double bisection_inverse_print(double (*fx)(double, c1 *), double y, c1 *info,
  double sol_l = 0, double sol_u = 10, 
  bool larger = false, bool smaller = false, bool exact = false,
  double prec = 1e-8) {
    Rcpp::Rcout << "# bisection_inverse # START" << std::endl;

  // Search target range //
  double itv = sol_u - sol_l;

    //Rcpp::Rcout << "# bisection_inverse # fx(lower_lim) = fx(" << sol_l << ")" << std::endl;
  double fx_l = (*fx)(sol_l, info);
    //Rcpp::Rcout << "# bisection_inverse # fx(lower_lim) = fx(" << sol_l << ") = " << fx_l << std::endl;
    //Rcpp::Rcout << "# bisection_inverse # fx(upper_lim) = fx(" << sol_u << ")" << std::endl;
  double fx_u = (*fx)(sol_u, info);
    //Rcpp::Rcout << "# bisection_inverse # fx(upper_lim) = fx(" << sol_u << ") = " << fx_u << std::endl;
  int sgn = (fx_l < fx_u) - (fx_l > fx_u);
    //glb_sgn = sgn;
    Rcpp::Rcout << "# bisection_inverse # ..... sol_l = " << sol_l << "; fx_l = " << fx_l 
      << "; sol_u = " << sol_u << "; fx_u = " << fx_u << "; sign = " << sgn << "." << std::endl;

  fx_l = sgn * fx_l;
  fx_u = sgn * fx_u;
  y = sgn * y;

  // Sub Settings //
  if ( exact ) { prec = 0; smaller = false; larger = false; }
  if ( prec == 0 ) { exact = TRUE; smaller = FALSE; larger = FALSE; }
  if ( larger ) { smaller = FALSE; exact = FALSE; }
  if ( smaller ) { larger = FALSE; exact = FALSE; }
  if ( sgn == -1 && !exact ) { larger = !larger; smaller = !smaller; }

    Rcpp::Rcout << "# bisection_inverse # ..... Passed initial check" << std::endl;

  // Search //
  int ii = 0;
  while( 1 ){
    R_CheckUserInterrupt();
    Rcpp::Rcout << "# bisection_inverse # ..... Search solution range [" << sol_l << ", " << sol_u << "]" << std::endl;
    Rcpp::Rcout << "# bisection_inverse # ..... sol_l = " << sol_l << "; fx_l = " << fx_l 
      << "; sol_u = " << sol_u << "; fx_u = " << fx_u << "; sign = " << sgn << "." << std::endl;
    if ( fx_l <= y && fx_u >= y ) {
      // Found //
      break;
    } else if ( fx_l < y && fx_u < y ) {
      // x > sol_u //
      sol_u = sol_u + itv;
      sol_l = sol_l + itv;
      fx_l = fx_u;
      fx_u = sgn * (*fx)(sol_u, info);
    } else if ( fx_l > y && fx_u > y ) {
      // x < sol_l //
      sol_u = sol_u - itv;
      sol_l = sol_l - itv;
      fx_u = fx_l;
      fx_l = sgn * (*fx)(sol_l, info);
    }
    ii += 1;
    Rcpp::Rcout << "ii: " << ii << std::endl;
  }
  Rcpp::Rcout << "# bisection_inverse # Solution lies between: [" << sol_l << ", " << sol_u << "]" << std::endl;

  // Search solution //
  Rcpp::Rcout << "# bisection_inverse # Is exact solution on the bound?" << std::endl;
  double sol;
  // Exact solution is on range limit //
  if ( fx_l == y ){
    sol = sol_l;
    //Rcpp::Rcout << "# bisection_inverse # Lower bound is exact." << std::endl;
    //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
    return sol;
  } else if ( fx_u == y ){
    sol = sol_u;
    //Rcpp::Rcout << "# bisection_inverse # Upper bound is exact." << std::endl;
    //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
    return sol;
  }

  // Search //
  Rcpp::Rcout << "# bisection_inverse # Search solution." << std::endl;
  double step = itv / 2., fx_m, fx_mw;
  sol = sol_l + step;
  while(1){
    R_CheckUserInterrupt();
    fx_m = (*fx)(sol, info);
    fx_mw = sgn * fx_m;
    Rcpp::Rcout << "# bisection_inverse # ..... x: " << sol << "; fx: " << fx_m << "." << std::endl;  //%%%%%%%%
    if ( fx_mw == y ) {
      // Exact solution //
        //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
      return sol;
    } else if ( fx_mw < y ) {
      // fx_l < x < fx_m //
      if ( smaller && step <= prec ) {
        //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
        return sol;
      } else if ( larger && step  <= prec ) {
        //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
        return sol_u;
      } else {
        fx_l = fx_mw;
        sol_l = sol;
        step = step / 2.;
        sol = sol + step;
      }
    } else if ( fx_mw > y ) {
      // fx_m < x < fx_u //
      if ( larger && step <= prec ) {
        //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
        return sol;
      } else if ( smaller && step  <= prec ) {
        //Rcpp::Rcout << "# bisection_inverse # RETURN" << std::endl;
        return sol_l;
      } else {
        fx_u = fx_mw;
        sol_u = sol;
        step = step / 2.;
        sol = sol - step;
      }
    }
  }
}
*/

// my struct base;
struct base {
  double* hypotheses;
  double* prior;
  int num_of_hypotheses;
  double effect_size;
  double cost0;
  double cost1;
  int analysis;
  double max_time;
  int work_KK;
  double* U_k;
  double fss_w;
  double seq_w;
  double stat;
  double* xdev;
  int xdev_l;
  double* up_wing_temp;
  double* up_wing_buffer;
  double* div_unit;
  int* up_wing_units;
  int* wing_l;
  int cc_k_add;
  double tol_boundary;
} ;

// my struct base_time;
struct base_time {
  struct base str_base;
  int stage;
  double time; // t_k
} ;

// my struct current_next;
struct current_next {
  struct base_time str_base_time;
  double time_1; // next time; "_1" indicates "next";
  double* gg_k_1;
  int gg_k_1_l;
  double* value_1;
  double* dummy;
  int* xdev_k;
  double* cc_k_1;
  int* cc_odd_n_k_1;
} ;

// my struct result
struct result {
  int cc_k_add;
  std::vector<double>* vss0;
  int* gg_odd_l;
  std::vector<double>* vgg_odd;
  int* gg_l;
  std::vector<double>* vgg;
  std::vector<double>* veta;
  std::vector<double>* vpr_rej_H0;
  std::vector<double>* vexp_time;
  int* cc_odd_n;
  int* cc_n;
  double* cc;
  double* val_k;
} ;

// my struct ground
struct ground {
  struct base str_base;
  struct result str_result;
} ;

// my function posterior
static std::vector<double> vposterior(double xx, struct base_time* pinfo) {
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# vposterior # START" << std::endl;
  struct base* str_base;
  str_base = &(pinfo->str_base);
  int num_of_hypotheses;
  num_of_hypotheses = str_base->num_of_hypotheses;

  std::vector<double> vdummyArr(num_of_hypotheses);
  //double* dummyArr = vdummyArr.data();
  double dummyVar[2] = {};

  if ( pinfo->time == 0 ) {
    for ( int ii = 0; ii < num_of_hypotheses; ii++ ) { 
      vdummyArr[ii] = (str_base->prior)[ii]; // vdummyArr = post
      //Rcpp::Rcout << "# vposterior # prior[" << ii << "] = " << vdummyArr[ii] << std::endl;
    }
  } else {

    for ( int ii = 0; ii < num_of_hypotheses; ii++ ) { // vdummyArr = log_lik, dummyVar[0] = mean(log_lik)
      vdummyArr[ii] = - pow((xx / pinfo->time) - str_base->hypotheses[ii], 2.) * pinfo->time / 2.;
      dummyVar[0] += vdummyArr[ii] / num_of_hypotheses;
        //if( glb_sgn == 1 ) Rcpp::Rcout << "# vposterior # log_lik(" << ii << ") = " << dummyVar[0] << std::endl;
    }
      //if( glb_sgn == 1 ) Rcpp::Rcout << "# vposterior # mean(log_lik): " << dummyVar[0] << std::endl;

    for ( int ii = 0; ii < num_of_hypotheses; ii++ ) { 
      vdummyArr[ii] = vdummyArr[ii] - dummyVar[0];  // vdummyArr = rel_log_lik_i
      vdummyArr[ii] = str_base->prior[ii] * exp(vdummyArr[ii]); // vdummyArr = rel_post_rel
      dummyVar[1] += vdummyArr[ii]; // dummyVar[1] = rel_pr_data
    }
      //if( glb_sgn == 1 ) Rcpp::Rcout << "# vposterior # rel_pr_data: " << dummyVar[1] << std::endl;

    for ( int ii = 0; ii < num_of_hypotheses; ii++ ) { 
      vdummyArr[ii] = vdummyArr[ii] / dummyVar[1]; // vdummyArr = post
      //if( glb_sgn == 1 ) Rcpp::Rcout << "# vposterior # posterior(" << ii << ") = " << dummyVar[1] << std::endl;
      //Rcpp::Rcout << "# vposterior # prior[" << ii << "] = " << vdummyArr[ii] << std::endl;
    }
  }

    //if( pinfo->stage == 0 ) {
    //  for ( int ii = 0; ii < num_of_hypotheses; ii++ ) { 
    //    Rcpp::Rcout << "post pr: " << vdummyArr[ii] << std::endl;
    //  }
    //}

    //for ( int ii = 0; ii < num_of_hypotheses; ii++ ) { 
    //  Rcpp::Rcout << "# vposterior # prior[" << ii << "] = vdummyArr = " << vdummyArr.at(ii) << std::endl;
    //}
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# vposterior # RETURN" << std::endl;
  return vdummyArr;
}

// my function posterior01
static std::vector<double> vposterior01(double xx, struct base_time* pinfo) {
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# vposterior01 # START" << std::endl;
  struct base* str_base;
  str_base = &(pinfo->str_base);
  int num_of_hypotheses;
  num_of_hypotheses = str_base->num_of_hypotheses;

  std::vector<double> vdummyArr(num_of_hypotheses);
  //double* dummyArr = vdummyArr.data();
  double dummyVar[2] = {};

  vdummyArr = vposterior(xx, pinfo);
    //Rcpp::Rcout << "# vposterior01 # num_of_hypotheses = " << num_of_hypotheses << "; vdummyArr.size() = " << vdummyArr.size() << std::endl;
    //Rcpp::Rcout << "# vposterior01 # dummyArr = " << dummyArr << "; vdummyArr.data() = " << vdummyArr.data() << std::endl;
    //for ( int ii = 0; ii < num_of_hypotheses; ii++ ) { 
    //  Rcpp::Rcout << "# vposterior01 # prior: vdummyArr.at(" << ii << ") = " << vdummyArr.at(ii) << "; vdummyArr[] = " << vdummyArr[ii] << std::endl;
    //}
  *dummyVar = 0;
  for ( int ii = 0; ii < num_of_hypotheses; ii++ ) {
    *dummyVar += ((double) (ii > 0)) * vdummyArr[ii];
  }
  double post_H_th0, post_H_th1;
  post_H_th0 = vdummyArr[0];
  post_H_th1 = *dummyVar;

  std::vector<double> vres(2);
  vres.at(0) = post_H_th0;
  vres.at(1) = post_H_th1;
    //Rcpp::Rcout << "# vposterior # posterior(H0) = " << post_H_th0 << "; posterior(H1) = " << post_H_th1 << std::endl;

    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# vposterior01 # RETURN" << std::endl;
  return vres;
}

// my function current_risk
static double current_risk(double xx, struct base_time* pinfo) {
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# current_risk # START" << std::endl;
  struct base* str_base;
  str_base = &(pinfo->str_base);
  int num_of_hypotheses;
  num_of_hypotheses = str_base->num_of_hypotheses;

  std::vector<double> vdummyVar(2);
  //double* dummyVar = vdummyVar.data();

  vdummyVar = vposterior01(xx, pinfo);
  double post_H_th0 = vdummyVar[0];
  double post_H_th1 = vdummyVar[1];
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# current_risk # posterior(H0) = " << post_H_th0 << "; posterior(H1) = " << post_H_th1 << std::endl;

  vdummyVar[0] = str_base->cost0 * post_H_th0;
  vdummyVar[1] = (str_base->fss_w * str_base->cost1 + str_base->seq_w *
                 (str_base->max_time - pinfo->time)) * post_H_th1;
  num_of_hypotheses = vdummyVar[0] > vdummyVar[1]; // num_of_hypotheses as dummy
  vdummyVar[0] = vdummyVar[num_of_hypotheses];

  return vdummyVar[0];
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# current_risk # RETURN" << std::endl;
}

// my function current_risk_balance
static double current_risk_balance(double xx, struct base_time* pinfo) {
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# current_risk_balance # START" << std::endl;
  struct base* str_base;
  str_base = &(pinfo->str_base);
  //int num_of_hypotheses;
  //num_of_hypotheses = str_base->num_of_hypotheses;

  std::vector<double> vdummyVar(2);
  //double* dummyVar = vdummyVar.data();

  vdummyVar = vposterior01(xx, pinfo);
  double post_H_th0 = vdummyVar[0];
  double post_H_th1 = vdummyVar[1];
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# current_risk_balance # posterior(H0) = " << post_H_th0 << "; posterior(H1) = " << post_H_th1 << std::endl;

  vdummyVar[0] = str_base->cost0 * post_H_th0;
  vdummyVar[1] = (str_base->fss_w * str_base->cost1 + str_base->seq_w *
                 (str_base->max_time - pinfo->time)) * post_H_th1;
    //Rcpp::Rcout << "# current_risk_balance # cost0 = " << str_base->cost0 << "; post_H_th0 = " << post_H_th0 << std::endl;
    //Rcpp::Rcout << "# current_risk_balance # cost(H0) = " << vdummyVar[0] << "; cost(H1) = " << vdummyVar[1] << "; cost(H0) - cost(H1) = " << vdummyVar[0] - vdummyVar[1] << std::endl;

    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# current_risk_balance # RETURN" << std::endl;
  return vdummyVar[0] - vdummyVar[1];
}



// my function future_risk
static double future_risk(double xx, struct current_next* pinfo) {
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # START" << std::endl;
  struct base_time* str_base_time;
  str_base_time = &(pinfo->str_base_time);
  struct base* str_base;
  str_base = &(str_base_time->str_base);

  int num_of_hypotheses;
  num_of_hypotheses = str_base->num_of_hypotheses;
  double* H_th = str_base->hypotheses;

  std::vector<double> vxx_posterior(num_of_hypotheses);
  //double* xx_posterior = vxx_posterior.data();
  vxx_posterior = vposterior(xx, str_base_time);

  //double post_H_th0 = *vxx_posterior;
  double post_H_th1 = 0;
  for ( int ii = 1; ii < num_of_hypotheses; ii++ ) { 
    post_H_th1 += vxx_posterior[ii];
  }
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # posterior(H0) = " << post_H_th0 << "; posterior(H1) = " << post_H_th1 << std::endl;

  double* xx_1 = pinfo->gg_k_1;
  //int xx_1_l = pinfo->gg_k_1_l;
  double* value_1 = pinfo->value_1;
  //double* doubleVar = pinfo->dummy;

    //for( int ii = 0; ii < xx_1_l; ii++ ) { 
    //  Rcpp::Rcout << "gg[" << ii << "]; " << xx_1[ii] << std::endl;
    //}

  int ll;
  double dif, pr, ww, xi1, xi2;
  double d_t = pinfo->time_1 - str_base_time->time;
  double sq_d_t = sqrt(d_t);

  int xdev_l = str_base->xdev_l;
  double* div_unit = str_base->div_unit;
  int* up_wing_units = str_base->up_wing_units;
  int* wing_l = str_base->wing_l;

  int kk = str_base_time->stage;
  int* xdev_k = pinfo->xdev_k;
  //double* cc_k_1 = pinfo->cc_k_1;
  int* cc_odd_n_k_1 = pinfo->cc_odd_n_k_1;
  int cc_k_add = str_base->cc_k_add;
    //Rcpp::Rcout << "# future_risk # cc_k_add = " << cc_k_add << std::endl;

  std::vector<int> vsimp_odd_1(xdev_l + cc_k_add + 1);  // 1 is added for overrun in 'loop for simp_odd_1'
  //vsimp_odd_1.resize(xdev_l + cc_k_add + 1);
  int* simp_odd_1 = vsimp_odd_1.data();
  std::vector<int> vsimp_1(2 * (xdev_l + cc_k_add) - 1);
  int* simp_1 = vsimp_1.data();

  // nodes for simpson's method //
  int gg_odd_k_1_i_xx = up_wing_units[kk + 1] - floor(xx / div_unit[kk + 1]);
  int gg_odd_k_1_i_cc_odd_n_k_1 = 0;
  int simp_odd_1_l = 0;
  int xdev_odd_k_1_base, xdev_odd_k_1, intVar, xdev_odd_k_1_proper;
    //Rcpp::Rcout << "gg_odd_k_1_i_xx" << ": " << gg_odd_k_1_i_xx << std::endl;
    //Rcpp::Rcout << "gg_odd_k_1_i_cc_odd_n_k_1" << ": " << gg_odd_k_1_i_cc_odd_n_k_1 << std::endl;
    //Rcpp::Rcout << "simp_odd_1_l" << ": " << simp_odd_1_l << std::endl;
  for ( int ii = 0; ii < xdev_l; ii++ ) { // loop for simp_odd_1
    xdev_odd_k_1_base = gg_odd_k_1_i_xx + xdev_k[ii];
      //Rcpp::Rcout << "xdev_odd_k_1_base[" << ii << "]: " << xdev_odd_k_1_base << std::endl;
    intVar = xdev_odd_k_1_base >= *cc_odd_n_k_1;
    xdev_odd_k_1 = xdev_odd_k_1_base + cc_k_add * intVar;
    simp_odd_1[simp_odd_1_l + cc_k_add * intVar] = xdev_odd_k_1;
    xdev_odd_k_1_proper = (xdev_odd_k_1 >= 0) && (xdev_odd_k_1 < (wing_l[kk + 1] + cc_k_add));
    simp_odd_1_l += xdev_odd_k_1_proper;
    gg_odd_k_1_i_cc_odd_n_k_1 += (1 - intVar) * xdev_odd_k_1_proper;
  }
  simp_odd_1_l += cc_k_add;
  for ( int ii = 0; ii < cc_k_add; ii++ ) {
    simp_odd_1[gg_odd_k_1_i_cc_odd_n_k_1 + ii] = *cc_odd_n_k_1 + ii;
  }
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # ... (xx_1_l + 1) / 2 = " << (xx_1_l + 1) / 2. << " ; " << wing_l[kk + 1] + cc_k_add << " = wing_l[" << kk << "] + cc_k_add." << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # simp_odd_1_l: " << simp_odd_1_l << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # gg_odd_k_1_i_cc_odd_n_k_1: " << gg_odd_k_1_i_cc_odd_n_k_1 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # cc_odd_n_k_1: " << *cc_odd_n_k_1 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # vsimp_odd_1.data(): " << vsimp_odd_1.data() << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # simp_odd_1*: " << simp_odd_1 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # vsimp_odd_1.size(): " << vsimp_odd_1.size() << std::endl;
    //for ( int ii = 0; ii < (xdev_l + cc_k_add + 1); ii++ ) {
    //  if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # simp_odd_1[" << ii << "]: " << simp_odd_1[ii] << std::endl;
    //}

  // add half points for nodes //
  *simp_1 = *simp_odd_1 * 2;
    //Rcpp::Rcout << "simp_1[" << 0 << "]: " << simp_1[0] << std::endl;
  for ( int ii = 1; ii < simp_odd_1_l; ii++ ) {
    simp_1[2 * ii] = simp_odd_1[ii] * 2;
    simp_1[2 * ii - 1] = simp_odd_1[ii - 1] + simp_odd_1[ii];
      //Rcpp::Rcout << "simp_1[" << 2 * ii - 1 << "]: " << simp_1[2 * ii - 1] << std::endl;
      //Rcpp::Rcout << "simp_1[" << 2 * ii << "]: " << simp_1[2 * ii] << std::endl;
  }

  xi1 = 0;
  xi2 = 0;
    //Rcpp::Rcout << "xx: " << xx << "; simp_odd_1_l: " << simp_odd_1_l << std::endl;
    //Rcpp::Rcout << "simp_index[" << "]*: " << simp_index << std::endl;
  int simp_1_l = simp_odd_1_l * 2 - 1;
  for ( int ii = 0; ii < simp_1_l; ii++ ) {
    for ( int jj = 0; jj < num_of_hypotheses; jj++ ) {  // for num_of_hypothesis;
      ll = simp_1[ii];
      dif = xx_1[ll] - xx;
      pr = R::dnorm(dif, H_th[jj] * d_t, sq_d_t, 0);
      ww = (xx_1[simp_1[ii - 1 + (ii == 0)]] - xx_1[simp_1[ii + 1 - (ii == (simp_1_l - 1))]]) * (1 + (ii % 2)) / 3.;
      xi2 += ww * vxx_posterior[jj] * pr * value_1[ll];
        //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # xx_1[" << ll << "]: " << xx_1[ll] << "; dif: " << dif << "; pr: " << pr << "; ww: " << ww << "; val: " << value_1[ll] << std::endl;
        //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # xx_1[" << ll << "]: " << xx_1[ll] << "; xi2: " << xi2 << "; xi2+: " << ww * vxx_posterior[jj] * pr * value_1[ll] << std::endl;
    }
  }
  xi1 = str_base->seq_w * post_H_th1 * d_t;
    //Rcpp::Rcout << "xi: " << xi1 << "; xi2: " << xi2 << std::endl;

    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # xi1 = " << xi1 << "; xi2 = " << xi2 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # RETURN" << std::endl;
  return xi1 + xi2;
}


// my function future_risk0
static double future_risk0(double xx, struct current_next* pinfo) {
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # START" << std::endl;
  struct base_time* str_base_time;
  str_base_time = &(pinfo->str_base_time);
  struct base* str_base;
  str_base = &(str_base_time->str_base);

  double effect_size = str_base->effect_size;

  double* xx_1 = pinfo->gg_k_1;
  //int xx_1_l = pinfo->gg_k_1_l;
  double* value_1 = pinfo->value_1;
  //double* doubleVar = pinfo->dummy;

  //for( int ii = 0; ii < *xx_1_l; ii++ ) { 
  //  Rcpp::Rcout << "gg[" << ii << "]; " << xx_1[ii] << std::endl;
  //}

  int ll;
  double dif, pr, ww, xi2;
  double d_t = pinfo->time_1 - str_base_time->time;
  double sq_d_t = sqrt(d_t);

  int xdev_l = str_base->xdev_l;
  double* div_unit = str_base->div_unit;
  int* up_wing_units = str_base->up_wing_units;
  int* wing_l = str_base->wing_l;

  int kk = str_base_time->stage;
  int* xdev_k = pinfo->xdev_k;
  //double* cc_k_1 = pinfo->cc_k_1;
  int* cc_odd_n_k_1 = pinfo->cc_odd_n_k_1;
  int cc_k_add = str_base->cc_k_add;
    //Rcpp::Rcout << "# future_risk # cc_k_add = " << cc_k_add << std::endl;

  int cc_k_add_ = 1; // only 1 cc_k is included in simp_odd_1 in this function
  std::vector<int> vsimp_odd_1(xdev_l + cc_k_add_ + 1);  // 1 is added for overrun in 'loop for simp_odd_1'
  //vsimp_odd_1.resize(xdev_l + cc_k_add_ + 1);
  int* simp_odd_1 = vsimp_odd_1.data();
  std::vector<int> vsimp_1(2 * (xdev_l + cc_k_add_) - 1);
  int* simp_1 = vsimp_1.data();

  // nodes for simpson's method //
  int gg_odd_k_1_i_xx = up_wing_units[kk + 1] - floor(xx / div_unit[kk + 1]);
  int gg_odd_k_1_i_cc_odd_n_k_1 = 0;
  int simp_odd_1_l = 0;
  int xdev_odd_k_1_base, xdev_odd_k_1, intVar, xdev_odd_k_1_proper;
    //Rcpp::Rcout << "gg_odd_k_1_i_xx" << ": " << gg_odd_k_1_i_xx << std::endl;
    //Rcpp::Rcout << "gg_odd_k_1_i_cc_odd_n_k_1" << ": " << gg_odd_k_1_i_cc_odd_n_k_1 << std::endl;
    //Rcpp::Rcout << "simp_odd_1_l" << ": " << simp_odd_1_l << std::endl;
  for ( int ii = 0; ii < xdev_l; ii++ ) {
    xdev_odd_k_1_base = gg_odd_k_1_i_xx + xdev_k[ii];
      //Rcpp::Rcout << "xdev_odd_k_1_base[" << ii << "]: " << xdev_odd_k_1_base << std::endl;
    intVar = xdev_odd_k_1_base >= *cc_odd_n_k_1;
    xdev_odd_k_1 = xdev_odd_k_1_base + cc_k_add * intVar;
    simp_odd_1[simp_odd_1_l + cc_k_add_ * intVar] = xdev_odd_k_1;
    xdev_odd_k_1_proper = (xdev_odd_k_1 >= (*cc_odd_n_k_1 + cc_k_add)) && (xdev_odd_k_1 < (wing_l[kk + 1] + cc_k_add)); // 0 -> *cc_odd_n_k_1 + cc_k_add from future_risk function
    simp_odd_1_l += xdev_odd_k_1_proper;
    gg_odd_k_1_i_cc_odd_n_k_1 += (1 - intVar) * xdev_odd_k_1_proper;
  }
  simp_odd_1_l += cc_k_add_;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # simp_odd_1_l: " << simp_odd_1_l << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # gg_odd_k_1_i_cc_odd_n_k_1: " << gg_odd_k_1_i_cc_odd_n_k_1 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # cc_odd_n_k_1: " << *cc_odd_n_k_1 << std::endl;
  for ( int ii = 0; ii < cc_k_add_; ii++ ) {
    simp_odd_1[gg_odd_k_1_i_cc_odd_n_k_1 + ii] = *cc_odd_n_k_1 + (cc_k_add - cc_k_add_) + ii;
  }
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # ... (xx_1_l + 1) / 2 = " << (xx_1_l + 1) / 2. << " ; " << wing_l[kk + 1] + cc_k_add << " = wing_l[" << kk << "] + cc_k_add." << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # simp_odd_1_l: " << simp_odd_1_l << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # gg_odd_k_1_i_cc_odd_n_k_1: " << gg_odd_k_1_i_cc_odd_n_k_1 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # cc_odd_n_k_1: " << *cc_odd_n_k_1 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # vsimp_odd_1.data(): " << vsimp_odd_1.data() << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # simp_odd_1*: " << simp_odd_1 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # vsimp_odd_1.size(): " << vsimp_odd_1.size() << std::endl;
    //for ( int ii = 0; ii < (xdev_l + cc_k_add_ + 1); ii++ ) {
    //  if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk # simp_odd_1[" << ii << "]: " << simp_odd_1[ii] << std::endl;
    //}

  // add half points for nodes //
  *simp_1 = *simp_odd_1 * 2;
    //Rcpp::Rcout << "simp_1[" << 0 << "]: " << simp_1[0] << std::endl;
  for ( int ii = 1; ii < simp_odd_1_l; ii++ ) {
    simp_1[2 * ii] = simp_odd_1[ii] * 2;
    simp_1[2 * ii - 1] = simp_odd_1[ii - 1] + simp_odd_1[ii];
      //Rcpp::Rcout << "simp_1[" << 2 * ii - 1 << "]: " << simp_1[2 * ii - 1] << std::endl;
      //Rcpp::Rcout << "simp_1[" << 2 * ii << "]: " << simp_1[2 * ii] << std::endl;
  }

  xi2 = 0;
    //Rcpp::Rcout << "xx: " << xx << "; simp_odd_1_l: " << simp_odd_1_l << std::endl;
    //Rcpp::Rcout << "simp_index[" << "]*: " << simp_index << std::endl;
  int simp_1_l = simp_odd_1_l * 2 - 1;
  for ( int ii = 0; ii < simp_1_l; ii++ ) {
    //int jj = 0; // for num_of_hypothesis;
    ll = simp_1[ii];
    dif = xx_1[ll] - xx;
    pr = R::dnorm(dif, effect_size * d_t, sq_d_t, 0);
    ww = (xx_1[simp_1[ii - 1 + (ii == 0)]] - xx_1[simp_1[ii + 1 - (ii == (simp_1_l - 1))]]) * (1 + (ii % 2)) / 3.;
    xi2 += ww * pr * value_1[ll];
      //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # xx_1[" << ll << "]: " << xx_1[ll] << "; dif: " << dif << "; pr: " << pr << "; ww: " << ww << "; val: " << value_1[ll] << std::endl;
      //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # xx_1[" << ll << "]: " << xx_1[ll] << "; xi2: " << xi2 << "; xi2+: " << ww * pr * value_1[ll] << std::endl;
  }

    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # xi2 = " << xi2 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # RETURN" << std::endl;
  return xi2;
}


// my function risk_balance
static double risk_balance(double xx, struct current_next* pinfo) {
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# risk_balance # START" << std::endl;
  struct base_time* str_base_time;
  str_base_time = &(pinfo->str_base_time);

  double f = future_risk(xx, pinfo);
  double c = current_risk(xx, str_base_time);

    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# risk_balance # current risk = " << c << "; future risk = " << f << std::endl;

    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# risk_balance # RETURN" << std::endl;
  return c - f;
}

// my function future_pr_rej_H0
// Remove the codes dealing with the critical value in the next stage (cc_n_k_1, cc_odd_n_k_1)
static double future_pr_rej_H0(double xx, struct current_next* pinfo) {
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_pr_rej_H0 # START" << std::endl;
  struct base_time* str_base_time;
  str_base_time = &(pinfo->str_base_time);
  struct base* str_base;
  str_base = &(str_base_time->str_base);

  double effect_size = str_base->effect_size;

  double* xx_1 = pinfo->gg_k_1;
  //int xx_1_l = pinfo->gg_k_1_l;
  double* value_1 = pinfo->value_1;
  //double* doubleVar = pinfo->dummy;

  //for( int ii = 0; ii < *xx_1_l; ii++ ) { 
  //  Rcpp::Rcout << "gg[" << ii << "]; " << xx_1[ii] << std::endl;
  //}

  int ll;
  double dif, pr, ww, xi2;
  double d_t = pinfo->time_1 - str_base_time->time;
  double sq_d_t = sqrt(d_t);

  int xdev_l = str_base->xdev_l;
  double* div_unit = str_base->div_unit;
  int* up_wing_units = str_base->up_wing_units;
  int* wing_l = str_base->wing_l;

  int kk = str_base_time->stage;
  int* xdev_k = pinfo->xdev_k;
  //double* cc_k_1 = pinfo->cc_k_1;
  //int* cc_odd_n_k_1 = pinfo->cc_odd_n_k_1;

  std::vector<int> vsimp_odd_1(xdev_l + 1);  // 1 is added for overrun in 'loop for simp_odd_1'
  vsimp_odd_1.resize(xdev_l);
  int* simp_odd_1 = vsimp_odd_1.data();
  std::vector<int> vsimp_1(2 * xdev_l - 1);
  int* simp_1 = vsimp_1.data();

  // nodes for simpson's method //
  int gg_odd_k_1_i_xx = up_wing_units[kk + 1] - floor(xx / div_unit[kk + 1]);
  int simp_odd_1_l = 0;
  int xdev_odd_k_1_base, xdev_odd_k_1, xdev_odd_k_1_proper;
    //Rcpp::Rcout << "gg_odd_k_1_i_xx" << ": " << gg_odd_k_1_i_xx << std::endl;
    //Rcpp::Rcout << "simp_odd_1_l" << ": " << simp_odd_1_l << std::endl;
  for ( int ii = 0; ii < xdev_l; ii++ ) {
    xdev_odd_k_1_base = gg_odd_k_1_i_xx + xdev_k[ii];
      //Rcpp::Rcout << "xdev_odd_k_1_base[" << ii << "]: " << xdev_odd_k_1_base << std::endl;
    xdev_odd_k_1 = xdev_odd_k_1_base;
    simp_odd_1[simp_odd_1_l] = xdev_odd_k_1;
      //Rcpp::Rcout << "# future_pr_rej_H0 # *cc_odd_n_k_1 + 1: " << *cc_odd_n_k_1 + 1 << std::endl;
    //xdev_odd_k_1_proper = (xdev_odd_k_1 >= (*cc_odd_n_k_1 + 1)) && (xdev_odd_k_1 < wing_l[kk + 1]); // *cc_odd_n_k_1 + 2 -> *cc_odd_n_k_1 + 1 from future_risk0 function
    xdev_odd_k_1_proper = (xdev_odd_k_1 >= 0) && (xdev_odd_k_1 < wing_l[kk + 1]); // *cc_odd_n_k_1 + 2 -> *cc_odd_n_k_1 + 1 from future_risk0 function
    simp_odd_1_l += xdev_odd_k_1_proper;
  }
  //simp_odd_1_l -= 1;  // += -> -= (or += 1 - 2) from future_risk0 function
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_pr_rej_H0 # ... (xx_1_l + 1) / 2 = " << (xx_1_l + 1) / 2. << " ; " << wing_l[kk + 1] << " = wing_l[" << kk << "]." << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_pr_rej_H0 # simp_odd_1_l: " << simp_odd_1_l << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_pr_rej_H0 # cc_odd_n_k_1: " << *cc_odd_n_k_1 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_pr_rej_H0 # vsimp_odd_1.data(): " << vsimp_odd_1.data() << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_pr_rej_H0 # simp_odd_1*: " << simp_odd_1 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_pr_rej_H0 # vsimp_odd_1.size(): " << vsimp_odd_1.size() << std::endl;
    //for ( int ii = 0; ii < (xdev_l + 1); ii++ ) {
    //  if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_pr_rej_H0 # simp_odd_1[" << ii << "]: " << simp_odd_1[ii] << std::endl;
    //}

  // add half points for nodes //
  *simp_1 = *simp_odd_1 * 2;
    //Rcpp::Rcout << "simp_1[" << 0 << "]: " << simp_1[0] << std::endl;
  for ( int ii = 1; ii < simp_odd_1_l; ii++ ) {
    simp_1[2 * ii] = simp_odd_1[ii] * 2;
    simp_1[2 * ii - 1] = simp_odd_1[ii - 1] + simp_odd_1[ii];
      //Rcpp::Rcout << "simp_1[" << 2 * ii - 1 << "]: " << simp_1[2 * ii - 1] << std::endl;
      //Rcpp::Rcout << "simp_1[" << 2 * ii << "]: " << simp_1[2 * ii] << std::endl;
  }

  xi2 = 0;
    //Rcpp::Rcout << "xx: " << xx << "; simp_odd_1_l: " << simp_odd_1_l << std::endl;
    //Rcpp::Rcout << "simp_index[" << "]*: " << simp_index << std::endl;
  int simp_1_l = simp_odd_1_l * 2 - 1;
  for ( int ii = 0; ii < simp_1_l; ii++ ) {
    //int jj = 0; // for num_of_hypothesis;
    ll = simp_1[ii];
    dif = xx_1[ll] - xx;
    pr = R::dnorm(dif, effect_size * d_t, sq_d_t, 0);
    ww = (xx_1[simp_1[ii - 1 + (ii == 0)]] - xx_1[simp_1[ii + 1 - (ii == (simp_1_l - 1))]]) * (1 + (ii % 2)) / 3.;
    xi2 += ww * pr * value_1[ll];
      //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # xx_1[" << ll << "]: " << xx_1[ll] << "; dif: " << dif << "; pr: " << pr << "; ww: " << ww << "; val: " << value_1[ll] << std::endl;
      //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # xx_1[" << ll << "]: " << xx_1[ll] << "; xi2: " << xi2 << "; xi2+: " << ww * pr * value_1[ll] << std::endl;
  }

    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # xi2 = " << xi2 << std::endl;
    //if ( glb_sgn == 1 ) Rcpp::Rcout << "# future_risk0 # RETURN" << std::endl;
  return xi2;
}


static double construct(double cost0, struct ground* ground) {
  struct base* str_base = &(ground->str_base);
  struct result* str_result = &(ground->str_result);

  ////Rcpp::Rcout << "# work_test_norm_c # # construct # START" << std::endl;
  ////Rcpp::Rcout << "# work_test_norm_c # # construct # cost for (rej H0 | H0): " << cost0 << std::endl;

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Recovery: Settings                                                        //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  int analysis = str_base->analysis;
  //double max_U0 = str_base->max_time;
  const int work_KK = str_base->work_KK;
  double* U_k = str_base->U_k;
  //const int prior_l = str_base->num_of_hypotheses;
  //double* H_th = str_base->hypotheses;
  //double* prior = str_base->prior;
  //double fss_weight = str_base->fss_w;
  //double seq_weight = str_base->seq_w;
  //double cost0 = str_base->cost0;
  double effect_size = str_base->effect_size;
  //double cost1 = str_base->cost1;
  double stat = str_base->stat;

  int xdev_l = str_base->xdev_l;
  double* xdev = str_base->xdev;
  double* up_wing_temp = str_base->up_wing_temp;
  double* up_wing_buffer = str_base->up_wing_buffer;
  double* div_unit = str_base->div_unit;
  //int* up_wing_units = str_base->up_wing_units;
  int* wing_l = str_base->wing_l;

  double tol_boundary = str_base->tol_boundary;

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Recovery: Memories for results                                            //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  int cc_k_add = str_result->cc_k_add;
  std::vector<double>* vss0 = str_result->vss0;
  int* gg_odd_l = str_result->gg_odd_l;
  std::vector<double>* vgg_odd = str_result->vgg_odd;
  int* gg_l = str_result->gg_l;
  std::vector<double>* vgg = str_result->vgg;
  std::vector<double>* veta = str_result->veta;
  std::vector<double>* vpr_rej_H0 = str_result->vpr_rej_H0;
  std::vector<double>* vexp_time = str_result->vexp_time;
  int* cc_odd_n = str_result->cc_odd_n;
  int* cc_n = str_result->cc_n;
  double* cc = str_result->cc;
  double* val_k = str_result->val_k;

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Main                                                                      //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Declarations (dummy box) ===//
  double doubleVar[8] = {};
  int intVar[8] = {};
  double doubleInf = std::numeric_limits<double>::infinity();

  //=== Declarations (status of the (k)th stage) ===//
  struct base_time str_base_time;
  struct current_next str_current_next;

  //=== Declarations (pointers for the (k)th stage) ===//
  int kk;
  double* ss_k;
  int gg_odd_k_l;
  double* gg_odd_k;
  int gg_k_l;
  double* gg_k;
  double* eta_k;
  double* pr_rej_H0_k;
  double* exp_time_k;
  int* cc_odd_n_k;
  int* cc_n_k;
  double* cc_k;
  double t_k;

  //=== Declarations (pointers for the (k + 1)th stage) ===//
  double t_k_1;
  double d_t;
  double sq_d_t;
  int* cc_odd_n_k_1;
  //int* cc_n_k_1;
  double* cc_k_1;
  std::vector<int> vxdev_k(xdev_l);
  int* xdev_k = vxdev_k.data();


  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // The Final Stage                                                           //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Initialize (variables for the stage) ===//
  kk = work_KK;
    //Rcpp::Rcout << "# work_test_norm_c # # construct # ..... Stage: " << kk << std::endl;

  ss_k = vss0[kk].data();
  gg_odd_k_l = gg_odd_l[kk];
  gg_odd_k = vgg_odd[kk].data();
  gg_k_l = gg_l[kk];
  gg_k = vgg[kk].data();
  eta_k = veta[kk].data();
  pr_rej_H0_k = vpr_rej_H0[kk].data();
  exp_time_k = vexp_time[kk].data();
  cc_odd_n_k = &(cc_odd_n[kk]);
  cc_n_k = &(cc_n[kk]);
  cc_k = &(cc[kk]);
  t_k = U_k[kk];
  // structure //
  (*str_base).cost0 = cost0;
  str_base_time.str_base = *str_base;
  str_base_time.stage = kk;
  str_base_time.time = t_k;

  //double res1 = obj_func(0, &str_base_time);
  //Rcpp::Rcout << res1 << std::endl;

 // *cc_k = Brent_fmin(ss_k[wing_l[kk] - 1], ss_k[0], (double (*)(double, void*)) obj_func, &str_base_time, 1e-4);
  //  Rcpp::Rcout << "# work_test_norm_c # # construct # upper boundary: " << *cc_k << std::endl;

  *cc_k = 
    bisection_inverse(current_risk_balance,
      0, &str_base_time, ss_k[wing_l[kk] - 1], ss_k[0],
      false, true, false, tol_boundary);
    //Rcpp::Rcout << "# work_test_norm_c # # construct # upper boundary: " << *cc_k << std::endl;

  //for( int ii = 0; ii < 100; ii++ ) {
  //  *doubleVar = current_risk_balance(9.10 + ii / 100., &str_base_time);
  //  Rcpp::Rcout << 9.10 + ii << ": " << *doubleVar << std::endl;
  //}
  //doubleVar[0] = 0;
  //doubleVar[1] = 0;
  //for( int ii = 0; ii < 100; ii++ ) {
  //  doubleVar[0] = obj_func(9.10 + ii / 100., &str_base_time);
  //  doubleVar[1] = doubleVar[0] - doubleVar[1];
  //  Rcpp::Rcout << 9.10 + ii << ": " << doubleVar[1] << std::endl;
  //  doubleVar[1] = doubleVar[0];
  //}

  //... Routine 2 ...//
  // gg_odd_k //
  *intVar = 0;
  *cc_odd_n_k = 0;
  for ( int ii = 0; ii < (gg_odd_k_l - cc_k_add); ii++ ) {
    *intVar = *cc_k > ss_k[ii];
    *cc_odd_n_k += *intVar;
    gg_odd_k[ii + *intVar * cc_k_add] = ss_k[ii];
  }
  *cc_odd_n_k = gg_odd_k_l - cc_k_add - *cc_odd_n_k;
    //Rcpp::Rcout << "# work_test_norm_c # # construct # upper boundary in gg_odd: " << "cc_odd_n_k = " << *cc_odd_n_k << std::endl;
  for ( int ii = 0; ii < cc_k_add; ii++ ) {
    gg_odd_k[*cc_odd_n_k + ii] = *cc_k;
  }
  *cc_n_k = *cc_odd_n_k * 2 + 1;
    //Rcpp::Rcout << "# work_test_norm_c # # construct # upper boundary in gg: " << "cc_n_k = " << *cc_n_k << std::endl;
    //for( int ii = 0; ii < gg_odd_k_l; ii++ ) {
    //  Rcpp::Rcout << "ii = " << ii << ": gg_odd = " << gg_odd_k[ii] << std::endl;
    //}
  // gg_k //
  *gg_k = *gg_odd_k;
  for ( int ii = 1; ii < gg_odd_k_l; ii++ ) {
    gg_k[2 * ii] = gg_odd_k[ii];
    gg_k[2 * ii - 1] = (gg_odd_k[ii - 1] + gg_odd_k[ii]) / 2.;
  }
    //for( int ii = 0; ii < gg_k_l; ii++ ) {
    //  Rcpp::Rcout << "ii = " << ii << ": gg = " << gg_k[ii] << std::endl;
    //}
  //... Routine 2 end ...//

  for ( int ii = 0; ii < gg_k_l; ii++ ) {
    eta_k[ii] = current_risk(gg_k[ii], &str_base_time);
    pr_rej_H0_k[ii] = (ii <= *cc_n_k);
    //pr_acc_H0_k[ii] = !pr_rej_H0_k[ii];
    //Rcpp::Rcout << "ii = " << ii << ": gg = " << gg_k[ii] << ": pr_rej_H0 = " << pr_rej_H0_k[ii] << std::endl;
    exp_time_k[ii] = 0;
  }

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // k th Stage                                                                //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  kk = kk - 1;
  for ( ; kk > 0; kk--) {
  //for( int kk = work_KK - 1; kk > (work_KK - 3); kk--) {
    //=== Initialize (variables for the stage) ===//
      //Rcpp::Rcout << "# work_test_norm_c # # construct # ..... Stage: " << kk << std::endl;

    ss_k = vss0[kk].data();
    gg_odd_k_l = gg_odd_l[kk];
    gg_odd_k = vgg_odd[kk].data();
    gg_k_l = gg_l[kk];
    gg_k = vgg[kk].data();
    eta_k = veta[kk].data();
    pr_rej_H0_k = vpr_rej_H0[kk].data();
    //pr_acc_H0_k = vpr_acc_H0[kk].data();
    exp_time_k = vexp_time[kk].data();
    cc_odd_n_k = &(cc_odd_n[kk]);
    cc_n_k = &(cc_n[kk]);
    cc_k = &(cc[kk]);
    t_k = U_k[kk];
    str_base_time.stage = kk;
    str_base_time.time = t_k;

    // *cc_k = Brent_fmin(ss_k[wing_l[kk] - 1], ss_k[0], (double (*)(double, void*)) obj_func, &str_base_time, 1e-5);
      //Rcpp::Rcout << *cc_k << std::endl;

    // Lower bound of the stopping boundary //
    *cc_k = 
      bisection_inverse(current_risk_balance,
        0, &str_base_time, ss_k[wing_l[kk] - 1], ss_k[0],
        false, true, false, tol_boundary);
      //Rcpp::Rcout << "# work_test_norm_c # # construct # lower bound for upper boundary: " << *cc_k << std::endl;

    cc_odd_n_k_1 = &(cc_odd_n[kk + 1]);
    //cc_n_k_1 = &(cc_n[kk + 1]);
    cc_k_1 = &(cc[kk + 1]);
    str_current_next.str_base_time = str_base_time;
    str_current_next.time_1 = U_k[kk + 1];
    str_current_next.gg_k_1 = vgg[kk + 1].data();
    str_current_next.gg_k_1_l = gg_l[kk + 1];
    str_current_next.value_1 = veta[kk + 1].data();
    str_current_next.dummy = val_k;
    str_current_next.cc_k_1 = cc_k_1;
    str_current_next.cc_odd_n_k_1 = cc_odd_n_k_1;
    str_current_next.xdev_k = xdev_k;

    t_k_1 = U_k[kk + 1];
    d_t = t_k_1 - t_k;
    sq_d_t = sqrt(d_t);
    for ( int ii = 0; ii < xdev_l; ii++ ) {
      *intVar = round(-(xdev[ii] * sqrt(d_t) + effect_size * d_t) / div_unit[kk + 1]);
      xdev_k[ii] = *intVar;
      //Rcpp::Rcout << "xdev_k[" << ii << "]: " << *intVar << std::endl;
      //xdev_k[ii] = round(-xdev[ii] * sqrt(d_t) / div_unit[kk + 1]);
    }
      //Rcpp::Rcout << "cc_odd_n_k_1" << *cc_odd_n_k_1 << std::endl;

      //doubleVar[0] = *cc_k; // cc_lw
      //doubleVar[1] = up_wing_temp[kk] + up_wing_buffer[kk]; // cc_up
      //Rcpp::Rcout << "cc_lw: " << doubleVar[0] << std::endl;
      //Rcpp::Rcout << "cc_up: " << doubleVar[1] << std::endl;
      //Rcpp::Rcout << "cc_up( " << up_wing_temp[kk] << "; " << up_wing_buffer[kk] << " )" << std::endl;
      //for( int ii = 0; ii < wing_l[kk]; ii++ ) {
      //  val_k[ii] = risk_balance(ss_k[ii], &str_current_next);
      //  Rcpp::Rcout << "ss[" << ii << "]: " << ss_k[ii] << "; " << current_risk(ss_k[ii], &str_base_time) << "; " << future_risk(ss_k[ii], &str_current_next) << "; " << risk_balance(ss_k[ii], &str_current_next) << std::endl;
      //}

      //Rcpp::Rcout << "xl: " << doubleVar[0] << "; f(xl): " << risk_balance(doubleVar[0], &str_current_next)
      //  << "; cur risk: " << current_risk(doubleVar[0], &str_base_time)
      //  << "; fut risk: " << future_risk(doubleVar[0], &str_current_next) << std::endl;
      //Rcpp::Rcout << "xu: " << doubleVar[1] << "; f(xu): " << risk_balance(doubleVar[1], &str_current_next)
      //  << "; cur risk: " << current_risk(doubleVar[1], &str_base_time)
      //  << "; fut risk: " << future_risk(doubleVar[1], &str_current_next) << std::endl;

    // Stopping boundary //
    *cc_k = 
      bisection_inverse(risk_balance,
        0, &str_current_next, *cc_k, up_wing_temp[kk] + up_wing_buffer[kk],
        false, true, false, tol_boundary);
      //Rcpp::Rcout << "# work_test_norm_c # # construct # upper boundary: " << *cc_k << std::endl;

    //... Routine 2 ...//
    // gg_odd_k //
    *intVar = 0;
    *cc_odd_n_k = 0;

    for ( int ii = 0; ii < (gg_odd_k_l - cc_k_add); ii++ ) {
      *intVar = *cc_k > ss_k[ii];
      *cc_odd_n_k += *intVar;
      gg_odd_k[ii + *intVar * cc_k_add] = ss_k[ii];
    }
    *cc_odd_n_k = gg_odd_k_l - cc_k_add - *cc_odd_n_k;
    for ( int ii = 0; ii < cc_k_add; ii++ ) {
      gg_odd_k[*cc_odd_n_k + ii] = *cc_k;
    }
    *cc_n_k = *cc_odd_n_k * 2 + 1;
    // gg_k //
    *gg_k = *gg_odd_k;
    for ( int ii = 1; ii < gg_odd_k_l; ii++ ) {
      gg_k[2 * ii] = gg_odd_k[ii];
      gg_k[2 * ii - 1] = (gg_odd_k[ii - 1] + gg_odd_k[ii]) / 2.;
    }
    //... Routine 2 end ...//
     // Rcpp::Rcout << (gg_odd_k_l - cc_k_add) << std::endl;

    // rej H0 //
      //Rcpp::Rcout << "cc_n_k: " << *cc_n_k << std::endl;
    for ( int ii = 0; ii < (*cc_n_k + 1); ii++ ) {
      eta_k[ii] = current_risk(gg_k[ii], &str_base_time);
      pr_rej_H0_k[ii] = 1;
      //pr_acc_H0_k[ii] = 0;
      exp_time_k[ii] = 0;
      //Rcpp::Rcout << "gg[" << ii << "]: " << gg_k[ii] << "; Pr(rej H0): " << pr_rej_H0_k[ii] << "; " << *doubleVar <<  std::endl;
    }
    // continue //
    for ( int ii = (*cc_n_k + 1); ii < gg_k_l; ii++ ) {
      eta_k[ii] = future_risk(gg_k[ii], &str_current_next);
    }
    //str_current_next.gg_k_1 = vgg[kk + 1].data() + (*cc_n_k_1 + 1); // for old future_risk0
    //str_current_next.gg_k_1_l = gg_k_l - (*cc_n_k_1 + 1); // for old future_risk0
    //str_current_next.value_1 = vpr_rej_H0[kk + 1].data() + (*cc_n_k_1 + 1); // for old future_risk0
    str_current_next.value_1 = vpr_rej_H0[kk + 1].data();
      //Rcpp::Rcout << "xl: " << doubleVar[0] << "; Pr(rej_H0 | xl): " << future_risk0(doubleVar[0], &str_current_next) << std::endl;
      //Rcpp::Rcout << "xu: " << doubleVar[1] << "; Pr(rej_H0 | xu): " << future_risk0(doubleVar[1], &str_current_next) << std::endl;
      //Rcpp::Rcout << "cc_n_k: " << *cc_n_k + 1 << std::endl;
    for ( int ii = (*cc_n_k + 1); ii < gg_k_l; ii++ ) {
      pr_rej_H0_k[ii] = future_risk0(gg_k[ii], &str_current_next)
                        + R::pnorm(-(*cc_k_1 - gg_k[ii]), - effect_size * d_t, sq_d_t, 1, 0);
      //pr_acc_H0_k[ii] = 1 - pr_rej_H0_k[ii];
      //exp_time_k[ii] = 0;
      //Rcpp::Rcout << "gg[" << ii << "]: " << gg_k[ii] << "; Pr(rej H0): " << pr_rej_H0_k[ii] <<  std::endl;
    }
  }

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // 1 st Stage                                                                //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Initialize (variables for the stage) ===//
  //kk = 0;
  cc_k_add = 0;
    //Rcpp::Rcout << "# work_test_norm_c # # construct # ..... Stage: " << kk << std::endl;

  ss_k = vss0[kk].data();
  gg_odd_k_l = gg_odd_l[kk];
  gg_odd_k = vgg_odd[kk].data();
  gg_k_l = gg_l[kk];
  gg_k = vgg[kk].data();
  eta_k = veta[kk].data();
  pr_rej_H0_k = vpr_rej_H0[kk].data();
  //pr_acc_H0_k = vpr_acc_H0[kk].data();
  exp_time_k = vexp_time[kk].data();
  cc_odd_n_k = &(cc_odd_n[kk]);
  cc_n_k = &(cc_n[kk]);
  cc_k = &(cc[kk]);
  t_k = U_k[kk];
  str_base_time.stage = kk;
  str_base_time.time = t_k;

  // *cc_k = Brent_fmin(ss_k[wing_l[kk] - 1], ss_k[0], (double (*)(double, void*)) obj_func, &str_base_time, 1e-5);
    //Rcpp::Rcout << *cc_k << std::endl;

  // Lower bound of the stopping boundary //
  //Rcpp::Rcout << ss_k[wing_l[kk] - 1] << "; " << ss_k[0] << std::endl;
  *cc_k = 
    bisection_inverse(current_risk_balance,
      0, &str_base_time, ss_k[wing_l[kk] - 1], ss_k[0],
      false, true, false, tol_boundary);
  //*cc_k = 
  //  bisection_inverse(current_risk_balance,
  //    0, &str_base_time, vss0[kk + 1].at(wing_l[kk + 1] - 1), vss0[kk + 1].at(0),
  //    false, true, false, tol_boundary);
    //Rcpp::Rcout << "# work_test_norm_c # # construct # lower bound for upper boundary: " << *cc_k << std::endl;

  cc_odd_n_k_1 = &(cc_odd_n[kk + 1]);
  //cc_n_k_1 = &(cc_n[kk + 1]);
  cc_k_1 = &(cc[kk + 1]);
  str_current_next.str_base_time = str_base_time;
  str_current_next.time_1 = U_k[kk + 1];
  str_current_next.gg_k_1 = vgg[kk + 1].data();
  str_current_next.gg_k_1_l = gg_l[kk + 1];
  str_current_next.value_1 = veta[kk + 1].data();
  str_current_next.dummy = val_k;
  str_current_next.cc_k_1 = cc_k_1;
  str_current_next.cc_odd_n_k_1 = cc_odd_n_k_1;
  str_current_next.xdev_k = xdev_k;

  t_k_1 = U_k[kk + 1];
  d_t = t_k_1 - t_k;
  sq_d_t = sqrt(d_t);
  for ( int ii = 0; ii < xdev_l; ii++ ) {
    *intVar = round(-(xdev[ii] * sqrt(d_t) + effect_size * d_t) / div_unit[kk + 1]);
    xdev_k[ii] = *intVar;
    //Rcpp::Rcout << "xdev_k[" << ii << "]: " << *intVar << std::endl;
    //xdev_k[ii] = round(-xdev[ii] * sqrt(d_t) / div_unit[kk + 1]);
  }
    //Rcpp::Rcout << "cc_odd_n_k_1" << *cc_odd_n_k_1 << std::endl;

    //doubleVar[0] = *cc_k; // cc_lw
    //doubleVar[1] = up_wing_temp[kk] + up_wing_buffer[kk]; // cc_up
    //Rcpp::Rcout << "cc_lw: " << doubleVar[0] << std::endl;
    //Rcpp::Rcout << "cc_up: " << doubleVar[1] << std::endl;
    //Rcpp::Rcout << "cc_up( " << up_wing_temp[kk] << "; " << up_wing_buffer[kk] << " )" << std::endl;
    //for( int ii = 0; ii < wing_l[kk]; ii++ ) {
    //  val_k[ii] = risk_balance(ss_k[ii], &str_current_next);
    //  Rcpp::Rcout << "ss[" << ii << "]: " << ss_k[ii] << "; " << current_risk(ss_k[ii], &str_base_time) << "; " << future_risk(ss_k[ii], &str_current_next) << "; " << risk_balance(ss_k[ii], &str_current_next) << std::endl;
    //}

    //Rcpp::Rcout << "xl: " << doubleVar[0] << "; f(xl): " << risk_balance(doubleVar[0], &str_current_next)
    //  << "; cur risk: " << current_risk(doubleVar[0], &str_base_time)
    //  << "; fut risk: " << future_risk(doubleVar[0], &str_current_next) << std::endl;
    //Rcpp::Rcout << "xu: " << doubleVar[1] << "; f(xu): " << risk_balance(doubleVar[1], &str_current_next)
    //  << "; cur risk: " << current_risk(doubleVar[1], &str_base_time)
    //  << "; fut risk: " << future_risk(doubleVar[1], &str_current_next) << std::endl;

  // Stopping boundary (old cc_k_add will be used) //
  *doubleVar = doubleInf;  // for t_k == 0;
  if ( t_k > 0 ) {
    *doubleVar = 
      bisection_inverse(risk_balance,
        0, &str_current_next, *cc_k, up_wing_temp[kk] + up_wing_buffer[kk],
        false, true, false, tol_boundary);
      //Rcpp::Rcout << "# work_test_norm_c # # construct # upper boundary: " << *doubleVar << std::endl;
  }
  // *cc_k = 
  //  bisection_inverse(risk_balance,
  //    0, &str_current_next, doubleVar[0], doubleVar[1],
  //    false, true, false, tol_boundary);
  *cc_k = -1;
    //Rcpp::Rcout << *cc_k << std::endl;

  //... Routine 2 ...//
  // gg_odd_k //
  *cc_odd_n_k = 0;
  *gg_odd_k = ss_k[1] + stat;
  // gg_k //
  *cc_n_k = 0;
  *gg_k = *gg_odd_k;
  //... Routine 2 end ...//

  // continue //
  *eta_k = future_risk(*gg_k, &str_current_next);
  str_current_next.value_1 = vpr_rej_H0[kk + 1].data();
  //double* gg_k_1 = vgg[kk + 1].data();
  //double* pr_rej_H0_kk_1 = vpr_rej_H0[kk + 1].data();
    //for( int ii = 0; ii < gg_l[kk + 1]; ii++ ) {
    //  Rcpp::Rcout << "gg_k_1[" << ii << "]: " << gg_k_1[ii] << "; rej0: " << pr_rej_H0_kk_1[ii] << std::endl;
    //}
  *pr_rej_H0_k = future_risk0(*gg_k, &str_current_next)
                 + R::pnorm(-(*cc_k_1 - *gg_k), - effect_size * d_t, sq_d_t, 1, 0);

  // for output //
  doubleVar[1] = std::numeric_limits<double>::infinity();
  *cc_k = doubleVar[1 - analysis];
  *pr_rej_H0_k = (*gg_k >= *cc_k) + (*gg_k < *cc_k) * *pr_rej_H0_k;
  //*pr_rej_H0_k += (stat >= *doubleVar)
  //Rcpp::Rcout << "# work_test_norm_c # # construct # Type I Err Pr: " << *pr_rej_H0_k << std::endl;
  ////Rcpp::Rcout << "# work_test_norm_c # # construct # END" << std::endl;

  return *pr_rej_H0_k;

}


// [[Rcpp::export]]
Rcpp::List work_test_norm_c(
const double overall_sig_level = 0.025,
const double work_beta = 0.05,
const double cond_alpha = 0.025,
const double cost_type_1_err = 0,
const double cost_type_2_err = 0,
const double prev_cost = 0,
const double min_effect_size = 1,
const double effect_size = 0,
const int basic_schedule_num = 50,
const int basic_schedule_power = 2,
Rcpp::NumericVector basic_schedule = Rcpp::NumericVector::create(0),
Rcpp::NumericVector prior_dist = Rcpp::NumericVector::create(0),
const double prev_time = 0,
const double time = 0,
const double next_time = 0,
const double stat = 0,
const bool input_check = true,
const bool out_process = false,
const int simpson_div = 6,
const double tol_boundary = 1e-8,
const double tol_cost = 1e-8) {
//, const int out = 0) {
  //glb_sgn = out;

  if ( input_check ) {
    int intVar[2];
    if ( (overall_sig_level < 0) || (overall_sig_level > 1) ) {
      Rcpp::stop("'overall_sig_level' should be a value in (0, 1).");
    }
    if ( (cond_alpha < 0) || (cond_alpha > 1) ) {
      Rcpp::stop("'cond_alpha' should be a value in [0, 1].");
    }
    if ( cost_type_1_err < 0 ) {
      Rcpp::stop("'cost_type_1_err' should be non-negative. When 0 is specified, its value will be calculated. For details, see help file.");
    }
    if ( cost_type_2_err < 0 ) {
      Rcpp::stop("'cost_type_2_err' should be non-negative. When 0 is specified, its value will be calculated. For details, see help file.");
    }
    if ( ((prev_time == 0) && (next_time > 0)) && ((cond_alpha == 0) || (cond_alpha == 1)) ) {
      Rcpp::stop("The valud of cond_alpha indicates that the trial has been completed.");
    }
    *intVar = 0;
    for ( int kk = 0; kk < basic_schedule.size(); kk++ ) {
      *intVar = *intVar + (basic_schedule.at(kk) == 0);
    }
    if ( *intVar > 1 ) {
      Rcpp::stop("'basic_schedule' is allowed to have only one '0' or no '0' among its elements.");
    }
    *intVar = (*intVar == 1);
    if ( (  ((basic_schedule.size() - *intVar) < 2) && (basic_schedule.at(0) == 0) && (basic_schedule_num > 100)) ||
         (!(((basic_schedule.size() - *intVar) < 2) && (basic_schedule.at(0) == 0)) && ((basic_schedule.size() - *intVar) > 100)) ) {
      Rcpp::stop("The length of 'basic_schedule' should be less than 100.");
    }
    if ( ((basic_schedule.size() - *intVar) < 2) && (basic_schedule.at(0) == 0) && (work_beta < 0) && (work_beta > 1) ) {
      Rcpp::stop("When 'basic_schedule' is unspecified, 'work_beta' should be a value between 0 and 1.");
    }
    if ( ((basic_schedule.size() - *intVar) < 2) && (basic_schedule.at(0) == 0) && (basic_schedule_num < 1) ) {
      Rcpp::stop("When 'basic_schedule' is unspecified (basic_schedule == c(0)), 'basic_schedule_num' should be equal to or greater than 1.");
    }
    intVar[1] = 0;
    for ( int kk = 0; kk < basic_schedule.size(); kk++ ) {
      intVar[1] = intVar[1] + (basic_schedule.at(kk) < 0);
    }
    if ( ((basic_schedule.size() - *intVar) > 0) && (intVar[1] > 0) ) {
      Rcpp::stop("All elements of 'basic_schedule' should be non-negative.");
    }
    intVar[1] = 0;
    for ( int kk = 1; kk < basic_schedule.size(); kk++ ) {
      intVar[1] = intVar[1] + (basic_schedule.at(kk - 1) >= basic_schedule.at(kk));
    }
    if ( ((basic_schedule.size() - *intVar) > 1) && (intVar[1] > 0) ) {
      Rcpp::stop("'basic_schedule' should be a strictly increasing sequence.");
    }
    if ( (!((prior_dist.size() == 1) && (prior_dist.at(0) == 0))) && (prior_dist.size() != 11)  ) {
      Rcpp::stop("When designating 'prior_dist', its length should be 11.");
    }
    intVar[1] = 0;
    for ( int kk = 0; kk < prior_dist.size(); kk++ ) {
      intVar[1] = intVar[1] + (prior_dist.at(kk) <= 0);
    }
    if ( (prior_dist.size() != 1) && (intVar[1] > 0)  ) {
      Rcpp::stop("All elements of 'prior_dist' should be positive.");
    }
    if ( prev_time < 0 ) {
      Rcpp::stop("'prev_time' should be non-negative.");
    }
    if ( time < 0 ) {
      Rcpp::stop("'time' should be non-negative.");
    }
    if ( next_time < 0 ) {
      Rcpp::stop("'next_time' should be non-negative.");
    }
    if ( (prev_time * next_time) > 0 ) {
      Rcpp::stop("Either 'prev_time' or 'next_time' should be 0.");
    }
    if ( ((prev_time == 0) && (next_time > 0)) && (time >= next_time) ) { // (update mode) && (invalid order) : error
      Rcpp::stop("To update the working test, 'time' < 'next_time' is necessary .");
    }
    if ( (next_time == 0) && (time > 0) && (prev_time >= time) ) { // (analysis mode) && (invalid order) : error
      Rcpp::stop("To conduct analysis (to calculate the conditinal Type I error probability, 'prev_time' < 'time' is necessary .");
    }
    if ( (prev_time == 0) && (time == 0) && (next_time >= 0) && (stat != 0) ) {
      Rcpp::stop("To construct the first working test at 'time' == 0, 'stat' should be 0.");
    }
    if ( ((prev_time == 0) && ((next_time > 0) || ((next_time + time) == 0))) && (effect_size != 0) ) {
      Rcpp::stop("To construct a working test , 'effect_size' should be 0.");
    }
    if ( (tol_boundary <= 0) ) {
      Rcpp::stop("'tol_boundary' should be positive.");
    }
    if ( (tol_cost <= 0) ) {
      Rcpp::stop("'tol_cost' should be positive.");
    }
  }


  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Settings                                                                  //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//

  //=== Dummy Box ===//
  double doubleVar[8] = {};
  int intVar[8] = {};

  //=== Initialize Inputs ===//
  // for iterative search of cost0 satisfying overall_sig_level //
  double cost_upper_lim;
  double theta = effect_size;
  if ( cost_type_1_err == 0 ) {
    theta = 0;
  }

  //Rcpp::Rcout << "# optim_design# " << "#======== Sig Lev : " << alpha << "; # of analyses : " << MM << "========#" << std::endl;

  //=== DESIGN DESIGNATION ===//
  //--- Basic Analysis Schedule ---//
  // Maximum of basic analysis schedule //
  double za = -R::qnorm(overall_sig_level, 0, 1, 1, 0);
  double work_zb = -R::qnorm(work_beta, 0, 1, 1, 0);
  double max_U0 = 0;
  // Number of analyses of the working test //
  *intVar = basic_schedule.size();
  intVar[1] = (basic_schedule.at(0) == 0);
  int work_KK = 0;
  if ( (*intVar - intVar[1]) > 0 ) {
    work_KK = *intVar - intVar[1];
  } else if ( basic_schedule_num > 0 ) {
    work_KK = basic_schedule_num;
  }
  int work_KK_1 = work_KK + 1;
    //Rcpp::Rcout << *intVar << "; " << intVar[1] << "; " << std::endl;
    //Rcpp::Rcout << "work_KK = " << work_KK << std::endl;
  std::vector<double> vU0(work_KK_1);
  vU0.resize(work_KK_1);
  double* U0 = vU0.data();
  // Basic analysis schedule //
  //U0[0] = 0;
  vU0.at(0) = 0; //..
  if ( (*intVar - intVar[1]) > 0 ) { // Use user-specified basic analysis schedule
    max_U0 = basic_schedule.at(basic_schedule.size() - 1);
    for ( int kk = 0; kk < *intVar; kk++ ) {
      //U0[kk + 1 - intVar[1]] = basic_schedule.at(kk);
      vU0.at(kk + 1 - intVar[1]) = basic_schedule.at(kk); //..
    }
  } else if ( basic_schedule_num > 0 ) {
    max_U0 = pow(- za - work_zb, 2.) / pow(min_effect_size, 2.);
    for ( int kk = 0; kk < work_KK_1; kk++ ) { // Use program-generating basic analysis schedule
      //U0[kk] = pow((double) kk / (double) work_KK * pow(max_U0, (double) 1 / (double) basic_schedule_power), (double) basic_schedule_power);
      vU0.at(kk) = pow((double) kk / (double) work_KK * pow(max_U0, (double) 1 / (double) basic_schedule_power), (double) basic_schedule_power); //..
    }
  }
    //Rcpp::Rcout << "U0[work_KK = " << work_KK << "] = " << U0[work_KK] << std::endl;
  if ( input_check ) {
    if ( (prev_time > U0[work_KK]) && (time > U0[work_KK]) ) Rcpp::stop("'prev_time' should be less than the maximum of the basic analysis schedule.");
    if ( (time > U0[work_KK]) && (next_time > U0[work_KK]) ) Rcpp::stop("'time' should be less than the maximum of the basic analysis schedule.");
  }
    //Rcpp::Rcout << *intVar << "; " << intVar[1] << "; " << std::endl;
    //Rcpp::Rcout << vU0.size() << "; " << std::endl;
    //for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //  Rcpp::Rcout << vU0.at(kk) << "; " << std::endl;
    //}

  // Current stage //
  int analysis = (time > next_time);
  double prev_time_ = analysis * prev_time + (1 - analysis) * time; // When (1 - analysis), work test will be updated.
  double time_ = analysis * time + (1 - analysis) * next_time;
  int time_kk = 0;
  int prev_time_kk = 0;
  for ( int kk = 1; kk < work_KK_1; kk++ ) {
    //*doubleVar = (U0[kk - 1] + U0[kk]) / 2.; // + 0.99 * (U0[kk + 1] - U0[kk]);
    *doubleVar = (vU0.at(kk - 1) + vU0.at(kk)) / 2.; //..
    time_kk += (*doubleVar < time_);
    prev_time_kk += (*doubleVar < prev_time_);
  }
    //Rcpp::Rcout << "time: " << time << "; time_kk: " << time_kk << std::endl;
    //Rcpp::Rcout << "prev time: " << prev_time << "; prev time_kk: " << prev_time_kk << std::endl;
  *intVar = (time_ > prev_time_);
  //intVar[1] = (time_kk == work_KK) && (time_ < U0[work_KK]);
  intVar[1] = (time_kk == work_KK) && (time_ < vU0.at(work_KK));
  //int too_close_to_maxU0 = intVar[1] && (((time_ - U0[work_KK - 1]) / (U0[work_KK] - U0[work_KK - 1])) > 0.99);
  int too_close_to_maxU0 = intVar[1] && (((time_ - vU0.at(work_KK - 1)) / (vU0.at(work_KK) - vU0.at(work_KK - 1))) > 0.99);
  if ( too_close_to_maxU0 && !analysis ) Rcpp::Rcout << "'next_time' is too close to the maximum of the basic analysis schedule. The calculated probability of rejecting the null hypothesis may be inaccurate." << std::endl;
  if ( too_close_to_maxU0 &&  analysis ) Rcpp::Rcout << "'time' is too close to the maximum of the basic analysis schedule. The calculated probability of rejecting the null hypothesis may be inaccurate." << std::endl;
    //Rcpp::Rcout << "U0[work_KK = " << work_KK << "] = " << U0[work_KK] << "; time_ = " << time_ << "; time_ < U0[work_KK] = " << (time_ < U0[work_KK]) << std::endl;
    //Rcpp::Rcout << "work_KK = " << work_KK << " - " << time_kk << " + " << *intVar << " + " << intVar[1] << std::endl;
  work_KK = work_KK - time_kk + *intVar + intVar[1];
    //Rcpp::Rcout << "work_KK = " << work_KK << std::endl;
  work_KK_1 = work_KK + 1;
  std::vector<double> vU_k(work_KK_1);
  vU_k.resize(work_KK_1);
  double* U_k = vU_k.data();
    //Rcpp::Rcout << U_k << std::endl;
  //U_k[0] = prev_time_;
  vU_k.at(0) = prev_time_; //..
  //U_k[*intVar] = time_;
  vU_k.at(*intVar) = time_; //..
  for ( int kk = 0; kk < (work_KK - *intVar); kk++ ) {  //!!!ok
    //U_k[kk + 1 + *intVar] = U0[kk + 1 + time_kk - intVar[1]];
    vU_k.at(kk + 1 + *intVar) = vU0.at(kk + 1 + time_kk - intVar[1]); //..
    //Rcpp::Rcout << "U0[" << kk - *intVar << "] = " << U0[kk - *intVar] << std::endl;
  }
  U_k += analysis;
  std::vector<double> vU_k0 = vU_k;
  if ( analysis ) { //.. // [deletion of the prev_time]
    vU_k.erase(vU_k.begin());
  }
  work_KK -= analysis;
  work_KK_1 -= analysis;
    //Rcpp::Rcout << vU_k.data() << std::endl;
    //Rcpp::Rcout << U_k << std::endl;
    //Rcpp::Rcout << work_KK << std::endl;
    //Rcpp::Rcout << work_KK_1 << std::endl;
    //Rcpp::Rcout << vU_k.size() << "; " << std::endl;
    //for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //  Rcpp::Rcout << vU_k.at(kk) << "; " << std::endl;
    //}
    //for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //  Rcpp::Rcout << U_k[kk] << "; " << std::endl;
    //}

  //--- Hypothesis ---//
  const int pri_0_l = 1;
  const int pri_1_l = 10;
  const int prior_l = pri_0_l + pri_1_l;
  // Null hypothesis //
  std::vector<double> vH_th(prior_l);
  //double* H_th = vH_th.data();
  for ( int ii = 0; ii < prior_l; ii++ ) {
    //H_th[ii] = (double)ii / 2. * min_effect_size;
    vH_th.at(ii) = (double)ii / 2. * min_effect_size;
  }

  //double* H_th0 = H_th;
  //double* H_th1 = H_th + pri_0_l;
  //for( int ii = 0; ii < pri_1_l; ii++ ) {
  //  Rcpp::Rcout << H_th1[ii] << std::endl;
  //}

  //--- Prior weights ---//
  std::vector<double> vprior(prior_l);
  double* prior = vprior.data();
    //Rcpp::Rcout << "# work_test_norm_c # prior = " << prior << std::endl;
  if ( prior_dist.size() == 11 ) {
    std::copy(prior_dist.begin(), prior_dist.end(), vprior.begin());
  } else {
    std::fill(vprior.begin(), vprior.end(), 1);
  }
    //Rcpp::Rcout << "# work_test_norm_c # prior = " << prior << std::endl;
    //for( int ii = 0; ii < vprior.size(); ii++ ) {
    //  Rcpp::Rcout << vprior.at(ii) << std::endl;
    //}
    //Rcpp::Rcout << "# work_test_norm_c # vprior.size = " << vprior.size() << std::endl;
  // Null hypothesis //
  //double* pri_0 = prior;
  // Alternative hypotheses //
  //double* pri_1 = prior + 1;
  //--- Loss Function ---//
  double fss_weight = 1;
  double seq_weight = 1;
  // Working Cost (will be updated appropriately) //
  double cost00 = max_U0 / overall_sig_level;
  double cost0 = cost_type_1_err;
  doubleVar[0] = max_U0;
  doubleVar[1] = cost_type_2_err;
  double cost1 = doubleVar[(int) (cost_type_2_err > 0)];

  //=== PRELIMINARIES ===//
  //--- Initialize ---//
  *doubleVar = std::accumulate(vprior.begin(), vprior.end(), 0.);
    //Rcpp::Rcout << "# work_test_norm_c # sum of prior = " << *doubleVar << std::endl;
  for ( int ii = 0; ii < prior_l; ii++ ) {
    //prior[ii] = prior[ii] / *doubleVar;
    vprior.at(ii) = vprior.at(ii) / *doubleVar; //..
  }
    //for( int ii = 0; ii < vprior.size(); ii++ ) {
    //  Rcpp::Rcout << vprior.at(ii) << std::endl;
    //}

  //--- Grid Points (Jennison & Turnbull (2000), chap.19) ---//
  // simpson_div
  int xdev_l = 6 * simpson_div - 1;
  std::vector<double> vxdev(xdev_l);
  double* xdev = vxdev.data();
  for ( int ii = 1; ii < simpson_div; ii++ ) {
    //xdev[ii - 1] = 3 + 4 * log( (double) simpson_div / (double) ii );
    vxdev.at(ii - 1) = 3 + 4 * log( (double) simpson_div / (double) ii ); //..
  }
  for ( int ii = 0; ii < (2 * simpson_div + 1); ii++ ) {
    //xdev[simpson_div - 1 + ii] = 3 - (double) 3 * ii / (double) (2 * simpson_div);
    vxdev.at(simpson_div - 1 + ii) = 3 - (double) 3 * ii / (double) (2 * simpson_div); //..
  }
  for ( int ii = 0; ii < (3 * simpson_div - 1); ii++ ) {
    //xdev[ii + 3 * simpson_div] = -xdev[3 * simpson_div - 2 - ii];
    vxdev.at(ii + 3 * simpson_div) = -vxdev.at(3 * simpson_div - 2 - ii); //..
  }
    //for( int ii = 0; ii < (xdev_l); ii++ ) {
    //  Rcpp::Rcout << xdev[ii] << std::endl;
    //}

  //--- Grid Points for Each Time (ss0) ---//
  // SPRT-based continuation region //
  std::vector<double> vdummy_prior(pri_1_l);
  std::vector<double> vupper_lim(pri_1_l);
  std::vector<double> vup_wing_temp(work_KK_1);
  std::vector<double> vup_wing_buffer(work_KK_1);
  std::vector<double> vdiv_unit(work_KK_1);
  std::vector<int> vup_wing_units(work_KK_1);
  std::vector<int> vlw_wing_units(work_KK_1);
  std::vector<int> vwing_l(work_KK_1);
  //int wing_l_max;

  //double* dummy_prior = vdummy_prior.data();
  //double* upper_lim = vupper_lim.data();
  //double* up_wing_temp = vup_wing_temp.data();
  //double* up_wing_buffer = vup_wing_buffer.data();
  //double* div_unit = vdiv_unit.data();
  //int* up_wing_units = vup_wing_units.data();
  //int* lw_wing_units = vlw_wing_units.data();
  //int* wing_l = vwing_l.data();

  // SPRT-based upper bound of grids
  //     Sum[ log{exp(-(x-mu1)^2/2) / exp(-(x-mu0)^2/2) * prior1 / prior0} ] >= -log(alp)
  // <=> Sum[d(x - m)] >= -log(alp) where (d = mu1 - mu0, m = (mu1 + mu0) / 2)
  for ( int ii = 0; ii < pri_1_l; ii++ ) {
    //dummy_prior[ii] = log(cost00 / cost1 * prior[0] / prior[ii + 1]) / H_th1[ii] + H_th1[ii] * max_U0 / 2.;
    vdummy_prior.at(ii) = log(cost00 / cost1 * vprior.at(0) / vprior.at(ii + 1)) / vH_th.at(ii + 1) + vH_th.at(ii + 1) * max_U0 / 2.;  //..
  }
  doubleVar[0] = *std::min_element(vdummy_prior.begin(), vdummy_prior.end());

  doubleVar[1] = R::pnorm(-doubleVar[0], vH_th.at(0), sqrt(max_U0), 1, 0); // alp_inf
  for ( int ii = 0; ii < pri_1_l; ii++ ) {
    //upper_lim[ii] = -log(doubleVar[1]) / H_th1[ii];
    vupper_lim.at(ii) = -log(doubleVar[1]) / vH_th.at(ii + 1); //..
  }
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    for ( int ii = 0; ii < pri_1_l; ii++ ) {
      //dummy_prior[ii] = upper_lim[ii] + U_k[kk] * H_th1[ii] / 2.;
      vdummy_prior.at(ii) = vupper_lim.at(ii) + vU_k.at(kk) * vH_th.at(ii + 1) / 2.; //..
    }
    doubleVar[0] = (1 - ((1 - analysis) * (kk == 0))) * *std::min_element(vdummy_prior.begin(), vdummy_prior.end());
    // sd_wing = *xdev;
    //up_wing_temp[kk] = ceil(doubleVar[0] * 10) / 10.;
    vup_wing_temp.at(kk) = ceil(doubleVar[0] * 10) / 10.; //..
    //div_unit[kk] = sqrt(U_k[kk] - U_k[kk - 1 + (kk == 0) - (kk == 0) * analysis]); // this div_unit is sqrt(diff(U_k));
    vdiv_unit.at(kk) = sqrt(vU_k.at(kk) - vU_k0.at(kk + ((kk == 0) - 1) * (1 - analysis))); //..
    //up_wing_buffer[kk] = div_unit[kk] * *xdev;
    vup_wing_buffer.at(kk) = vdiv_unit.at(kk) * *xdev; //..
      //Rcpp::Rcout << "up_buf: " << up_wing_buffer[kk] << std::endl;
    //div_unit[kk] = div_unit[kk] * 3 / (double) (2 * simpson_div);
    vdiv_unit.at(kk) = vdiv_unit.at(kk) * 3 / (double) (2 * simpson_div); //..
      //Rcpp::Rcout << "div_unit: " << div_unit[kk] << std::endl;
    doubleVar[1] = (kk == 0) * (1 - analysis); // for avoidance of division by 0 [2018.9.18 new]
    //up_wing_units[kk] = ceil((up_wing_temp[kk] + up_wing_buffer[kk] * 2) / div_unit[kk]);
    vup_wing_units.at(kk) = ceil((vup_wing_temp.at(kk) + vup_wing_buffer.at(kk) * 2) / (vdiv_unit.at(kk) + doubleVar[1])); //.. [2018.9.18 add doubleVar[1]]
    //lw_wing_units[kk] = ceil((-sqrt(U_k[kk]) * *xdev) / div_unit[kk]);
    vlw_wing_units.at(kk) = ceil((-sqrt(vU_k.at(kk)) * *xdev) / (vdiv_unit.at(kk) + doubleVar[1])); //.. [2018.9.18 add doubleVar[1]]
    //wing_l[kk] = up_wing_units[kk] - lw_wing_units[kk] + 1;
    vwing_l.at(kk) = vup_wing_units.at(kk) - vlw_wing_units.at(kk) + 1; //..
      //Rcpp::Rcout << "up_wing_units: " << up_wing_units[kk] << "; lw_wing_units: " << lw_wing_units[kk] << "; wing_l: " << wing_l[kk] << std::endl;
  }
  //wing_l[0] = 3;
  vwing_l.at(0) = 3; //..

  std::vector<double> vss0[105];  // num of basic_schedule <= 100
  vss0[0].reserve(3);
  vss0[0].resize(3);
  for ( int kk = 1; kk < work_KK_1; kk++ ) {
    //vss0[kk].reserve(wing_l[kk]);
    vss0[kk].reserve(vwing_l.at(kk)); //..
    //vss0[kk].resize(wing_l[kk]);
    vss0[kk].resize(vwing_l.at(kk)); //..
    //Rcpp::Rcout << "capacity (kk): " << vss0[kk].capacity() << std::endl;
  }

  if ( analysis ) {
    //vss0[0].at(0) = up_wing_units[0] * div_unit[0];
    vss0[0].at(0) = vup_wing_units.at(0) * vdiv_unit.at(0); //..
    vss0[0].at(1) = 0;
    //vss0[0].at(2) = lw_wing_units[0] * div_unit[0];
    vss0[0].at(2) = vlw_wing_units.at(0) * vdiv_unit.at(0);
  } else {
    std::fill(vss0[0].begin(), vss0[0].end(), 0);
  }
  //up_wing_units[0] = 1;
  vup_wing_units[0] = 1; //..
  //lw_wing_units[0] = -1;
  vlw_wing_units[0] = -1; //..
    //for ( int ii = 0; ii < wing_l[0]; ii++ ) {
    //  Rcpp::Rcout << "vss0[0].at(" << ii << "): " << vss0[0].at(ii) << std::endl;
    //}
  for ( int kk = 1; kk < work_KK_1; kk++ ) {
    double* ss0 = vss0[kk].data();
    for ( int ii = 0; ii < vwing_l.at(kk); ii++ ) {
      //ss0[ii] = (up_wing_units[kk] - ii) * div_unit[kk];
      ss0[ii] = (vup_wing_units.at(kk) - ii) * vdiv_unit.at(kk); //..
      //ss0[ii] = (up_wing_units[kk] - ii) * div_unit[kk];
      ss0[ii] = (vup_wing_units.at(kk) - ii) * vdiv_unit.at(kk); //..
    }
    //Rcpp::Rcout << "stage: " << kk << "; min(ss0): " << ss0[wing_l[kk] - 1] << "; max(ss0): " << ss0[0] << std::endl;
  }

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Memories for results                                                      //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Initialize ===//
  int cc_k_add = 2;

  std::vector<int> vgg_odd_l(work_KK_1);
  int* gg_odd_l = vgg_odd_l.data();
  std::vector<double> vgg_odd[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //gg_odd_l[kk] = wing_l[kk] + (kk > 0) * cc_k_add;
    vgg_odd_l.at(kk) = vwing_l.at(kk) + (kk > 0) * cc_k_add; //..
    vgg_odd[kk].reserve(gg_odd_l[kk]);
    vgg_odd[kk].resize(gg_odd_l[kk]);
      //Rcpp::Rcout << "gg_odd_l[" << kk << "]: " << gg_odd_l[kk] << std::endl;
  }
  std::vector<int> vgg_l(work_KK_1);
  int* gg_l = vgg_l.data();
  std::vector<double> vgg[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //gg_l[kk] = gg_odd_l[kk] * 2 - 1;
    vgg_l.at(kk) = vgg_odd_l.at(kk) * 2 - 1; //..
    vgg[kk].reserve(gg_l[kk]);
    vgg[kk].resize(gg_l[kk]);
      //Rcpp::Rcout << "gg_l[" << kk << "]: " << gg_l[kk] << std::endl;
  }
  std::vector<double> veta[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //veta[kk].reserve(gg_l[kk]);
    veta[kk].reserve(vgg_l.at(kk)); //..
    //veta[kk].resize(gg_l[kk]);
    veta[kk].resize(vgg_l.at(kk)); //..
  }
  std::vector<double> vpr_rej_H0[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    vpr_rej_H0[kk].reserve(gg_l[kk]);
    vpr_rej_H0[kk].resize(gg_l[kk]);
  }
  //std::vector<double> vpr_acc_H0[work_KK_1];
  //for( int kk = 0; kk < work_KK_1; kk++ ) {
  //  vpr_acc_H0[kk].reserve(gg_l[kk]);
  //  vpr_acc_H0[kk].resize(gg_l[kk]);
  //}
  std::vector<double> vexp_time[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    vexp_time[kk].reserve(gg_l[kk]);
    vexp_time[kk].resize(gg_l[kk]);
  }
  std::vector<int> vcc_odd_n(work_KK_1);
  int* cc_odd_n = vcc_odd_n.data();
  std::vector<int> vcc_n(work_KK_1);
  int* cc_n = vcc_n.data();
  std::vector<double> vcc(work_KK_1);
  double* cc = vcc.data();

  // Initialize (dummy variables) //
  int gg_l_max = *std::max_element(vgg_l.begin(), vgg_l.end());
  std::vector<double> vval_k(gg_l_max);
  double* val_k = vval_k.data();

  //=== Initialize (status variables in the structure) ===//
  struct base str_base;
  str_base.hypotheses = vH_th.data();
  str_base.prior = prior;
  str_base.num_of_hypotheses = prior_l;
  str_base.effect_size = theta;
  //str_base.cost0 = cost0;
  str_base.cost1 = cost1;
  str_base.analysis = analysis;
  str_base.max_time = max_U0;
  str_base.work_KK = work_KK; //
  str_base.U_k = vU_k.data(); //
  str_base.fss_w = fss_weight;
  str_base.seq_w = seq_weight;
  str_base.stat = stat;
  str_base.xdev = vxdev.data(); //
  str_base.xdev_l = xdev_l;
  str_base.up_wing_temp = vup_wing_temp.data(); //
  str_base.up_wing_buffer = vup_wing_buffer.data(); //
  str_base.div_unit = vdiv_unit.data();
  str_base.up_wing_units = vup_wing_units.data();
  str_base.wing_l = vwing_l.data();
  str_base.cc_k_add = cc_k_add;
  str_base.tol_boundary = tol_boundary;

  struct result str_result;
  str_result.cc_k_add = cc_k_add;
  str_result.vss0 = vss0;
  str_result.gg_odd_l = gg_odd_l;
  str_result.vgg_odd = vgg_odd;
  str_result.gg_l = gg_l;
  str_result.vgg = vgg;
  str_result.veta = veta;
  str_result.vpr_rej_H0 = vpr_rej_H0;
  str_result.vexp_time = vexp_time;
  str_result.cc_odd_n = cc_odd_n;
  str_result.cc_n = cc_n;
  str_result.cc = cc;
  str_result.val_k = val_k;

  struct ground str_ground;
  str_ground.str_base = str_base;
  str_ground.str_result = str_result;

// vectorizeここまで
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Design Construction                                                       //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Cost for erroneous rejection of H0 ===//
  //cost_upper_lim = (prev_cost != 0) * prev_cost + (prev_cost == 0) * (max_U0 * prior_l * 10);
  doubleVar[1] = sqrt(max_U0) * (-R::qnorm(cond_alpha, 0, 1, 1, 0));
  doubleVar[2] = doubleVar[1] / max_U0;  // MLE at critical value of z-test
    //Rcpp::Rcout << "# work_test_norm_c # MLE at critical value of z-test = " << doubleVar[2] << std::endl;
  *doubleVar = 1 / vprior.at(0) * exp((doubleVar[1] - max_U0 * doubleVar[2] / 2) * doubleVar[2] + log((1 - vprior.at(0)) * cost1));
    //Rcpp::Rcout << "# work_test_norm_c # cost for z-test = " << *doubleVar << std::endl;
  cost_upper_lim = (prev_cost != 0) * prev_cost + (prev_cost == 0) * (*doubleVar);

  if ( cost0 == 0) {
    ////Rcpp::Rcout << "# work_test_norm_c # cost search start." << std::endl;
    cost0 = bisection_inverse(construct,
            cond_alpha, &str_ground, 1e-6, cost_upper_lim,
            false, true, false, tol_cost);
  }
  double pr_type_1_err = construct(cost0, &str_ground);
  ////Rcpp::Rcout << "# work_test_norm_c # cost for (rej H0 | H0): " << cost0 << std::endl;

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Operating Characteristic                                                  //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //std::vector<double> 

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Output                                                                    //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//

    //Rcpp::Rcout << "# TEST # " << 1 << std::endl;
  //if( analysis ) { vU_k.erase(vU_k.begin() + 0); } // this has been done at the place labeled [deletion of the prev_time]. comment out at Aug 14, 2018
  Rcpp::List d_par = Rcpp::List::create(
    Rcpp::Named("overall_sig_level") = overall_sig_level,
    Rcpp::Named("work_beta") = work_beta,
    Rcpp::Named("min_effect_size") = min_effect_size,
    Rcpp::Named("prev_time") = prev_time,
    Rcpp::Named("time") = time,
    Rcpp::Named("next_time") = next_time,
    Rcpp::Named("stat") = stat,
    Rcpp::Named("U_0") = vU0,
    Rcpp::Named("U_k") = vU_k,
    Rcpp::Named("H_th") = vH_th,
    Rcpp::Named("prior") = vprior,
    Rcpp::Named("cost0") = cost0,
    Rcpp::Named("cost1") = cost1,
    Rcpp::Named("fss_weight") = fss_weight,
    Rcpp::Named("seq_weight") = seq_weight,
    Rcpp::Named("simpson_div") = simpson_div,
    Rcpp::Named("tol_boundary") = tol_boundary,
    Rcpp::Named("tol_cost") = tol_cost
    );

    //Rcpp::Rcout << "# TEST # " << 2 << std::endl;
  Rcpp::List rss0(work_KK_1);
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    rss0[kk] = vss0[kk];
  }
    //Rcpp::Rcout << "# TEST # " << 3 << std::endl;
  Rcpp::List rgg_odd(work_KK_1);
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    rgg_odd[kk] = vgg_odd[kk];
  }
    //Rcpp::Rcout << "# TEST # " << 4 << std::endl;
  Rcpp::List rgg(work_KK_1);
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    rgg[kk] = vgg[kk];
  }
    //Rcpp::Rcout << "# TEST # " << 5 << std::endl;
  Rcpp::List rpr_rej_H0(work_KK_1);
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    rpr_rej_H0[kk] = vpr_rej_H0[kk];
  }
    //Rcpp::Rcout << "# TEST # " << 6 << std::endl;
  Rcpp::List reta(work_KK_1);
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    reta[kk] = veta[kk];
  }

    //Rcpp::Rcout << "# TEST # " << 7 << std::endl;
  Rcpp::List d_char;
  if ( out_process ) {
    //Rcpp::Rcout << "# TEST # " << 8 << std::endl;
    d_char = Rcpp::List::create(
      Rcpp::Named("boundary") = vcc,
      Rcpp::Named("effect_size") = theta,
      Rcpp::Named("pr_rej_H0") = pr_type_1_err,
      //expected_loss
      //exp_time
      Rcpp::Named("grid_base_k") = rss0,
      Rcpp::Named("grid_odd_k") = rgg_odd,
      Rcpp::Named("grid_k") = rgg,
      Rcpp::Named("boundary_odd_n") = vcc_odd_n,
      Rcpp::Named("boundary_n") = vcc_n,
      Rcpp::Named("pr_rej_H0_k") = rpr_rej_H0,
      Rcpp::Named("expected_loss_k") = reta
      //exp_time_k
    );
  } else {
    //Rcpp::Rcout << "# TEST # " << 9 << std::endl;
    d_char = Rcpp::List::create(
      Rcpp::Named("boundary") = vcc,
      Rcpp::Named("effect_size") = theta,
      Rcpp::Named("pr_rej_H0") = pr_type_1_err
      //expected_loss
      //exp_time
    );
  }

    //Rcpp::Rcout << "# TEST # " << 10 << std::endl;
  Rcpp::List d_out = Rcpp::List::create(
    Rcpp::Named("par") = d_par,
    Rcpp::Named("char") = d_char
    );

    //Rcpp::Rcout << "# TEST # " << 11 << std::endl;
  return d_out;
    //Rcpp::Rcout << "# TEST # " << 12 << std::endl;

// 残り：仮説ごとの棄却確率（ステージ、全体）、期待中止時間（（ベイズ、仮説ごと）×（ステージ、全体））、time statで条件付け
}


static double pr_rej_H0_lower_bound(double final_analysis, struct ground* ground) {
  struct base* str_base = &(ground->str_base);
  struct result* str_result = &(ground->str_result);

  //Rcpp::Rcout << "# pr_rej_H0_lower_bound # START" << std::endl;
  //Rcpp::Rcout << "# pr_rej_H0_lower_bound # final_analysis: " << final_analysis << std::endl;

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Recovery: Settings                                                        //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  const int work_KK = str_base->work_KK;
  double* U_k = str_base->U_k;
  double effect_size = str_base->effect_size;
  double stat = str_base->stat;

  int xdev_l = str_base->xdev_l;
  double* xdev = str_base->xdev;
  double* div_unit = str_base->div_unit;
  //int* up_wing_units = str_base->up_wing_units;
  //int* wing_l = str_base->wing_l;

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Recovery: Memories for results                                            //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  std::vector<double>* vss0 = str_result->vss0;
  int* gg_odd_l = str_result->gg_odd_l;
  std::vector<double>* vgg_odd = str_result->vgg_odd;
  int* gg_l = str_result->gg_l;
  std::vector<double>* vgg = str_result->vgg;
  std::vector<double>* vpr_rej_H0 = str_result->vpr_rej_H0;
  int* cc_odd_n = str_result->cc_odd_n;
  int* cc_n = str_result->cc_n;
  double* cc = str_result->cc;
  double* val_k = str_result->val_k;

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Main                                                                      //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Declarations (dummy box) ===//
  double doubleVar[8] = {};
  int intVar[8] = {};

  //=== Declarations (status of the (k)th stage) ===//
  struct base_time str_base_time;
  struct current_next str_current_next;

  //=== Declarations (pointers for the (k)th stage) ===//
  int kk;
  double* ss_k;
  int gg_odd_k_l;
  double* gg_odd_k;
  int gg_k_l;
  double* gg_k;
  double* pr_rej_H0_k;
  int* cc_odd_n_k;
  int* cc_n_k;
  double* cc_k;
  double t_k;

  //=== Declarations (pointers for the (k + 1)th stage) ===//
  double t_k_1;
  double d_t;
  double sq_d_t;
  int* cc_odd_n_k_1;
  //int* cc_n_k_1;
  double* cc_k_1;
  std::vector<int> vxdev_k(xdev_l);
  int* xdev_k = vxdev_k.data();

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // The Final Stage                                                           //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Initialize (variables for the stage) ===//
  kk = work_KK;
    //Rcpp::Rcout << "# pr_rej_H0_lower_bound # ..... Stage: " << kk << std::endl;

  ss_k = vss0[kk].data();
  gg_odd_k_l = gg_odd_l[kk];
  gg_odd_k = vgg_odd[kk].data();
  gg_k_l = gg_l[kk];
  gg_k = vgg[kk].data();
  pr_rej_H0_k = vpr_rej_H0[kk].data();
  cc_odd_n_k = &(cc_odd_n[kk]);
  cc_n_k = &(cc_n[kk]);
  cc_k = &(cc[kk]);
  t_k = U_k[kk];
  // structure //
  str_base_time.str_base = *str_base;
  str_base_time.stage = kk;
  str_base_time.time = t_k;


  //... Routine 2 ...//
  // gg_odd_k //
  for ( int ii = 0; ii < gg_odd_k_l; ii++ ) {
    gg_odd_k[ii] = ss_k[ii];
  }
  // gg_k //
  *gg_k = *gg_odd_k;
  for ( int ii = 1; ii < gg_odd_k_l; ii++ ) {
    gg_k[2 * ii] = gg_odd_k[ii];
    gg_k[2 * ii - 1] = (gg_odd_k[ii - 1] + gg_odd_k[ii]) / 2.;
  }
    //for( int ii = 0; ii < gg_k_l; ii++ ) {
    //  Rcpp::Rcout << "ii = " << ii << ": gg = " << gg_k[ii] << std::endl;
    //}
  //... Routine 2 end ...//

  // Conditional Power //
  d_t = final_analysis - t_k;
  sq_d_t = sqrt(d_t);
  if ( d_t > 1e-8 ) {
    for ( int ii = 0; ii < gg_k_l; ii++ ) {
      *doubleVar = -R::qnorm(val_k[ii], 0, 1, 1, 0) * sq_d_t; // Critial Value
      pr_rej_H0_k[ii] = R::pnorm(-*doubleVar, - effect_size * d_t, sq_d_t, 1, 0);
        //Rcpp::Rcout << "ii = " << ii << ": gg = " << gg_k[ii] << ": cond_err = " << val_k[ii] << "; pr_rej_H0 = " << pr_rej_H0_k[ii] << std::endl;
    }
  } else {
    for ( int ii = 0; ii < gg_k_l; ii++ ) {
      pr_rej_H0_k[ii] = 0;
        //Rcpp::Rcout << "ii = " << ii << ": gg = " << gg_k[ii] << ": cond_err = " << val_k[ii] << "; pr_rej_H0 = " << pr_rej_H0_k[ii] << std::endl;
    }
  }

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // k th Stage                                                                //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  kk = kk - 1;
  for ( ; kk > 0; kk--) {
  //for( int kk = work_KK - 1; kk > (work_KK - 3); kk--) {
    //=== Initialize (variables for the stage) ===//
      //Rcpp::Rcout << "# pr_rej_H0_lower_bound # ..... Stage: " << kk << std::endl;

    ss_k = vss0[kk].data();
    gg_odd_k_l = gg_odd_l[kk];
    gg_odd_k = vgg_odd[kk].data();
    gg_k_l = gg_l[kk];
    gg_k = vgg[kk].data();
    pr_rej_H0_k = vpr_rej_H0[kk].data();
    cc_odd_n_k = &(cc_odd_n[kk]);
    cc_n_k = &(cc_n[kk]);
    cc_k = &(cc[kk]);
    t_k = U_k[kk];
    str_base_time.stage = kk;
    str_base_time.time = t_k;

    cc_odd_n_k_1 = &(cc_odd_n[kk + 1]);
    //cc_n_k_1 = &(cc_n[kk + 1]);
    cc_k_1 = &(cc[kk + 1]);
    str_current_next.str_base_time = str_base_time;
    str_current_next.time_1 = U_k[kk + 1];
    str_current_next.gg_k_1 = vgg[kk + 1].data();
    str_current_next.gg_k_1_l = gg_l[kk + 1];
    str_current_next.value_1 = vpr_rej_H0[kk + 1].data();
    str_current_next.dummy = val_k;
    str_current_next.cc_k_1 = cc_k_1;
    str_current_next.cc_odd_n_k_1 = cc_odd_n_k_1;
    str_current_next.xdev_k = xdev_k;

    t_k_1 = U_k[kk + 1];
    d_t = t_k_1 - t_k;
    sq_d_t = sqrt(d_t);
    for ( int ii = 0; ii < xdev_l; ii++ ) {
      *intVar = round(-(xdev[ii] * sqrt(d_t) + effect_size * d_t) / div_unit[kk + 1]);
      xdev_k[ii] = *intVar;
      //Rcpp::Rcout << "xdev_k[" << ii << "]: " << *intVar << std::endl;
      //xdev_k[ii] = round(-xdev[ii] * sqrt(d_t) / div_unit[kk + 1]);
    }

    //... Routine 2 ...//
    // gg_odd_k //
    for ( int ii = 0; ii < gg_odd_k_l; ii++ ) {
      gg_odd_k[ii] = ss_k[ii];
    }
    // gg_k //
    *gg_k = *gg_odd_k;
    for ( int ii = 1; ii < gg_odd_k_l; ii++ ) {
      gg_k[2 * ii] = gg_odd_k[ii];
      gg_k[2 * ii - 1] = (gg_odd_k[ii - 1] + gg_odd_k[ii]) / 2.;
    }
      //for( int ii = 0; ii < gg_k_l; ii++ ) {
      //  Rcpp::Rcout << "ii = " << ii << ": gg = " << gg_k[ii] << std::endl;
      //}
    //... Routine 2 end ...//

    // continue //
    for ( int ii = 0; ii < gg_k_l; ii++ ) {
      pr_rej_H0_k[ii] = future_pr_rej_H0(gg_k[ii], &str_current_next)
                        + R::pnorm(-(*cc_k_1 - gg_k[ii]), - effect_size * d_t, sq_d_t, 1, 0);
        //Rcpp::Rcout << "gg[" << ii << "]: " << gg_k[ii] << "; Pr(rej H0): " << pr_rej_H0_k[ii] <<  std::endl;
    }
  }


  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // 1 st Stage                                                                //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Initialize (variables for the stage) ===//
  //kk = 0;
    //Rcpp::Rcout << "# pr_rej_H0_lower_bound # ..... Stage: " << kk << std::endl;

  ss_k = vss0[kk].data();
  gg_odd_k_l = gg_odd_l[kk];
  gg_odd_k = vgg_odd[kk].data();
  gg_k_l = gg_l[kk];
  gg_k = vgg[kk].data();
  pr_rej_H0_k = vpr_rej_H0[kk].data();
  //pr_acc_H0_k = vpr_acc_H0[kk].data();
  cc_odd_n_k = &(cc_odd_n[kk]);
  cc_n_k = &(cc_n[kk]);
  cc_k = &(cc[kk]);
  t_k = U_k[kk];
  str_base_time.stage = kk;
  str_base_time.time = t_k;

  cc_odd_n_k_1 = &(cc_odd_n[kk + 1]);
  //cc_n_k_1 = &(cc_n[kk + 1]);
  cc_k_1 = &(cc[kk + 1]);
  str_current_next.str_base_time = str_base_time;
  str_current_next.time_1 = U_k[kk + 1];
  str_current_next.gg_k_1 = vgg[kk + 1].data();
  str_current_next.gg_k_1_l = gg_l[kk + 1];
  str_current_next.value_1 = vpr_rej_H0[kk + 1].data();
  str_current_next.dummy = val_k;
  str_current_next.cc_k_1 = cc_k_1;
  str_current_next.cc_odd_n_k_1 = cc_odd_n_k_1;
  str_current_next.xdev_k = xdev_k;

  t_k_1 = U_k[kk + 1];
  d_t = t_k_1 - t_k;
  sq_d_t = sqrt(d_t);
  for ( int ii = 0; ii < xdev_l; ii++ ) {
    *intVar = round(-(xdev[ii] * sqrt(d_t) + effect_size * d_t) / div_unit[kk + 1]);
    xdev_k[ii] = *intVar;
    //Rcpp::Rcout << "xdev_k[" << ii << "]: " << *intVar << std::endl;
    //xdev_k[ii] = round(-xdev[ii] * sqrt(d_t) / div_unit[kk + 1]);
  }
  *doubleVar = *cc_k; // for later recovery
  *cc_k = -1;

  //... Routine 2 ...//
  // gg_odd_k //
  *cc_odd_n_k = 0;
  *gg_odd_k = *ss_k + stat;
  // gg_k //
  *cc_n_k = 0;
  *gg_k = *gg_odd_k;
  //... Routine 2 end ...//

  // continue //
  //double* gg_k_1 = vgg[kk + 1].data();
  //double* pr_rej_H0_kk_1 = vpr_rej_H0[kk + 1].data();
    //for( int ii = 0; ii < gg_l[kk + 1]; ii++ ) {
    //  Rcpp::Rcout << "gg_k_1[" << ii << "]: " << gg_k_1[ii] << "; rej0: " << pr_rej_H0_kk_1[ii] << std::endl;
    //}
  *pr_rej_H0_k = future_pr_rej_H0(*gg_k, &str_current_next)
                 + R::pnorm(-(*cc_k_1 - *gg_k), - effect_size * d_t, sq_d_t, 1, 0);

  //Rcpp::Rcout << "# pr_rej_H0_lower_bound # Pr(rej H0 | " << effect_size << "): " << *pr_rej_H0_k << std::endl;
  //Rcpp::Rcout << "# pr_rej_H0_lower_bound # END" << std::endl;

  //for ( int kk = work_KK; kk >= 0; kk-- ) {
  //  Rcpp::Rcout << "!!! Stage: " << kk << std::endl;
  //for ( int ii = 0; ii < wing_l[kk]; ii++ ) {
  //  Rcpp::Rcout << "!!! " << vpr_rej_H0[kk][ii] << std::endl;
  //}
  //}

  *cc_k = *doubleVar;
  return *pr_rej_H0_k;
}


// [[Rcpp::export]]
double sample_size_norm_c(
  Rcpp::List initial_test = 0,
  const bool sample_size = true,
  const double effect_size = 0,
  const double time = 0,
  const double target_power = 0.8,
  const double final_time = 0,
  const double tol_sample_size = 1e-5,
  const bool input_check = true
  ) {
  double final_analysis = final_time;

  if ( input_check ) {
    if ( time < 0 ) {
      Rcpp::stop("'time' should be non-negative.");
    }
    if ( final_time < 0 ) {
      Rcpp::stop("'final_time' should be non-negative.");
    }
    if ( (target_power < 0) || (target_power > 1) ) {
      Rcpp::stop("'target_power' should be a value in (0, 1).");
    }
    if ( tol_sample_size <= 0 ) {
      Rcpp::stop("'tol_sample_size' should be positive.");
    }

    if ( (!sample_size) && (time > final_time) ) {
      Rcpp::Rcout << "NOTE: Because 'sample_size' is FLASE but 'final_time' is omitted or less than 'time', 'time' was substituted into 'final_time'." << std::endl;
    }
    if ( sample_size && (effect_size <= 0) ) {
      Rcpp::stop("When 'sample_size' is TRUE, 'effect_size' should be positive.");
    }
  }

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Settings                                                                  //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Dummy Box ===//
  double doubleVar[8] = {};
  int intVar[8] = {};
  Rcpp::NumericVector rvec;

  //=== Input Working Test ===//
  Rcpp::List d_par;
  Rcpp::List d_char;
  d_par = initial_test["par"];
  d_char = initial_test["char"];

  //double time = d_par["time"];
  //double stat = 0; //d_par["stat"];
  std::vector<double> vU0 = d_par["U_0"];
  int work_KK_1 = vU0.size();
  int work_KK = work_KK_1 - 1;
  //double* U0 = vU0.data();

  if ( (!sample_size) && (final_analysis < time) ) {
    final_analysis = time;
  }

  // Current stage //
  int time_kk = 0;
  for ( int kk = 1; kk < work_KK_1; kk++ ) {
    //*doubleVar = (U0[kk - 1] + U0[kk]) / 2.; // + 0.99 * (U0[kk + 1] - U0[kk]);
    *doubleVar = (vU0.at(kk - 1) + vU0.at(kk)) / 2.; //..
    time_kk += (*doubleVar < time);
  }
    //Rcpp::Rcout << "# sample_size_norm_c # time: " << time << "; time_kk: " << time_kk << std::endl;
  *intVar = (time_kk == 0) && (time > 0);
  //intVar[1] = (time_kk == work_KK) && (time < U0[work_KK]);
  intVar[1] = (time_kk == work_KK) && (time < vU0.at(work_KK));
  if ( (*intVar + intVar[1]) > 0 ) {
    vU0.insert(vU0.begin() + *intVar + intVar[1] * work_KK, 0);
  }
  //if ( *intVar ) {
  //  vU0.insert(vU0.begin() + 1, 0);
  //}
  time_kk += *intVar;
  vU0.at(time_kk) = time;
  work_KK += *intVar + intVar[1];
  work_KK_1 = work_KK + 1;
    //for ( int ii = 0; ii < vU0.size(); ii++ ) {
    //  Rcpp::Rcout << "# sample_size_norm_c # vU0.at(" << ii << ") = " << vU0.at(ii) << std::endl;
    //}

  Rcpp::NumericVector nvU0(work_KK_1); // std::vector to Rcpp::NumericVector
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    nvU0.at(kk) = vU0.at(kk);
      //Rcpp::Rcout << "nvU0.at(" << kk << ") = " << nvU0.at(kk) << std::endl;
  }
  Rcpp::NumericVector nvprior(11);
  nvprior = d_par["prior"];

  Rcpp::List initial_test_detail;
  initial_test_detail = work_test_norm_c(
    d_par["overall_sig_level"],  //overall_sig_level
    d_par["work_beta"],  //work_beta
    d_par["overall_sig_level"],  //cond_alpha
    d_par["cost0"],  //cost_type_1_err
    d_par["cost1"],  //cost_type_2_err
    0,  //prev_cost
    d_par["min_effect_size"],  //min_effect_size
    0,  //effect_size
    0,  //basic_schedule_num
    0,  //basic_schedule_power
    nvU0,  //basic_schedule
    nvprior,  //prior_dist
    0,  //prev_time
    0,  //time
    0,  //next_time
    0,  //stat
    false,  //input_check
    true,  //out_process
    d_par["simpson_div"],  //simpson_div
    d_par["tol_boundary"],  //tol_boundary
    d_par["tol_cost"]);  //tol_cost
  
  //d_par = initial_test_detail["par"];
  d_char = initial_test_detail["char"];


  //=== Initialize (status variables in the structure) ===//
  //struct base str_base;
  //str_base.effect_size = effect_size;
  //str_base.work_KK = work_KK;
  //str_base.U_k = U_k; //
  //str_base.stat = stat;
  //str_base.div_unit = div_unit;
  //str_base.up_wing_units = up_wing_units;
  //str_base.wing_l = wing_l;


  //U0 = vU0.data();
  std::vector<double> vU_k = vU0; //..
  //double* U_k = U0;
  work_KK = time_kk;  //******************************************************************************
  work_KK_1 = work_KK + 1;   //******************************************************************************
  int simpson_div = d_par["simpson_div"];
  std::vector<double> vcc = d_char["boundary"];
  //double* cc = vcc.data();
  std::vector<double> vdiv_unit(work_KK_1);
  //double* div_unit = vdiv_unit.data();
  std::vector<int> vup_wing_units(work_KK_1);
  //int* up_wing_units = vup_wing_units.data();
  std::vector<int> vlw_wing_units(work_KK_1);
  //int* lw_wing_units = vlw_wing_units.data();
  std::vector<int> vwing_l(work_KK_1);
  //int* wing_l = vwing_l.data();

    //Rcpp::Rcout << "time: " << time << std::endl;
    //Rcpp::Rcout << "stat: " << stat << std::endl;
    //Rcpp::Rcout << "vU0.size: " << vU0.size() << std::endl;
    //Rcpp::Rcout << "work_KK: " << work_KK << std::endl;
    //Rcpp::Rcout << "simpson_div: " << simpson_div << std::endl;
    //for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //  Rcpp::Rcout << "vU0.at(" << kk << ") = " << vU0.at(kk) << std::endl;
    //}
    //for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //  Rcpp::Rcout << "vcc.at(" << kk << ") = " << vcc.at(kk) << std::endl;
    //}

  //--- Grid Points (Jennison & Turnbull (2000), chap.19) ---//
  // simpson_div
  int xdev_l = 6 * simpson_div - 1;
  std::vector<double> vxdev(xdev_l);
  double* xdev = vxdev.data();
  for ( int ii = 1; ii < simpson_div; ii++ ) {
    //xdev[ii - 1] = 3 + 4 * log( (double) simpson_div / (double) ii );
    vxdev.at(ii - 1) = 3 + 4 * log( (double) simpson_div / (double) ii ); //..
  }
  for ( int ii = 0; ii < (2 * simpson_div + 1); ii++ ) {
    //xdev[simpson_div - 1 + ii] = 3 - (double) 3 * ii / (double) (2 * simpson_div);
    vxdev.at(simpson_div - 1 + ii) = 3 - (double) 3 * ii / (double) (2 * simpson_div); //..
  }
  for ( int ii = 0; ii < (3 * simpson_div - 1); ii++ ) {
    //xdev[ii + 3 * simpson_div] = -xdev[3 * simpson_div - 2 - ii];
    vxdev.at(ii + 3 * simpson_div) = -vxdev.at(3 * simpson_div - 2 - ii); //..
  }

  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //div_unit[kk] = sqrt(U_k[kk] - U_k[kk - 1 + (kk == 0)]); // this div_unit is sqrt(diff(U_k));
    vdiv_unit.at(kk) = sqrt(vU_k.at(kk) - vU_k.at(kk - 1 + (kk == 0))); //..
    //div_unit[kk] = div_unit[kk] * 3 / (double) (2 * simpson_div);
    vdiv_unit.at(kk) = vdiv_unit.at(kk) * 3 / (double) (2 * simpson_div); //..
      //Rcpp::Rcout << "div_unit: " << div_unit[kk] << std::endl;
    doubleVar[1] = (kk == 0); // for avoidance of division by 0 [2018.9.18 new]
    //up_wing_units[kk] = ceil(cc[kk] / div_unit[kk]);
    if ( kk > 0 ) { vup_wing_units.at(kk) = ceil(vcc.at(kk) / (vdiv_unit.at(kk) + doubleVar[1])); } //.. because vcc.at(0) == Inf for kk == 0
    //lw_wing_units[kk] = ceil((-sqrt(U_k[kk]) * *xdev) / div_unit[kk]);
    vlw_wing_units.at(kk) = ceil((-sqrt(vU_k.at(kk)) * *xdev) / (vdiv_unit.at(kk) + doubleVar[1])); //.. [2018.9.18 add doubleVar[1]]
      //up_wing_units[0] = 0;
      vup_wing_units.at(0) = 0; //..
      //lw_wing_units[0] = 0;
      vlw_wing_units.at(0) = 0; //..
    //wing_l[kk] = up_wing_units[kk] - lw_wing_units[kk] + 1;
    vwing_l.at(kk) = vup_wing_units.at(kk) - vlw_wing_units.at(kk) + 1; //..
      //Rcpp::Rcout << "up_wing_units: " << up_wing_units[kk] << "; lw_wing_units: " << lw_wing_units[kk] << "; wing_l: " << wing_l[kk] << std::endl;
  }


  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Memories for results                                                      //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  std::vector<double> vss0[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //vss0[kk].reserve(wing_l[kk]);
    vss0[kk].reserve(vwing_l.at(kk)); //..
    //vss0[kk].resize(wing_l[kk]);
    vss0[kk].resize(vwing_l.at(kk)); //..
    //Rcpp::Rcout << "capacity (kk): " << vss0[kk].capacity() << std::endl;
  }
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    double* ss0 = vss0[kk].data();
    for ( int ii = 0; ii < vwing_l.at(kk); ii++ ) {
      //ss0[ii] = (kk > 0) * (up_wing_units[kk] - ii) * div_unit[kk];
      ss0[ii] = (kk > 0) * (vup_wing_units.at(kk) - ii) * vdiv_unit.at(kk);
    }
    //ss0[0] = cc[kk];
    ss0[0] = vcc.at(kk); //..
      //Rcpp::Rcout << "stage: " << kk << "; min(ss0): " << ss0[wing_l[kk] - 1] << "; max(ss0): " << ss0[0] << std::endl;
  }
  vss0[0].at(0) = 0;

  std::vector<int> vgg_odd_l(work_KK_1);
  int* gg_odd_l = vgg_odd_l.data();
  std::vector<double> vgg_odd[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //gg_odd_l[kk] = wing_l[kk];
    vgg_odd_l.at(kk) = vwing_l.at(kk); //..
    vgg_odd[kk].reserve(gg_odd_l[kk]);
    vgg_odd[kk].resize(gg_odd_l[kk]);
      //Rcpp::Rcout << "gg_odd_l[" << kk << "]: " << gg_odd_l[kk] << std::endl;
  }
  std::vector<int> vgg_l(work_KK_1);
  int* gg_l = vgg_l.data();
  std::vector<double> vgg[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //gg_l[kk] = gg_odd_l[kk] * 2 - 1;
    vgg_l.at(kk) = vgg_odd_l.at(kk) * 2 - 1; //..
    vgg[kk].reserve(gg_l[kk]);
    vgg[kk].resize(gg_l[kk]);
      //Rcpp::Rcout << "gg_l[" << kk << "]: " << gg_l[kk] << std::endl;
  }

  std::vector<double> vpr_rej_H0[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    vpr_rej_H0[kk].reserve(gg_l[kk]);
    vpr_rej_H0[kk].resize(gg_l[kk]);
  }

  std::vector<int> vcc_odd_n(work_KK_1);
  //int* cc_odd_n = vcc_odd_n.data();
  std::vector<int> vcc_n(work_KK_1);
  //int* cc_n = vcc_n.data();
  //std::vector<double> vcc(work_KK_1);
  //double* cc = vcc.data();

  // Initialize (dummy variables) //
  int gg_l_max = *std::max_element(vgg_l.begin(), vgg_l.end());
  std::vector<double> vval_k(gg_l_max);
  //double* val_k = vval_k.data();
  // Store conditional error at the final interim analysis into "val_k" //
  //vcc_n = d_char["boundary_n"];
  rvec = d_char["boundary_n"];
  vcc_n.resize(rvec.size());
  std::copy(rvec.begin(), rvec.end(), vcc_n.begin());
    //Rcpp::Rcout << "vval_k.data = " << vval_k.data() << std::endl;
  vval_k.resize(gg_l[work_KK]);
    //Rcpp::Rcout << "vval_k.data = " << vval_k.data() << std::endl;
  d_par = d_char["pr_rej_H0_k"];  // confirmation of copy region
  //d_par = d_char["grid_k"];
  std::vector<double> dummy = d_par[work_KK];
  std::copy(dummy.begin() + vcc_n.at(work_KK) + 1, dummy.end(), vval_k.begin());
    //for ( int ii = 0; ii < vval_k.size(); ii++ ) { //gg_l[work_KK]
    //  Rcpp::Rcout << "vval_k.at(" << ii << ") = " << vval_k.at(ii) << std::endl;
    //}
  // Set all cc_odd_n and cc_n to -1 //
  std::fill(vcc_odd_n.begin(), vcc_odd_n.end(), -1);
  std::fill(vcc_n.begin(), vcc_n.end(), -1);


  //=== Initialize (status variables in the structure) ===//
  struct base str_base;
  str_base.effect_size = effect_size;
  //str_base.max_time = max_U0;
  str_base.work_KK = work_KK; //
  str_base.U_k = vU_k.data(); //
  str_base.stat = 0; //stat;
  str_base.xdev = vxdev.data(); //
  str_base.xdev_l = xdev_l;
  str_base.div_unit = vdiv_unit.data();
  str_base.up_wing_units = vup_wing_units.data();
  str_base.wing_l = vwing_l.data();

  struct result str_result;
  str_result.vss0 = vss0;
  str_result.gg_odd_l = vgg_odd_l.data();
  str_result.vgg_odd = vgg_odd;
  str_result.gg_l = vgg_l.data();
  str_result.vgg = vgg;
  str_result.vpr_rej_H0 = vpr_rej_H0;
  str_result.cc_odd_n = vcc_odd_n.data();
  str_result.cc_n = vcc_n.data();
  str_result.cc = vcc.data();
  str_result.val_k = vval_k.data();

  struct ground str_ground;
  str_ground.str_base = str_base;
  str_ground.str_result = str_result;



    //int kkk = 2;
    //Rcpp::Rcout << "gg_l = " << vgg_l.at(kkk) << std::endl;
    //for ( int ii = 0; ii < wing_l[kkk]; ii++ ) {
    //  Rcpp::Rcout << "ss0.at(" << ii << ") = " << vss0[kkk].at(ii) << std::endl;
    //}
    //Rcpp::List rgg_ = d_char["grid_k"];
    //std::vector<double> vgg_ = rgg_[kkk];
    //std::vector<int> vcc_n_ = d_char["boundary_n"];
    //std::vector<double> vcc_ = d_char["boundary"];
    //Rcpp::Rcout << "vgg_.size = " << vgg_.size() << std::endl;
    //Rcpp::Rcout << "vgg_.capacity = " << vgg_.capacity() << std::endl;
    //Rcpp::Rcout << "vcc_n_.size = " << vcc_n_.size() << std::endl;
    //Rcpp::Rcout << "vcc_n_.capacity = " << vcc_n_.capacity() << std::endl;
    //Rcpp::Rcout << "cc_k = " << vcc_.at(kkk) << std::endl;
    //Rcpp::Rcout << "cc_k_n = " << vcc_n_.at(kkk) << std::endl;
    //for ( int ii = 0; ii < vgg_.size(); ii++ ) {
    //  Rcpp::Rcout << "gg_.at(" << ii << ") = " << vgg_.at(ii) << std::endl;
    //}

    //std::vector<double> vgg_;
    //Rcpp::List rgg_ = d_char["grid_k"];
    //std::vector<int> vcc_n_ = d_char["boundary_n"];
    //Rcpp::Rcout << "work_KK_1 = " << work_KK_1 << std::endl;
    //for ( int kkk = 0; kkk < work_KK_1; kkk++ ) {
    //  vgg_ = rgg_[kkk];
    //  Rcpp::Rcout << "Old gg_l[" << kkk << "] = " << vgg_.size() - 1 - vcc_n_.at(kkk) << "; New gg_l[" << kkk << "] = " << vgg_l.at(kkk) << std::endl;
    //}

  if ( sample_size ) {
    doubleVar[1] = pr_rej_H0_lower_bound(time, &str_ground);  // Cumulative power until the final interim analysis
    if ( doubleVar[1] >= target_power ) {
      *doubleVar = time;
    } else {
      *doubleVar = bisection_inverse(pr_rej_H0_lower_bound,
        target_power, &str_ground, time, time * 10,
        true, false, false, tol_sample_size);
    }
    ////Rcpp::Rcout << "# sample_size_norm_c # Sample Size = " << *doubleVar << std::endl;
  } else {
    *doubleVar = pr_rej_H0_lower_bound(final_analysis, &str_ground);
    //Rcpp::Rcout << "# sample_size_norm_c # pr_rej_H0 = " << *doubleVar << std::endl;
    ////Rcpp::Rcout << "# sample_size_norm_c # Pr(rej H0 | " << effect_size << "): " << *doubleVar << std::endl;
  }

  return doubleVar[0];
}


// my struct base_test;
struct base_test {
  Rcpp::List* test;
  std::vector<double>* times;
  std::vector<double>* stats;
  bool sig_level_exhaustive;
  double ci_coef;
  double* cost0;
  Rcpp::List* U_k;
  Rcpp::List* cc;
} ;

static double pr_rej_H0(
std::vector<double>* pvU_k = 0,
std::vector<double>* pvcc = 0,
double stat = 0,
const double effect_size = 0,
const int simpson_div = 6) {

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Settings                                                                  //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Dummy Box ===//
  double doubleVar[8] = {};
  //int intVar[8] = {};

  double* U_k = (*pvU_k).data();
  std::vector<double> vU_k = (*pvU_k);
  int work_KK = (*pvU_k).size() - 1;
  int work_KK_1 = work_KK + 1;
  std::vector<double> vcc = (*pvcc);
  //double* cc = (*pvcc).data();

    //for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //  Rcpp::Rcout << "# pr_rej_H0 # U_k[ " << kk << " ] = " << U_k[kk] << "; boundary = " << cc[kk] << std::endl;
    //}

  std::vector<double> vdiv_unit(work_KK_1);
  //double* div_unit = vdiv_unit.data(); //..
  std::vector<int> vup_wing_units(work_KK_1);
  //int* up_wing_units = vup_wing_units.data(); //..
  std::vector<int> vlw_wing_units(work_KK_1);
  //int* lw_wing_units = vlw_wing_units.data(); //..
  std::vector<int> vwing_l(work_KK_1);
  //int* wing_l = vwing_l.data(); //..

    //Rcpp::Rcout << "time: " << time << std::endl;
    //Rcpp::Rcout << "stat: " << stat << std::endl;
    //Rcpp::Rcout << "vU0.size: " << vU0.size() << std::endl;
    //Rcpp::Rcout << "work_KK: " << work_KK << std::endl;
    //Rcpp::Rcout << "simpson_div: " << simpson_div << std::endl;
    //for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //  Rcpp::Rcout << "vU0.at(" << kk << ") = " << vU0.at(kk) << std::endl;
    //}
    //for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //  Rcpp::Rcout << "vcc.at(" << kk << ") = " << vcc.at(kk) << std::endl;
    //}

  //--- Grid Points (Jennison & Turnbull (2000), chap.19) ---//
  // simpson_div
  int xdev_l = 6 * simpson_div - 1;
  std::vector<double> vxdev(xdev_l);
  //double* xdev = vxdev.data(); //..
  for ( int ii = 1; ii < simpson_div; ii++ ) {
    //xdev[ii - 1] = 3 + 4 * log( (double) simpson_div / (double) ii );
    vxdev.at(ii - 1) = 3 + 4 * log( (double) simpson_div / (double) ii ); //..
  }
  for ( int ii = 0; ii < (2 * simpson_div + 1); ii++ ) {
    //vxdev[simpson_div - 1 + ii] = 3 - (double) 3 * ii / (double) (2 * simpson_div);
    vxdev.at(simpson_div - 1 + ii) = 3 - (double) 3 * ii / (double) (2 * simpson_div); //..
  }
  for ( int ii = 0; ii < (3 * simpson_div - 1); ii++ ) {
    //xdev[ii + 3 * simpson_div] = -xdev[3 * simpson_div - 2 - ii];
    vxdev.at(ii + 3 * simpson_div) = -vxdev.at(3 * simpson_div - 2 - ii); //..
  }

  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //div_unit[kk] = sqrt(U_k[kk] - U_k[kk - 1 + (kk == 0)]); // this div_unit is sqrt(diff(U_k));
    vdiv_unit.at(kk) = sqrt(vU_k.at(kk) - vU_k.at(kk - 1 + (kk == 0))); // this div_unit is sqrt(diff(U_k)); //..
    //div_unit[kk] = div_unit[kk] * 3 / (double) (2 * simpson_div);
    vdiv_unit.at(kk) = vdiv_unit.at(kk) * 3 / (double) (2 * simpson_div); //..
      //Rcpp::Rcout << "div_unit: " << div_unit[kk] << std::endl;
    //up_wing_units[kk] = ceil(cc[kk] / div_unit[kk]);
    vup_wing_units.at(kk) = ceil(vcc.at(kk) / vdiv_unit.at(kk)); //..
    //lw_wing_units[kk] = ceil((-sqrt(U_k[kk]) * *xdev) / div_unit[kk]);
    //lw_wing_units[kk] = ceil((-sqrt(U_k[kk]) * *xdev + effect_size * U_k[kk]) / div_unit[kk]);
    vlw_wing_units.at(kk) = ceil((-sqrt(vU_k.at(kk)) * vxdev.at(0) + effect_size * vU_k.at(kk)) / vdiv_unit.at(kk)); //..
    //doubleVar[0] = ceil((-sqrt(U_k[kk]) * *xdev) / div_unit[kk]);
    //doubleVar[0] = -abs(up_wing_units[kk]);
    doubleVar[0] = -abs(vup_wing_units.at(kk)); //..
    //doubleVar[1] = up_wing_units[kk] - 1;
    doubleVar[1] = vup_wing_units.at(kk) - 1; //..
    //lw_wing_units[kk] = doubleVar[(int) (doubleVar[0] > doubleVar[1])];
    vlw_wing_units.at(kk) = doubleVar[(int) (doubleVar[0] > doubleVar[1])]; //..
    //lw_wing_units[kk] = ceil((-sqrt(U_k[kk]) * *xdev) / div_unit[kk]);
      //up_wing_units[0] = 0;
      vup_wing_units.at(0) = 0; //..
      //lw_wing_units[0] = 0;
      vlw_wing_units.at(0) = 0; //..
    //wing_l[kk] = up_wing_units[kk] - lw_wing_units[kk] + 1;
    vwing_l.at(kk) = vup_wing_units.at(kk) - vlw_wing_units.at(kk) + 1; //..
      //Rcpp::Rcout << "up_wing_units: " << up_wing_units[kk] << "; lw_wing_units: " << lw_wing_units[kk] << "; wing_l: " << wing_l[kk] << std::endl;
  }

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Memories for results                                                      //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  std::vector<double> vss0[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    vss0[kk].reserve(vwing_l[kk]);
    vss0[kk].resize(vwing_l[kk]);
    //Rcpp::Rcout << "capacity (kk): " << vss0[kk].capacity() << std::endl;
  }
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    double* ss0 = vss0[kk].data();
    //for ( int ii = 0; ii < wing_l[kk]; ii++ ) {
    for ( int ii = 0; ii < vwing_l.at(kk); ii++ ) { //..
      //ss0[ii] = (kk > 0) * (up_wing_units[kk] - ii) * div_unit[kk];
      ss0[ii] = (kk > 0) * (vup_wing_units.at(kk) - ii) * vdiv_unit.at(kk); //..
    }
    //ss0[0] = cc[kk];
    ss0[0] = vcc.at(kk); //..
      //Rcpp::Rcout << "stage: " << kk << "; min(ss0): " << ss0[wing_l[kk] - 1] << "; max(ss0): " << ss0[0] << std::endl;
  }
  vss0[0].at(0) = 0;

  std::vector<int> vgg_odd_l(work_KK_1);
  //int* gg_odd_l = vgg_odd_l.data(); //..
  std::vector<double> vgg_odd[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //gg_odd_l[kk] = wing_l[kk];
    vgg_odd_l.at(kk) = vwing_l.at(kk); //..
    //vgg_odd[kk].reserve(gg_odd_l[kk]);
    vgg_odd[kk].reserve(vgg_odd_l.at(kk)); //..
    //vgg_odd[kk].resize(gg_odd_l[kk]);
    vgg_odd[kk].resize(vgg_odd_l.at(kk));
      //Rcpp::Rcout << "gg_odd_l[" << kk << "]: " << gg_odd_l[kk] << std::endl;
  }
  std::vector<int> vgg_l(work_KK_1);
  //int* gg_l = vgg_l.data(); //..
  std::vector<double> vgg[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //gg_l[kk] = gg_odd_l[kk] * 2 - 1;
    vgg_l.at(kk) = vgg_odd_l.at(kk) * 2 - 1; //..
    //vgg[kk].reserve(gg_l[kk]);
    vgg[kk].reserve(vgg_l.at(kk));
    //vgg[kk].resize(gg_l[kk]);
    vgg[kk].resize(vgg_l.at(kk));
      //Rcpp::Rcout << "gg_l[" << kk << "]: " << gg_l[kk] << std::endl;
  }

  std::vector<double> vpr_rej_H0[105];  // num of basic_schedule <= 100
  for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //vpr_rej_H0[kk].reserve(gg_l[kk]);
    vpr_rej_H0[kk].reserve(vgg_l.at(kk));
    //vpr_rej_H0[kk].resize(gg_l[kk]);
    vpr_rej_H0[kk].resize(vgg_l.at(kk));
  }
  std::vector<int> vcc_odd_n(work_KK_1);
  //int* cc_odd_n = vcc_odd_n.data(); //..
  std::vector<int> vcc_n(work_KK_1);
  //int* cc_n = vcc_n.data(); //..

  // Initialize (dummy variables) //
  int gg_l_max = *std::max_element(vgg_l.begin(), vgg_l.end());
  std::vector<double> vval_k(gg_l_max);
  //double* val_k = vval_k.data(); //..
  // Set all cc_odd_n and cc_n to -1 //
  std::fill(vcc_odd_n.begin(), vcc_odd_n.end(), -1);
  std::fill(vcc_n.begin(), vcc_n.end(), -1);

  //=== Initialize (status variables in the structure) ===//
  struct base str_base;
  str_base.effect_size = effect_size;
  //str_base.max_time = max_U0;
  str_base.work_KK = work_KK; //
  str_base.U_k = U_k; //
  str_base.stat = stat;
  str_base.xdev = vxdev.data(); //
  str_base.xdev_l = xdev_l;
  //str_base.div_unit = div_unit;
  str_base.div_unit = vdiv_unit.data(); //..
  //str_base.up_wing_units = up_wing_units;
  str_base.up_wing_units = vup_wing_units.data();
  //str_base.wing_l = wing_l;
  str_base.wing_l = vwing_l.data(); //..

  struct result str_result;
  str_result.vss0 = vss0;
  //str_result.gg_odd_l = gg_odd_l;
  str_result.gg_odd_l = vgg_odd_l.data(); //..
  str_result.vgg_odd = vgg_odd;
  //str_result.gg_l = gg_l;
  str_result.gg_l = vgg_l.data(); //..
  str_result.vgg = vgg;
  str_result.vpr_rej_H0 = vpr_rej_H0;
  //str_result.cc_odd_n = cc_odd_n;
  str_result.cc_odd_n = vcc_odd_n.data(); //..
  //str_result.cc_n = cc_n;
  str_result.cc_n = vcc_n.data(); //..
  //str_result.cc = cc;
  str_result.cc = vcc.data(); //..
  //str_result.val_k = val_k;
  str_result.val_k = vval_k.data(); //..

  struct ground str_ground;
  str_ground.str_base = str_base;
  str_ground.str_result = str_result;


    //Rcpp::Rcout << "# pr_rej_H0 # time = " << U_k[work_KK] << std::endl;
  doubleVar[1] = 0;
  doubleVar[2] = pr_rej_H0_lower_bound(U_k[work_KK], &str_ground);
  doubleVar[3] = 1;
  // max([1], min([2], [3])) //
  *doubleVar = ((doubleVar[1] <= doubleVar[2]) && (doubleVar[2] <= doubleVar[3])) * doubleVar[2] + (doubleVar[2] > doubleVar[3]);

    //Rcpp::Rcout << "# pr_rej_H0 # Pr(rej H0 | " << effect_size << ") = " << *doubleVar << std::endl;
  return *doubleVar;
}

struct arg_pr_rej_H0 {
  std::vector<double>* vU_k;
  std::vector<double>* vcc;
  double stat;
  double effect_size;
  double simpson_div;
  double* cc_fin_kk;
} ;

static double pr_rej_H0_sol1(
double fin_cc,
struct arg_pr_rej_H0* parg) {
  double* cc_fin_kk = parg->cc_fin_kk;
  *cc_fin_kk = fin_cc;
  //Rcpp::Rcout << "# pr_rej_H0_sol1 # cc[fin_kk] = " << *cc_fin_kk << std::endl;
  double cond_power = pr_rej_H0(parg->vU_k, parg->vcc, parg->stat, parg->effect_size, parg->simpson_div);
  return cond_power;
}

static double project_power(
const double effect_size,
struct base_test* ptest
) {

  //Rcpp::Rcout << "# project_power # START" << std::endl;
  Rcpp::List initial_test = *(ptest->test);
  std::vector<double> vtimes = *(ptest->times);
  std::vector<double> vstats = *(ptest->stats);
  //bool sig_level_exhaustive = ptest->sig_level_exhaustive;
  //double ci_coef = ptest->ci_coef;
  Rcpp::List* prU_k = ptest->U_k;
  Rcpp::List* prcc = ptest->cc;

  int max_kk = vtimes.size();
  // after insertion of 0, time length should be greater than 3

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Settings                                                                  //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Dummy Box ===//
  double doubleVar[8] = {};
  //int intVar[8] = {};
  Rcpp::NumericVector rvec;

  //=== Initial Working Test ===//
  Rcpp::List init_par = initial_test["par"];
  Rcpp::List init_char = initial_test["char"];
  Rcpp::List d_par;
  Rcpp::List d_char;

  std::vector<double> vU0 = init_par["U_0"];
  int work_KK_1 = vU0.size();
  int work_KK = work_KK_1 - 1;
  int simpson_div = init_par["simpson_div"];
  //double* U0 = vU0.data();
  Rcpp::NumericVector nvU0(work_KK_1); // std::vector to Rcpp::NumericVector
  std::copy(vU0.begin(), vU0.end(), nvU0.begin());
    //for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //  Rcpp::Rcout << "nvU0.at(" << kk << ") = " << nvU0.at(kk) << std::endl;
    //}
  std::vector<double> vU_k(work_KK_1);
  std::vector<double> vcc(work_KK_1);

  //=== Variables ===//
  int kk;
  double cond_power;
  double pre_cond_power;
  double prev_time;
  double time;
  //double next_time;
  double prev_stat;
  double stat;
  double d_s, d_t;
  int fin_kk;
  double fin_cc;
  //double* cost0 = ptest->cost0;

  kk = max_kk - 1;
  prev_stat = vstats.at(kk - 1);
  stat = vstats.at(kk);
  prev_time = vtimes.at(kk - 1);
  time = vtimes.at(kk);
  d_s = stat - prev_stat;
  d_t = time - prev_time;
  cond_power = R::pnorm(-d_s, -effect_size * d_t, sqrt(d_t), 1, 0);
    //Rcpp::Rcout << "# project_power # Stage " << kk << "; cond power (effect_size = " << effect_size << ") = " << cond_power << std::endl;

  kk -= 1;
  fin_kk = 1;
  for ( ; kk > 0; kk-- ) { // final projection (kk == 0) is unnecessary because projection preserves cond power.
    //Rcpp::Rcout << "# project_power # Stage " << kk << std::endl;

    prev_stat = vstats.at(kk - 1 + (kk == 0));
    stat = vstats.at(kk);
    prev_time = vtimes.at(kk - 1 + (kk == 0));
    time = vtimes.at(kk);
      //Rcpp::Rcout << "# project_power # Stage " << kk << "; prev_time = " << prev_time << "; time = " << time << "; prev_stat = " << prev_stat << std::endl;
    //Rcpp::List work_test;
    //work_test = work_test_norm_c(
    //  init_par["overall_sig_level"],  //overall_sig_level
    //  init_par["work_beta"],  //work_beta
    //  0,  //cond_alpha
    //  cost0[kk],  //cost_type_1_err
    //  0,  //cost_type_2_err
    //  0,  //prev_cost
    //  init_par["min_effect_size"],  //min_effect_size
    //  0,  //effect_size
    //  0,  //basic_schedule_num
    //  0,  //basic_schedule_power
    //  nvU0,  //basic_schedule
    //  0,  //prior_dist
    //  0,  //prev_time
    //  prev_time,  //time
    //  time,  //next_time
    //  prev_stat,  //stat
    //  false,  //input_check
    //  false,  //out_process
    //  simpson_div,  //simpson_div
    //  init_par["tol_boundary"],  //tol_boundary
    //  init_par["tol_cost"]);  //tol_cost
    //d_par = work_test["par"];
    //d_char = work_test["char"];
    //vU_k = d_par["U_k"];
    //vcc = d_char["boundary"];
    //vU_k.erase(vU_k.begin());
    //vcc.erase(vcc.begin());
    //vU_k = (*prU_k)[kk];
    //vcc = (*prcc)[kk];
    rvec = (*prU_k)[kk];
    vU_k.resize(rvec.size());
    std::copy(rvec.begin(), rvec.end(), vU_k.begin());
    rvec = (*prcc)[kk];
    vcc.resize(rvec.size());
    std::copy(rvec.begin(), rvec.end(), vcc.begin());
      //for ( int kk = 0; kk < vU_k.size(); kk++ ) {
      //  Rcpp::Rcout << "# project_power # U_k[ " << kk << " ] = " << vU_k.at(kk) << "; boundary = " << vcc.at(kk) << std::endl;
      //}

    work_KK_1 = vU_k.size();
    work_KK = work_KK_1 - 1;
    std::vector<double> vU_k_(work_KK_1);
    std::vector<double> vcc_(work_KK_1);

    fin_kk = 1;
    d_s = vcc.at(1) - stat;
    d_t = vU_k.at(1) - time;
      //Rcpp::Rcout << "# project_power # next_stat = " << vcc.at(1) << "; stat = " << stat << "; d_s = " << d_s << std::endl;
      //Rcpp::Rcout << "# project_power # next_time = " << vU_k.at(1) << "; time = " << time << "; d_t = " << d_t << std::endl;
    *doubleVar = R::pnorm(-d_s, -effect_size * d_t, sqrt(d_t), 1, 0);
    doubleVar[2] = -R::qnorm(cond_power, 0, 1, 1, 0) * sqrt(d_t) + effect_size * d_t + stat;
    d_t = d_t + time - prev_time;
    doubleVar[3] = (((double) (cond_power < 0.5)) - 0.5) * 2; // sign 
    doubleVar[1] = doubleVar[3] * -R::qnorm(1e-8 , 0, 1, 1, 0) * sqrt(d_t) + effect_size * d_t + prev_stat;
      //Rcpp::Rcout << "sign = " << doubleVar[3] << "; sel = " << ((doubleVar[3] * (doubleVar[2] - doubleVar[1])) < 0) + 1 << std::endl;
    fin_cc = doubleVar[((doubleVar[3] * (doubleVar[2] - doubleVar[1])) < 0) + 1]; // prevent infinity
      //Rcpp::Rcout << "fin_cc = " << fin_cc << "; new = " << doubleVar[1] << "; old = " << doubleVar[2] << std::endl;
    //fin_cc = -R::qnorm(cond_power, 0, 1, 1, 0) * sqrt(d_t) + effect_size * d_t + stat;
      //Rcpp::Rcout << "# project_power # Stage " << kk << "; cc1 = " << vcc.at(1) << "; stat = " << stat << "; U_k[1] = " << vU_k.at(1) << "; U_k[0] = " << time << std::endl;
    pre_cond_power = 0;  // preserve
     //Rcpp::Rcout << "# project_power # Stage " << kk << "; temp. conditional power = " << *doubleVar << "; fin_cc = " << fin_cc << std::endl;
    if ( *doubleVar < cond_power ) {
      pre_cond_power = *doubleVar;  // preserve
      fin_kk += 1;
      for ( ; fin_kk < work_KK; fin_kk++ ) {
        // Check for fin_kk == work_KK is unnecessary because, even if *doubleVar < cond_power at fin_kk == work_KK,
        // then fin_kk = work_KK will be used.
        vU_k_ = vU_k;
        vU_k_.resize(fin_kk + 1);
        vcc_ = vcc;
        vcc_.resize(fin_kk + 1);
        *doubleVar = pr_rej_H0(&vU_k_, &vcc_, stat, effect_size, simpson_div);
         //Rcpp::Rcout << "# project_power # Stage " << kk << "; temp. conditional power = " << *doubleVar << std::endl;
        if ( *doubleVar >= cond_power ) break;
      pre_cond_power = *doubleVar;  // preserve
      }

      vU_k_ = vU_k;
      vU_k_.resize(fin_kk + 1);
      vcc_ = vcc;
      vcc_.resize(fin_kk + 1);
      struct arg_pr_rej_H0 str_arg_pr_rej_H0;
      str_arg_pr_rej_H0.vU_k = &(vU_k_);
      str_arg_pr_rej_H0.vcc = &(vcc_);
      str_arg_pr_rej_H0.stat = stat;
      str_arg_pr_rej_H0.effect_size = effect_size;
      str_arg_pr_rej_H0.simpson_div = simpson_div;
      str_arg_pr_rej_H0.cc_fin_kk = vcc_.data() + fin_kk;
        //Rcpp::Rcout << "# project_power # &(cc[ " << 0 << "]) = " << vcc_.data() << std::endl;
        //Rcpp::Rcout << "# project_power # &(cc[ " << fin_kk << "]) = " << vcc_.data() + fin_kk << std::endl;
       //Rcpp::Rcout << "# project_power # Stage " << kk << "; temp. conditional power = " << *doubleVar << std::endl;
      //*doubleVar = pr_rej_H0_sol1(15, &str_arg_pr_rej_H0);
      //*doubleVar = vcc.at(fin_kk);
      d_t = vU_k.at(fin_kk) - time;
        // prevent cond_power = 1;
        doubleVar[2] = 1 - 1e-8;
        doubleVar[3] = pr_rej_H0_sol1(R::qnorm(1e-8, 0, 1, 1, 0) * sqrt(d_t) + effect_size * d_t + stat, &str_arg_pr_rej_H0);
        doubleVar[2] = doubleVar[(doubleVar[3] <= doubleVar[2]) + 2];  // first selection: select smaller one
//Rcpp::Rcout << "doubleVar[2] = " << doubleVar[2] << std::endl;
        doubleVar[3] = cond_power;
        doubleVar[2] = doubleVar[(doubleVar[3] <= doubleVar[2]) + 2];  // second selection: select smaller one
//Rcpp::Rcout << "doubleVar[2] = " << doubleVar[2] << std::endl;
        doubleVar[3] = pr_rej_H0_sol1(-R::qnorm(1e-8, 0, 1, 1, 0) * sqrt(d_t) + effect_size * d_t + stat, &str_arg_pr_rej_H0);
        doubleVar[2] = doubleVar[(doubleVar[2] < doubleVar[3]) + 2];  // third selection: select larger one
//Rcpp::Rcout << "doubleVar[2] = " << doubleVar[2] << std::endl;
        doubleVar[3] = 1e-8;
        doubleVar[2] = doubleVar[(doubleVar[2] < doubleVar[3]) + 2];  // third selection: select larger one
        cond_power = doubleVar[2];
//Rcpp::Rcout << "modified cond_power = " << cond_power << std::endl;
        *doubleVar = vcc.at(fin_kk) * 2;
        doubleVar[1] = -R::qnorm(cond_power - pre_cond_power, 0, 1, 1, 0) * sqrt(d_t) + effect_size * d_t + stat;
        *doubleVar = doubleVar[(fin_kk == work_KK) + 0];
//Rcpp::Rcout << "cond_power - pre_cond_power = " << cond_power - pre_cond_power << std::endl;
//Rcpp::Rcout << "work_KK = " << work_KK << "; fin_kk = " << fin_kk << std::endl;
//Rcpp::Rcout << "cond_power = " << cond_power << std::endl;
//Rcpp::Rcout << "pr_rej_H0_sol1 (" << -100 << ") = " << pr_rej_H0_sol1(-100, &str_arg_pr_rej_H0) << std::endl;
//Rcpp::stop("complete");
      if ( (cond_power - pre_cond_power) <= 1e-8 ) {
        fin_kk -= 1;
        fin_cc = vcc.at(fin_kk);
        vU_k_ = vU_k;
        vU_k_.resize(fin_kk + 1);
        vcc_ = vcc;
        vcc_.resize(fin_kk + 1);
        str_arg_pr_rej_H0.vU_k = &(vU_k_);
        str_arg_pr_rej_H0.vcc = &(vcc_);
      }  else {
//Rcpp::Rcout << "pr_rej_H0_sol1 sol_l = " << *doubleVar / 2. << std::endl;
//Rcpp::Rcout << "pr_rej_H0_sol1 sol_u = " << *doubleVar << "; vcc.at(" << fin_kk << ") = " << vcc.at(fin_kk) << std::endl;
//Rcpp::Rcout << "pr_rej_H0_sol1 (" << *doubleVar / 2. << ") = "; Rcpp::Rcout << pr_rej_H0_sol1(*doubleVar / 2., &str_arg_pr_rej_H0) << std::endl;
//Rcpp::Rcout << "pr_rej_H0_sol1 (" << *doubleVar << ") = "; Rcpp::Rcout << pr_rej_H0_sol1(*doubleVar, &str_arg_pr_rej_H0) << std::endl;
        fin_cc = bisection_inverse(pr_rej_H0_sol1,
          cond_power, &str_arg_pr_rej_H0, *doubleVar / 2., *doubleVar,
          false, true, false, init_par["tol_boundary"]); // smaller
      }
      *doubleVar = pr_rej_H0_sol1(fin_cc, &str_arg_pr_rej_H0);
      //Rcpp::Rcout << "# project_power # Stage " << kk << "; temp. conditional power = " << *doubleVar << "; fin_cc = " << fin_cc << std::endl;
    } // ここまでOK

      // Projected working test;
      //vU_k = d_par["U_k"];
      //vcc = d_char["boundary"];
      vU_k.insert(vU_k.begin(), prev_time);
      vcc.insert(vcc.begin(), std::numeric_limits<double>::infinity());
        //for ( int kk = 0; kk < vU_k.size(); kk++ ) {
        //  Rcpp::Rcout << "# project_power # U_k[ " << kk << " ] = " << vU_k.at(kk) << "; boundary = " << vcc.at(kk) << std::endl;
        //}
      vU_k.resize(fin_kk + 1 + 1);
      vcc.resize(fin_kk + 1 + 1);
      vcc.at(fin_kk + 1) = fin_cc;
        //for ( int kk = 0; kk < vU_k.size(); kk++ ) {
        //  Rcpp::Rcout << "# project_power # U_k[ " << kk << " ] = " << vU_k.at(kk) << "; boundary = " << vcc.at(kk) << std::endl;
        //}

        //Rcpp::Rcout << "# project_power # prev_stat = " << prev_stat << std::endl;
        //Rcpp::Rcout << "# project_power # prev_time = " << prev_time << std::endl;
        //Rcpp::Rcout << "# project_power # fin_kk = " << fin_kk << std::endl;
        //Rcpp::Rcout << "# project_power # fin_cc = " << fin_cc << std::endl;
//if ( kk == 1 ) {
//  for ( int ii = 0; ii < vU_k.size(); ii++ ) {
//    Rcpp::Rcout << "vU_k(" << ii << ") = " << vU_k.at(ii) << "; boundary = " << vcc.at(ii) << std::endl;
//  }
//  Rcpp::Rcout << "prev_stat = " << prev_stat << "; effect_size = " << effect_size << "; simpson_div = " << simpson_div << std::endl;
//}
      cond_power = pr_rej_H0(&vU_k, &vcc, prev_stat, effect_size, simpson_div);
      //Rcpp::Rcout << "# project_power # Stage " << kk << "; cond power (effect_size = " << effect_size << ") = " << cond_power << std::endl;
  }
    //Rcpp::Rcout << "# project_power # effect_size = " << effect_size << "; conditional power = " << cond_power << "; fin_kk = " << fin_kk << "; fin_time = " << vU_k.at(1) << "; fin_cc = " << fin_cc << std::endl;

    //Rcpp::Rcout << "# project_power # effect_size = " << effect_size << "; conditional power = " << cond_power << std::endl;
    //Rcpp::Rcout << "# project_power # END" << std::endl;
  return cond_power;
}


static double pr_rej_H0_sol2(
const double effect_size,
struct arg_pr_rej_H0* parg
) {
  double prob = pr_rej_H0(parg->vU_k, parg->vcc, parg->stat, effect_size, parg->simpson_div);
  return prob;
}



// [[Rcpp::export]]
Rcpp::List exact_est_norm_c(
  Rcpp::List initial_test = 0,
  Rcpp::NumericVector times = 0,
  Rcpp::NumericVector stats = 0,
  Rcpp::NumericVector costs = 0,
  const bool final_analysis = true,
  const bool estimate = true,
  const double ci_coef = 0.95,
  const double tol_est = 1e-8,
  const bool input_check = true
  ) {
  bool sig_level_exhaustive = final_analysis;

  if ( input_check ) {
    if ( stats.size() != times.size() ) {
      Rcpp::stop("'times' and 'stats' should have the same length.");
    }
    if ( costs.size() > times.size() ) {
      Rcpp::stop("When 'costs' is specified, its length should be equal to or less than that of 'times'.");
    }
    if ( (ci_coef < 0) || (ci_coef > 1) ) {
      Rcpp::stop("'ci_coef' should be a value in (0, 1).");
    }
    if ( tol_est <= 0 ) {
      Rcpp::stop("'tol_est' should be positive.");
    }
  }

  if ( times.at(0) > 0 ) {
    times.insert(times.begin(), 0);
    stats.insert(stats.begin(), 0);
    costs.insert(costs.begin(), 0);
  }
  int max_kk = times.size();
  std::vector<double> vtimes(max_kk);
  std::vector<double> vstats(max_kk);
  std::vector<double> vcost0(max_kk);
  std::copy(times.begin(), times.end(), vtimes.begin());
  std::copy(stats.begin(), stats.end(), vstats.begin());
  std::copy(costs.begin(), costs.end(), vcost0.begin());
  std::vector<double> vcond_alpha(max_kk);
  std::vector<double> vboundary(max_kk);

  if ( input_check ) {
    if ( vtimes.at(0) < 0 ) {
      Rcpp::stop("All values of 'times' should be non-negative.");
    }
    for ( int ii = 1; ii < (int) vtimes.size(); ii++ ) {
      if ( vtimes.at(ii) < 0 ) {
        Rcpp::stop("All values of 'times' should be non-negative.");
      }
      if ( (vtimes.at(ii) - vtimes.at(ii - 1)) < 0 ) {
        Rcpp::stop("All intervals of 'times' should be positive.");
      }
    }
    if ( *std::min_element(vcost0.begin(), vcost0.end()) < 0 ) {
      Rcpp::stop("All values of 'costs' should be non-negative.");
    }
  }

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Settings                                                                  //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  //=== Dummy Box ===//
  double doubleVar[8] = {};
  int intVar[8] = {};
  Rcpp::NumericVector rvec;

  //=== Initial Working Test ===//
  Rcpp::List init_par = initial_test["par"];
  Rcpp::List init_char = initial_test["char"];
  Rcpp::List d_par;
  Rcpp::List d_char;
  Rcpp::List rU_k(max_kk);
  Rcpp::List rcc(max_kk);
  rU_k[0] = init_par["U_k"];
  rcc[0] = init_char["boundary"];
  vcost0.at(0) = init_par["cost0"];

  // init_test should be prev_time = time = next_time = 0;
  doubleVar[0] = init_par["prev_time"];
  doubleVar[1] = init_par["time"];
  doubleVar[2] = init_par["next_time"];
  if( !((doubleVar[0] == 0) && (doubleVar[1] == 0) && (doubleVar[2] == 0)) ) {
    Rcpp::stop("Initial working test should be loaded through the argument 'inital_test'.");
  }

  std::vector<double> vU00 = init_par["U_0"];
  std::vector<double> vU0 = vU00;
  int work_KK_1 = vU0.size();
  //int work_KK = work_KK_1 - 1;
  //double* U0 = vU0.data();
  Rcpp::NumericVector nvU0(work_KK_1); // std::vector to Rcpp::NumericVector
  std::copy(vU0.begin(), vU0.end(), nvU0.begin());
    //Rcpp::Rcout << "# exact_est_norm_c # Check 1 [START]" << std::endl;
    //for ( int kkk = 0; kkk < (vU0.end() - vU0.begin()); kkk++ ) {
    //  *doubleVar = nvU0.at(kkk);
    //}
    //Rcpp::Rcout << "# exact_est_norm_c # Check 1 [END]" << std::endl;
    //for ( int kk = 0; kk < work_KK_1; kk++ ) {
    //  Rcpp::Rcout << "nvU0.at(" << kk << ") = " << nvU0.at(kk) << std::endl;
    //}
  Rcpp::NumericVector nvprior(11);
  nvprior = init_par["prior"];

  //=== Variables ===//
  int kk;
  double cond_alpha;
  double prev_time;
  double time;
  double next_time;
  double stat;
  //double* cost0 = vcost0.data();
  //*cost0 = init_par["cost0"];
    //Rcpp::Rcout << "# exact_est_norm_c # Stage " << 0 << "; cost0 = " << vcost0.at(0) << std::endl;

  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Analysis                                                                  //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  cond_alpha = init_par["overall_sig_level"];
  vcond_alpha.at(0) = cond_alpha;
  vboundary.at(0) = std::numeric_limits<double>::infinity();
  kk = 1;
  time = 0;
  stat = 0;
    //Rcpp::Rcout << "# exact_est_norm_c # Stage " << 0 << "; Cond err = " << cond_alpha << std::endl;
  for ( ; kk < (max_kk - sig_level_exhaustive); kk++ ) {
    Rcpp::List work_test;
    if ( vcost0.at(kk) == 0 ) {
      //test1 <- work_test_norm_c(time = 0, next_time=9.2, cost_type_1_err = 0, min_effect_size = -log(0.65), simpson_div = 4)
      prev_time = 0;
      time = vtimes.at(kk - 1);
      next_time = vtimes.at(kk);
      stat = vstats.at(kk - 1);
        //Rcpp::Rcout << "# exact_est_norm_c # Stage " << kk << "; prev_time = " << prev_time << "; time = " << time << "; next_time = " << next_time << "; stat = " << stat << std::endl;
      work_test = work_test_norm_c(
        init_par["overall_sig_level"],  //overall_sig_level
        init_par["work_beta"],  //work_beta
        cond_alpha,  //cond_alpha
        0,  //cost_type_1_err
        0,  //cost_type_2_err
        0,  //prev_cost
        init_par["min_effect_size"],  //min_effect_size
        0,  //effect_size
        0,  //basic_schedule_num
        0,  //basic_schedule_power
        nvU0,  //basic_schedule
        nvprior,  //prior_dist
        prev_time,  //prev_time
        time,  //time
        next_time,  //next_time
        stat,  //stat
        false,  //input_check
        false,  //out_process
        init_par["simpson_div"],  //simpson_div
        init_par["tol_boundary"],  //tol_boundary
        init_par["tol_cost"]);  //tol_cost
      d_par = work_test["par"];
      d_char = work_test["char"];
      vcost0.at(kk) = d_par["cost0"];
        //Rcpp::Rcout << "# exact_est_norm_c # Stage " << kk << "; cost0 = " << vcost0.at(kk) << std::endl;
    }

    //result1 <- work_test_norm_c(prev_time = 0, time=9.2, stat=4.35, cost_type_1_err = test1$par$cost0, min_effect_size = -log(0.65), simpson_div = 4)
    prev_time = vtimes.at(kk - 1);
    time = vtimes.at(kk);
    next_time = 0;
    stat = vstats.at(kk);
      //Rcpp::Rcout << "# exact_est_norm_c # Stage " << kk << "; prev_time = " << prev_time << "; time = " << time << "; next_time = " << next_time << "; stat = " << stat << std::endl;
    work_test = work_test_norm_c(
      init_par["overall_sig_level"],  //overall_sig_level
      init_par["work_beta"],  //work_beta
      0,  //cond_alpha
      vcost0.at(kk),  //cost_type_1_err
      0,  //cost_type_2_err
      0,  //prev_cost
      init_par["min_effect_size"],  //min_effect_size
      0,  //effect_size
      0,  //basic_schedule_num
      0,  //basic_schedule_power
      nvU0,  //basic_schedule
      nvprior,  //prior_dist
      prev_time,  //prev_time
      time,  //time
      next_time,  //next_time
      stat,  //stat
      false,  //input_check
      false,  //out_process
      init_par["simpson_div"],  //simpson_div
      init_par["tol_boundary"],  //tol_boundary
      init_par["tol_cost"]);  //tol_cost
    d_par = work_test["par"];
    d_char = work_test["char"];
    rU_k[kk] = d_par["U_k"];
    rcc[kk] = d_char["boundary"];

    cond_alpha = d_char["pr_rej_H0"];
      //Rcpp::Rcout << "# exact_est_norm_c # Stage " << kk << "; Cond err = " << cond_alpha << std::endl;
    vcond_alpha.at(kk) = cond_alpha;
    rvec = d_char["boundary"]; // U0 as an dummy;
    vU0.resize(rvec.size());
    std::copy(rvec.begin(), rvec.end(), vU0.begin());
      //Rcpp::Rcout << "# exact_est_norm_c # Check 2 [START]" << std::endl;
      //for ( int kkk = 0; kkk < (rvec.end() - rvec.begin()); kkk++ ) {
      //  *doubleVar = vU0.at(kkk);
      //}
      //Rcpp::Rcout << "# exact_est_norm_c # Check 2 [END]" << std::endl;
    vboundary.at(kk) = vU0.at(0);

    if ( (cond_alpha >= 1) || (kk == (max_kk - sig_level_exhaustive - 1)) ) break;
  }
  if ( (cond_alpha < 1) && sig_level_exhaustive ) {
    kk += (max_kk > 2);
    prev_time = vtimes.at(kk - 1);
    time = vtimes.at(kk);
    next_time = 0;
    stat = vstats.at(kk);
      //Rcpp::Rcout << "# exact_est_norm_c # Stage " << kk << "; prev_time = " << prev_time << "; time = " << time << "; next_time = " << next_time << "; stat = " << stat << std::endl;
    //*doubleVar = (stat - vstats.at(kk - 1)) / sqrt(time - prev_time);  // z-statistic
    //doubleVar[1] = -R::qnorm(cond_alpha, 0, 1, 1, 0);  // conditinal critical value
    //cond_alpha = *doubleVar >= doubleVar[1];
      //Rcpp::Rcout << "# exact_est_norm_c # z-stat = " << *doubleVar << "; adaptive critical value = " << doubleVar[1] << std::endl;
    *doubleVar = -R::qnorm(cond_alpha, 0, 1, 1, 0) * sqrt(time - prev_time) + vstats.at(kk - 1);
    cond_alpha = stat >= *doubleVar;
      //Rcpp::Rcout << "# exact_est_norm_c # Stage " << kk << "; Cond err = " << cond_alpha << std::endl;
    vcond_alpha.at(kk) = cond_alpha;
    vboundary.at(kk) = *doubleVar;
  }
  vtimes.resize(kk + 1);
  vstats.resize(kk + 1);
  vcost0.resize(kk + 1);
  vcond_alpha.resize(kk + 1);
  vboundary.resize(kk + 1);
  std::vector<bool> vrej_H0(kk + 1);
  std::fill(vrej_H0.begin(), vrej_H0.end(), false);
  vrej_H0.at(kk) = (cond_alpha >= 1);

  int stp_kk = kk;
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  // Projection                                                                //
  //-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-//
  double pval, mue; // p-value, median unbiased est
  std::vector<double> vest(2); // , lower lim, upper lim
  double lbb_pval, lbb_mue; // lower bound-based est
  std::vector<double> vlbb_est(2); // lower bound-based est
  struct base_test str_test;
  if ( estimate ) {
      //Rcpp::Rcout << "# exact_est_norm_c # Stage " << stp_kk << "; prev_time = " << prev_time << "; time = " << time << "; next_time = " << next_time << "; stat = " << stat << std::endl;

    str_test.test = &initial_test;
    str_test.times = &vtimes;
    str_test.stats = &vstats;
    str_test.sig_level_exhaustive = sig_level_exhaustive;
    str_test.ci_coef = ci_coef;
    str_test.cost0 = vcost0.data();
    str_test.U_k = &rU_k;
    str_test.cc = &rcc;

      //for( int ii = 0; ii < vtimes.size(); ii++ ) {
      //  Rcpp::Rcout << "vtimes.at(" << ii << ") = " << vtimes.at(ii) << std::endl;
      //}
      //for( int ii = 0; ii < vstats.size(); ii++ ) {
      //  Rcpp::Rcout << "vstats.at(" << ii << ") = " << vstats.at(ii) << std::endl;
      //}
      //  Rcpp::Rcout << "sig_level_exhaustive = " << sig_level_exhaustive << std::endl;
      //  Rcpp::Rcout << "ci_coef = " << ci_coef << std::endl;
      //for( int ii = 0; ii < vcost0.size(); ii++ ) {
      //  Rcpp::Rcout << "vcost0.at(" << ii << ") = " << vcost0.at(ii) << std::endl;
      //}
      //for( int kkk = 0; kkk < stp_kk; kkk++ ) {
      //Rcpp::NumericVector vvv = rU_k[kkk];
      //for( int ii = 0; ii < vvv.size(); ii++ ) {
      //  Rcpp::Rcout << "rU_k[" << kkk << "].at(" << ii << ") = " << vvv.at(ii) << std::endl;
      //}}
      //for( int kkk = 0; kkk < stp_kk; kkk++ ) {
      //Rcpp::NumericVector vvv = rcc[kkk];
      //for( int ii = 0; ii < rcc.size(); ii++ ) {
      //  Rcpp::Rcout << "rcc[" << kkk << "].at(" << ii << ") = " << vvv.at(ii) << std::endl;
      //}}

    double ciVar[8];
    *ciVar = (1. - ci_coef) / 2.;
      //Rcpp::Rcout << "# exact_est_norm_c # ci_coef = " << ci_coef << std::endl;

    //double* est = vest.data();
    mue = stat / time; // temporarily MLE
    ciVar[2] = sqrt(1. / time); // naive SE for MLE
    ciVar[3] = -R::qnorm(*ciVar, 0, 1, 1, 0);
    ciVar[5] = ciVar[3] * ciVar[2]; // CI wing
      //Rcpp::Rcout << "# exact_est_norm_c # ciVar[0] = " << ciVar[0] << std::endl;
      //Rcpp::Rcout << "# exact_est_norm_c # ciVar[2] = " << ciVar[2] << std::endl;
      //Rcpp::Rcout << "# exact_est_norm_c # ciVar[3] = " << ciVar[3] << std::endl;
      //Rcpp::Rcout << "# exact_est_norm_c # ciVar[5] = " << ciVar[5] << std::endl;

    if ( stp_kk == 1 ) {  // "max_kk == 2" -> "stp_kk == 1"
      pval = R::pnorm(-stat / sqrt(time), 0, 1, 1, 0);
      // mue = MLE;
      vest.at(0) = mue - ciVar[5];
      vest.at(1) = mue + ciVar[5];

      lbb_pval = pval;
      lbb_mue = mue;
      vlbb_est = vest;
    } else {
      ciVar[1] = project_power(mue, &str_test); // p-value at the MLE
      ciVar[4] = mue + -R::qnorm(ciVar[1], 0, 1, 1, 0) * ciVar[2];
        //Rcpp::Rcout << "# exact_est_norm_c # Temp. Lower Conf Lim = " << ciVar[4] - ciVar[5] << std::endl;
        //Rcpp::Rcout << "# exact_est_norm_c # Temp. Med Unbias Est = " << ciVar[4] << std::endl;
        //Rcpp::Rcout << "# exact_est_norm_c # Temp. Upper Conf Lim = " << ciVar[4] + ciVar[5] << std::endl;

  // Lower
  //Rcpp::Rcout << "effect_size = " << ciVar[4] - ciVar[5] * 3/2 << std::endl;
  //Rcpp::Rcout << "pval = " << project_power(ciVar[4] - ciVar[5] * 3/2, &str_test) << std::endl;
  //Rcpp::Rcout << "effect_size = " << ciVar[4] - ciVar[5] / 2. << std::endl;
  //Rcpp::Rcout << "pval = " << project_power(ciVar[4] - ciVar[5] / 2., &str_test) << std::endl;
  // Upper
  //Rcpp::Rcout << "effect_size = " << ciVar[4] + ciVar[5] * 3 / 2. << std::endl;
  //Rcpp::Rcout << "pval = " << project_power(ciVar[4] + ciVar[5] * 3 / 2., &str_test) << std::endl;
  //Rcpp::stop("complete");
   Rcpp::Rcout << "Computing exact estimates." << std::endl;
     //Rcpp::Rcout << "# exact_est_norm_c # Exact inference [START]" << std::endl;
     //Rcpp::Rcout << "# exact_est_norm_c # Computing exact P-value" << std::endl;
      pval = project_power(0, &str_test); // p-value at the null
      // Median Unbiased Est

     //Rcpp::Rcout << "# exact_est_norm_c # Computing exact median unbiased estimator" << std::endl;
      mue = bisection_inverse(project_power,
        0.5, &str_test, ciVar[4] - ciVar[5] / 2., ciVar[4] + ciVar[5] / 2.,
        false, true, false, tol_est); // smaller
      // Lower Limit
     //Rcpp::Rcout << "# exact_est_norm_c # Computing exact lower confidence limit" << std::endl;
      vest.at(0) = bisection_inverse(project_power,
        *ciVar, &str_test, ciVar[4] - ciVar[5] * 3/2, mue,
        true, false, false, tol_est); // larger
      // Upper Limit
     //Rcpp::Rcout << "# exact_est_norm_c # Computing exact upper confidence limit" << std::endl;
      vest.at(1) = bisection_inverse(project_power,
        1 - *ciVar, &str_test, mue, ciVar[4] + ciVar[5] * 3 / 2.,
        false, true, false, tol_est); // smaller
        //Rcpp::Rcout << "# exact_est_norm_c # P-value        = " << pval   << std::endl;
        //Rcpp::Rcout << "# exact_est_norm_c # Lower Conf Lim = " << est[0] << std::endl;
        //Rcpp::Rcout << "# exact_est_norm_c # Med Unbias Est = " << mue    << std::endl;
        //Rcpp::Rcout << "# exact_est_norm_c # Upper Conf Lim = " << est[1] << std::endl;
     //Rcpp::Rcout << "# exact_est_norm_c # Exact inference [END]" << std::endl;

      {
   Rcpp::Rcout << "Computing approximate estimates." << std::endl;
     //Rcpp::Rcout << "# exact_est_norm_c # Approx inference [START]" << std::endl;
          //Rcpp::Rcout << "# exact_est_norm_c # stp_kk(max stage excluding t=0) = " << stp_kk << std::endl;
        time = vtimes.at(stp_kk - 1);

        vU0 = vU00;
        int work_KK_1 = vU0.size();
        int work_KK = work_KK_1 - 1;

        // Current stage //
        int time_kk = 0;
        for ( int kk = 1; kk < work_KK_1; kk++ ) {
          //*doubleVar = (U0[kk - 1] + U0[kk]) / 2.; // + 0.99 * (U0[kk + 1] - U0[kk]);
          *doubleVar = (vU0.at(kk - 1) + vU0.at(kk)) / 2.; //..
          time_kk += (*doubleVar < time);
        }
          //Rcpp::Rcout << "# sample_size_norm_c # time: " << time << "; time_kk: " << time_kk << std::endl;
        *intVar = (time_kk == 0) && (time > 0);
        intVar[1] = (time_kk == work_KK) && (time < vU0.at(work_KK));
        if ( (*intVar + intVar[1]) > 0 ) {
          vU0.insert(vU0.begin() + *intVar + intVar[1] * work_KK, 0);
        }
        work_KK_1 = vU0.size(); // added at Aug 20, 2018
        time_kk += *intVar;
        vU0.at(time_kk) = time;
        Rcpp::NumericVector nvU0(work_KK_1); // std::vector to Rcpp::NumericVector
        std::copy(vU0.begin(), vU0.end(), nvU0.begin());
          //Rcpp::Rcout << "# exact_est_norm_c # Check 3 [START]" << std::endl;
          //for ( int kkk = 0; kkk < (vU0.end() - vU0.begin()); kkk++ ) {
          //  *doubleVar = nvU0.at(kkk);
          //}
          //Rcpp::Rcout << "# exact_est_norm_c # Check 3 [END]" << std::endl;

        std::vector<double> vU_k(time_kk + 1);
        std::copy(vU0.begin(), vU0.begin() + time_kk + 1, vU_k.begin());
          //Rcpp::Rcout << "# exact_est_norm_c # Check 4 [START]" << std::endl;
          //for ( int kkk = 0; kkk < (time_kk + 1); kkk++ ) {
          //  *doubleVar = vU_k.at(kkk);
          //}
          //Rcpp::Rcout << "# exact_est_norm_c # Check 4 [END]" << std::endl;

        vU_k.push_back(vtimes.at(stp_kk));
        work_KK = vU_k.size() - 1;
        work_KK_1 = work_KK + 1;
          //for ( int ii = 0; ii < vU_k.size(); ii++ ) {
          //  Rcpp::Rcout << "# sample_size_norm_c # vU_k.at(" << ii << ") = " << vU_k.at(ii) << std::endl;
          //}

        Rcpp::List work_test;
        work_test = work_test_norm_c(
          init_par["overall_sig_level"],  //overall_sig_level
          init_par["work_beta"],  //work_beta
          init_par["overall_sig_level"],  //cond_alpha
          0,  //cost_type_1_err
          0,  //cost_type_2_err
          0,  //prev_cost
          init_par["min_effect_size"],  //min_effect_size
          0,  //effect_size
          0,  //basic_schedule_num
          0,  //basic_schedule_power
          nvU0,  //basic_schedule
          nvprior,  //prior_dist
          0,  //prev_time
          0,  //time
          0,  //next_time
          0,  //stat
          false,  //input_check
          false,  //out_process
          init_par["simpson_div"],  //simpson_div
          init_par["tol_boundary"],  //tol_boundary
          init_par["tol_cost"]);  //tol_cost
        d_par = work_test["par"];
        d_char = work_test["char"];
        std::vector<double> vwcc = d_char["boundary"];
        std::vector<double> vcc(time_kk + 2);
        std::copy(vwcc.begin(), vwcc.begin() + time_kk + 1, vcc.begin());
          //Rcpp::Rcout << "# exact_est_norm_c # Check 5 [START]" << std::endl;
          //for ( int kkk = 0; kkk < (time_kk + 1); kkk++ ) {
          //  *doubleVar = vcc.at(kkk);
          //}
          //Rcpp::Rcout << "# exact_est_norm_c # Check 5 [END]" << std::endl;
        vcc.at(time_kk + 1) = vstats.at(stp_kk);
          //for ( int ii = 0; ii < vcc.size(); ii++ ) {
          //  Rcpp::Rcout << "# sample_size_norm_c # vU_k.at(" << ii << ") = " << vU_k.at(ii) << "; vcc.at() = " << vcc.at(ii) << std::endl;
          //}

        struct arg_pr_rej_H0 str_arg_pr_rej_H0;
        str_arg_pr_rej_H0.vU_k = &vU_k;
        str_arg_pr_rej_H0.vcc = &vcc;
        str_arg_pr_rej_H0.stat = 0;
        str_arg_pr_rej_H0.simpson_div = init_par["simpson_div"];

         // P-value
          //Rcpp::Rcout << "# exact_est_norm_c # Computing lower bound-based P-value" << std::endl;
        lbb_pval = pr_rej_H0_sol2(0, &str_arg_pr_rej_H0); // p-value at the null
          //Rcpp::Rcout << "# exact_est_norm_c # Computing lower bound-based P-value: " << lbb_pval << std::endl;

         // Median Unbiased Est
          //Rcpp::Rcout << "# exact_est_norm_c # Computing lower bound-based median unbiased estimator" << std::endl;
         lbb_mue = bisection_inverse(pr_rej_H0_sol2,
           0.5, &str_arg_pr_rej_H0, ciVar[4] - ciVar[5] / 2., ciVar[4] + ciVar[5] / 2.,
           false, true, false, tol_est); // smaller
          //Rcpp::Rcout << "# exact_est_norm_c # Computing lower bound-based median unbiased estimator: " << lbb_mue << std::endl;
         // Lower Limit
          //Rcpp::Rcout << "# exact_est_norm_c # Computing lower bound-based lower confidence limit" << std::endl;
         vlbb_est.at(0) = bisection_inverse(pr_rej_H0_sol2,
           *ciVar, &str_arg_pr_rej_H0, ciVar[4] - ciVar[5] * 3/2, lbb_mue,
           true, false, false, tol_est); // larger
          //Rcpp::Rcout << "# exact_est_norm_c # Computing lower bound-based lower confidence limit: " << vlbb_est.at(0) << std::endl;
         // Upper Limit
          //Rcpp::Rcout << "# exact_est_norm_c # Computing lower bound-based upper confidence limit" << std::endl;
         vlbb_est.at(1) = bisection_inverse(pr_rej_H0_sol2,
           1 - *ciVar, &str_arg_pr_rej_H0, lbb_mue, ciVar[4] + ciVar[5] * 3 / 2.,
           false, true, false, tol_est); // smaller
          //Rcpp::Rcout << "# exact_est_norm_c # Computing lower bound-based upper confidence limit: " << vlbb_est.at(1) << std::endl;

          //Rcpp::Rcout << "# exact_est_norm_c # P-value        = " << lbb_pval   << std::endl;
          //Rcpp::Rcout << "# exact_est_norm_c # Lower Conf Lim = " << vlbb_est.at(0) << std::endl;
          //Rcpp::Rcout << "# exact_est_norm_c # Med Unbias Est = " << lbb_mue    << std::endl;
          //Rcpp::Rcout << "# exact_est_norm_c # Upper Conf Lim = " << vlbb_est.at(1) << std::endl;
          //Rcpp::Rcout << "# exact_est_norm_c # Approx inference [END]" << std::endl;
      }
    }
  }

  *doubleVar = init_par["overall_sig_level"];
  doubleVar[1] = init_par["min_effect_size"];

  d_par = Rcpp::List::create(
    Rcpp::Named("overall_sig_level") = *doubleVar,
    Rcpp::Named("min_effect_size") = doubleVar[1],
    Rcpp::Named("analyses") = vtimes.size() - 1,
    Rcpp::Named("times") = vtimes,
    Rcpp::Named("stats") = vstats,
    Rcpp::Named("final_analysis") = final_analysis
  );
  Rcpp::NumericVector rcost0;
  rcost0 = vcost0;
  if ( final_analysis && (vtimes.size() == max_kk) ) rcost0.at(vtimes.size() - 1) = NA_REAL; // '&& (vtimes.size() == max_kk)' was added at Aug 21, 2018.
  d_char = Rcpp::List::create(
    Rcpp::Named("cost0") = rcost0,
    Rcpp::Named("boundary") = vboundary,
    Rcpp::Named("cond_type_I_err") = vcond_alpha,
    Rcpp::Named("rej_H0") = vrej_H0
  );

  Rcpp::List d_out;
  if ( estimate ) {
    Rcpp::List d_est = Rcpp::List::create(
      Rcpp::Named("p_value") = pval,
      Rcpp::Named("ci_coef") = ci_coef,
      Rcpp::Named("median_unbiased") = mue,
      Rcpp::Named("conf_limits") = vest
      );
    Rcpp::List d_lbb_est = Rcpp::List::create(
      Rcpp::Named("p_value") = lbb_pval,
      Rcpp::Named("ci_coef") = ci_coef,
      Rcpp::Named("median_unbiased") = lbb_mue,
      Rcpp::Named("conf_limits") = vlbb_est
      );

    d_out = Rcpp::List::create(
      Rcpp::Named("par") = d_par,
      Rcpp::Named("char") = d_char,
      Rcpp::Named("est") = d_est,
      Rcpp::Named("lbb_est") = d_lbb_est
      );
  } else {
    d_out = Rcpp::List::create(
      Rcpp::Named("par") = d_par,
      Rcpp::Named("char") = d_char
      );
  }

  return d_out;

}
