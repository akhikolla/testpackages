#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <vector>
using namespace std;
using namespace Rcpp;


NumericVector func (NumericVector para, Function fun) {
  NumericVector loss_i_temp = fun(para);
  return loss_i_temp;
}

NumericVector var_funcc (NumericVector para_0, int fun_length, NumericVector rf) {
  NumericVector ret_var_func(fun_length);
  for(int k = 0; k < (fun_length); k++) {
    ret_var_func[k] = para_0[k] + (round(R::runif(0, rf[k]) * 100) / 100) * ((R::rbinom(1, 0.5) * -2) + 1);
  }
  return ret_var_func;
}

// [[Rcpp::export]]
List main_loop (double temp, double t_min, double r, int fun_length, int nlimit, NumericVector para_0, NumericVector para_i, Function var_func, bool vf_user,
                bool trace, NumericVector rf, NumericVector lower, NumericVector upper, Function fun, double loss_0, double k, double loss_opt, NumericVector para_opt,
                bool dyn_rf, double maxgood, double ac_acc, int stopac, bool maximization) {
  // Initializating variables
  IntegerVector n_oob(fun_length);
  int n_outer = 0;
  int savei = 0;
  double savet = 0;
  int ac = 0;
  double loss_i = 0;
  bool firsttemp = true;
  // Init trace
  vector<int> trace_n_outer;
  vector<double> trace_loss;
  vector<double> row;
  vector< vector<double> > trace_para;
  vector<int> trace_n_inner;
  vector<double> trace_temp;
  vector<int> trace_goodcounter;
  vector<double> row_rf;
  vector< vector<double> > trace_rf;

  // The outer while loop: Number of repetitions depends on cooling function and the temp. limit.
  while (temp > t_min) {
    // Initializing and resetting variables inside the while loop
    int goodcounter = 0;
    int n_inner = 0;
    n_outer++;
    std::fill(n_oob.begin(), n_oob.end(), 0);

    // If the tzemperature is < 2 for the first time, the temp. optimum is overwritten by the global optimum. Since traps cannot be left (practically) for
    if(temp <= 2 && firsttemp) {
      // tbd. r could be reduced at this point
      firsttemp = false;
      loss_0 = loss_opt;
      para_0 = para_opt;
    }

    for (int i = 0; i < nlimit; i++) { // Inner loop, no. of repetitions depends on the break criteria or on nlimit if no break criterion stops the loop.
      // Changing the variables
      n_inner++;
      if(!vf_user){ // Variation of the variables...
        para_i = var_funcc(para_0, fun_length, rf); // ...by the default function
      } else {
        para_i = var_func(para_0, fun_length, rf, temp); // ...by a user declared function. This is an SEXP. The algorithm is therefore much slower with it.
      }

      // Counting the variables which are out of bounds and change them.
      for(int j = 0; j < fun_length; j++) {
        if(para_i[j] < lower[j] || para_i[j] > upper[j]) {
          n_oob[j]++;
          // Generate new values for the variable until it is within the boundaries.
          int emergency_stop = 0;
          while (para_i[j] < lower[j] || para_i[j] > upper[j]) {
            emergency_stop++;
            NumericVector temp_para_i(1);

            if(!vf_user) { // Variation of the variables.
              NumericVector para_0_j(1);
              para_0_j = para_0[j];
              NumericVector rf_j(1);
              rf_j = rf[j];
              temp_para_i = var_funcc(para_0_j, 1, rf_j); // By the default function
              //temp_para_i = var_func(para_0[j], 1, rf[j]);
            } else {
              temp_para_i = var_func(para_0[j], 1, rf[j], temp); // By a user declared function. This is an SEXP. The algorithm is therefore much slower with it.
            }
            // NumericVector temp_para_i = var_func(para_0[i], 1, rf[i]); // MUST BE UPDATED: C FUN NEEDED

            para_i[j] = temp_para_i[0];
            if (emergency_stop > 1000){stop("The restrictions cannot be hold. Try different combination of starting values, boundaries or random factor.");}
          }

        }

      }

      //NumericVector loss_i_temp = fun(para_i);
      NumericVector loss_i_temp = func(para_i, fun);

      loss_i = loss_i_temp[0];
      double delta = loss_i - loss_0;
      // Check, if the loss has improved
      if(maximization) {delta = delta * -1;}
      if (delta < 0) {
        loss_0 = loss_i;
        para_0 = para_i;
      } else{ // This is the difference between Sim. Ann. and other Algorithms. It ist the prob. of accepting the worse loss.
        // If a loss_i is not defined (e. g. due to restrictions of the loss function [NA in ther R function]), the if cannot be true
        // loss_0 and para_0 are thus never updated with undefined values.

        if (R::runif(0, 1) < exp (- fabs (delta) / (k * temp) )) {
          loss_0 = loss_i;
          para_0 = para_i;
        }
      }
      if(maximization) {
        if (loss_0 > loss_opt) {
          goodcounter++;
          loss_opt = loss_0;
          para_opt = para_0;
          savei = n_outer;
          savet = temp;
        }
      } else {
        if (loss_0 < loss_opt) {
          goodcounter++;
          loss_opt = loss_0;
          para_opt = para_0;
          savei = n_outer;
          savet = temp;
        }
      }

      // Check for break criterions.
      if (goodcounter > maxgood) {break;}
      if (fabs(loss_0 - loss_opt) < ac_acc){
        ac++;

      }else{
        ac = 0;
      }
    if (ac > stopac){break;}
    } // End of the inner loop.

    if(trace) {
      trace_n_outer.push_back(n_outer);
      trace_loss.push_back(loss_0);

      for(k = 0; k < fun_length; k++) {
        row.push_back((double)para_0[k]);
        row_rf.push_back((double)rf[k]);
      }
      trace_para.push_back(row);
      row.erase (row.begin(), row.end());
      trace_n_inner.push_back(n_inner);
      trace_temp.push_back(temp);
      trace_goodcounter.push_back(goodcounter);
      trace_rf.push_back(row_rf);
      row_rf.erase (row_rf.begin(), row_rf.end());
    }

    temp = temp * r; // Temperature reduction.

    // Calculation of rf for the next iteration step according to the ratio of random values out of bounds (Corana et al. 1987).
    if (dyn_rf == true) {
      NumericVector ratio_noob (n_oob.size());
      for(int j = 0; j < n_oob.size(); j++) {
        ratio_noob[j] = ( (double) n_inner - (double) n_oob[j]) / (double) n_inner;
        if (ratio_noob[j] < 0.4 || ratio_noob[j] > 0.6) {
          if (ratio_noob[j] < 0.4) {
            rf[j] = rf[j] * (1.0 / (1.0 + 2.0 * ((0.4 - (double) ratio_noob[j]) / 0.4)));
          } else {
            rf[j] = rf[j] * (1.0 + (2.0 * (( (double) ratio_noob[j] - 0.6) / 0.4)));
          }
        }


      }

      // Downscaling of the rf makes the algorithm more efficient (Pronzato et al. 1984)
      // Could be relative instead of 5
      NumericVector ds (n_oob.size());
      for(int j = 0; j < n_oob.size(); j++){
        if (n_outer <= 5) {
          ds[j] = 1.0;
       }else {
          ds[j] = 1.0 / ( (double) n_outer  / 5.0);
        }
        if (rf[j] * ds[j] <= 0.1) {
          rf[j] = 0.1;
        }else{
          rf[j] = rf[j] * ds[j];
        }
     }
    }
  }


  List ret;
  ret["savei"] = savei;
  ret["savet"] = savet;
  ret["loss_opt"] = loss_opt;
  ret["para_opt"] = para_opt;
  ret["n_outer"] = trace_n_outer;
  ret["loss"] = trace_loss;
  ret["para"] = trace_para;
  ret["n_inner"] = trace_n_inner;
  ret["temp"] = trace_temp;
  ret["goodcounter"] = trace_goodcounter;
  ret["rf"] = trace_rf;

  return ret;
}
