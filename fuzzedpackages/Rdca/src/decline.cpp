
#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// ********** Prediction Models **********

// [[Rcpp::export]]
arma::mat exponential(Rcpp::List lst, arma::vec time) {

    //  Di: nominal decline rate
    int lr = time.size();
    int ll = lst.length();
    double qi = as<double>(lst["qi"]);
    double Di = as<double>(lst["Di"]);
    double b = as<double>(lst["b"]);
    double q_abnd = 0.0;
    int lc;
    if (ll == 6) {
        lc = 5;
    } else {
        q_abnd = as<double>(lst["q_abnd"]);
        lc = 7;
    }
    arma::mat results_table(lr,lc);
    if (lc == 5) {
        for (int i = 0; i < lr; i++) {
            results_table(i,0) = time(i);
            results_table(i,1) = qi * exp(-Di * time(i));
            results_table(i,2) = (qi -  results_table(i,1)) / Di;
            results_table(i,3) = Di;
            results_table(i,4) = b;
        }
    }
    if (lc == 7) {
        for (int i = 0; i < lr; i++) {
            results_table(i,0) = time(i);
            results_table(i,1) = qi * exp(-Di * time(i));
            results_table(i,2) = (qi -  results_table(i,1)) / Di;
            results_table(i,3) = Di;
            results_table(i,4) = b;
            results_table(i,5) = -(log(q_abnd / qi)) / Di;
            results_table(i,6) = (qi - q_abnd) / Di;
        }
    }
    return(results_table);
}


// [[Rcpp::export]]
arma::mat harmonic(Rcpp::List lst, arma::vec time) {

    //  Di: nominal decline rate
    int lr = time.size();
    int ll = lst.length();
    double qi = as<double>(lst["qi"]);
    double Di = as<double>(lst["Di"]);
    double b = as<double>(lst["b"]);
    double q_abnd = 0.0;
    int lc;
    if (ll == 6) {
        lc = 5;
    } else {
        q_abnd = as<double>(lst["q_abnd"]);
        lc = 7;
    }
    arma::mat results_table(lr,lc);
    if (lc == 5) {
        for (int i = 0; i < lr; i++) {
            results_table(i,0) = time(i);
            results_table(i,1) = qi / (1 + Di * time(i));
            results_table(i,2) = (qi / Di) * log(qi / results_table(i,1));
            results_table(i,3) = Di / (1 + b * Di * time(i));
            results_table(i,4) = b;
        }
    }
    if (lc == 7) {
        for (int i = 0; i < lr; i++) {
            results_table(i,0) = time(i);
            results_table(i,1) = qi / (1 + Di * time(i));
            results_table(i,2) = (qi / Di) * log(qi / results_table(i,1));
            results_table(i,3) = Di / (1 + b * Di * time(i));
            results_table(i,4) = b;
            results_table(i,5) = (qi / q_abnd - 1) / b / Di;
            results_table(i,6) = (qi / Di) * log(qi / q_abnd);
        }
    }
    return(results_table);
}


// [[Rcpp::export]]
arma::mat hyperbolic(Rcpp::List lst, arma::vec time) {

    //  Di: nominal decline rate
    int lr = time.size();
    int ll = lst.length();
    double qi = as<double>(lst["qi"]);
    double Di = as<double>(lst["Di"]);
    double b = as<double>(lst["b"]);
    double q_abnd = 0.0;
    int lc;
    if (ll == 6) {
        lc = 5;
    } else {
        q_abnd = as<double>(lst["q_abnd"]);
        lc = 7;
    }
    arma::mat results_table(lr,lc);
    if (lc == 5) {
        for (int i = 0; i < lr; i++) {
            results_table(i,0) = time(i);
            results_table(i,1) = qi / pow((1. + b * Di * time(i)),(1. / b));
            results_table(i,2) = (pow(qi,b) / Di / (1. - b)) * (pow(qi,(1. - b)) - pow(results_table(i,1), (1. - b)));
            results_table(i,3) = Di / (1. + b * Di * time(i));
            results_table(i,4) = b;
        }
    }
    if (lc == 7) {
        for (int i = 0; i < lr; i++) {
            results_table(i,0) = time(i);
            results_table(i,1) = qi / pow((1. + b * Di * time(i)),(1. / b));
            results_table(i,2) = (pow(qi,b) / Di / (1. - b)) * (pow(qi,(1. - b)) - pow(results_table(i,1), (1. - b)));
            results_table(i,3) = Di / (1. + b * Di * time(i));
            results_table(i,4) = b;
            results_table(i,5) = (pow((qi / q_abnd),b) - 1.) / Di / b;
            results_table(i,6) = (pow(qi,b) / Di / (1. - b)) * (pow(qi,(1 - b)) - pow(q_abnd,(1. - b)));
        }
    }
    return(results_table);
}


// [[Rcpp::export]]
arma::mat modified_hyperbolic(Rcpp::List lst, arma::vec time) {

    //  Di: nominal decline rate
    int lr = time.size();
    int ll = lst.length();
    double qi = as<double>(lst["qi"]);
    double Di = as<double>(lst["Di"]);
    double b = as<double>(lst["b"]);
    double Dt = as<double>(lst["Dt"]);
    double q_abnd = 0.0;
    double D_star; double t_star; double q_star;
    int lc;
    if (ll == 7) {
        lc = 5;
    } else {
        q_abnd = as<double>(lst["q_abnd"]);
        lc = 7;
    }
    arma::mat results_table(lr,lc);
    if (lc == 5) {
        if (b == 1) {
            if (ll == 7) {
                Rcpp::List lst_1 = List::create(Named("input_unit") = lst[0], Named("output_unit") = lst[1], Named("fluid") = lst[2], Named("qi") = qi, Named("Di") = Di, Named("b") = b);
                results_table = harmonic(lst_1, time);
            } else {
                Rcpp::List lst_1 = List::create(Named("input_unit") = lst[0], Named("output_unit") = lst[1], Named("fluid") = lst[2], Named("qi") = qi, Named("Di") = Di, Named("b") = b, Named("q_abnd") = q_abnd);
                results_table = harmonic(lst_1, time);
            }
        } else if (Dt == 0.) {
            if (ll == 7) {
                Rcpp::List lst_1 = List::create(Named("input_unit") = lst[0], Named("output_unit") = lst[1], Named("fluid") = lst[2], Named("qi") = qi, Named("Di") = Di, Named("b") = b);
                results_table = hyperbolic(lst_1, time);
            } else {
                Rcpp::List lst_1 = List::create(Named("input_unit") = lst[0], Named("output_unit") = lst[1], Named("fluid") = lst[2], Named("qi") = qi, Named("Di") = Di, Named("b") = b, Named("q_abnd") = q_abnd);
                results_table = hyperbolic(lst_1, time);
            }
            // results_table = hyperbolic(lst(idx), time);
        } else {
            D_star = Dt;
            t_star = (Di / D_star - 1.) / (b * Di);
            q_star = (qi / pow((1. + b * Di * t_star),(1. / b))) * (1. / exp(-D_star * t_star));
            for (int i = 0; i < lr; i++) {
                results_table(i,0) = time(i);
                if (time(i) <= t_star) {
                    results_table(i,1) = qi / pow((1. + b * Di * time(i)),(1. / b));
                    results_table(i,2) = (pow(qi,b) / Di / (1. - b)) * (pow(qi,(1. - b)) - pow(results_table(i,1), (1. - b)));
                    results_table(i,3) = Di / (1. + b * Di * time(i));
                    results_table(i,4) = b;
                } else {
                    results_table(i,1) = q_star * exp(-D_star * time(i));
                    results_table(i,2) = (qi / (1. - b) / Di) * (1. - pow((1. + b * Di * t_star),(1. - (1. / b)))) + (q_star * (-exp(-D_star * time(i)) + exp(-D_star * t_star)) / D_star);
                    results_table(i,3) = D_star;
                    results_table(i,4) = 0.;
                }
            }
        }
    }
    if (lc == 7) {
        if (b == 1) {
            results_table = harmonic(lst, time);
        } else if (Dt == 0.) {
            results_table = hyperbolic(lst, time);
        } else {
            D_star = Dt;
            t_star = (Di / D_star - 1.) / (b * Di);
            q_star = (qi / pow((1. + b * Di * t_star),(1. / b))) * (1. / exp(-D_star * t_star));
            double time_abnd;
            if (q_abnd < q_star * exp(-D_star * t_star)) {
                time_abnd = -log(q_abnd / q_star) / D_star;
            } else {
                time_abnd = (pow((qi / q_abnd),b) - 1.) / Di / b;
            }
            for (int i = 0; i < lr; i++) {
                results_table(i,0) = time(i);
                if (time(i) <= t_star) {
                    results_table(i,1) = qi / pow((1. + b * Di * time(i)),(1. / b));
                    results_table(i,2) = (pow(qi,b) / Di / (1. - b)) * (pow(qi,(1. - b)) - pow(results_table(i,1), (1. - b)));
                    results_table(i,3) = Di / (1. + b * Di * time(i));
                    results_table(i,4) = b;
                } else {
                    results_table(i,1) = q_star * exp(-D_star * time(i));
                    results_table(i,2) = (qi / (1. - b) / Di) * (1. - pow((1. + b * Di * t_star),(1. - (1. / b)))) + (q_star * (-exp(-D_star * time(i)) + exp(-D_star * t_star)) / D_star);
                    results_table(i,3) = D_star;
                    results_table(i,4) = 0.;
                }
                results_table(i,5) = -log(q_abnd / q_star) / D_star;
                if (time_abnd <= t_star) {
                    results_table(i,6) = (pow(qi,b) / Di / (1. - b)) * (pow(qi,(1. - b)) - pow(q_abnd,(1. - b)));
                } else {
                    results_table(i,6) = (qi / (1. - b) / Di) * (1. - pow((1. + b * Di * t_star),(1. - (1. / b)))) + (q_star * (-exp(-D_star * time_abnd) + exp(-D_star * t_star)) / D_star);
                }
            }
        }
    }
    return(results_table);
}

