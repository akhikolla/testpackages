//
//  MLModelSelection.hpp
//  
//
//  Created by kuojung on 2019/12/11.
//

#ifndef MLModelSelection_hpp
#define MLModelSelection_hpp

//#include <math.h>

#include <stdio.h>
//#include <algorithm>
//#include <assert.h>
#include <cmath>
#include <ctime>    // For time()
//#include <cstdlib>  // For srand() and rand()
//#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iomanip>
//#include <list>
//#include <limits>
#include <vector>
#include <string>
//#include <sstream>
#include <algorithm>


#include<iostream>
//#include<chrono> //for sleeping
#include<thread> // --do--
#include<cstdlib>//for random increments


//#include <RcppArmadillo.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace std;

using namespace Rcpp;
using namespace arma;


class MLModelSelection{
private:
    int Num_of_iterations, AR_Order, ISD_Model;
    int Num_of_obs, Num_of_covariates, Num_of_attributes;
    int Num_of_alphas;
    cube Y, X;
    field <mat>X_all;
    vec TimePointsAvailable;//group_indices,
    List Data, InitialValues, HyperPara, UpdatePara, TuningPara;
    cube beta_samples, lambda_samples;
    mat nu_samples;
    field<cube> alpha_samples, delta_samples;
    double sigma2_beta, sigma2_alpha, sigma2_lambda, sigma2_nu;
    double tuning_alpha, tuning_lambda, tuning_nu;
    
    bool updatebeta, updatealpha, updatedelta, updatelambda, updatenu;
    
    mat K_diag;
    int lambda_size;
    vec nu_mean;
    mat beta_mean, lambda_mean;
    cube alpha_mean, delta_mean;
    double AIC, BIC, DIC, MPL, MSPE;
    cube pred_y;
    

public:
    MLModelSelection(int iNum_of_iterations, List list_Data, List list_InitialValues, List list_HyperPara, List list_UpdatePara, List list_TuningPara);
    
    //mat Dit(int i, int t, vec nu, mat lambda);
    mat Dit(int i, int t, vec nu, mat lambda);
    mat Di(int i, vec nu, mat lambda, int tp);
    
    //mat Ti_Mat(cube alpha);
    mat Ti_Mat(cube alpha, int Subject);
    
    void Update_beta(int iter);
    void Update_alpha(int iter);
    void Update_delta(int iter);
    void Update_delta_alpha(int iter);
    void Update_lambda(int iter);
    void Update_nu(int iter);
    
    void ParameterEstimation(); 
    SEXP MCMC_Procedure();
    //SEXP MCMC_Procedure_Robust();
};


#endif /* MLModelSelection_hpp */
