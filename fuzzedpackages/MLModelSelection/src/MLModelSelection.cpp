//
//  MLModelSelection.cpp
//  
//
//  Created by kuojung on 2019/12/11.
//

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "MLModelSelection.h"
RNGScope scope;
//#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

MLModelSelection::MLModelSelection(int iNum_of_iterations, List list_Data, List list_InitialValues, List list_HyperPara, List list_UpdatePara, List list_TuningPara)
{
    Num_of_iterations = iNum_of_iterations;
    Data = list_Data;
    InitialValues = list_InitialValues;
    HyperPara = list_HyperPara;
    UpdatePara = list_UpdatePara;
    TuningPara = list_TuningPara;
    
    Rcout<< "Read Data" << endl;
    Y = as<cube>(Data["Y"]);
    X = as<cube>(Data["X"]);
    TimePointsAvailable = as<vec>(Data["TimePointsAvailable"]);
    //X_all = as<cube>(Data["X.all"]);
    //group_indices = as<vec>(Data["group.indices"]);
    AR_Order = as<int>(Data["AR.Order"]);
    //CovariateForSigma = as<int>(Data["CovariateForSigma"])-1;
    ISD_Model = as<int>(Data["ISD.Model"]);
    
    Num_of_obs = Y.n_rows;
    Num_of_covariates = X.n_slices;
    Num_of_attributes = Y.n_slices;
    //TimePoints = Y.n_cols;
    pred_y.set_size(size(Y));
    
    updatebeta = as<bool>(UpdatePara["UpdateBeta"]);
    updatedelta = as<bool>(UpdatePara["UpdateDelta"]);
    updatealpha = as<bool>(UpdatePara["UpdateAlpha"]);
    updatelambda = as<bool>(UpdatePara["UpdateLambda"]);
    updatenu = as<bool>(UpdatePara["UpdateNu"]);


    Rcout<< "Read Hyperparameters." << endl;
    // Hyperparameters
    sigma2_beta = as<double>(HyperPara["sigma2.beta"]);
    sigma2_alpha = as<double>(HyperPara["sigma2.alpha"]);
    sigma2_lambda = as<double>(HyperPara["sigma2.lambda"]);
    sigma2_nu = as<double>(HyperPara["sigma2.nu"]);

    Rcout<< "Read Tuning parameters." << endl;
    // Hyperparameters
    tuning_alpha = as<double>(TuningPara["TuningAlpha"]);
    tuning_lambda = as<double>(TuningPara["TuningLambda"]);
    tuning_nu = as<double>(TuningPara["TuningNu"]);

    
    beta_samples.set_size(Num_of_covariates, Num_of_attributes, Num_of_iterations);
    beta_samples.zeros();
    delta_samples.set_size(Num_of_iterations);
    alpha_samples.set_size(Num_of_iterations);
    
    
    // 2 and 3 would be modified later with input
    mat lambda_ini_mat = as<mat>(InitialValues["lambda"]);
    lambda_size = lambda_ini_mat.n_cols;
    
    
    lambda_samples.set_size(Num_of_attributes, lambda_size, Num_of_iterations);
    lambda_samples.zeros();
    nu_samples.set_size(3, Num_of_iterations);
    nu_samples.zeros();
    
    Rcout<< "Initial Values" << endl;
    beta_samples.slice(0) = as<mat>(InitialValues["beta"]);
    delta_samples(0) = as<cube>(InitialValues["delta"]);
    alpha_samples(0) = as<cube>(InitialValues["alpha"]);
    
    Num_of_alphas = alpha_samples(0).n_slices;
    
    lambda_samples.slice(0) = lambda_ini_mat; //as<mat>(InitialValues["lambda"]);
    nu_samples.col(0) = as<vec>(InitialValues["nu"]);

    
    
    //T_diag.eye(TimePoints, TimePoints);
    K_diag.eye(Num_of_attributes, Num_of_attributes);
    //T_diag.set_size(Num_of_obs);
    X_all.set_size(Num_of_obs);
    
    
    //X_all.set_size(TimePoints*Num_of_attributes, Num_of_attributes*Num_of_covariates, Num_of_obs);
    
    rowvec X_row;
    mat Xi_tmp;
    //Rcout << "X(0, 0, span::all)  = " << X(span(0), span(0), span::all) << endl << X_row << endl;
    Rcout << "Create X.all" << endl;
    for(int i=0; i<Num_of_obs; i++){
        for(int t=0; t<TimePointsAvailable(i); t++){
            //T_diag.set_size(TimePointsAvailable(i), TimePointsAvailable(i));
            //T_diag.eye();
            rowvec X_row(X(span(i), span(t), span::all));
            if(t==0){
                Xi_tmp = kron(K_diag, X_row);
                //Rcout << "t=0 " << Xi_tmp << endl;
            }
            else{
                Xi_tmp = join_cols(Xi_tmp, kron(K_diag, X_row));
                //Rcout << "t=1 " << Xi_tmp << endl;
            }
        }
        X_all(i) = Xi_tmp;
    }
    Rcout << "Create X.all Done" << endl;
    //Rcout << X_all.slice(0) << endl;
    //mat Ti_tmp = kron(T_diag,K_diag);
    //mat tmp_mat(TimePoints, TimePoints);
    //uvec aa = regspace<uvec>(1,  9), bb=regspace<uvec>(2,  9);
    //tmp_mat(1, 9, 2, 9) = 1;
    //Rcout << "tmp_mat =" << endl << tmp_mat << endl;
    
    // For parameter estimates
    beta_mean.set_size(Num_of_covariates, Num_of_attributes);
    beta_mean.zeros();
    delta_mean.set_size(Num_of_attributes, Num_of_attributes, Num_of_alphas);
    alpha_mean.set_size(Num_of_attributes, Num_of_attributes, Num_of_alphas);
    delta_mean.zeros();
    alpha_mean.zeros();
    
    lambda_mean.set_size(Num_of_attributes, 2);
    lambda_mean.zeros();
    nu_mean.set_size(3);
    nu_mean.zeros();

    AIC = 0.;
    BIC = 0.;
    DIC = 0.;
    MPL = 0.;
    MSPE = 0.;
    
    Rcout<< "All set" << endl;
}

/*
mat MLModelSelection::Dit(int i, int t, vec nu, mat lambda)
{
    mat F_tmp(Num_of_attributes, Num_of_attributes), F(Num_of_attributes, Num_of_attributes);
    F.zeros();
    
    //Rcout << "Dit 1" << endl;
    for(int j=0; j<Num_of_attributes; j++)
        for(int i=0; i<Num_of_attributes; i++)
            F_tmp(i, j) = nu(i) + nu(j);
    F_tmp = datum::pi*exp(F_tmp)/(1.+exp(F_tmp));
    F(0, 0) = 1;
    //Rcout << "Dit 2" << endl;
    for(int l=1; l<Num_of_attributes; l++)
        F(l, 0) = cos(F_tmp(l, 0));
    for(int m = 1; m<Num_of_attributes-1; m++)
        for(int l = m+1; l<Num_of_attributes; l++)
            F(l, m) = cos(F_tmp(l, m))*prod(sin(F_tmp(l, span(0, m-1) )));
    //Rcout << "Dit 3" << endl;
    for(int m=1; m<Num_of_attributes; m++)
        F(m, m) = prod(sin(F_tmp(m, span(0, m-1) )));
    //Rcout << "Dit 4" << endl;
    mat Ri = F * F.t();
    //Rcout << "Dit 5" << endl;
    vec Si_vec(lambda.n_rows);
    
    for(int ii=0; ii<lambda.n_rows; ii++)
        Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, 1));
    //Rcout << "Dit 6" << endl;
    mat Si = diagmat(Si_vec);
    return (Si*Ri*Si);
}
*/

mat MLModelSelection::Dit(int i, int t, vec nu, mat lambda)
{
    mat F_tmp(Num_of_attributes, Num_of_attributes), F(Num_of_attributes, Num_of_attributes);
    F.zeros();
    
    //Rcout << "Dit 1" << endl;
    for(int j=0; j<Num_of_attributes; j++)
        for(int i=0; i<Num_of_attributes; i++)
            F_tmp(i, j) = nu(i) + nu(j);
    F_tmp = datum::pi*exp(F_tmp)/(1.+exp(F_tmp));
    F(0, 0) = 1;
    //Rcout << "Dit 2" << endl;
    for(int l=1; l<Num_of_attributes; l++)
        F(l, 0) = cos(F_tmp(l, 0));
    for(int m = 1; m<Num_of_attributes-1; m++)
        for(int l = m+1; l<Num_of_attributes; l++)
            F(l, m) = cos(F_tmp(l, m))*prod(sin(F_tmp(l, span(0, m-1) )));
    //Rcout << "Dit 3" << endl;
    for(int m=1; m<Num_of_attributes; m++)
        F(m, m) = prod(sin(F_tmp(m, span(0, m-1) )));
    //Rcout << "Dit 4" << endl;
    mat Ri = F * F.t();
    //Rcout << "Dit 5" << endl;
    vec Si_vec(lambda.n_rows);
    for(int ii=0; ii<lambda.n_rows; ii++){
        //if(AR_Order == 1) //For simulation only
            //Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, p));
            //Rcout << "X(" << i << ", " <<  t << "," << p <<") = " << X(i, t, p) << endl;
        if(ISD_Model == 1)
            Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, 1)); //Arm
        if(ISD_Model == 2)
            Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, 4)); //Sex
        if(ISD_Model == 3)
            Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, 5)); //Meta
        if(ISD_Model == 4)
            Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, 1)+lambda(ii, 2)*X(i, t, 4)); // Arm, Sex
        if(ISD_Model == 5)
            Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, 1)+lambda(ii, 2)*X(i, t, 5)); // Arm, Meta
        if(ISD_Model == 6)
            Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, 4)+lambda(ii, 2)*X(i, t, 5)); // Arm, Meta
                //Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, p) + lambda(ii, 2)*pow(X(i, t, p), 2.));
        if(ISD_Model == 7)
            Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, 1)+lambda(ii, 2)*X(i, t, 4) + lambda(ii, 3)*X(i, t, 5));
        if(ISD_Model == 8)
            Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, 3));
        if(ISD_Model == 9)
            Si_vec(ii) = exp(lambda(ii, 0) + lambda(ii, 1)*X(i, t, 3) + lambda(ii, 2)*pow(X(i, t, 3), 2.));
        if(ISD_Model == 10)
            Si_vec(ii) = exp(lambda(ii, 0));
    }
    //Rcout << "Dit 6" << endl;
    mat Si = diagmat(Si_vec);
    return (Si*Ri*Si);
}

mat MLModelSelection::Di(int i, vec nu, mat lambda, int tp)
{
    //Rcout << "Di" << endl;
    mat Di_mat(Num_of_attributes*tp, Num_of_attributes*tp, fill::zeros);
    //cube Di_tmp(Num_of_attributes)
    for(int t=0; t<tp; t++)
        Di_mat.submat(Num_of_attributes*t, Num_of_attributes*t, Num_of_attributes*(t+1)-1, Num_of_attributes*(t+1)-1) = Dit(i, t, nu, lambda);
    return (Di_mat);
}

/*
mat MLModelSelection::Ti_Mat(cube alpha)
{
    mat Ti = kron(T_diag,K_diag);
    mat phi_mat(Num_of_attributes, Num_of_attributes, fill::zeros);
    int K = Num_of_attributes;
    //Rcout << "Ti_Mat" << endl;
    //Rcout << size(alpha) << endl;
    for(int t=1; t<TimePoints; t++)
        for(int j=0; j<t; j++){
            for(int k=0; k<Num_of_attributes; k++)
                for(int g=0; g<Num_of_attributes; g++){
                    //Rcout << "(k, g) = " << "(" << k << "," << g << ")="<<  alpha(k, g, 0) + alpha(k, g, 1)*abs(t-j) + alpha(k, g, 2)*pow(t-j, 2) << endl;
                    phi_mat(k, g) = alpha(k, g, 0) + alpha(k, g, 1)*abs(t-j) + alpha(k, g, 2)*pow(t-j, 2);
                }
            //Rcout << "phi_mat = " << endl << phi_mat << endl;
            
            Ti.submat(t*K,  j*K, (t+1)*K-1, (j+1)*K-1) = -phi_mat;
            
            //Rcout << "Ti.submat(t*K,  j*K, (t+1)*K-1, (j+1)*K-1) =" << endl << Ti.submat(t*K,  j*K, (t+1)*K-1, (j+1)*K-1) << endl;
        }
    //Rcout << "Ti_Mat end" << endl;
    return(Ti);
}
*/



mat MLModelSelection::Ti_Mat(cube alpha, int Subject)
{
    int tp = TimePointsAvailable(Subject);
    mat Ti = kron(eye(tp,tp),K_diag);
    mat phi_mat(Num_of_attributes, Num_of_attributes, fill::zeros);
    int K = Num_of_attributes;
    //Rcout << "Ti_Mat" << endl;
    //Rcout << size(alpha) << endl;
    for(int t=1; t<tp; t++)
        for(int j=0; j<t; j++){
            for(int k=0; k<Num_of_attributes; k++)
                for(int g=0; g<Num_of_attributes; g++){
                    //Rcout << "(k, g) = " << "(" << k << "," << g << ")="<<  alpha(k, g, 0) + alpha(k, g, 1)*abs(t-j) + alpha(k, g, 2)*pow(t-j, 2) << endl;
                    if(AR_Order == 1)
                        phi_mat(k, g) = alpha(k, g, 0)*(abs(t-j)==1);
                        //phi_mat(k, g) = alpha(k, g, 0) + alpha(k, g, 1)*abs(t-j) + alpha(k, g, 2)*pow(t-j, 2); //alpha(k, g, 0)*(abs(t-j)==1) + alpha(k, g, 1)*(abs(t-j)==2)+ alpha(k, g, 2)*(abs(t-j)==3);//
                    if(AR_Order == 2)
                        phi_mat(k, g) = alpha(k, g, 0)*(abs(t-j)==1) + alpha(k, g, 1)*(abs(t-j)==2);
                    if(AR_Order == 3)
                        phi_mat(k, g) = alpha(k, g, 0)*(abs(t-j)==1) + alpha(k, g, 1)*(abs(t-j)==2)+ alpha(k, g, 2)*(abs(t-j)==3);
                    if(AR_Order == 4)
                        phi_mat(k, g) = alpha(k, g, 0)*(abs(t-j)==1) + alpha(k, g, 1)*(abs(t-j)==2)+ alpha(k, g, 2)*(abs(t-j)==3)+ alpha(k, g, 3)*(abs(t-j)==4);
                    if(AR_Order == 5)
                        phi_mat(k, g) = alpha(k, g, 0)*(abs(t-j)==1) + alpha(k, g, 1)*(abs(t-j)==2)+ alpha(k, g, 2)*(abs(t-j)==3)+ alpha(k, g, 3)*(abs(t-j)==4)+ alpha(k, g, 4)*(abs(t-j)==5);
                    if(AR_Order == 6)
                        phi_mat(k, g) = alpha(k, g, 0) + alpha(k, g, 1)*t/10. + alpha(k, g, 2)*pow(t/10., 2.);
                }
            //Rcout << "phi_mat = " << endl << phi_mat << endl;
            
            Ti.submat(t*K,  j*K, (t+1)*K-1, (j+1)*K-1) = -phi_mat;
            //Rcout << "(t, j, K) = " << "(" << t <<"," << j << "," << K <<")"<<endl;
            //Rcout << "Ti.submat(t*K,  j*K, (t+1)*K-1, (j+1)*K-1) =" << endl << Ti.submat(t*K,  j*K, (t+1)*K-1, (j+1)*K-1) << endl;
        }
    //Rcout << "Ti_Mat end" << endl;
    
    return(Ti);
}




void MLModelSelection::Update_beta(int iter)
{
    //Rcout << "Generate beta" << endl;
    mat Sigma_tmp(Num_of_covariates*Num_of_attributes, Num_of_covariates*Num_of_attributes, fill::zeros), Sigma_tmp_i;
    Sigma_tmp.zeros();
    vec mu_tmp(Num_of_covariates*Num_of_attributes);
    mu_tmp.zeros();
    mat Ti, Di_mat_inv;
    vec yi;
    //cube alpha_tmp_cube = alpha_samples(iter);
    //Ti = Ti_Mat(alpha_samples(iter));
    //Rcout << " Ti " << endl << Ti.submat(0, 0, 8, 8) << endl;
    //Rcout << X_all.subcube(1, span::all, span::all) << endl;
    
    //Rcout <<" X.all.row(1) = " << endl << X_all.row(1) << endl;
    
    //Rcout << "Dit = "<< Dit(1, 1, nu_samples.col(iter), lambda_samples.slice(iter)) << endl;
   
    //Rcout << vectorise((mat)Y.row(0), 1) << endl;
    int tp;
    for(int i=0; i<Num_of_obs; i++){
        tp = TimePointsAvailable(i);
        Ti = Ti_Mat(alpha_samples(iter), i);
        //Rcout << "Sigma_tmp_i" << endl;
        //Rcout << "X.all = " << size(X_all(i))  << endl;
        //Rcout << "i = " << i <<"\ttp = "<< tp << endl;
        //Rcout << "Ti = " << size(Ti) << endl;
        
        Di_mat_inv = Di(i, nu_samples.col(iter), lambda_samples.slice(iter), tp).i();
        
        //Rcout << "Di_mat_inv = " << size(Di_mat_inv) << endl;
        //A.raw_print(cout, "A:");
        
        //Ti.submat(0, 0, 14, 14).raw_print(std::cout);
        //Rcout.precision(2) << Ti.submat(0, 0, 14, 14) << endl;
        
        Sigma_tmp_i = X_all(i).t() * (Ti.t() * Di_mat_inv * Ti);
        
        //Rcout << "Sigma_tmp_i = " << size(Sigma_tmp_i) << endl;
        
        yi =  vectorise((mat)Y(span(i), span(0, tp-1), span::all), 1).t();
        
        //Rcout << "i = " << i << "\tYi = " << endl << yi << endl;
        
        //yi_tmp = vectorise((mat)Y.row(i), 1).t();
        
        
        
        //yi_tmp.head(TimePointsAvailable(i));
        //Rcout << "yi = " << yi << endl;
        mu_tmp += Sigma_tmp_i *yi;
        //Rcout << "check 1" << endl;
        Sigma_tmp += Sigma_tmp_i * X_all(i);
    }
    //Rcout << "Sigma_tmp_i done" << endl;
    //Rcout << "Check 2" << endl;
    mat diag_sigma2_beta(Num_of_covariates*Num_of_attributes, Num_of_covariates*Num_of_attributes, fill::zeros);
    diag_sigma2_beta.diag() += 1./sigma2_beta;
    //Rcout << "Check 2" << endl;
    mat Sigma = (Sigma_tmp + diag_sigma2_beta).i();
    vec mu = Sigma * mu_tmp;
    //Rcout << "Check 3" << endl;
    mat beta_sim = rmvnorm(1, mu, Sigma);
    //Rcout << "Check 4" << endl;
    for(int k=0; k<Num_of_attributes; k++)
        for(int p = 0; p<Num_of_covariates; p++)
            beta_samples(p, k, iter+1) = beta_sim(Num_of_covariates*k + p);
     
    //Rcout << "done beta" << endl;
}

void MLModelSelection::Update_delta(int iter)
{
    //Rcout << "Update delta" << endl;
    
    //Rcout << vectorise(beta_samples.slice(iter+1)) << endl;
    
    //int M = nu_samples.n_rows;
    cube alpha_cube;
    cube delta_cube = delta_samples(iter);
    mat Ti_0, Ti_1;
    double delta_1 = 0., delta_0 = 0.;
    mat Di_mat_inv;
    vec yi, yi_tmp;
    double prob_delta = 0.;
    double double_alpha_tmp_old;
    int tp;
    for(int m=0; m<Num_of_alphas; m++)
        for(int k=0; k<Num_of_attributes; k++)
            for(int g=0; g<Num_of_attributes; g++){
                alpha_cube = alpha_samples(iter);
                if(alpha_cube(k, g, m)==0){
                    delta_1 = -0.5*log(sigma2_alpha*2*datum::pi); //Rf_dnorm4(0, 0., sqrt(sigma2_alpha), 1);
                    prob_delta = 1./(1.+exp(-delta_1));
                    delta_cube(k, g, m) = Rf_rbinom(1.,  prob_delta);
                }
                else{
                    double_alpha_tmp_old = alpha_cube(k, g, m);
                    alpha_cube(k, g, m) = 0;
                    delta_1 = delta_0 = 0.;
                    for(int i=0; i<Num_of_obs; i++){
                        tp = TimePointsAvailable(i);
                        Ti_0 = Ti_Mat(alpha_cube, i);
                        Ti_1 = Ti_Mat(alpha_samples(iter), i);
                        
                        yi =  vectorise((mat)Y(span(i), span(0, tp-1), span::all), 1).t();
                        //yi = vectorise((mat)Y.row(i), 1).t();
                        //Dit_tmp = Dit(i, 1, nu_samples.col(iter), lambda_samples.slice(iter), 3);
                        
                        Di_mat_inv = Di(i, nu_samples.col(iter), lambda_samples.slice(iter), tp).i();
                        
                        yi_tmp = yi-X_all(i)*vectorise(beta_samples.slice(iter+1));
                        
                        //Rcout << "check 2" << endl;
                        
                        delta_1 += -0.5*as_scalar( yi_tmp.t()* Ti_1.t()*Di_mat_inv*Ti_1*yi_tmp );
                        
                        delta_0 += -0.5*as_scalar( yi_tmp.t()* Ti_0.t()*Di_mat_inv*Ti_0*yi_tmp );
                    }
                    
                    
                    delta_1 += Rf_dnorm4(double_alpha_tmp_old, 0., sqrt(sigma2_alpha), 1);
                    prob_delta = 1./(1.+exp(delta_0-delta_1));
                    delta_cube(k, g, m) = Rf_rbinom(1.,  prob_delta);
                }
            }
    delta_samples(iter+1) = delta_cube;
}

void MLModelSelection::Update_alpha(int iter)
{
    //Rcout << "Update alpha" << endl;
    cube alpha_cube_old = alpha_samples(iter);
    cube alpha_cube_new =  alpha_samples(iter);
    cube delta_cube = delta_samples(iter+1);
    double alpha_den = 0., alpha_num = 0.;
    mat Ti_den, Ti_num;
    vec yi, yi_tmp;
    int tp;
    mat Di_mat_inv;
    
    for(int m=0; m<Num_of_alphas; m++)
        for(int k=0; k<Num_of_attributes; k++)
            for(int g=0; g<Num_of_attributes; g++){
                if(delta_cube(k, g, m)==0)
                    alpha_cube_new(k, g, m) = 0;
                else{
                    alpha_den = alpha_num =0.;
                    alpha_cube_old = alpha_cube_new;
                    alpha_cube_new(k, g, m) =  Rf_rnorm(alpha_cube_old(k, g, m), tuning_alpha);
                    
                    for(int i=0; i<Num_of_obs; i++){
                        tp = TimePointsAvailable(i);
                        yi = vectorise((mat)Y(span(i), span(0, tp-1), span::all), 1).t(); //vectorise((mat)Y.row(i), 1).t();
                        //Dit_mat = Di(i, nu_samples.col(iter), lambda_samples.slice(iter), 3, tp);
                        Di_mat_inv = Di(i, nu_samples.col(iter), lambda_samples.slice(iter), tp).i();
                        //Dit_tmp = Dit(i, 1, nu_samples.col(iter), lambda_samples.slice(iter), 3);
                        yi_tmp = yi-X_all(i)*vectorise(beta_samples.slice(iter+1));
                        
                        Ti_den = Ti_Mat(alpha_cube_old, i);
                        alpha_den += -0.5*as_scalar( yi_tmp.t()* Ti_den.t()*Di_mat_inv*Ti_den*yi_tmp );
                        
                        Ti_num = Ti_Mat(alpha_cube_new, i);
                        alpha_num += -0.5*as_scalar( yi_tmp.t()* Ti_num.t()*Di_mat_inv*Ti_num*yi_tmp );
                    }
                    
                    //Ti_den.submat(0, 0, 11, 11).raw_print(cout, "Ti_den = ");
                    
                    alpha_den += -0.5*pow(alpha_cube_old(k, g, m), 2.)/sigma2_alpha;
                    alpha_num += -0.5*pow(alpha_cube_new(k, g, m), 2.)/sigma2_alpha;
                    //Rcout << " alpha_den = " << alpha_den << "\t alpha_num = " << alpha_num << "\talpha_num-alpha_den =" << alpha_num-alpha_den << endl;
                    if(alpha_num-alpha_den < log(Rf_runif(0., 1.)))
                        alpha_cube_new(k, g, m) = alpha_cube_old(k, g, m);
                }
            }
    alpha_samples(iter+1) = alpha_cube_new;
}


void MLModelSelection::Update_delta_alpha(int iter)
{
    Rcout << "Generate delta and alpha" << endl;
    cube alpha_cube_old = alpha_samples(iter);
    cube alpha_cube_new =  alpha_samples(iter);
    cube delta_cube_new = delta_samples(iter);
    cube delta_cube_old = delta_samples(iter);
    double alpha_den = 0., alpha_num = 0.;
    mat Ti_den, Ti_num;
    vec yi, yi_tmp;
    int tp;
    mat Di_mat_inv;
    
    for(int m=0; m<Num_of_alphas; m++)
        for(int k=0; k<Num_of_attributes; k++)
            for(int g=0; g<Num_of_attributes; g++){
                alpha_cube_old = alpha_cube_new;
                delta_cube_old = delta_cube_new;
                
                delta_cube_new(k, g, m) = Rf_rbinom(1.,  0.5);
                
                if(delta_cube_new(k, g, m)==0)
                    alpha_cube_new(k, g, m) = 0;
                else
                    alpha_cube_new(k, g, m) =  Rf_rnorm(alpha_cube_old(k, g, m), tuning_alpha);
                
                alpha_den = alpha_num =0.;
                for(int i=0; i<Num_of_obs; i++){
                    //Rcout << "(k, g, m) = (" << k << ", "<< g << ", " << m << ")" << endl;
                    tp = TimePointsAvailable(i);
                    yi = vectorise((mat)Y(span(i), span(0, tp-1), span::all), 1).t(); //vectorise((mat)Y.row(i), 1).t();
                    //Dit_mat = Di(i, nu_samples.col(iter), lambda_samples.slice(iter), 3, tp);
                    Di_mat_inv = Di(i, nu_samples.col(iter), lambda_samples.slice(iter), tp).i();
                    //Dit_tmp = Dit(i, 1, nu_samples.col(iter), lambda_samples.slice(iter), 3);
                    yi_tmp = yi-X_all(i)*vectorise(beta_samples.slice(iter+1));
                    
                    Ti_den = Ti_Mat(alpha_cube_old, i);
                    //Ti_den.submat(0, 0, 11, 11).raw_print(cout, "Ti_den = ");
                    //Rcout << "Ti_den = " << endl << Ti_den.submat(0, 0, 8, 8).raw_print() << endl;
                    alpha_den += -0.5*as_scalar( yi_tmp.t()* Ti_den.t()*Di_mat_inv*Ti_den*yi_tmp );
                    
                    Ti_num = Ti_Mat(alpha_cube_new, i);
                    
                    //Rcout << "Ti_num = " << endl << Ti_num.submat(0, 0, 8, 8) << endl;
                    //Ti_num.submat(0, 0, 11, 11).raw_print(cout, "Ti_num = ");
                    alpha_num += -0.5*as_scalar( yi_tmp.t()* Ti_num.t()*Di_mat_inv*Ti_num*yi_tmp );
                }
                
                //alpha_den += -0.5*pow(alpha_cube_old(k, g, m), 2.)/sigma2_alpha;
                //alpha_num += -0.5*pow(alpha_cube_new(k, g, m), 2.)/sigma2_alpha;
                //Rcout << " alpha_den = " << alpha_den << "\t alpha_num = " << alpha_num << "\talpha_num-alpha_den =" << alpha_num-alpha_den << endl;
                if(alpha_num-alpha_den < log(Rf_runif(0., 1.))){
                    alpha_cube_new(k, g, m) = alpha_cube_old(k, g, m);
                    delta_cube_new(k, g, m) = delta_cube_old(k, g, m);
                }
            }
    alpha_samples(iter+1) = alpha_cube_new;
    delta_samples(iter+1) = delta_cube_new;
}


void MLModelSelection::Update_lambda(int iter)
{
    //Rcout << "Generate lambda" << endl;
    double lambda_den = 0, lambda_num = 0.;
    mat Ti, Di_mat_inv;
    vec lambda_vec = vectorise(lambda_samples.slice(iter));
    mat lambda_Sigma(lambda_vec.n_elem, lambda_vec.n_elem, fill::eye);
    lambda_samples.slice(iter+1) = reshape(rmvnorm(1, lambda_vec, tuning_lambda*lambda_Sigma), Num_of_attributes, lambda_size);
    vec yi, yi_tmp;
    int tp;
    for(int i=0; i<Num_of_obs; i++){
        tp = TimePointsAvailable(i);
        yi = vectorise((mat)Y(span(i), span(0, tp-1), span::all), 1).t();
        Ti = Ti_Mat(alpha_samples(iter+1), i);
        yi_tmp = yi-X_all(i)*vectorise(beta_samples.slice(iter+1));
        Di_mat_inv = Di(i, nu_samples.col(iter), lambda_samples.slice(iter), tp).i();
        lambda_den += 0.5*log(det(Di_mat_inv)) - 0.5*as_scalar(yi_tmp.t()* Ti.t()*Di_mat_inv*Ti*yi_tmp);
        Di_mat_inv = Di(i, nu_samples.col(iter), lambda_samples.slice(iter+1), tp).i();
        lambda_num += 0.5*log(det(Di_mat_inv)) - 0.5*as_scalar(yi_tmp.t()* Ti.t()*Di_mat_inv*Ti*yi_tmp);
    }
    lambda_den = lambda_den-0.5*accu(square(lambda_samples.slice(iter)))/sigma2_lambda;
    lambda_num = lambda_num-0.5*accu(square(lambda_samples.slice(iter+1)))/sigma2_lambda;
    
    if(log(Rf_runif(0., 1.)) > lambda_num - lambda_den )
        lambda_samples.slice(iter+1) = lambda_samples.slice(iter);
    
}

void MLModelSelection::Update_nu(int iter)
{
    //Rcout << "Generate nu" << endl;
    double nu_den = 0., nu_num = 0.;
    vec nu_cand = vectorise(rmvnorm(1, nu_samples.col(iter), tuning_nu*K_diag));
    vec yi_tmp, yi;
    mat Ti, Di_mat_inv;
    int tp;
    for(int i=0; i<Num_of_obs; i++){
        tp = TimePointsAvailable(i);
        yi = vectorise((mat)Y(span(i), span(0, tp-1), span::all), 1).t();
        Ti = Ti_Mat(alpha_samples(iter+1), i);
        yi_tmp = yi-X_all(i)*vectorise(beta_samples.slice(iter+1));
        Di_mat_inv = Di(i, nu_samples.col(iter), lambda_samples.slice(iter+1), tp).i();
        nu_den += 0.5*log(det(Di_mat_inv)) - 0.5*as_scalar(yi_tmp.t()* Ti.t()*Di_mat_inv*Ti*yi_tmp);
        Di_mat_inv = Di(i, nu_cand, lambda_samples.slice(iter+1), tp).i(); //Dit(i, 1, (nu_cand), lambda_samples.slice(iter+1), 3);
        nu_num += 0.5*log(det(Di_mat_inv)) - 0.5*as_scalar(yi_tmp.t()* Ti.t()*Di_mat_inv*Ti*yi_tmp);
     }
    
    nu_den = nu_den - 0.5*accu(nu_samples.col(iter))/sigma2_nu;
    nu_num = nu_num - 0.5*accu(nu_samples.col(iter+1))/sigma2_nu;
    
    if(log(Rf_runif(0., 1.)) < nu_num - nu_den )
        nu_samples.col(iter+1) = nu_cand;
    else
        nu_samples.col(iter+1) = nu_samples.col(iter);
}


void MLModelSelection::ParameterEstimation()
{
    for(int iter=0; iter<Num_of_iterations; iter++){
        alpha_mean += alpha_samples(iter)/Num_of_iterations;
        delta_mean += delta_samples(iter)/Num_of_iterations;
    }
    lambda_mean = mean(lambda_samples, 2);
    beta_mean = mean(beta_samples, 2);
    nu_mean = mean(nu_samples, 1);
    double res, phi_est;
    int tp;
    rowvec x;

    for(int i=0; i<Num_of_obs; i++){
        tp = TimePointsAvailable(i);
        for(int t=0; t<tp; t++)
            for(int k=0; k<Num_of_attributes; k++)
                if(t==0){
                    //conv_to< rowvec >::from(X.tube(i, t))
                    //reshape(X.tube(i, t), Num_of_covariates, 1)
                    x = X.tube(i, t);
                    pred_y(i, t, k) = as_scalar( x*beta_mean.col(k) );
                }
                else{
                    res = 0;
                    for(int j=0; j<t; j++)
                        for(int g=0; g<Num_of_attributes; g++){
                            if(AR_Order == 1)
                                phi_est = alpha_mean(k, g, 0)*(abs(t-j)==1);
                                //phi_est = alpha_mean(k, g, 0) + alpha_mean(k, g, 1)*abs(t-j) + alpha_mean(k, g, 2)*pow(t-j, 2); //alpha(k, g, 0)*(abs(t-j)==1) + alpha(k, g, 1)*(abs(t-j)==2)+ alpha(k, g, 2)*(abs(t-j)==3);//
                            if(AR_Order == 2)
                                phi_est = alpha_mean(k, g, 0)*(abs(t-j)==1) + alpha_mean(k, g, 1)*(abs(t-j)==2);
                            if(AR_Order == 3)
                                phi_est = alpha_mean(k, g, 0)*(abs(t-j)==1) + alpha_mean(k, g, 1)*(abs(t-j)==2)+ alpha_mean(k, g, 2)*(abs(t-j)==3);
                            if(AR_Order == 4)
                                phi_est = alpha_mean(k, g, 0)*(abs(t-j)==1) + alpha_mean(k, g, 1)*(abs(t-j)==2)+ alpha_mean(k, g, 2)*(abs(t-j)==3)+ alpha_mean(k, g, 3)*(abs(t-j)==4);
                            if(AR_Order == 5)
                                phi_est = alpha_mean(k, g, 0)*(abs(t-j)==1) + alpha_mean(k, g, 1)*(abs(t-j)==2)+ alpha_mean(k, g, 2)*(abs(t-j)==3)+ alpha_mean(k, g, 3)*(abs(t-j)==4)+ alpha_mean(k, g, 4)*(abs(t-j)==5);
                            if(AR_Order == 6)
                                phi_est = alpha_mean(k, g, 0) + alpha_mean(k, g, 1)*t/10 + alpha_mean(k, g, 2)*pow(t/10., 2.);
                            
                            x = X.tube(i, j);
                            res += phi_est * ( pred_y(i, j, g) - as_scalar(x*beta_mean.col(g)) );
                        }
                    x = X.tube(i, t);
                    pred_y(i, t, k) =as_scalar(x*beta_mean.col(k))+ res;
                }
    }
    
    
    
    // Calculate Measurements
    
    double logL=0;
    //Rcout << "nu = "  << nu_cand << endl;
    vec yi_tmp, yi;
    //after check, change iter to iter+1
    mat Ti, Di_mat_inv;
    //int tp;
    //Rcout << "check 1" << endl;
    
    for(int i=0; i<Num_of_obs; i++){
        tp = TimePointsAvailable(i);
        yi = vectorise((mat)Y(span(i), span(0, tp-1), span::all), 1).t();
        Ti = Ti_Mat(alpha_mean, i);
        
        //yi = vectorise((mat)Y.row(i), 1).t();
        
        yi_tmp = yi-X_all(i)*vectorise(beta_mean);
        //Dit_tmp = Dit(i, 1, nu_samples.col(iter), lambda_samples.slice(iter+1), 3);
        Di_mat_inv = Di(i, nu_mean, lambda_mean, tp).i();
        //Rcout << Di_mat << endl;
        //Rcout << "check 2" << endl;
        
        logL += -0.5*Di_mat_inv.n_cols*log(2*datum::pi)+ 0.5*log(det(Di_mat_inv)) - 0.5*as_scalar(yi_tmp.t()* Ti.t()*Di_mat_inv*Ti*yi_tmp);
        //Rcout << "check 3" << endl;
        
    }
    //Rcout << size(beta_mean) << endl;
    double dM = beta_mean.n_elem + alpha_mean.n_elem + lambda_mean.n_elem + nu_mean.n_elem;
    AIC = -2*logL + 2*dM;
    BIC = -2*logL + dM*log(Num_of_obs);

    double pD = 0.;

    mat SigmaMat_inv;
    vec mu_tmp; //, vec_tmp;
    //double num_of_data_points = accu(TimePointsAvailable);
    double logL_tmp;
    vec  CPO = zeros<vec>(Num_of_obs);
    
    for(int iter=500; iter<Num_of_iterations; iter++){
        for(int i=0; i<Num_of_obs; i++){
            tp = TimePointsAvailable(i);
            yi = vectorise((mat)Y(span(i), span(0, tp-1), span::all), 1).t();
            Ti = Ti_Mat(alpha_samples(iter), i);
            
            //yi = vectorise((mat)Y.row(i), 1).t();
            mu_tmp = X_all(i)*vectorise(beta_samples.slice(iter));
            yi_tmp = yi-mu_tmp;
            //Dit_tmp = Dit(i, 1, nu_samples.col(iter), lambda_samples.slice(iter+1), 3);
            Di_mat_inv = Di(i, nu_samples.col(iter), lambda_samples.slice(iter), tp).i();
            //Rcout << Di_mat << endl;
            //Rcout << "check 2" << endl;
            logL_tmp = -0.5*Di_mat_inv.n_cols*log(2*datum::pi)+0.5*log(det(Di_mat_inv)) - 0.5*as_scalar(yi_tmp.t()* Ti.t()*Di_mat_inv*Ti*yi_tmp);
            pD += logL_tmp/Num_of_iterations;
            //vec_tmp = dmvnorm(yi, mu_tmp, SigmaMat_inv.i());
            //MPL += //dmvnorm(yi, mu_tmp, SigmaMat_inv.i())(0)/(num_of_data_points*Num_of_iterations);
            

            CPO(i) += exp(-logL_tmp);
            //Rcout << "check 3" << endl;
        }
    }
    
    //Rcout << CPO(0) << "\t" << CPO(1) << "\t" << CPO(2) << endl;
    ///CPO = 1./CPO;
    //Rcout << CPO(0) << "\t" << CPO(1) << "\t" << CPO(2) << endl;
    
    MPL = Num_of_iterations*accu(CPO);
    //Rcout << "pD = " << pD << endl;
    //Rcout << "-4.*logL = " << -4.*logL << endl;
    DIC = 2*logL + -4*pD;
    
    //Rcout << "DIC = " << DIC << endl;
    
}


SEXP MLModelSelection::MCMC_Procedure()
{
    Rcout << "============= FMR: MCMC Starts=============="<< endl;
    List PosteriorSamples;
    List PosteriorEstimates;
    List Posterior;
    
    time_t start = time(NULL);
    
    int iter = 0;
    
    while(iter < Num_of_iterations-1){
        if(updatebeta)
            Update_beta(iter);
        else
            beta_samples.slice(iter+1) = beta_samples.slice(iter);
        
        
        if(updatedelta)
            Update_delta(iter);
        else
            delta_samples(iter+1) = delta_samples(iter);
        
        if(updatealpha)
            Update_alpha(iter);
        else
            alpha_samples(iter+1) = alpha_samples(iter);
        /*
        if(updatedelta & updatealpha)
            Update_delta_alpha(iter);
        else{
            alpha_samples(iter+1) = alpha_samples(iter);
            delta_samples(iter+1) = delta_samples(iter);
        }
        */
        if(updatelambda)
            Update_lambda(iter);
        else
            lambda_samples.slice(iter+1) = lambda_samples.slice(iter);
    
        if(updatenu)
            Update_nu(iter);
        else
            nu_samples.col(iter+1) = nu_samples.col(iter);

        iter++;
        if((iter+1)%100==0)
            Rcout << iter+1 << endl;
        //Rcout<<"\t"<< round((iter+1.)/Num_of_iterations*100)<<" %"<<'\r';
    }
    
    ParameterEstimation();
    
    Rcout << "============= FMR: MCMC: Done =============="<< endl;
    time_t end = time(NULL);
    Rcout<<"Execution Time: "<< (double)(end-start)<<" Seconds"<<std::endl;
    
    PosteriorSamples["beta.samples"] = beta_samples;
    PosteriorSamples["delta.samples"] = delta_samples;
    PosteriorSamples["alpha.samples"] = alpha_samples;
    PosteriorSamples["lambda.samples"] = lambda_samples;
    PosteriorSamples["nu.samples"] = nu_samples;
    
    PosteriorEstimates["beta.mean"] = beta_mean;
    PosteriorEstimates["delta.mean"] = delta_mean;
    PosteriorEstimates["alpha.mean"] = alpha_mean;
    PosteriorEstimates["lambda.mean"] = lambda_mean;
    PosteriorEstimates["nu.mean"] = nu_mean;
    PosteriorEstimates["AIC"] = AIC;
    PosteriorEstimates["BIC"] = BIC;
    PosteriorEstimates["DIC"] = DIC;
    PosteriorEstimates["MPL"] = MPL;
    PosteriorEstimates["pred.y"] = pred_y;
    
    Posterior["PosteriorSamples"] = PosteriorSamples;
    Posterior["PosteriorEstimates"] = PosteriorEstimates;
    
    return (Posterior);
}

