#include <Rcpp.h>
using namespace Rcpp;
#include "IncludeLib.h"

// [[Rcpp::export]]

double HMMll (Rcpp::List tuningparameters){

    std::vector<double> obs = Rcpp::as<std::vector<double> >(tuningparameters["Y"]);
    std::vector<double> mu_init = Rcpp::as<std::vector<double> >(tuningparameters["Mu"]);
    std::vector<double> sigma2_init = Rcpp::as<std::vector<double> >(tuningparameters["Sigma2"]);
    std::vector<double> pi_init = Rcpp::as<std::vector<double> >(tuningparameters["Pi"]);
    std::vector<std::vector<double> > A_init;
    std::vector<double> Atemp = Rcpp::as<std::vector<double> >(tuningparameters["A"]);
    int DimensionA = (int)sqrt( (double) Atemp.size());
    A_init.resize(DimensionA);
    int j, k, i = 0;
    for (j = 0; j < DimensionA; j ++){
        A_init[j].resize(DimensionA);
        for (k = 0; k < DimensionA; k ++){
            A_init[j][k] = Atemp[i];
            i ++;
        }
    }
    int T = (int)obs.size();
    HMM phmm;
    phmm.initialize(obs, DimensionA);
    phmm.set_HMM_num(T, DimensionA, mu_init, sigma2_init, pi_init, A_init);
    phmm.Forwardlog();
    return(phmm.get_logprobf());

}


// Input requires: Method, Kfit, Nstart, Ngibbs, Burnin, Thin, Y, MuInit, Sigma2Init, PiInit, nu, s2, alphaPi, MuPriorMean, MuPriorVar, AInit, AlphaPrior, updates2EM, verbose
// Output: resultHMM (EM result or posterior samples)

// [[Rcpp::export]]


Rcpp::NumericMatrix HMMfitting(Rcpp::List tuningparameters){

    int Method = Rcpp::as<int>(tuningparameters["Method"]);
    int K = Rcpp::as<int>(tuningparameters["Kfit"]);
    int N_starting_points = Rcpp::as<int>(tuningparameters["Nstart"]);
    int num_gibbs = Rcpp::as<int>(tuningparameters["Ngibbs"]);
    int burnin = Rcpp::as<int>(tuningparameters["Burnin"]);
    int thin = Rcpp::as<int>(tuningparameters["Thin"]);

    std::vector<double> obs = Rcpp::as<std::vector<double> >(tuningparameters["Y"]);
    int T = (int)obs.size();
    std::vector<double> mu_init = Rcpp::as<std::vector<double> >(tuningparameters["MuInit"]);
    std::vector<double> sigma2_init = Rcpp::as<std::vector<double> >(tuningparameters["Sigma2Init"]);
    std::vector<double> pi_init = Rcpp::as<std::vector<double> >(tuningparameters["PiInit"]);
    std::vector<double> nu = Rcpp::as<std::vector<double> >(tuningparameters["nu"]);
    std::vector<double> s2 = Rcpp::as<std::vector<double> >(tuningparameters["s2"]);
    std::vector<double> alpha_pi = Rcpp::as<std::vector<double> >(tuningparameters["alphaPi"]);
    std::vector<double> mu_prior_mean = Rcpp::as<std::vector<double> >(tuningparameters["MuPriorMean"]);
    std::vector<double> mu_prior_var = Rcpp::as<std::vector<double> >(tuningparameters["MuPriorVar"]);

    std::vector<std::vector<double> > A_init, alpha_prior;
    std::vector<double> Atemp = Rcpp::as<std::vector<double> >(tuningparameters["AInit"]);
    int DimensionA = (int)sqrt( (double) Atemp.size());
    A_init.resize(DimensionA);
    int j, k, i = 0;
    for (j = 0; j < DimensionA; j ++){
        A_init[j].resize(DimensionA);
        for (k = 0; k < DimensionA; k ++){
            A_init[j][k] = Atemp[i];
            i ++;
        }
    }
    Atemp = Rcpp::as<std::vector<double> >(tuningparameters["AlphaPrior"]);
    DimensionA = (int)sqrt( (double) Atemp.size());
    alpha_prior.resize(DimensionA);
    i = 0;
    for (j = 0; j < DimensionA; j ++){
        alpha_prior[j].resize(DimensionA);
        for (k = 0; k < DimensionA; k ++){
            alpha_prior[j][k] = Atemp[i];
            i ++;
        }
    }

    bool update_s2_EM = Rcpp::as<int>(tuningparameters["updates2EM"]);
    bool bool_verbose = Rcpp::as<int>(tuningparameters["verbose"]);

    std::vector<std::vector<double> > result;
    HMM phmm;
    phmm.initialize(obs, K);
    phmm.set_vec_nu(nu);
    phmm.set_vec_s2(s2);
    phmm.set_bool_update_s2(update_s2_EM);
    phmm.set_alpha_prior(alpha_prior);
    phmm.set_alpha_pi(alpha_pi);
    for(k = 0; k < K; k ++){
        phmm.set_muprior(k, mu_prior_mean[k]);
        phmm.set_mupriorvar(k, mu_prior_var[k]);
    }
    phmm.set_HMM_num(T, K, mu_init, sigma2_init, pi_init, A_init);
    phmm.set_mu_init(mu_init);
    phmm.set_pi_init(pi_init);
    phmm.set_sigma2_init(sigma2_init);
    phmm.set_A_init(A_init);
    result = phmm.fit_hmm_uneqvar(Method, obs, K, bool_verbose, N_starting_points, num_gibbs, burnin, thin);

    int M = result.size(), N = result[0].size();
    Rcpp::NumericMatrix resultHMM(M, N);
    for (k = 0; k < M; k ++)
        for(j = 0; j < N; j ++)
            resultHMM(k, j) = result[k][j];

    return resultHMM;
}


// Input list require: K, T, seed, Nrep, Nnoisetype
// Return list contains: obs, hidden, TransmatType, NoiseType

// [[Rcpp::export]]


Rcpp::List HMMrepsim(Rcpp::List tuningparameters){

    int K = Rcpp::as<int>(tuningparameters["K"]);
    int T = Rcpp::as<int>(tuningparameters["T"]);
    int SEED = Rcpp::as<int>(tuningparameters["seed"]);
    int Nrep = Rcpp::as<int>(tuningparameters["Nrep"]);
    int Nnoisetype = Rcpp::as<int>(tuningparameters["Nnoisetype"]);
    int para_type = ceil(1.0 * SEED / Nrep);
    int transmat_type = ceil(1.0 * para_type / Nnoisetype);
    int noise_type = para_type - Nnoisetype * (transmat_type - 1);
    SEED = SEED - Nrep * (para_type - 1);

    std::vector<double> mu_true, sigma2_true, pi_true;
    std::vector<std::vector<double> > A_true;

    int i, j;

    for(j = 0; j < K; j ++){
        mu_true.push_back(1.0 * (j + 1));
        pi_true.push_back(1.0 / K);
        sigma2_true.push_back(0.01 * (noise_type + 1) * (noise_type + 1));
    }

    if(transmat_type == 1){
        for(j = 0; j < K; j ++){
            A_true.push_back(std::vector<double>());
            for(i = 0; i < K; i ++){
                A_true[j].push_back(1.0/K);
            }
        }
    }
    else if(transmat_type == 2){
        for(j = 0; j < K; j ++){
            A_true.push_back(std::vector<double>());
            for(i = 0; i < K; i ++){
                if(i == j){
                    A_true[j].push_back(0.8);
                }
                else{
                    A_true[j].push_back(0.2 / (K-1));
                }
            }
        }
    }
    else if(transmat_type == 3){
        for(j = 0; j < K; j ++){
            A_true.push_back(std::vector<double>());
            for(i = 0; i < K; i ++){
                if(i == j){
                    A_true[j].push_back(0.95);
                }
                else{
                    A_true[j].push_back(0.05 / (K-1));
                }
            }
        }
    }
    else if(transmat_type == 4){
        for(j = 0; j < K; j ++){
            A_true.push_back(std::vector<double>());
            for(i = 0; i < K; i ++){
                if(i == j){
                    A_true[j].push_back(0.1);
                }
                else{
                    A_true[j].push_back(0.9 / (K-1));
                }
            }
        }
    }
    else{
        Rprintf("Wrong transition matrix type.\n");
        stop(0);
    }


    /* Simulate HMM */

    std::vector<double> Y;
    std::vector<int> X;
    HMM hmm_simulated;
    Y.resize(T, 0.0);
    hmm_simulated.initialize(Y, K);
    hmm_simulated.set_HMM_num(T, K, mu_true, sigma2_true, pi_true, A_true);
    hmm_simulated.hmm_sim (Y, X);


    return Rcpp::List::create(Rcpp::Named("obs") = Y, Rcpp::Named("hidden") = X, Rcpp::Named("TransmatType") = transmat_type, Rcpp::Named("NoiseType") = noise_type, Rcpp::Named("seed") = SEED);


}



// Input list require: T, Mu, Sigma2, Pi, A
// Return list contains: obs, hidden

// [[Rcpp::export]]

Rcpp::List HMMsimulate(Rcpp::List tuningparameters){

    int T = Rcpp::as<int>(tuningparameters["T"]);

    std::vector<double> mu_true = Rcpp::as<std::vector<double> >(tuningparameters["Mu"]);
    std::vector<double> sigma2_true = Rcpp::as<std::vector<double> >(tuningparameters["Sigma2"]);
    std::vector<double> pi_true = Rcpp::as<std::vector<double> >(tuningparameters["Pi"]);

    std::vector<std::vector<double> > A_true;
    std::vector<double> Atemp = Rcpp::as<std::vector<double> >(tuningparameters["A"]);
    int K = (int)sqrt( (double) Atemp.size());
    A_true.resize(K);
    int j, k, i = 0;
    for (j = 0; j < K; j ++){
        A_true[j].resize(K);
        for (k = 0; k < K; k ++){
            A_true[j][k] = Atemp[i];
            i ++;
        }
    }

    /* Simulate HMM */

    std::vector<double> Y;
    std::vector<int> X;
    HMM hmm_simulated;
    Y.resize(T, 0.0);
    hmm_simulated.initialize(Y, K);
    hmm_simulated.set_HMM_num(T, K, mu_true, sigma2_true, pi_true, A_true);
    hmm_simulated.hmm_sim (Y, X);


    return Rcpp::List::create(Rcpp::Named("obs") = Y, Rcpp::Named("hidden") = X);


}



