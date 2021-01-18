

class HMM {
public:
    HMM(){};
    ~HMM(){};
    void initialize(std::vector<double>, int);
    void set_HMM_num(int, int, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<std::vector<double> >);
    void hmm_sim(std::vector<double>&, std::vector<int>&);
    void calculateBmatrix();
    void Backwardlog();
    void Forwardlog();
    void ComputeGamma();
    void ComputeXi();
    void BaumWelch(bool);
    void BaumWelch_multi_starting_point_fitting(int, std::vector<double>, std::vector<double>, std::vector<std::vector<double> >, bool);
    void BaumWelch_after_gibbs (std::vector<std::vector<double> >, bool);
    void samplenewparameters(bool);
    void samplehiddenstates(bool);
    void samplenewparameters_gm(bool);
    void samplehiddenstates_gm(bool);
    double calculate_full_likelihood();
    double calculate_marginal_posterior();
    double calculate_marginal_posterior_gm();
    std::vector<std::vector<double> > gibbs_sampling(int, int, int, bool);
    std::vector<std::vector<double> > gibbs_sampling_gm(int, int, int, bool);
    std::vector<std::vector<double> > fit_hmm_uneqvar(int, std::vector<double>, int, bool, int, int, int, int);

    int get_N(){ return(N); }
    void set_N(int Nset){ N = Nset; }
    int get_M(){ return(M); }
    void set_M(int Mset){ M = Mset; }
    void resize_vec_s2(int a){ s2.resize(a); }
    double get_vec_s2(int a){ return(s2[a]); }
    void set_s2(int a, double set_s2_val){ s2[a] = set_s2_val; }
    void set_vec_s2(std::vector<double> set_s2_val){
        s2.resize(set_s2_val.size());
        for(int k = 0; k < (int)s2.size(); k ++)
            s2[k] = set_s2_val[k];
    }
    void resize_vec_nu(int a){ nu.resize(a); }
    double get_vec_nu(int a){ return(nu[a]); }
    void set_nu(int a, double nu_init){ nu[a] = nu_init; }
    void set_vec_nu(std::vector<double> nu_init){
        nu.resize(nu_init.size());
        for(int k = 0; k < (int)nu.size(); k ++)
            nu[k] = nu_init[k];
    }
    double get_A(int i, int j){ return(A[i][j]); }
    void set_A(int i, int j, double aset){ A[i][j] = aset; }
    double get_B(int i, int j){ return(B[i][j]); }
    void set_B(int i, int j, double bset){ B[i][j] = bset; }
    double get_logbeta(int i, int j){ return(logbeta[i][j]); }
    void set_logbeta(int i, int j, double betaset){ logbeta[i][j] = betaset; }
    double get_logalpha(int i, int j){ return(logalpha[i][j]); }
    void set_logalpha(int i, int j, double logalphaset){ logalpha[i][j] = logalphaset; }
    double get_loggamma(int i, int j){ return(loggamma[i][j]); }
    void set_loggamma(int i, int j, double loggammaset){ loggamma[i][j] = loggammaset; }
    double get_logxi(int i, int j, int k){ return(logxi[i][j][k]); }
    void set_logxi(int i, int j, int k, double logxiset){ logxi[i][j][k] = logxiset; }
    double get_muprior(int i){ return(muprior[i]); }
    void set_muprior(int i, double mupriorset){ muprior[i] = mupriorset; }
    double get_mupriorvar(int i){ return(mupriorvar[i]); }
    void set_mupriorvar(int i, double mupriorvarset){ mupriorvar[i] = mupriorvarset; }
    double get_pi(int i){ return(pi[i]); }
    void set_pi(int i, double piset){ pi[i] = piset; }
    double get_mu(int i){ return(mu[i]); }
    void set_mu(int i, double muset){ mu[i] = muset; }
    double get_sigma2(int i){ return(sigma2[i]); }
    void set_sigma2(int i, double sigma2set){ sigma2[i] = sigma2set; }
    double get_logprobf() {return(logprobf); }
    void set_logprobf(double logprobfset){ logprobf = logprobfset; }
    int get_z(int t){ return(z[t]); }
    void set_z(int t, int i){ z[t] = i; }
    double get_O(int i){ return(O[i]); }
    void set_Obs(std::vector<double> O_obs){
        O.clear();
        for(int i = 0; i < (int)O_obs.size(); i ++)
            O.push_back(O_obs[i]);
        Rprintf("set\t observations\n");
        return;
    }
    void resize_mu(int a){ mu.resize(a); }
    void resize_sigma2(int a){ sigma2.resize(a); }
    void resize_pi(int a){ pi.resize(a); }
    void resize_A(int a){
        A.resize(a);
        for(int i = 0; i < a; i ++)
            A[i].resize(a);
        return;
    }
    void resize_muprior(int a){ muprior.resize(a); }
    void resize_mupriorvar(int a){ mupriorvar.resize(a); }
    void resize_B(int a, int b){
        B.resize(a);
        for(int i = 0; i < a; i ++)
            B[i].resize(b);
        return;
    }
    void resize_logbeta(int a, int b){
        logbeta.resize(a);
        for(int i = 0; i < a; i ++)
            logbeta[i].resize(b);
        return;
    }
    void resize_logalpha(int a, int b){
        logalpha.resize(a);
        for(int i = 0; i < a; i ++)
            logalpha[i].resize(b);
        return;
    }
    void resize_loggamma(int a, int b){
        loggamma.resize(a);
        for(int i = 0; i < a; i ++)
            loggamma[i].resize(b);
        return;
    }
    void resize_logxi(int a, int b, int c){
        logxi.resize(a);
        for(int i = 0; i < a; i ++){
            logxi[i].resize(b);
            for(int j = 0; j < b; j ++){
                logxi[i][j].resize(c);
            }
        }
        return;
    }
    void resize_alpha_prior(int K){
        alpha_prior.resize(K);
        for(int k = 0; k < K; k ++){
            alpha_prior[k].resize(K);
        }
    }
    void resize_alpha_pi(int K){
        alpha_prior.resize(K);
    }
    void set_alpha_prior(std::vector<std::vector<double> > alpha_prior_init){
        int K = (int)alpha_prior_init.size();
        alpha_prior.resize(K);
        for(int k = 0; k < K; k ++){
            alpha_prior[k].resize(K);
            for(int j = 0; j < K; j ++){
                alpha_prior[k][j] = alpha_prior_init[k][j];
            }
        }
    }
    void set_alpha_pi(std::vector<double> alpha_pi_init){
        int K = (int)alpha_pi_init.size();
        alpha_pi.resize(K);
        for(int k = 0; k < K; k ++){
            alpha_pi[k] = alpha_pi_init[k];
        }
    }
    void set_mu_init(std::vector<double> temp){
        mu_init = temp;
    }
    void set_sigma2_init(std::vector<double> temp){
        sigma2_init = temp;
    }
    void set_pi_init(std::vector<double> temp){
        pi_init = temp;
    }
    void set_A_init(std::vector<std::vector<double> > temp){
        A_init = temp;
    }
    double get_o(int i){ return(O[i]); }
    bool get_bool_update_s2(){ return(bool_update_s2); }
    void set_bool_update_s2(bool update){ bool_update_s2 = update; }


private:
    /* N: number of states; M: number of observations */
    /* A: transition matrix */
    /* B: LOG probability density of observed values at different states*/
    /* pi: starting point distribution */
    /* logprobf: marginal log likelihood p(observation|parameters), from forward probability */
    int M, N;
    double logprobf, MAXob, MINob;
    std::vector<std::vector<double> > A, B, logbeta, logalpha, loggamma, alpha_prior;
    std::vector<std::vector<std::vector<double> > > logxi;
    std::vector<double> O, muprior, mupriorvar, mu, sigma2, pi, nu, s2, alpha_pi;
    std::vector<double> mu_init, sigma2_init, pi_init;
    std::vector<std::vector<double> > A_init;
    std::vector<int> z;
    bool bool_update_s2;
};

/* set starting values in a deterministic or random way */

void set_starting_value(std::vector<double> &mu_init, int K, std::vector<double> Obs, bool bool_equal_space){

    int i;

    std::vector<double> O = Obs;
    mu_init.resize(K);
    std::sort(O.begin(), O.end());
    if(bool_equal_space){
        double f = 0.05;
        for(i = 0; i < K; i ++){
            mu_init[i] = quantileinline (O, f);
            f += 0.9 / (K - 1);
        }
    }
    else{
        std::vector<double> f_vec;
        GetRNGstate();
        for(i = 0; i < K; i ++){
            f_vec.push_back(runiforminline());
        }
        std::sort(f_vec.begin(), f_vec.end());
        for(i = 0; i < K; i ++){
            mu_init[i] = quantileinline (O, f_vec[i]);
        }
    }
    return;

}

/* write a 2 dimension array to file */

void output2dimarray(std::string filename, std::vector<std::vector<double> > results){

    int i, j;
    FILE *output;
    const char * c;
    c = filename.c_str();
    output = fopen(c, "a");
    for(i = 0; i < (int)results.size(); i ++){
        for(j = 0; j < (int)results[i].size(); j ++)
            fprintf(output, "%f\t", results[i][j]);
        fprintf(output, "\n");
    }
    fclose(output);

    return;
}

/* initialization of HMM with observations O_init and number of hidden states K */
/* set Minob and Maxob (the lower and upper limits of observations) at initialization */
/* MINob and MAXob will be used for the lower and upper bounds of the means of hidden states */

void HMM::initialize(std::vector<double> O_init, int K){

    M = O_init.size();
    O.clear();
    for(int i = 0; i < M; i ++){
        O.push_back(O_init[i]);
        if (i == 0){
            MAXob = O[i];
            MINob = O[i];
        }
        else{
            MAXob = MAX(MAXob, O[i]);
            MINob = MIN(MINob, O[i]);
        }
    }
    N = K;
    mu.resize(N);
    sigma2.resize(N);
    pi.resize(N);
    A.resize(N);
    for (int j = 0; j < N; j ++)
        A[j].resize(N);
    muprior.resize(N);
    mupriorvar.resize(N);
    nu.resize(K);
    s2.resize(K);
    return;

}


/* set numerical values of HMM */

void HMM::set_HMM_num(int T, int K, std::vector<double> mu_num, std::vector<double> sigma2_num, std::vector<double> pi_num, std::vector<std::vector<double> > A_num){

    M = T;
    N = K;
    mu.resize(K);
    pi.resize(K);
    sigma2.resize(K);
    A.resize(K);
    for(int k = 0; k < K; k ++){
        A[k].resize(K);
        mu[k] = mu_num[k];
        sigma2[k] = sigma2_num[k];
        pi[k] = pi_num[k];
        for(int j = 0; j < K; j ++){
            A[k][j] = A_num[k][j];
        }
    }
    return;
}


/* simulate from an HMM */
/* sim_hidden: M-vector of hidden states 0 to K - 1 (states 1 to K) */
/* need pi, A, mu, sigma */

void HMM::hmm_sim(std::vector<double> &sim_data, std::vector<int> &sim_hidden){

    double r;
    int i, j;
    double sum;
    sim_data.resize(M);
    sim_hidden.resize(M);

    GetRNGstate();

    /* Simulate Hidden States */

    r = runiforminline();
    sim_hidden[0] = 0;
    sum = pi[0];
    for(i = 1; i < N; i ++){
        if(r > sum){
            sim_hidden[0] ++;
            sum += pi[i];
        }
    }

    for(i = 1; i < M; i ++){
        r = runiforminline();
        sim_hidden[i] = 0;
        sum = A[sim_hidden[i - 1]][0];
        for(j = 1; j < N; j ++){
            if(r > sum){
                sim_hidden[i] ++;
                sum += A[sim_hidden[i - 1]][j];
            }
        }
    }

    /* Sample Observations */

    for(i = 0; i <  M; i ++){
        sim_data[i] = mu[sim_hidden[i]] + sqrt(sigma2[sim_hidden[i]]) * rnorminline();
    }
    return;
}

/* calculate the Gaussian densities at each hidden state */
/* B is N*M matrix with the (i, j) entry the log-Gaussian-density */
/* of the jth observation evaluated at the ith hidden state */
/* require O, mu, sigma */

void HMM::calculateBmatrix()
{
    int i, t;
    B.resize(N);
    for(i = 0; i < N; i ++){
        B[i].resize(M);
        for(t = 0; t < M; t ++)
            B[i][t] = ldnorminline(O[t], mu[i], sqrt(sigma2[i]));
    }
    return;
}

/* Backward algorithm */
/* requires the B matrix */
/* logbeta is a M*N matrix */

void HMM::Backwardlog()
{
    int i, j, t;
    if((int)B.size() != N || (int)B[0].size() != M){
        Rprintf("Backwardlog:\trun\t calculateBmatrix\t first\n");
        stop(0);
    }
    double sum;
    double maxtemp=0.0;
    std::vector<double> temp;
    temp.resize(N);
    logbeta.resize(M);
    for(i = 0; i < M; i ++)
        logbeta[i].resize(N);

    for(i = 0; i < N; i++)
        logbeta[M - 1][i] = 0.0;
    for(t = M - 1; t > 0; t--){
        for(i = 0; i < N; i ++){
            sum = 0.0;
            for(j = 0; j < N; j ++){
                temp[j] = log(A[i][j]) + B[j][t] + logbeta[t][j];
                if(j == 0) maxtemp = temp[j];
                else maxtemp = MAX(maxtemp, temp[j]);
            }
            for(j = 0; j < N; j ++)
                sum += exp(temp[j] - maxtemp);
            logbeta[t - 1][i] = log(sum) + maxtemp;
        }
    }
    return;
}

/* forward algorithm */
/* requires the B matrix */
/* logalpha is a M*N matrix */

void HMM::Forwardlog()
{
    int i, j, k, t;
    if((int)B.size() != N || (int)B[0].size() != M){
        Rprintf("Forwardlog:\trun\t calculateBmatrix\t first\n");
        stop(0);
    }
    std::vector<double> temp;
    temp.resize(N);
    logalpha.resize(M);
    for (i = 0; i < M; i ++){
        logalpha[i].resize(N);
    }
    double sum;
    double maxtemp = 0.0;
    for(i = 0; i < N; i ++){
        if(pi[i] > 0.0000000000001){
            logalpha[0][i] = log(pi[i]) + B[i][0];
        }
        else{
            logalpha[0][i] = log(pi[i] + 0.0000000000001) + B[i][0];
        }
    }
    for(t = 1; t < M; t ++){
        for(i = 0; i < N; i ++){
            sum = 0.0;
            for(j = 0; j < N; j ++){
                temp[j] = log(A[j][i]) + logalpha[t - 1][j];
                if(j == 0) maxtemp = temp[j];
                else maxtemp = MAX(maxtemp,temp[j]);
            }
            for(j = 0; j < N; j ++)
                sum += exp(temp[j] - maxtemp);
            logalpha[t][i] = log(sum) + maxtemp + B[i][t];
        }
    }
    double logprobftemp = 0.0;
    for(j = 0; j < N; j ++){
        if(j == 0) maxtemp = logalpha[M - 1][j];
        else maxtemp = MAX(maxtemp, logalpha[M - 1][j]);
    }
    for(k = 0; k < N; k ++)
        logprobftemp += exp(logalpha[M - 1][k] - maxtemp);
    logprobf = log(logprobftemp) + maxtemp;
    return;
}

/* Compute gamma (require forward and backward first) */
/* loggamma is a M*N matrix */

void HMM::ComputeGamma()
{
    int i, j, k, t;
    if((int)logalpha.size() != M || (int)logalpha[0].size() != N){
        Rprintf("ComputeGamma:\t run\t Forwardlog\t first\n");
        stop(0);
    }
    if((int)logbeta.size() != M || (int)logbeta[0].size() != N){
        Rprintf("ComputeGamma:\t run\t Backwardlog\t first\n");
        stop(0);
    }
    loggamma.resize(M);
    for(i = 0; i < M; i ++)
        loggamma[i].resize(N);
    double denominator;
    double maxtemp = 0.0;
    for(t = 0; t < M; t ++){
        denominator = 0.0;
        for(j = 0; j < N; j ++){
            loggamma[t][j] = logalpha[t][j] + logbeta[t][j];
            if(j == 0) maxtemp = loggamma[t][j];
            else maxtemp = MAX(maxtemp, loggamma[t][j]);
        }
        for(k = 0; k < N; k ++)
            denominator += exp(loggamma[t][k] - maxtemp);
        for(i = 0; i < N; i ++)
            loggamma[t][i] = loggamma[t][i] - log(denominator) - maxtemp;
    }
    return;
}

/* compute xi, require forward, backward and B matrix */
/* logxi is (M-1)*N*N array */

void HMM::ComputeXi()
{
    int i, j, t;
    if((int)B.size() != N || (int)B[0].size() != M){
        Rprintf("ComputeXi:\trun\t calculateBmatrix\t first\n");
        stop(0);
    }
    if((int)logalpha.size() != M || (int)logalpha[0].size() != N){
        Rprintf("ComputeXi:\t run\t Forwardlog\t first\n");
        stop(0);
    }
    if((int)logbeta.size() != M || (int)logbeta[0].size() != N){
        Rprintf("ComputeXi:\t run\t Backwardlog\t first\n");
        stop(0);
    }
    double sum, maxtemp;
    logxi.resize(M - 1);
    for(i = 0; i < M - 1; i ++){
        logxi[i].resize(N);
        for (j = 0; j < N; j ++){
            logxi[i][j].resize(N);
        }
    }
    maxtemp = logalpha[0][0] + logbeta[1][0] + log(A[0][0]) + B[0][1];
    for(t = 0; t < M - 1; t ++){
        sum = 0.0;
        for(i = 0; i < N; i ++){
            for(j = 0; j < N; j ++){
                logxi[t][i][j] = logalpha[t][i] + logbeta[t + 1][j] +
                log(A[i][j]) + B[j][t + 1];
                maxtemp = MAX(maxtemp, logxi[t][i][j]);
            }
        }
        for(i = 0; i < N; i ++)
            for(j = 0; j < N; j ++)
                sum += exp(logxi[t][i][j] - maxtemp);
        for(i = 0; i < N; i ++)
            for(j = 0; j < N; j ++)
                logxi[t][i][j] = logxi[t][i][j] - log(sum) - maxtemp;
    }
    return;
}

/* Baum Welch algorithm */

#define DELTA 0.000001

void HMM::BaumWelch(bool verbose)
{
    int i, j, t, niter = 0;
    double tempsum, numeratorA, denominatorA, numeratorB, denominatorB, numeratorC, delta, logprobprev;
    std::vector<double> temp;
    temp.resize(N);

    if(verbose){
        Rprintf("Initializing Baum Welch\n");}

    calculateBmatrix();
    Forwardlog();
    Backwardlog();
    ComputeGamma();
    ComputeXi();
    logprobprev = logprobf;

    Rprintf("begin Baum Welch\n");

    do{
        tempsum = 0.0;
        for(i = 0; i < N; i ++){
            pi[i] = exp(loggamma[0][i]);
            tempsum += pi[i];
        }
        for(i = 0; i < N; i ++)
            pi[i] = pi[i] / tempsum;
        for(i = 0; i < N; i ++){
            denominatorA = 0.0;
            for(t = 0; t < M - 1; t++)
                denominatorA += exp(loggamma[t][i]);
            for(j = 0; j < N; j ++){
                numeratorA = 0.0;
                for(t = 0; t < M - 1; t++)
                    numeratorA += exp(logxi[t][i][j]);
                if(numeratorA/denominatorA < 10e-7){
                    A[i][j] = 0.000001 + 0.999999 * numeratorA/denominatorA;
                }
                else{
                    A[i][j] = numeratorA/denominatorA;
                }
            }
            denominatorB = denominatorA + exp(loggamma[M - 1][i]);
            numeratorB = 0.0;
            numeratorC = 0.0;
            for(t = 0; t < M; t ++){
                numeratorB += O[t] * exp(loggamma[t][i]);
                numeratorC += pow(O[t] - mu[i], 2.0) * exp(loggamma[t][i]);
            }
            mu[i] = numeratorB/denominatorB;
            sigma2[i] = (numeratorC + nu[i] * s2[i]) / (denominatorB + nu[i] + 2.0);

        }
        std::sort(mu.begin(), mu.end());

        calculateBmatrix();
        Forwardlog();
        Backwardlog();
        ComputeGamma();
        ComputeXi();
        delta = logprobf - logprobprev;
        logprobprev = logprobf;
        niter ++;
        if(verbose){ Rprintf("%d\t %f\n", niter, delta);}
    }
    while (delta > DELTA || niter < 10);

    Rprintf("Baum-Welch\t Results:\n");
    for(i = 0; i < N; i ++){
        Rprintf("%f\t", mu[i]);
    }
    Rprintf("\n");
    for(i = 0; i < N; i ++){
        Rprintf("%f\t", sigma2[i]);
    }
    Rprintf("\n");
    for(i = 0; i < N; i ++){
        for(j = 0; j < N; j ++){
            Rprintf("%f\t", A[i][j]);
        }
        Rprintf("\n");
    }
    Rprintf("LL %f;\t s2:\n", logprobf);
    for(i = 0; i < N; i ++)
        Rprintf("%f\t", s2[i]);
    Rprintf("\n");

    return;
}

/* Multiple starting point -- Baum Welch */

void HMM::BaumWelch_multi_starting_point_fitting(int N_starting_points, std::vector<double> sigma2_init_input, std::vector<double> pi_init_input, std::vector<std::vector<double> > A_init_input, bool verbose){

    int i, j, k;
    double current_max_ll = 0.0;
    std::vector<double> current_max_mu, current_max_sigma2, current_max_pi, current_max_s2, mu_init_temp;
    std::vector<std::vector<double> > current_max_A;

    if(N_starting_points > 1){
        for(i = 0; i < N_starting_points; i ++){
            set_starting_value(mu_init_temp, N, O, 0);
            set_HMM_num((int)O.size(), N, mu_init_temp, sigma2_init_input, pi_init_input, A_init_input);
            BaumWelch(verbose);
            if(i == 0){
                current_max_ll = logprobf;
                current_max_mu.resize(N);
                current_max_pi.resize(N);
                current_max_s2.resize(N);
                current_max_sigma2.resize(N);
                current_max_A.resize(N);
                for(j = 0; j < N; j ++)
                    current_max_A[j].resize(N);

            }
            if(i == 0 || current_max_ll < logprobf){
                for(k = 0; k < N; k ++){
                    current_max_mu[k] = mu[k];
                    current_max_sigma2[k] = sigma2[k];
                    current_max_s2[k] = s2[k];
                    current_max_pi[k] = pi[k];
                    for(j = 0; j < N; j ++)
                        current_max_A[k][j] = A[k][j];
                }
                current_max_ll = logprobf;
            }
        }
        logprobf = current_max_ll;
        for(k = 0; k < N; k ++){
            mu[k] = current_max_mu[k];
            pi[k] = current_max_pi[k];
            sigma2[k] = current_max_sigma2[k];
            s2[k] = current_max_s2[k];
            for(j = 0; j < N; j ++)
                A[k][j] = current_max_A[k][j];
        }
    }
    else{
        set_starting_value(mu_init_temp, N, O, 1);
        set_HMM_num((int)O.size(), N, mu_init_temp, sigma2_init_input, pi_init_input, A_init_input);
        BaumWelch(verbose);
    }
    return;

}

/* Taking Gibbs samples result as starting point of Baum Welch */

void HMM::BaumWelch_after_gibbs (std::vector<std::vector<double> > hmm_gibbs_result, bool verbose){

    int Ngibbs = (int)hmm_gibbs_result.size(), i, j;

    mu_init.resize(N, 0.0);
    sigma2_init.resize(N, 0.0);
    pi_init.resize(N, 1.0 / N);
    A_init.resize(N);
    for(i = 0; i < N; i ++)
        A_init[i].resize(N, 1.0 / N);
    for(i = 0; i < Ngibbs; i ++){
        for(j = 0; j < N; j ++){
            mu_init[j] += hmm_gibbs_result[i][j];
            sigma2_init[j] += hmm_gibbs_result[i][j + N];
        }
    }
    for(j = 0; j < N; j ++){
        mu_init[j] /= Ngibbs;
        sigma2_init[j] /= Ngibbs;
    }

    set_HMM_num(M, N, mu_init, sigma2_init, pi_init, A_init);
    BaumWelch(verbose);

    return;

}


/* sample hidden states */
/* INPUT: B matrix calculated and parameter values set */


void HMM::samplehiddenstates (bool verbose){

    GetRNGstate();

    double temp_max, r, sum;
    int j, k;
    std::vector<double> temp;
    temp.resize(N);
    z.resize(M);

    calculateBmatrix();
    Backwardlog();
    Forwardlog();
    ComputeGamma();

    /* Update z0 */

    sum = 0.0;
    for(j = 0; j < N; j ++){
        temp[j] = exp(loggamma[0][j]);
        sum += temp[j];
    }
    if(sum < 0.999 || sum > 1.001){
        Rprintf("samplehiddenstates:gamma\t is\t incorrect.\n");
        stop(0);
    }
    r = runiforminline();
    z[0] = 0;
    sum = temp[0];
    for(j = 1; j < N; j ++){
        if(r > sum){
            z[0] ++;
            sum += temp[j];
        }
        else{
            break;
        }
    }

    /* Update z1,...,zn sequentially */

    for(k = 1; k < M; k++){

        z[k] = 0;
        for(j = 0; j < N; j++){
            temp[j] = B[j][k] + logbeta[k][j] + log(A[z[k - 1]][j]);
            if(j == 0)
                temp_max = temp[j];
            else
                temp_max = MAX (temp_max, temp[j]);
        }
        sum = 0.0;
        for(j = 0; j < N; j++){
            temp[j] = exp(temp[j] - temp_max);
            sum += temp[j];
        }
        for(j = 0; j < N; j++){
            temp[j] = temp[j] / sum;
        }
        r = runiforminline();
        z[k] = 0;
        sum = temp[0];
        for(j = 1; j< N; j++){
            if(r > sum){
                z[k] ++;
                sum += temp[j];
            }
            else{
                break;
            }
        }
    }
    if(verbose){
        for(j = 0; j < M; j ++)
            Rprintf("%d\t", z[j]);
        Rprintf("\n");
    }
    return;
}

/* sample parameters */
/* INPUT: need to have hidden states and old parameters set */

void HMM::samplenewparameters (bool verbose)
{

    int i, j, k;
    GetRNGstate();
    std::vector<double> tempA, sum;
    std::vector<int> num;
    num.resize(N);
    std::vector<std::vector<double> > count;
    count.resize(N);
    for(i = 0; i < N; i++){
        count[i].resize(N);
    }

    /* Sample pi */
    //  keep pi fixed

    /* Sample transition matrix */

    for(j = 0; j < N; j++){
        for(k = 0; k < N; k++){
            count[j][k] = alpha_prior[j][k];
        }
    }
    for(i = 1; i < M; i++)
        count[z[i - 1]][z[i]] += 1.0;
    for(i = 0; i < N; i++){
        tempA = randomdirichlet (count[i]);
        for(j = 0;  j < N; j++){
            A[i][j] = tempA[j];
        }
    }

    /* Sample mu, sigma2 */
    sum.resize(N);
    for(i = 0; i < N; i++){
        sum[i] = 0.0;
        num[i] = 0;
    }
    for(i = 0; i < M; i++){
        sum[z[i]] += O[i];
        num[z[i]] ++;
    }
    for(j = 0; j < N; j++){
        mu[j] = (sum[j]/ sigma2[j] + muprior[j] / mupriorvar[j]) / (num[j] / sigma2[j] + 1.0 / mupriorvar[j]) + rnorminline() / sqrt(1.0 * num[j] / sigma2[j] + 1.0 / mupriorvar[j]);
        i = 0;
        while (mu[j] > MAXob || mu[j] < MINob){
            if (i < 20){
                mu[j] = (sum[j]/ sigma2[j] + muprior[j] / mupriorvar[j]) / (num[j] / sigma2[j] + 1.0 / mupriorvar[j]) + rnorminline() / sqrt(1.0 * num[j] / sigma2[j] + 1.0 / mupriorvar[j]);
                i ++;
            }
            else{
                mu[j] = muprior[j];
            }
        }
    }
    std::sort(mu.begin(), mu.end());
    for(i = 0; i < N; i++)
        sum[i] = 0.0;
    for(i = 0; i < M; i++)
        sum[z[i]] += pow(O[i] - mu[z[i]], 2.0);
    for(j = 0; j < N; j++){
        sigma2[j] = (sum[j] + nu[j] * s2[j]) / randomchisq(nu[j] + num[j]);
    }

    if(verbose){
        for(i = 0; i < N ; i ++)
            Rprintf("%f\t", mu[i]);
        Rprintf("\n");
        for(i = 0; i < N; i ++)
            Rprintf("%f\t", sigma2[i]);
        Rprintf("\n");
        for(i = 0; i < N; i ++){
            for(j = 0; j < N; j ++){
                Rprintf("%f\t", A[i][j]);
            }
            Rprintf("\n");
        }
    }

    return;

}

/* calculate joint likelihood */
/* need model parameters and hidden states */

double HMM::calculate_full_likelihood(){

    double LL = 0.0, temp;

    int j, t;

    /* likelihood part */

    for(t = 0; t < M; t ++){
        LL += ldnorminline(O[t], mu[z[t]], sqrt(sigma2[z[t]]));
        if(t > 0){
            LL += log(A[z[t - 1]][z[t]]);
        }
    }

    /* prior part */

    for(j = 0; j < N; j++){
        // mu
        LL += ldnorminline(mu[j], muprior[j], sqrt(mupriorvar[j]));
        // sigma2
        LL += ldscaleinvchisq(sigma2[j], nu[j], s2[j]);
    }
    // transition matrix
    for(j = 0; j < N; j ++){
        temp = 0.0;
        for(t = 0; t < N; t ++){
            LL += (alpha_prior[j][t] - 1.0) * log(A[j][t]) - std::lgamma(alpha_prior[j][t]);
            temp += alpha_prior[j][t];
        }
        LL += std::lgamma(temp);
    }

    return(LL);
}

/* calculate marginal posterior */
/* need model parameters */

double HMM::calculate_marginal_posterior(){

    calculateBmatrix();
    Forwardlog();
    int j, t;
    double LL = logprobf, temp;
    for(j = 0; j < N; j++){
        // mu
        LL += ldnorminline(mu[j], muprior[j], sqrt(mupriorvar[j]));
        // sigma2
        LL += ldscaleinvchisq(sigma2[j], nu[j], s2[j]);
    }
    // transition matrix
    for(j = 0; j < N; j ++){
        temp = 0.0;
        for(t = 0; t < N; t ++){
            LL += (alpha_prior[j][t] - 1.0) * log(A[j][t]) - std::lgamma(alpha_prior[j][t]);
            temp += alpha_prior[j][t];
        }
        LL += std::lgamma(temp);
    }
    return(LL);

}

/* Gibbs sampling function */
/* output is 2dim array, each row corresponds to a parameter */
/* the transition matrix only contains off-diagonal elements */
/* last column records the log posterior (marginal) */

std::vector<std::vector<double> > HMM::gibbs_sampling(int num_gibbs, int burnin, int thin, bool verbose){

    int j, k;
    std::vector<std::vector<double> > temp_output;

    for(int iter = 0; iter < burnin + num_gibbs * thin; iter ++){
        if(verbose){
            if(iter >= burnin && (iter - burnin) % thin == 0){
                samplehiddenstates(1);
                samplenewparameters(1);
            }
            else{
                samplehiddenstates(0);
                samplenewparameters(0);
            }
        }
        else{
            samplehiddenstates(0);
            samplenewparameters(0);
        }
        if(iter >= burnin && (iter - burnin) % thin == 0){
            temp_output.push_back(std::vector<double>());
            for(k = 0; k < N; k ++)
                temp_output[(int)temp_output.size() - 1].push_back(mu[k]);
            for(k = 0; k < N; k ++)
                temp_output[(int)temp_output.size() - 1].push_back(sigma2[k]);
            for(k = 0; k < N; k ++){
                for(j = 0; j < N; j ++){
                    if(j != k)
                        temp_output[(int)temp_output.size() - 1].push_back(A[k][j]);
                }
            }
            temp_output[(int)temp_output.size() - 1].push_back(calculate_marginal_posterior());
        }
    }
    return(temp_output);
}


/* Gaussian mixture sample hidden states */

void HMM::samplehiddenstates_gm (bool verbose){

    GetRNGstate();
    double temp_max, r, sum;
    int j, k;
    std::vector<double> temp;
    temp.resize(N);
    z.resize(M);

    for(k = 0; k < M; k ++){
        for(j = 0; j < N; j ++){
            temp[j] = log(pi[j]) + ldnorminline(O[k], mu[j], sqrt(sigma2[j]));
            if(j == 0){
                temp_max = temp[j];
            }
            else{
                temp_max = MAX(temp_max, temp[j]);
            }
        }
        sum = 0.0;
        for(j = 0; j < N; j ++){
            temp[j] = exp(temp[j] - temp_max);
            sum += temp[j];
        }
        for(j = 0; j < N; j ++){
            temp[j] = temp[j] / sum;
        }
        r = runiforminline();
        z[k] = 0;
        sum = temp[0];
        for(j = 1; j< N; j++){
            if(r > sum){
                z[k]++;
                sum += temp[j];
            }
        }
    }
    if(verbose){
        for(j = 0; j < M; j ++)
            Rprintf("%d\t", z[j]);
        Rprintf("\n");
    }
    return;
}

/* Gaussian mixture sample new parameters */


void HMM::samplenewparameters_gm (bool verbose)
{

    int i, j;
    GetRNGstate();
    std::vector<double> count, sum, tempA;
    std::vector<int> num;
    num.resize(N);
    if (fabs((double) alpha_pi.size() - N) > 0.1){
        Rprintf("Error: prior for gaussian mixture is not correct.\n");
        stop(0);
    }
    sum.resize(N);

    /* Sample pi */

    count = alpha_pi;
    for(i = 0; i < M; i++)
        count[z[i]] += 1.0;
    pi = randomdirichlet (count);


    /* Sample mu, sigma2 */

    for(i = 0; i < N; i++){
        sum[i] = 0.0;
        num[i] = 0;
    }
    for(i = 0; i < M; i++){
        sum[z[i]] += O[i];
        num[z[i]] ++;
    }
    for(j = 0; j < N; j++){
        mu[j] = (sum[j]/ sigma2[j] + muprior[j] / mupriorvar[j]) / (num[j] / sigma2[j] + 1.0 / mupriorvar[j]) + rnorminline() / sqrt(1.0 * num[j]/ sigma2[j] + 1.0 / mupriorvar[j]);
        i = 0;
        while(mu[j] < MINob || mu[j] > MAXob){
            if (i < 20){
                mu[j] = (sum[j]/ sigma2[j] + muprior[j] / mupriorvar[j]) / (num[j] / sigma2[j] + 1.0 / mupriorvar[j]) + rnorminline() / sqrt(1.0 * num[j]/ sigma2[j] + 1.0 / mupriorvar[j]);
                i ++;
            }
            else{
                mu[j] = muprior[j];
            }
        }
    }
    std::sort(mu.begin(), mu.end());
    for(i = 0; i < N; i++)
        sum[i] = 0.0;
    for(i = 0; i < M; i++)
        sum[z[i]] += pow(O[i] - mu[z[i]], 2.0);
    for(j = 0; j < N; j++){
        sigma2[j] = (sum[j] + nu[j] * s2[j]) / randomchisq(nu[j] + num[j]);
    }

    if(verbose){
        for(i = 0; i < N ; i ++)
            Rprintf("%f\t", mu[i]);
        Rprintf("\n");
        for(i = 0; i < N; i ++)
            Rprintf("%f\t", sigma2[i]);
        Rprintf("\n");
        for(j = 0; j < N; j ++){
            Rprintf("%f\t", pi[j]);
        }
        Rprintf("\n");
    }

    return;

}

/* Marginal likelihood of gaussian mixture model */

double HMM::calculate_marginal_posterior_gm(){
    double LL = 0.0, tempsum = 0.0;
    int i, j;
    for(j = 0; j < N; j++){
        LL += ldnorminline(mu[j], muprior[j], sqrt(mupriorvar[j]));
        LL += ldscaleinvchisq(sigma2[j], nu[j], s2[j]);
        LL += (alpha_pi[j] - 1.0) * log(pi[j] + 0.000000000001) - std::lgamma(alpha_pi[j]);
        tempsum += alpha_pi[j];
    }
    LL += std::lgamma(tempsum);
    for(i = 0; i < M; i ++){
        LL += ldnormmixinline(O[i], pi, mu, sigma2);
    }
    return(LL);
}

/* Gibbs sampling of Gaussian mixture model */

std::vector<std::vector<double> > HMM::gibbs_sampling_gm(int num_gibbs, int burnin, int thin, bool verbose){

    int k;
    std::vector<std::vector<double> > temp_output;

    for(int iter = 0; iter < burnin + num_gibbs * thin; iter ++){
        if(verbose){
            if(iter >= burnin && (iter - burnin) % thin == 0){
                samplehiddenstates_gm(1);
                samplenewparameters_gm(1);
            }
            else{
                samplehiddenstates_gm(0);
                samplenewparameters_gm(0);
            }
        }
        else{
            samplehiddenstates_gm(0);
            samplenewparameters_gm(0);
        }
        if(iter >= burnin && (iter - burnin) % thin == 0){
            temp_output.push_back(std::vector<double>());
            for(k = 0; k < N; k ++)
                temp_output[(int)temp_output.size() - 1].push_back(mu[k]);
            for(k = 0; k < N; k ++)
                temp_output[(int)temp_output.size() - 1].push_back(sigma2[k]);
            for(k = 1; k < N; k ++)
                temp_output[(int)temp_output.size() - 1].push_back(pi[k]);
            temp_output[(int)temp_output.size() - 1].push_back(calculate_marginal_posterior_gm());
        }
    }
    return(temp_output);
}

/* Model fitting: EM, Gibbs (HMM), Gibbs (GM) */

std::vector<std::vector<double> > HMM::fit_hmm_uneqvar(int Method, std::vector<double> obs, int K, bool bool_verbose, int N_starting_points, int num_gibbs, int burnin, int thin){

    int i, j;

    std::vector<std::vector<double> > hmm_gibbs_result, gm_gibbs_result, HMMembefore, resultstemp;

    if (Method == 1){

        /* EM fitting */

        BaumWelch_multi_starting_point_fitting(N_starting_points, sigma2_init, pi_init, A_init, bool_verbose);

        HMMembefore.resize(K + 4);
        for(i = 0; i < K; i ++){
            HMMembefore[0].push_back(pi[i]);
            HMMembefore[1].push_back(mu[i]);
            HMMembefore[2].push_back(sigma2[i]);
            for(j = 0; j < K; j ++){
                HMMembefore[i + 3].push_back(A[i][j]);
            }
            HMMembefore[K + 3].push_back(logprobf);
        }
        resultstemp = HMMembefore;
    }

    if (Method == 2){

        /* Gibbs sampling */

        hmm_gibbs_result = gibbs_sampling(num_gibbs, burnin, thin, bool_verbose);

        resultstemp = hmm_gibbs_result;
    }

    if (Method == 3) {

        gm_gibbs_result = gibbs_sampling_gm(num_gibbs, burnin, thin, bool_verbose);

        resultstemp = gm_gibbs_result;
    }

    return(resultstemp);

}




































































































































































































































































