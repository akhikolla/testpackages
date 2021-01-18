// losses.h

//#include "RcppArmadillo.h"
#include<string> // for string class 

#define LOSS_LINEAR 0
#define LOSS_LOGIT 1



double R_linear(const arma::vec & beta, 
                const arma::mat & X,
                const arma::mat & y){
    arma::vec r = y - X*beta;
    return arma::sum(arma::square(r))/(2*y.n_elem);
}
double R_linear2(const arma::vec & beta, 
                const arma::mat & X,
                double y){
    arma::vec r = y - X*beta;
    return arma::sum(arma::square(r))/(2*X.n_rows);
}

arma::vec grad_R_linear(const arma::vec & beta, 
                        const arma::mat & X,
                        const arma::mat & y){
    return -X.t()*(y - X*beta)/y.n_elem;
}

arma::vec grad_R_linear2(const arma::vec & beta, 
                        const arma::mat & X,
                        double y){
    return -X.t()*(y - X*beta)/X.n_elem;
}

// Logit
double log1pexp(double eta){
    if(eta <= 18 && eta >= -37){
        return std::log(1 + std::exp(eta));
    }
    if(eta > 18 && eta <= 33.3){
        return eta + std::exp(-eta);
    }
    if(eta > 33.3){
        return eta;
    }
    return exp(eta);
}

double R_logit(const arma::vec & beta, 
               const arma::mat & X,
               const arma::mat & y){
    arma::vec eta = X*beta;
    // return arma::mean(arma::log(1 + arma::exp(eta)) - y % eta);

    // There is a more efficient and stable computation
    double rlogit = 0;
    for (unsigned int i = 0; i < eta.n_elem; i++)
    {
        rlogit += eta[i] = log1pexp(eta[i]) - y[i] * eta[i];
    }
    return rlogit/eta.n_elem;
}
double R_logit2(const arma::vec & beta, 
               const arma::mat & X,
               double y){
    arma::vec eta = X*beta;
    // There is a more efficient and stable computation
    double rlogit = 0;
    for (unsigned int i = 0; i < eta.n_elem; i++)
    {
        rlogit += eta[i] = log1pexp(eta[i]) - y * eta[i];
    }
    return rlogit/eta.n_elem;
}

double inv1pexpm(double eta){
    if(eta >= -30){
        return 1/(1 + std::exp(-eta));
    }else{
        return std::exp(eta);
    }
}
arma::vec grad_R_logit(const arma::vec & beta, 
                       const arma::mat & X,
                       const arma::mat & y){
    // arma::vec eta = X*beta;
    // arma::vec d = 1/(1 + arma::exp(-eta)) - y;
    // return X.t() * d/y.n_elem;

    // More stable
    arma::vec tmp = X*beta;
    for (unsigned int i = 0; i < tmp.n_elem; i++)
    {
        tmp[i] = inv1pexpm(tmp[i]) - y[i];
    }
    return X.t()*tmp/y.n_elem;
}
arma::vec grad_R_logit2(const arma::vec & beta, 
                       const arma::mat & X,
                       double y){
    // More stable
    arma::vec tmp = X*beta;
    for (unsigned int i = 0; i < tmp.n_elem; i++)
    {
        tmp[i] = inv1pexpm(tmp[i]) - y;
    }
    return X.t()*tmp/tmp.n_elem;
}

//===============================
// Class: Loss

class Loss
{
public:
    double (*R)(const arma::vec & beta, const arma::mat & X, const arma::mat & y);
    arma::vec (*gradR)(const arma::vec & beta, const arma::mat & X, const arma::mat & y);
    double (*R2)(const arma::vec & beta, const arma::mat & X, double y);
    arma::vec (*gradR2)(const arma::vec & beta, const arma::mat & X, double y);
    
    std::string name;
    
    Loss(){
        name = "linear";
        R = R_linear;
        gradR = grad_R_linear;
        R2 = R_linear2;
        gradR2 = grad_R_linear2;
    }

    Loss(int type){
        switch (type)
        {
        case LOSS_LOGIT:
            this->name = "logit";
            this->R = R_logit;
            this->gradR = grad_R_logit;
            this->R2 = R_logit2;
            this->gradR2 = grad_R_logit2;
            break;    
        default:
            this->name = "linear";
            this->R = R_linear;
            this->gradR = grad_R_linear;
            this->R2 = R_linear2;
            this->gradR2 = grad_R_linear2;
            break;
            } 
    }
};

// This class is not exposed...
