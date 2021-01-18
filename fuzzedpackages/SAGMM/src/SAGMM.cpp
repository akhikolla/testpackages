// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]   

//'@importFrom Rcpp sourceCpp
//'@useDynLib SAGMM

#include "RcppArmadillo.h"

Rcpp::NumericVector export_vec(arma::vec y)
{
    Rcpp::NumericVector tmp = Rcpp::wrap(y);
    tmp.attr("dim") = R_NilValue;
    return tmp;
}


// [[Rcpp::export]]
double mahalanobis_HD(arma::rowvec y, arma::rowvec mu, arma::mat sigma)
{
    arma::rowvec d = y-mu;
    arma::mat sigma1 =  sigma;
    sigma1.diag() =  1/sigma.diag();
    double delta = (as_scalar(d*sigma1*d.t()));
    return delta;
}

// [[Rcpp::export]]
double norm_HD(arma::rowvec y, arma::rowvec mu, arma::mat sigma)
{
    double f1 = std::sqrt(arma::det(arma::datum::pi*2.0*sigma));
    
    double f2 = std::exp(-0.5*mahalanobis_HD(y, mu, sigma));
    
    double f = f2/f1;
    
    return f;
}

// [[Rcpp::export]]
Rcpp::List main_loop(arma::mat& X, int Dimensions, int Number, int Groups, arma::mat MU_O, arma::vec LAMBDA_O,arma::vec PISTAR_O,  
                     arma::cube SIGMA, arma::vec GAMMA){
    
    double LogLike =0.0;
    
    arma::mat MU = arma::zeros(Groups,Dimensions);
    
    arma::vec LAMBDA = arma::zeros(Groups);
    arma::vec LAMBDA_P = arma::zeros(Groups);
    arma::vec PISTAR = arma::zeros(Groups);
    arma::vec LAMBDASTAR = arma::zeros(Groups);
    arma::vec LAMBDASQ = arma::zeros(Groups);
    arma::vec PIvec = arma::zeros(Groups);

    arma::vec Tau = arma::zeros(Groups);
    arma::vec Tau2 = arma::zeros(Groups);
    arma::vec Comps = arma::zeros(Groups);
    arma::mat TauMAT = arma::zeros(Number,Groups);
    
    double tmp, tmp2;

    arma::mat EYE(Dimensions,Dimensions,arma::fill::eye);

    
    LAMBDA = LAMBDA_O;
    MU = MU_O;
    
    for(int ii=0; ii< Number; ii++) {

        PIvec = arma::exp(PISTAR_O)/arma::sum(arma::exp(PISTAR_O));
         
        // ### Tau # this is the component/overall_distribution that goes at the start of each gradient component
        LAMBDASQ = 0.5*arma::pow(LAMBDA,2.0);
   
        for(int gg=0; gg< Groups; gg++) {
            Tau(gg) = PIvec(gg)*norm_HD(X.row(ii), MU.row(gg), EYE*(LAMBDASQ(gg)));
        }
        
        Tau = Tau / sum(Tau);

        LAMBDA_P = arma::pow(LAMBDA_O,2.0);
        LAMBDASTAR = arma::log(LAMBDA);
        
        for(int gg=0; gg< Groups; gg++) {
            MU.row(gg) = MU_O.row(gg) + GAMMA(gg)* Tau(gg)*2/ LAMBDA_P(gg)*(X.row(ii)-MU_O.row(gg));
            LAMBDASTAR(gg) = LAMBDASTAR(gg) - GAMMA(ii)*Tau(gg)*(Dimensions-(2.0/LAMBDA_P(gg))*arma::dot(MU_O.row(gg) - X.row(ii), MU_O.row(gg) - X.row(ii)));
        }
     
        LAMBDA = arma::exp(LAMBDASTAR);
        LAMBDASQ = 0.5*arma::pow(LAMBDA,2.0);
    
        for(int gg=0; gg< Groups; gg++) {
            SIGMA.slice(gg) = EYE*LAMBDASQ(gg);
        }
    
        // # ### Pi
        for(int gg=0; gg< (Groups-1); gg++) {
            PISTAR(gg) = PISTAR_O(gg) + GAMMA(ii)*Tau(gg)*(2 - PIvec(gg));
        }
    
        for(int gg=0; gg< Groups; gg++) {
            tmp2 = std::log(norm_HD(X.row(ii), MU_O.row(gg), SIGMA.slice(gg)));
            Comps(gg) = std::log(PISTAR(gg))/arma::sum(PISTAR) + tmp2;
            Tau2(gg) = log(PIvec(gg)) + tmp2;
        }
        
        TauMAT.row(ii) =  Tau2.t();
        
        LogLike =  LogLike + arma::max(Comps) + std::log(arma::sum(arma::exp(Comps-arma::max(Comps))));
        LAMBDA_O = LAMBDA;
        MU_O = MU;
        PISTAR_O = PISTAR;
 
    }

    PIvec = arma::exp(PISTAR)/arma::sum(arma::exp(PISTAR));

    
    Rcpp::List retList = Rcpp::List::create(
        Rcpp::Named("PI")= export_vec(PIvec),
        Rcpp::Named("MU")= (MU),
        Rcpp::Named("LAMBDA")= export_vec(LAMBDA),
        //Rcpp::Named("Comps")= export_vec(Comps),
        Rcpp::Named("LogLike")= LogLike,
        Rcpp::Named("Tau")= TauMAT
        
    );

    
    return(retList);
}
