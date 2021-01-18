#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppTN.h>
// [[Rcpp::depends(RcppTN)]]
using namespace Rcpp;

const double log2pi = log(2.0 * M_PI);


double rhoObjFunc(double rho, arma::mat ystar, arma::mat X1, arma::mat X2, arma::mat betaDraw, double rho0, double v0){
  //' rho Objective Function.
  //'
  //' rhoObjFunc() maximized for rho draw. 
  //'
  //' @param double rho
  //' @param arma::mat ystar
  //' @param arma::mat X1
  //' @param arma::mat X2
  //' @param arma::betaDraw
  //' @param double rho0
  //' @param double v0
  //' @return If the previous rho value is greater than 1, or the new 
  //'         rho value is not finite, return 0. Else, return the new
  //'         rho value.
  
   // Bug: This condition happens often if used for numDeriv::Hessian
   // This is due to the numerical algorithm evaluating outside the 
   // limits of the domain [-1,1]
  if (std::fabs(rho) > 1)
    return 99999999999;
  
  int k1 = X1.n_cols;
  int k2 = X2.n_cols;
  int k = k1 + k2;
  int N = ystar.n_cols;
  double detSigma = 1-pow(rho,2);
  arma::mat SigmaInv = {{1, -rho}, {-rho, 1}};
  SigmaInv = SigmaInv / detSigma;
  
  double tempsum = 0;
  arma::rowvec k1_vector = arma::zeros<arma::rowvec>(k1);
  arma::rowvec k2_vector = arma::zeros<arma::rowvec>(k2);
  for (int i = 0; i < N; i++)  {
    arma::mat Xi(2, k);
    Xi.row(0) = join_rows(X1.row(i), k2_vector);
    Xi.row(1) = join_rows(k1_vector, X2.row(i));
    arma::mat ystarmXB = ystar.row(i).t() - Xi * betaDraw.t();
    arma::mat result = ystarmXB.t() * SigmaInv * ystarmXB;
    tempsum = tempsum + result(0,0);
  }
  
  double out = -N * log(detSigma) / 2 - pow(rho-rho0, 2) / (2*v0) - .5 * tempsum - log( 2 * M_PI );

  return out;
}



double optim_rcpp(arma::mat ystar, arma::mat X1, arma::mat X2, arma::mat betaDraw, double rho0, double v0){
  
  //' Optimization.
  //'
  //' optim_rcpp() calls rhoObjFunc() and maximizes rho. 
  //'
  //' @param arma::mat ystar
  //' @param arma::mat X1
  //' @param arma::mat X2
  //' @param arma::betaDraw
  //' @param double rho0
  //' @param double v0
  //' @return double rho
  
  // Extract R's optim function
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optimize = stats["optimize"];
 
  // Call the optim function from R in C++ 
  Rcpp::List opt_results = optimize(Rcpp::_["f"] = Rcpp::InternalFunction(rhoObjFunc),
                                    Rcpp::_["lower"] = -1,
                                    Rcpp::_["upper"] =  1,
                                    Rcpp::_["ystar"] = ystar,
                                    Rcpp::_["X1"] = X1,
                                    Rcpp::_["X2"] = X2,
                                    Rcpp::_["betaDraw"] = betaDraw,
                                    Rcpp::_["rho0"] = rho0,
                                    Rcpp::_["v0"] = v0,
                                    Rcpp::_["maximum"] = true);
  
  
  
  // Extract out the estimated parameter values
  arma::vec out = Rcpp::as<arma::vec>(opt_results[0]);
  
  // Return estimated values
  return out[0];
}


double hessian_rcpp(double rhoHat, arma::mat ystar, arma::mat X1, arma::mat X2, arma::mat betaDraw, double rho0, double v0){
  
  
  //' Hessian.
  //'
  //' hessian_rcpp() calls rhoObjFunc() and approximates the Hessian matrix. 
  //'
  //' @param double rhoHat
  //' @param arma::mat ystar
  //' @param arma::mat X1
  //' @param arma::mat X2
  //' @param arma::betaDraw
  //' @param double rho0
  //' @param double v0
  //' @return double hessian
  //' 
  //Extract R's optim function
  Rcpp::Environment numDeriv("package:numDeriv");
  Rcpp::Function hessian = numDeriv["hessian"];

  // Call the optim function from R in C++
  Rcpp::List hessian_results = hessian(Rcpp::_["func"] = Rcpp::InternalFunction(rhoObjFunc),
                                 Rcpp::_["x"] = rhoHat,
                                 Rcpp::_["ystar"] = ystar,
                                 Rcpp::_["X1"] = X1,
                                 Rcpp::_["X2"] = X2,
                                 Rcpp::_["betaDraw"] = betaDraw,
                                 Rcpp::_["rho0"] = rho0,
                                 Rcpp::_["v0"] = v0);

  arma::vec out = Rcpp::as<arma::vec>(hessian_results[0]);

  // Return estimated values
  return out[0];
}


void set_seed(double seed) {
  //' Set Seed.
  //'
  //' set_seed() sets the seed for replication. 
  //'
  //' @param double seed
  //' @return double void
  //' 

  // set seed
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}


arma::mat mvrnorm_arma(int n, arma::vec mu, arma::mat sigma) {
  //' Draw from multivariate normal
  //'
  //' mvrnorm_arma(n,mu,sigma) takes n draws from a multivariate 
  //' normal with mean mu and covariance sigma
  //'
  //' @param int n 
  //' @param arma::vec mu 
  //' @param arma::mat sigma
  //' @return arma::mat RandomDraws
  //' 

  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


arma::vec dmvnorm_arma(arma::mat x,  
                      arma::mat mu,  
                      arma::mat sigma, 
                      bool logd = false) { 
  
  //' Density of a multivariate normal
  //'
  //' dmvnorm_arma(x, mu, sigma, logd) evaluates the density of a 
  //' multivariate normal with mean mu and covariance sigma at x  
  //'
  //' @param arma::mat x value to be evaluated at 
  //' @param arma::vec mu mean
  //' @param arma::mat sigma covariance
  //' @param bool logd evaluate on the log scale
  //' @return arma::vec density
  //' 
  
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv_sympd(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mu.row(i)) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

//' MCMC algorithm to sample from bivariate probit with partial
//' observability
//'
//' MCMC1() produces MCMC draws from the posterior of the bivariate
//' probit with partial observability. It does not perform input validation. 
//' It is reccomended to use \code{BiProbitPartial} instead of this function.
//' \code{BiProbitPartial} performs input validation
//' and then calls this function if \code{philosophy == "bayesian"}.
//'
//' @param X1 a matrix of covariates for the first equation
//' @param X2 a matrix of covariates for the second equation
//' @param Z a matrix of response values
//' @param beta1 a matrix of starting values for beta1
//' @param beta2 a matrix of starting values for beta2
//' @param rho a numeric starting value for rho
//' @param fixrho a logical determining if rho is fixed
//' @param S a numeric for the number of MCMC iterations
//' @param beta0 a matrix of the beta prior mean parameter
//' @param B0inv a matrix of the inverse of beta prior covariance parameter
//' @param rho0 a numeric for the mu prior parameter for rho
//' @param v0 a numeric for the Sigma prior parameter for rho
//' @param nu a numeric for MCMC tuning parameter 1
//' @param P a numeric for MCMC tuning parameter 2
//' @param tauSq a numeric for MCMC tuning parameter 3
//' @param seed a numeric seed for determining the random draw sequence
//' @return a matrix of MCMC draws
// [[Rcpp::export]]
arma::mat MCMC1(const arma::mat X1, const arma::mat X2, const arma::mat Z, 
                NumericVector beta1,
                NumericVector beta2,
                double rho, bool fixrho, int S, 
                NumericVector beta0 ,
                NumericMatrix B0inv ,
                double rho0, double v0,
                double nu, double P, double tauSq, 
                int seed)
{
  
  // Get parameter space size and sample size
  int k1 = X1.n_cols;
  int k2 = X2.n_cols;
  int k = k1 + k2;
  int N = Z.n_rows;
  int NZ1 = accu(Z);
  int NZ0 = N - NZ1;
  
  arma::vec beta1_arma = beta1;
  arma::vec beta2_arma = beta2;
  arma::vec beta0_arma = beta0;
  arma::mat B0inv_arma = as<arma::mat>(B0inv); 


  // Initialize MCMC storage
  arma::mat betaStorage = arma::zeros(S, k);
  arma::mat betaStorage1;
  betaStorage1.insert_cols(betaStorage1.n_cols, beta1_arma.t());
  betaStorage1.insert_cols(betaStorage1.n_cols, beta2_arma.t());
  betaStorage.row(0) = betaStorage1; 

  arma::vec rhoStorage(S);
  rhoStorage(0) = rho;
  double detSigma = 1 - pow(rho,2);
  
  // Precomputations for loop efficiency
  arma::mat Sigma = {{1, rho}, {rho, 1}};
  arma::mat SigmaInv = {{1, -rho}, {-rho, 1}};
  SigmaInv /= detSigma;
  
  arma::vec B0invbeta0 = B0inv_arma * beta0_arma;
  
  arma::uvec Z1Dex = find(Z == 1);
  arma::uvec Z0Dex = find(Z == 0);
  arma::mat X1_Z1 = X1.rows(Z1Dex);
  arma::mat X2_Z1 = X2.rows(Z1Dex);
  arma::mat X1_Z0 = X1.rows(Z0Dex);
  arma::mat X2_Z0 = X2.rows(Z0Dex);
  arma::mat X12 = arma::join_rows(X1,X2);
  const arma::mat Mk1 = arma::ones(k1, k1);
  const arma::mat Mk2 = arma::ones(k2, k2);
  arma::mat Mrhok2k1 = arma::zeros(k2,k1);
  arma::mat Mrhok1k2 = arma::zeros(k1,k2);
  arma::rowvec k1_vector = arma::zeros<arma::rowvec>(k1);
  arma::rowvec k2_vector = arma::zeros<arma::rowvec>(k2);
  
  
  arma::mat tempsum1 = X12.t() * X12 / detSigma;
  tempsum1 = tempsum1 % arma::join_rows(arma::join_cols(Mk1,
                                                        Mrhok2k1.fill(-rho)),
                                                        arma::join_cols(Mrhok1k2.fill(-rho),
                                                                        Mk2));
  arma::mat B = B0inv_arma + tempsum1;
  arma::mat Binv = arma::inv_sympd(B);
  
  // Prep MCMC loops for first draw 
  arma::vec beta1Draw = beta1;
  arma::vec beta2Draw = beta2;
  
  arma::vec y1star = arma::ones(N);
  arma::vec y2star = arma::ones(N);
  y2star.rows(Z0Dex) = -arma::ones(NZ0);
  arma::mat ystar = join_cols(y1star.t(),y2star.t());
  arma::vec y1star_Z1 = y1star.rows(Z1Dex);
  arma::vec y2star_Z1 = y2star.rows(Z1Dex);


  arma::mat X1beta1 = X1 * beta1Draw;
  arma::mat X2beta2 = X2 * beta2Draw;
  arma::mat Xbeta;
  Xbeta.insert_cols(Xbeta.n_cols,X1beta1);
  Xbeta.insert_cols(Xbeta.n_cols,X2beta2);
  
  arma::mat X1beta1_Z1 = X1beta1.rows(Z1Dex);
  arma::mat X2beta2_Z1 = X2beta2.rows(Z1Dex);
  arma::mat X1beta1_Z0 = X1beta1.rows(Z0Dex);
  arma::mat X2beta2_Z0 = X2beta2.rows(Z0Dex);
  arma::mat Xbeta_Z0 = Xbeta.rows(Z0Dex);
  
  for (int s = 1; s < S; s++)
  {
    
    if (seed != 0) {set_seed(seed+s);}
    
    
    for (int i = 0; i < NZ1; i++){
      y1star_Z1[i] = RcppTN::rtn1(X1beta1_Z1[i] + rho * (y2star_Z1[i] - X2beta2_Z1[i]), sqrt(1-pow(rho,2)), 0,HUGE_VAL);
      y2star_Z1[i] = RcppTN::rtn1(X2beta2_Z1[i] + rho * (y1star_Z1[i] - X1beta1_Z1[i]), sqrt(1-pow(rho,2)), 0,HUGE_VAL);
    }

    // Z == 0, Accept Reject algorthm
    arma::mat ytemp_Z0(2,NZ0);
    
    for (int i = 0; i < NZ0; i++){
      arma::vec xbeta_z0_val = Xbeta_Z0.row(i).t() ;
      
      ytemp_Z0.col(i) = mvrnorm_arma(1, xbeta_z0_val, Sigma).t();
    }
    
    arma::vec Quad1Dex_temp = min(ytemp_Z0.t(), 1);
    arma::uvec Quad1Dex = find(Quad1Dex_temp > 0);
    
    int j = 0;
    while ((Quad1Dex.size() > 0)){
      for (unsigned i = 0; i < Quad1Dex.size(); i++){
        arma::vec xbeta_z0_val = Xbeta_Z0.row(Quad1Dex(i)).t() ;
        ytemp_Z0.col(Quad1Dex(i)) = mvrnorm_arma(1, xbeta_z0_val, Sigma).t();
      }
      Quad1Dex_temp = min(ytemp_Z0.t(), 1);
      Quad1Dex = find(Quad1Dex_temp > 0);
      j+= 1;
    }
    
    arma::mat ytemp_Z1 = arma::join_cols(y1star_Z1.t(),y2star_Z1.t());
    
    
    ystar.cols(Z1Dex) = ytemp_Z1;
    ystar.cols(Z0Dex) = ytemp_Z0;
      
    arma::vec onemrho = arma::ones(2);
    onemrho(1) = -rho;
    arma::vec mrhoone = arma::ones(2);
    mrhoone(0) = -rho;
    arma::mat tempysum1 = arma::zeros(N, k1);
    tempysum1.each_col()= ystar.t() * onemrho;
    arma::mat tempysum2 = arma::zeros(N, k2);
    tempysum2.each_col()= ystar.t() * mrhoone;
    arma::mat temp1 = X1 % tempysum1;
    arma::mat temp2 = X2 % tempysum2;
    arma::mat tempsum2 = sum(arma::join_rows(temp1, temp2),0).t();
    tempsum2 /= detSigma;
      
    arma::vec betahat = Binv * (B0invbeta0 + tempsum2);
      
    
    // Draw beta, Gibbs
    arma::mat betaDraw = mvrnorm_arma(1, betahat, Binv);
    
    beta1Draw = betaDraw.cols(0, k1-1).row(0).t();
    beta2Draw = betaDraw.cols(k1, k1+k2-1).row(0).t();

    betaStorage.row(s) = betaDraw;

    
    // Draw rho
    // Precomputation
    X1beta1 = X1 * beta1Draw;
    X2beta2 = X2 * beta2Draw;
    Xbeta = join_rows(X1beta1, X2beta2);
    
    X1beta1_Z1 = X1beta1.rows(Z1Dex);
    X2beta2_Z1 = X2beta2.rows(Z1Dex);
    X1beta1_Z0 = X1beta1.rows(Z0Dex);
    X2beta2_Z0 = X2beta2.rows(Z0Dex);
    Xbeta_Z0 = Xbeta.rows(Z0Dex);
     
    if ((fixrho == true) && (s == 1)){
      arma::vec rhoStorage (N);
      rhoStorage.fill(rho);
      }
    else if (fixrho == false){

      double rhoHat = optim_rcpp(ystar.t(), X1, X2, betaDraw, rho0, v0);

      double SHat = -1/hessian_rcpp(rhoHat, ystar.t(), X1, X2, betaDraw, rho0, v0)/N;

      double rhoStar;
      if (SHat < 0){ 
        rhoStar = rho;
      }
      else
        rhoStar = rt(1,nu)[0]*sqrt(SHat * tauSq) + (rhoHat + P*(rho - rhoHat));
      
      if(std::fabs(rhoStar)<1){

        double PriorRatio = log(RcppTN::dtn1(rhoStar,rho0,sqrt(v0),-1,1)) - log(RcppTN::dtn1(rho,rho0,sqrt(v0),-1,1));

        double LikelihoodRatio = sum(dmvnorm_arma(ystar.t(), Xbeta, {{1, rhoStar}, {rhoStar, 1}}, true)) - sum(dmvnorm_arma(ystar.t(), Xbeta, {{1, rhoStar}, {rhoStar, 1}},true));

        double ProposalRatio = R::dt((rhoStar-rhoHat)/sqrt(SHat * tauSq),nu, true) - R::dt((rho-rhoHat)/sqrt(SHat * tauSq),nu, true);

        double alpha = PriorRatio + LikelihoodRatio + ProposalRatio;
        
        if(runif(1)(0) <= exp(alpha)){
 
          rho = rhoStar;
          
          Sigma = {{1, rho}, {rho, 1}};
          detSigma = 1 - pow(rho,2);
          SigmaInv = {{1, -rho}, {-rho, 1}};
          SigmaInv /= detSigma;
          
          tempsum1 = X12.t() * X12 / detSigma;
          tempsum1 = tempsum1 % arma::join_rows(arma::join_cols(Mk1,
                                                                Mrhok2k1.fill(-rho)),
                                                                arma::join_cols(Mrhok1k2.fill(-rho),
                                                                                Mk2));
          
          B = B0inv_arma + tempsum1;
          Binv = arma::inv_sympd(B);
        }
      }
          
      rhoStorage(s) = rho;
    }

  }


  arma::mat results = betaStorage;
  results.insert_cols(results.n_cols, rhoStorage);
  return results; 
}



