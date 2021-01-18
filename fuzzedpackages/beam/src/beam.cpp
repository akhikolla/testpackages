//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(BH)]]
#include <R.h>
#include <tuple>
#include <boost/math/tools/minima.hpp>
#include <boost/math/special_functions/beta.hpp>
#undef NDEBUG
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;
using boost::math::tools::brent_find_minima;
using boost::math::ibetac;


//////////////////////////////////////////////
//----------- Internal functions -----------//
//////////////////////////////////////////////

arma::mat cov2cor(arma::mat c){
  colvec sds = 1/sqrt(c.diag());
  c.each_col() %= sds;
  c.each_row() %= sds.t();
  c.diag() = c.diag() / c.diag().max();
  return c;
}

arma::colvec get_p_cor(arma::mat T){
  T = cov2cor(T);
  T += 2;
  T = trimatu(T, 1);
  colvec Tvec = nonzeros(T);
  Tvec -= 2;
  return(-Tvec);
}

void standardize(arma::mat* x){
  arma::rowvec cm = mean(*x);
  x->each_row() -= cm;
  rowvec cs = 1/sqrt(sum(square(*x)) / x->n_rows);
  x->each_row() %= cs;
}

void center(arma::mat* x){
  arma::rowvec cm = mean(*x);
  x->each_row() -= cm;
}

void scale(arma::mat* x){
  rowvec cs = 1/sqrt(sum(square(*x)) / x->n_rows);
  x->each_row() %= cs;
}

double alphaToDelta(double alpha, int n, int p){
  return (alpha*n+(1-alpha)*p+(1-alpha))/(1-alpha);
}

double deltaToAlpha(double delta, int n, int p){
  return (delta-p-1)/(n+delta-p-1);
}

double lpvarGamma(const double x, const int p) {
  double ans = (p * (p - 1) * 0.25) * log(datum::pi);
  for(int j = 1; j < (p + 1); ++j){
    ans += std::lgamma(x - ((j - 1.0) * 0.5));
  }
  return ans;
}

double logML(const double delta, const int p, const int n, colvec eigs, const double logdetD){
  double out = -0.5*n*p*std::log(datum::pi);
  out += lpvarGamma((delta+n)*0.5, p);
  out -= lpvarGamma(delta*0.5, p);
  out += 0.5*delta*p*std::log(delta-p-1);
  out -= 0.5*(delta+n)*sum(arma::log((delta-p-1)+eigs));
  if(logdetD!=0){
    out -= 0.5*n*logdetD;
  }
  return(out);
}

double getDeltaOpt(const int n, const int p, colvec eigs, const double logdetD){
  const double lowerVal = alphaToDelta(0.001, n, p);
  const double upperVal = alphaToDelta(0.999, n, p);
  const auto obj = [p, n, eigs, logdetD](double x) { return -logML(x, p, n, eigs, logdetD); };
  boost::uintmax_t it = 1000;
  const auto result = brent_find_minima(obj, lowerVal, upperVal, 1000, it);
  //std::pair<double, double> result = brent_find_minima(obj, lowerVal, upperVal, 1000, it);
  auto deltaOpt = 0.0, valOpt = 0.0;
  std::tie(deltaOpt, valOpt) = result;
  return(deltaOpt);
}

arma::mat get_rzij2(arma::mat* TT, arma::mat* FF){
  arma::mat tiitjj = TT->diag() * TT->diag().t();
  arma::mat fiifjj = FF->diag() * FF->diag().t();
  arma::mat den1 = tiitjj - square(*TT);
  arma::mat den2 = fiifjj - square(*FF);
  arma::mat zij = - *TT/den1 + *FF/den2;
  arma::mat tiifjj = TT->diag() * FF->diag().t();
  tiifjj += tiifjj.t();
  arma::mat ziizjj = tiitjj/square(den1) + fiifjj/square(den2) - tiifjj/(den1%den2);
  tiitjj.clear();
  fiifjj.clear();
  tiifjj.clear();
  den1.clear();
  den2.clear();
  arma::mat rzij2 = zij/sqrt(ziizjj);
  rzij2 = square(rzij2);
  rzij2 = trimatu(rzij2, 1);
  return(nonzeros(rzij2));
}

arma::colvec get2_rzij2(arma::mat* TT, const double d){
  
  arma::mat den1 = - *TT;
  arma::mat rzij2 = den1;
  
  den1 %= *TT;
  arma::mat ziizjj = TT->diag() * TT->diag().t();
  den1 += ziizjj;
  den1 = 1/den1;
  ziizjj %= den1;
  ziizjj %= den1;
  ziizjj += 1/(d*d);
  arma::mat tiifjj( TT->n_rows , TT->n_cols , fill::zeros);
  tiifjj.each_col() = TT->diag()*d;
  tiifjj += tiifjj.t();
  tiifjj %= den1;
  tiifjj = tiifjj * (1/(d*d));
  ziizjj -= tiifjj;
  tiifjj.clear();
  
  rzij2 %= den1;
  den1.clear();
  rzij2 /= sqrt(ziizjj);
  ziizjj.clear();
  rzij2 = square(rzij2);
  rzij2 += 2;
  rzij2 = trimatu(rzij2, 1);
  rzij2 = nonzeros(rzij2);
  rzij2 -= 2;
  
  return(rzij2);
}

arma::colvec getTails(arma::colvec q, const double s1, const double s2){
  const auto fun = [s1, s2](double x) { return R::pbeta(x, s1, s2, false, false); };
  q.transform(fun);
  return(q);
}

arma::sp_mat getTails2(arma::sp_mat q, const double s1, const double s2){
  const auto fun = [s1, s2](double x) { return ibetac(s1, s2, x); };
  q.transform(fun);
  return(q);
}

arma::colvec get_p_BF(arma::colvec rqij, arma::colvec rgij, const double delta, const int n){
  double k1 = std::lgamma((delta+n)/2);
  k1 -= std::lgamma(delta/2);
  k1 += std::lgamma((delta+n-1)/2);
  k1 -= std::lgamma((delta-1)/2);
  k1 += 2*std::lgamma((delta+1)/2);
  k1 -= 2*std::lgamma((delta+n+1)/2);
  colvec logBF = ones(size(rqij));
  logBF -= square(rqij);
  logBF = log(logBF);
  logBF = - logBF*0.5*(delta+n);
  logBF += k1;
  if(!rgij.is_empty()){
    logBF += 0.5*delta*log(ones(size(rqij))-square(rgij));
  }
  return(logBF);
}

arma::colvec get_m_BF(arma::colvec rtij, arma::colvec rfij, const double delta, const int n, const int p){
  double k2 = lpvarGamma((delta+n-p+2)/2, 2);
  k2 -= lpvarGamma((delta-p+2)/2, 2);
  k2 += 2*std::lgamma((delta-p+3)/2);
  k2 -= 2*std::lgamma((delta+n-p+3)/2);
  arma::colvec logBF = ones(size(rtij));
  logBF -= square(rtij);
  logBF = log(logBF);
  logBF = - logBF*0.5*(delta+n-p+2);
  logBF += k2;
  if(!rfij.is_empty()){
    logBF += 0.5*(delta-p+2)*log(ones(size(rtij))-square(rfij));
  }
  return(logBF);
}


//////////////////////////////////////////////
//----------- External functions -----------//
//////////////////////////////////////////////

// [[Rcpp::export(".beam")]]
Rcpp::List beam(arma::mat X, std::string type, arma::colvec ronly, arma::mat D, bool verbose = true){
  
  // Dimension data
  const int n = X.n_rows;
  const int p = X.n_cols;
  
  // Center data
  center(&X);
  
  // Sample variances
  rowvec s = sum(square(X), 0)/n;
  
  // Scale data
  scale(&X);
  
  // Scatter matrix
  arma::mat XTX;
  if(type == "both" || type == "marginal"){
    XTX = X.t()*X;
  }
  
  arma::mat Dinv, cholD;
  double logdetD = 0.0;
  bool isD = !all(D.diag()==0);
  if(isD){
    cholD = chol(D);
    if(all(cholD.diag()>0)){
      logdetD = 2*sum(log(cholD.diag()));
      cholD = inv(cholD);
      Dinv = cholD*cholD.t();
    }else{
      cholD.clear();
      isD = false;
      Rcpp::Rcout << "Warning: D is not positive definite and has been set to the identity matrix" << std::endl;
    }
  }
  
  // Eigenvalue decomposition
  arma::colvec eigs = zeros(p);
  double deltaOpt;
  double d;
  arma::mat Tinv;
  
  if(n >= p){
    
    // Compute eigenvalues
    if(XTX.is_empty()){
      XTX = X.t()*X;
    }
    
    if(isD){
      
      arma::mat DinvXTX = Dinv*XTX;
      cx_vec eigsval = eig_gen(DinvXTX);
      DinvXTX.clear();
      eigs += sort(real(eigsval), "descend");
      eigsval.clear();
      
    }else{
      
      eigs += sort(eig_sym(XTX), "descend");
      
    }
    
    // Optimal shrinkage
    deltaOpt = getDeltaOpt(n, p, eigs, logdetD);
    
    // Compute Tinv
    d = 1/(deltaOpt-p-1);
    if(isD){
      
      XTX += (1/d)*D;
      Tinv = inv_sympd(XTX);
      XTX -= (1/d)*D;
      
    }else{
      
      XTX.diag() += (1/d);
      Tinv = inv_sympd(XTX);
      XTX.diag() -= (1/d);
      
    }
    
    if(type == "conditional"){
      
      XTX.clear();
      
    }		
    
  }else{
    
    if(isD){
      
      // Compute eigenvalues
      arma::mat XDinvXT = X*Dinv*X.t();
      eigs.subvec(0, n-1) += sort(eig_sym(XDinvXT), "descend");
      
      // Optimal shrinkage
      deltaOpt = getDeltaOpt(n, p, eigs, logdetD);
      
      // Posterior partial correlations
      d = 1/(deltaOpt-p-1);
      Tinv = d*Dinv;
      Tinv -= d*d*Dinv*X.t()*inv_sympd(eye(n, n) + d*XDinvXT)*X*Dinv;
      XDinvXT.clear();
      
    }else{
      
      // Compute eigenvalues
      arma::mat XXT = X*X.t();
      eigs.subvec(0, n-1) += sort(eig_sym(XXT), "descend");
      
      // Optimal shrinkage
      deltaOpt = getDeltaOpt(n, p, eigs, logdetD);
      
      // Posterior partial correlations
      d = 1/(deltaOpt-p-1);
      Tinv = d*eye(p, p);
      Tinv -= d*d*X.t()*inv_sympd(eye(n, n) + d*XXT)*X;
      XXT.clear();
      
    }
    
  }
  
  X.clear();
  const double valOpt = logML(deltaOpt, p, n, eigs, logdetD);
  const double alphaOpt = deltaToAlpha(deltaOpt, n, p);
  
  // Values of log-ML for a grid of alphas
  arma::mat gridAlpha(100, 3, fill::zeros);
  gridAlpha.col(1) += linspace(0.01, 0.99, 100);
  for(int j=0; j<100; j++){
    gridAlpha(j, 0) = alphaToDelta(gridAlpha(j, 1), n, p);
    gridAlpha(j, 2) = logML(gridAlpha(j, 0), p, n, eigs, logdetD);
  }
  eigs.clear();
  
  // Results
  int number_columns = sum(ronly);
  if(type == "both"){
    number_columns = number_columns*2;
  }
  arma::mat results(0.5*p*(p-1), number_columns, fill::zeros);
  int current_free_column = 0;
  
  //----------------------------------------//
  //          MARGINAL DEPENDENCIES         //
  //----------------------------------------//
  
  if(type == "both" || type == "marginal"){
    
    if(verbose){
      Rcpp::Rcout << "Marginal dependencies:" << std::endl;
      Rcpp::Rcout << "--> compute marginal correlation estimates...";
    }
 
    // Prior and posterior marginal correlations
    arma::colvec rsij =  nonzeros(trimatu(XTX/n, 1));
    XTX.clear();
    arma::colvec rtij = (1-alphaOpt)*rsij;
    arma::colvec rfij;

    if(isD){
      rfij = nonzeros(trimatu(D+1, 1))-1;
      rtij += alphaOpt*rfij;
    }
    if(ronly(0)==1){
      results.col(current_free_column) += rtij;
      current_free_column++;
    }
    if(verbose){
      Rcpp::Rcout << "DONE" << std::endl;
    }
    
    // Scaled Bayes factors
    if(ronly(1)==1){
      
      if(verbose){
        Rcpp::Rcout << "--> compute Bayes factors...";
      }
      
      results.col(current_free_column) += get_m_BF(rtij, rfij, deltaOpt, n, p);
      current_free_column++;
      
      if(verbose){
        Rcpp::Rcout << "DONE" << std::endl;
      }
      
    }
    rtij.clear();

    // Tail/Error probabilities
    if(ronly(2)==1){
      
      if(verbose){
        Rcpp::Rcout << "--> compute tail probabilities...";
      }
      
      results.col(current_free_column) = getTails(square(rsij), 0.5, (n-1)*0.5);
      current_free_column++;
      
      if(verbose){
        Rcpp::Rcout << "DONE" << std::endl;
      }
    }
    rsij.clear();
  }
  
  //----------------------------------------//
  //        CONDITIONAL DEPENDENCIES        //
  //----------------------------------------//
  
  arma::colvec rgij, rqij;
  
  if(type == "both" || type == "conditional"){
    
    if(verbose){
      Rcpp::Rcout << "Conditional dependencies:" << std::endl;
      Rcpp::Rcout << "--> compute partial correlation estimates...";
    }
    
    // Prior and posterior partial correlations
    rqij = get_p_cor(Tinv);
    if(ronly(0)==1){
      results.col(current_free_column) += rqij;
      current_free_column++;
    }
    
    if(verbose){
      Rcpp::Rcout << "DONE" << std::endl;
    }
    
    // Compute scaled log-Bayes factors
    if(ronly(1)==1){
      
      if(verbose){
        Rcpp::Rcout << "--> compute Bayes factors...";
      }
      
      if(isD){
        rgij = get_p_cor(d*Dinv);
      }
      results.col(current_free_column) += get_p_BF(rqij, rgij, deltaOpt, n);
      current_free_column++;
      
      if(verbose){
        Rcpp::Rcout << "DONE" << std::endl;
      }
    }
    
    // Compute tail/error probabilities
    if(ronly(2)==1){
      
      if(verbose){
        Rcpp::Rcout << "--> compute tail probabilities...";
      }
      
      arma::colvec rzij2;
      if(isD){
        Dinv = Dinv*d;
        rzij2 = get_rzij2(&Tinv, &Dinv);
        Dinv = Dinv*(1/d);
      }else{
        rzij2 = get2_rzij2(&Tinv, d);
      }
      
      results.col(current_free_column) = getTails(rzij2, 0.5, (n-1)*0.5);
      current_free_column++;
      
      if(verbose){
        Rcpp::Rcout << "DONE" << std::endl;
      }
    }
  }
  
  // Output object
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("table") = results,
                                      Rcpp::Named("deltaOpt") = deltaOpt,
                                      Rcpp::Named("alphaOpt") = alphaOpt,
                                      Rcpp::Named("gridAlpha") = gridAlpha,
                                      Rcpp::Named("valOpt") = valOpt,
                                      Rcpp::Named("TinvStdev") = sqrt(Tinv.diag()),
                                      Rcpp::Named("s") = s);
  
  return out;
}

// [[Rcpp::export(".lightbeam")]]
arma::sp_mat lightbeam(arma::mat X, const double thres, bool verbose = true){
  
  // Dimension data
  const int n = X.n_rows;
  const int p = X.n_cols;
  
  // Standardize data
  standardize(&X);

  if(verbose){
    Rcpp::Rcout << "DONE" << std::endl;
    Rcpp::Rcout << "--> compute eigenvalues and optimal shrinkage... ";
  }
  
  arma::colvec eigs = zeros(p);
  double deltaOpt;
  double d;
  arma::mat Tinv;
  if(n >= p){
    
    // Compute eigenvalues
    arma::mat XTX = X.t()*X;
    eigs += sort(eig_sym(XTX), "descend");
    
    // Optimal shrinkage
    deltaOpt = getDeltaOpt(n, p, eigs, 0);

    // Posterior partial correlations
    d = 1/(deltaOpt-p-1);
    Tinv = inv_sympd((1/d)*eye(p, p) + XTX);
    XTX.clear();
    
  }else{
    
    // Compute eigenvalues
    arma::mat XXT = X*X.t();
    eigs.subvec(0, n-1) += sort(eig_sym(XXT), "descend");
    
    // Optimal shrinkage
    deltaOpt = getDeltaOpt(n, p, eigs, 0);

    // Posterior partial correlations
    d = 1/(deltaOpt-p-1);
    Tinv = d*eye(p, p);
    Tinv -= d*d*X.t()*inv_sympd(eye(n, n) + d*XXT)*X;
    XXT.clear();
    
  }
  X.clear();
  eigs.clear();
  
  if(verbose){
    Rcpp::Rcout << "DONE" << std::endl;
    Rcpp::Rcout << "--> compute rzij2... ";
  }
  
  // Compute rzij2
  arma::mat den1 = - Tinv;
  arma::mat rzij2 = den1;
  den1 %= Tinv;
  arma::mat ziizjj = Tinv.diag() * Tinv.diag().t();
  den1 += ziizjj;
  den1 = 1/den1;
  ziizjj %= den1;
  ziizjj %= den1;
  ziizjj = ziizjj + 1/(d*d);
  arma::mat tiifjj(size(Tinv), fill::zeros);
  tiifjj.each_col() = Tinv.diag() * d;
  tiifjj += tiifjj.t();
  tiifjj %= den1;
  tiifjj = tiifjj * (1/(d*d));
  ziizjj -= tiifjj;
  tiifjj.clear();
  Tinv.clear();
  rzij2 %= den1;
  den1.clear();
  rzij2 /= sqrt(ziizjj);
  ziizjj.clear();
  rzij2 = square(rzij2);
  
  // Sparsify rzij2
  rzij2.elem(find(rzij2 <= thres)).zeros();
  
  if(verbose){
    Rcpp::Rcout << "DONE" << std::endl;
    Rcpp::Rcout << "--> compute tail probabilities... ";
  }
  
  // Compute tail/error probabilities
  arma::sp_mat tails = sp_mat(trimatu(rzij2, 1));
  rzij2.clear();
  tails = getTails2(tails, 0.5, (n-1)*0.5);
  
  if(verbose){
    Rcpp::Rcout << "DONE" << std::endl;
  }
  
  return tails;
}

