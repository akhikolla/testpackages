#include "EMalgo.h"
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/distributions/gamma.hpp>

using namespace std;
using namespace arma;
using namespace Rcpp;

using namespace boost::math;

//[[Rcpp::export]]
List  EMmhmmCPP(const S4 yiR, List paramR, const double tolR, const double nbKeepR, const double itersmallR){
  const Data * mydata_p =  new Data(yiR);
  EMalgo algo(mydata_p, paramR, tolR, nbKeepR, itersmallR);
  algo.Run();
  S4 paramOut = paramR[0];
  S4 * paramOut_p = & paramOut;
  algo.Output(paramOut_p);
  return List::create(Named("param") = paramOut,
                      Named("loglike") = algo.m_loglikeoutput,
                      Named("probaClusters") = algo.m_tau,
                      Named("probaStates") = algo.m_eta);
}

//code dgamma



//[[Rcpp::export]]
List oneEMgammaCPP(const NumericVector& my_xs, const NumericVector& my_ws, const int g, const NumericVector& val){
  Col<double> my_x = as< Col<double> >(my_xs), my_w = as<vec>(my_ws);
  my_w = my_w / sum(my_w);
  Col<double> prod_xw = my_x % my_w;
  int n = my_x.size();
  double tmpmean = sum(prod_xw) / sum(my_w), vartotal = sum( pow(my_x - tmpmean, 2) % my_w) / sum(my_w);
  Col<double> tmp = randu(g), tmp2 = cumsum(my_w);
  Col<double> b = ones(g), a = val + 0.1, norm = zeros(g); 
  Col<double> startb=b, starta=a;
  Mat<double> logprobs = zeros(n, g), tik = zeros(n, g);
  for (int k=0; k<g; k++) logprobs.col(k) = a(k) * log(b(k)) + (a(k) - 1) * log(my_x) - b(k) * my_x - lgamma(a(k)) - log(g);
  Col<double> maxprobs = max(logprobs, 1);
  Col<double> logdensity = maxprobs + log(sum(exp(logprobs.each_col() - maxprobs),1));
  double loglike = sum(my_w % logdensity), prec= log(0), lower_mean = 0, upper_mean = 0, upper_2=0, cst=0, tmp_prec_alpha = 0;   
  while ((loglike - prec)>0.0001){
      // Estep
      tik = exp(logprobs.each_col() - logdensity);
      //# Mstep
      norm = trans(trans(my_w) * tik);
      for (int k=0; k<g; k++){
        lower_mean = sum(tik.col(k) % my_w);
        upper_mean = sum(tik.col(k) % my_w % my_x);
        upper_2 = sum( log(my_x) % tik.col(k) % my_w);
        cst = (log(upper_mean / lower_mean) - upper_2 / lower_mean);
        tmp_prec_alpha = 0;
        a(k) =  0.1;
        int cp = 0;
        while ( (abs(tmp_prec_alpha - a(k)) > 0.00001) && (cp<50)){
          tmp_prec_alpha = a(k);
          a(k) = a(k) - (log(a(k)) - digamma(a(k)) - cst) / (1/a(k) - polygamma(1,a(k)) );
          cp++;
        }
        b(k) =  a(k)  * lower_mean / upper_mean;
      }
      //loglike
      for (int k=0; k<g; k++) logprobs.col(k) = a(k) * log(b(k)) + (a(k) - 1) * log(my_x) - b(k) * my_x - lgamma(a(k)) - log(g);
      maxprobs = max(logprobs, 1);
      logdensity = maxprobs + log(sum(exp(logprobs.each_col() - maxprobs),1));
      prec = loglike;
      loglike = sum(my_w % logdensity);
  }
  //list(a = a, b = b, loglike = loglike, error = loglike<prec)*/
  return List::create(Named("a") = a, Named("b") = b, Named("loglike") = loglike, Named("error") = loglike < prec, Named("starta") = starta, Named("startb") = startb);
}
