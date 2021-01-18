#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"SurvCD.h"
#include"QR.h"
#include"Robust.h"
#include"NonRobust.h"
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;
//using namespace R;

// [[Rcpp::export()]]
arma::vec RunNetSurv(arma::mat& xc, arma::mat& xg, arma::vec& y, double lamb1, double lamb2, arma::vec bc0, arma::vec bg0, 
            double r, arma::mat& a, int p, int pc, bool robust)
{
  int n = xc.n_rows, count = 0;
  arma::vec bnew(p, fill::zeros), yc, yg, bc = bc0, bg = bg0;
  double diff;
  
  while(count < 20){
    yc = y - xg * bg;
    //std::cout << "bc: "<< bc << std::endl;
    if(robust){
      bc = QRWMR(xc, yc, bc);
      yg = y - xc * bc;
      bnew = LadNet(xg, yg, lamb1, lamb2, bg, r, a, n, p);
    }else{
      bc = fastLm(yc, xc);
      yg = y - xc * bc;
      bnew = LSNet(xg, yg, lamb1, lamb2, bg, r, a, n, p);
    }
    
    diff = arma::accu(arma::abs(bg - bnew));
    if(diff < 0.001) break;
    else{
      bg = bnew;
      count++;
    }
  }
  //std::cout << "diff: " << diff <<std::endl;
  arma::vec b1(pc+p, fill::none);
  b1.subvec(0, pc-1) = bc;
  b1.subvec(pc, pc+p-1) = bnew;
  return(b1);
}


// [[Rcpp::export()]]
arma::vec RunMCPSurv(arma::mat xc, arma::mat xg, arma::vec y, double lamb1, arma::vec bc, arma::vec bg, double r, int p, int pc, bool robust)
{
  int n = xc.n_rows, count = 0;
  arma::vec bnew(p, fill::zeros), yc, yg;
  
  while(count < 20){
    yc = y - xg * bg;
    //std::cout << "mark_1 "<<  std::endl;
    if(robust){
      bc = QRWMR(xc, yc, bc);
      yg = y - xc * bc;
      bnew = LadMCP(xg, yg, lamb1, bg, r, n, p);
    }else{
      bc = fastLm(yc, xc);
      yg = y - xc * bc;
      bnew = LSMCP(xg, yg, lamb1, bg, r, n, p);
    }
    double diff = arma::accu(arma::abs(bg - bnew));
    //std::cout << "diff: " << diff <<std::endl;
    if(diff < 0.001) break;
    else{
      bg = bnew;
      count++;
    }
  }
  arma::vec b1(pc+p, fill::none);
  b1.subvec(0, pc-1) = bc;
  b1.subvec(pc, pc+p-1) = bnew;
  return(b1);
}


// [[Rcpp::export()]]
arma::vec RunLassoSurv(arma::mat xc, arma::mat xg, arma::vec y, double lamb1, arma::vec bc, arma::vec bg, int p, int pc, bool robust)
{
  int n = xc.n_rows, count = 0;
  arma::vec bnew(p, fill::zeros), yc, yg;
  
  while(count < 20){
    yc = y - xg * bg;
    if(robust){
      bc = QRWMR(xc, yc, bc);
      yg = y - xc * bc;
      bnew = LadLasso(xg, yg, lamb1, bg, n, p);
    }else{
      bc = fastLm(yc, xc);
      yg = y - xc * bc;
      bnew = LSLasso(xg, yg, lamb1, bg, n, p);
    }
    double diff = arma::accu(arma::abs(bg - bnew));
    //std::cout << "diff: " << diff <<std::endl;
    if(diff < 0.001) break;
    else{
      bg = bnew;
      count++;
    }
  }
  arma::vec b1(pc+p, fill::none);
  b1.subvec(0, pc-1) = bc;
  b1.subvec(pc, pc+p-1) = bnew;
  return(b1);
}
