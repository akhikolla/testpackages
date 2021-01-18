// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

Rcpp::List findpcs(arma::mat x, arma::mat B, arma::rowvec mu, int h, double s,
                   int N1, int N2, int Npc, double tol, bool fixed)
  {
  int iter = 0, it=0;
  int p = B.n_rows;
//  int h = hsub.n_rows;
  int n = x.n_rows;
  int q = B.n_cols;
  double sc0,sc,s0,delta1,delta2;
  arma::uvec x_index;
  arma::mat xcent;
  arma::mat hsub, hsubcent;
  arma::mat A(n,q); //
  arma::mat aprov(q,h); // Note that this is the transposed matrix for convenience
  arma::mat hsubnorm;
  arma::mat xnorm;

//  Rcpp::Rcout << "B is " << B << std::endl;
//  Rcpp::Rcout << "The center is " << mu << std::endl;
  xcent=x.each_row() - mu;

  if ((N1>0) || fixed)
    {
    A = xcent*B;
//    Rcpp::Rcout << "The scores are " << A << std::endl;
    }
  else
    {
    for (int kloop = 0; kloop < n; kloop++)
      {
      //          arma::vec aprov = solve( trans(B) * B, trans(B) * hsubcent.col(kloop));
      A.row(kloop) = trans(solve( trans(B) * B, trans(xcent.row(kloop)*B)));
      }
    }
  xnorm = xcent -  A*trans(B);
  colvec sqresid_dist = sum((xnorm % xnorm),1);
//  Rcpp::Rcout << "The distances are " << trans(sqresid_dist) << std::endl;
  delta1 = 1;
  s0 = s;

  while( ++it <= (N1+N2) && delta1>tol)
    {
//  Make h subset
    //    ws<-rep(1,n)
    x_index = sort_index(sqresid_dist);
    hsub = x.rows(x_index.head(h));
    mu = mean(hsub);
//    Rcpp::Rcout << "The subset is " << trans(sort(x_index.head(h)+1)) << std::endl;
//    Rcpp::Rcout << "The center is " << mu << std::endl;
//    Rcpp::Rcout << "The condition is " <<(!fixed & (it>N1)) << std::endl;
    if(!fixed & (it>N1))
      {
      iter = 0;
//      sc0 = 1e10;
//      sc0 = std::numeric_limits<double>::infinity();
      if (s == 0) {sc0 = sum(sqresid_dist(x_index.head(h))) / h;}
      else {sc0 = s;}
      delta2 =1; // delta = 1-sc/sc0;
                // sc0 = sc;

      while(++iter < Npc && delta2>tol )
        {
        hsubcent = hsub.each_row() - mu;
//      Rcpp::Rcout << "The matrix is " << hsubcent << std::endl;
        try
          {
          for (int kloop =0; kloop<h; kloop++)
            {
//          arma::vec aprov = solve( trans(B) * B, trans(B) * hsubcent.col(kloop));
            aprov.col(kloop) = solve( trans(B) * B, trans(hsubcent.row(kloop)*B));
            }
//          Rcpp::Rcout << "The matrix A is " << aprov << std::endl;
          for(int jloop =0; jloop<p; jloop++)
            {
            B.row(jloop) = trans(solve(aprov * trans(aprov), aprov*hsubcent.col(jloop)));
//          Rcpp::Rcout << "The product is " << mean(trans(B.row(jloop)*aprov)) << std::endl;
            mu(jloop) = mean( hsub.col(jloop) - trans(B.row(jloop)*aprov));
            }
//          Rcpp::Rcout << "The matrix is " << B  << std::endl;
//          Rcpp::Rcout << "The center in the loop is " << mu << std::endl;
          hsubcent = hsub.each_row() - mu;
          hsubnorm = hsubcent -  trans(B*aprov);
          colvec sqresid_dist0 = sum((hsubnorm % hsubnorm),1);
          sc = sum(sqresid_dist0) / h; // why not mean??
          delta2 = 1-sc/sc0;
          sc0 = sc;
//          Rcpp::Rcout << "The scale in the loop is " << sc << std::endl;
          }
        catch(const std::runtime_error&e)
        {
          arma::vec eigval;
          arma::mat eigvec;
          mu = mean(hsub);
          eig_sym(eigval, eigvec, cov(hsub));
          B = eigvec.cols(1,B.n_cols);
          iter=Npc;
          throw e;
          }
        }
      xcent=x.each_row() - mu;
      for (int kloop = 0; kloop < n; kloop++)
        {
        //          arma::vec aprov = solve( trans(B) * B, trans(B) * hsubcent.col(kloop));
        A.row(kloop) = trans(solve( trans(B) * B, trans(xcent.row(kloop)*B)));
        }
      }
    else
      {
      xcent=x.each_row() - mu;
      A = xcent*B;
//      Rcpp::Rcout << "xcent is " << xcent.row(0)  << std::endl;
//      Rcpp::Rcout << "xcent is " << xcent  << std::endl;
//      Rcpp::Rcout << "B is " << B  << std::endl;
//      Rcpp::Rcout << "A is " << A.row(1)  << std::endl;
      A = (x.each_row() - mu) * B;
//      Rcpp::Rcout << "The first row is " << A  << std::endl;
      }
    xnorm = xcent -  A * trans(B);
//    Rcpp::Rcout << "xnorm is " << xnorm.row(0)  << std::endl;
    sqresid_dist = sum((xnorm % xnorm),1);
//    Rcpp::Rcout << "dist is " << sqresid_dist  << std::endl;
    if (N1==0 && s!=0)
      {
//      Rcpp::Rcout << "B is " << B << std::endl;
      x_index = sort_index(sqresid_dist);
      s = sum(sqresid_dist(x_index.head(h))) / h;
      delta1 = 1-s/s0;
      s0 = s;
//      Rcpp::Rcout << "The scale is " << s << std::endl;
      }
    }
  if(s==0)
    {
    x_index = sort_index(sqresid_dist);
    s = sum(sqresid_dist(x_index.head(h))) / h;
    }
  x_index = x_index.head(h)+1;
//  return(B);
return Rcpp::List::create(Rcpp::Named("Bmat") = B, Rcpp::Named("mu") = mu, Rcpp::Named("scale") = s,
                          Rcpp::Named("index") = trans(x_index));
// return Rcpp::List::create(Rcpp::Named("Bmat") = mu);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//


