#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rcpp.h>
// [[Rcpp::depends("RcppArmadillo")]]


template< typename ARMA_VECTOR_TYPE >
ARMA_VECTOR_TYPE vintersection( ARMA_VECTOR_TYPE first, ARMA_VECTOR_TYPE second )
{
  std::vector< typename ARMA_VECTOR_TYPE::value_type > output ;
  std::set_intersection( first.begin(), first.end(), second.begin(), second.end(),
                         std::back_inserter(output) ) ;
  std::reverse( output.begin(), output.end() ) ;

  ARMA_VECTOR_TYPE result = arma::conv_to< ARMA_VECTOR_TYPE >::from(output);
  return result ;
}


// [[Rcpp::export]]
const arma::vec imapThetaFast(const arma::vec& theta0) {
  int n = theta0.size();
  arma::vec temp(n);
  temp(0) = theta0(0);
  if (n > 1)
    temp.subvec(1, n-1) = arma::log(theta0.subvec(1, n-1) - theta0.subvec(0, n-2));
  return temp;
}

// [[Rcpp::export]]
const arma::vec fscale_cutsFast(const arma::vec& par) {
  arma::vec temp(par.size());
  temp(0) = par(0);
  if (par.size() > 1)
    temp.subvec(1,par.size()-1) = arma::cumsum(arma::exp(par.subvec(1, par.size()-1)))+par(0);
  return temp;
}

// [[Rcpp::export]]
const arma::mat tableFast(const arma::vec& x, const arma::vec& y, const arma::vec& w) {
  const arma::vec xUni = arma::sort(arma::unique(x));
  const arma::vec yUni = arma::sort(arma::unique(y));
  int xsize = xUni.size();
  int ysize = yUni.size();
  arma::mat tab(xsize, ysize);
  for (int i = 0 ; i < xsize; i+=1) {
    const arma::uvec x_ind = arma::find(x==xUni(i));
    for(int j = 0; j < ysize; j+=1) {
      //const arma::vec y_val = y.elem(x_ind);
      //const arma::vec temp = arma::conv_to<arma::vec>::from(arma::find(y_val==yUni(j)));
      //int n = temp.size();
      const arma::uvec y_ind = arma::find(y==yUni(j));
      const arma::uvec intersect = vintersection(x_ind, y_ind);
      tab(i, j) = arma::sum(w.elem(intersect));
    }
  }
  return tab;
}

// [[Rcpp::export]]
int discord(const arma::mat& xytab) {
  //const arma::mat xytab = tableFast(x,y,w);
  int i = 0;
  int j = 0;
  bool foundConcord = false;
  bool foundDiscord = false;
  int ncols = xytab.n_cols-1;
  int nrows = xytab.n_rows-1;
  while(j < ncols) {
    if(i<nrows && j<ncols) {

      const arma::uvec ind1 = arma::linspace<arma::uvec>(i+1, nrows);
      const arma::uvec ind2 = arma::linspace<arma::uvec>(j+1,ncols);

      bool temp1 = xytab(i,j)>0;

      double temp2 = arma::accu(xytab.elem(ind1, ind2));

      bool temp3 = temp2 > 0;
      if(temp1 && temp3) {
        foundConcord = true;
        break;
      }
    }
    if(i>0 && j >0) {
      const arma::uvec ind1 = arma::linspace<arma::uvec>(0, i-1);
      const arma::uvec ind2 = arma::linspace<arma::uvec>(0,j-1);
      bool temp1 = xytab(i,j)>0;
      double temp2 = arma::accu(xytab.elem(ind1, ind2));
      bool temp3 = temp2 >0;
      if(temp1 && temp3 ) {
        foundConcord = true;
        break;
      }
    }
    i += 1;
    if(i>nrows) {
      i = 0;
      j = j + 1;
    }
  }


  i = 0;
  j = 0;
  while(j < ncols) {
    if(i>0 && j<ncols) {

      const arma::uvec ind1 = arma::linspace<arma::uvec>(0, i-1);
      const arma::uvec ind2 = arma::linspace<arma::uvec>(j+1,ncols);
      bool temp1 = xytab(i,j)>0;
      double temp2 = arma::accu(xytab.elem(ind1, ind2));
      bool temp3 = temp2 >0;
      if(temp1 && temp3 ) {

        foundDiscord = true;
        break;
      }
    }
    if(i<nrows && j > 0) {
      const arma::uvec ind1 = arma::linspace<arma::uvec>(i+1, nrows);
      const arma::uvec ind2 = arma::linspace<arma::uvec>(0,j-1);
      bool temp1 = xytab(i,j)>0;
      double temp2 = arma::accu(xytab.elem(ind1, ind2));
      bool temp3 = temp2 > 0;
      if(temp1 && temp3 ) {
        foundDiscord = true;
        break;
      }
    }
    i += 1;
    if(i>nrows) {
      i = 0;
      j = j + 1;
    }
  }
  if(foundDiscord == false)
    return 1;
  if(foundConcord == false)
    return -1;
  return 0;
}
// [[Rcpp::export]]

double lnlFast(const arma::mat& xytab, const arma::mat& pm) {
  arma::mat lpm = arma::log(pm);
  //lpm.elem(arma::find(lpm, std::numeric_limits<double>::infinity())).fill(arma::datum::log_min);
  //Rcpp::Rcout << lpm;
  lpm.elem(arma::find_nonfinite(lpm)).fill(arma::datum::log_min);
  double sum = arma::accu(xytab % lpm);

  return sum;
}

// double bNormaL(double prob,  arma::vec& lower,  arma::vec& upper, const arma::vec& mean, const arma::mat& S) {
//   int nu = 0;
//   const arma::vec sd = arma::sqrt(S.diag());
//   double rho = S(0,1)/pow(S(0,0)*S(1,1), 0.5);
//   lower = (lower - mean)/sd;
//   upper = (upper-mean)/sd;
//   arma::vec infin(2);
//   infin.fill(2);
//   const arma::vec a = Rcpp::.Fortran("smvbvt", prob, nu, lower, upper, infin, rho,
//                 PACKAGE = "mnormt");
//   return a(0);
// }
