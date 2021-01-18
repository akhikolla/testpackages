# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

////// Type Conversions //////

int sign (double x){
  return((x>0)-(x<0));
  }

// [[Rcpp::export()]]
List strat_cpp (arma::vec y, arma::vec r, arma::vec w){

  int n=y.n_rows;
  double nume=0;
  double deno=0;

  for (int i=0; i<n; i++){
    for (int j=i+1; j<n; j++){
      if(r(j)<=r(i) || y(j)==y(i)) continue;
      double wij=w(i)*w(j);
      deno=deno+wij;
      nume=nume+sign(y(j)-y(i))*wij;
    }
  }

  return List::create(
    _["deno"] = deno,
    _["strat"] = nume/deno
  );
}


// [[Rcpp::export()]]
List strat_cpp_by (arma::vec y, arma::vec r,
                    arma::vec w, arma::vec c){

  int n=y.n_rows;

  arma::vec uniq_c=unique(c);
  int m=uniq_c.n_rows;

  arma::vec nume_by=arma::zeros(m);
  arma::vec deno_by=arma::zeros(m);
  double nume_between=0;
  double deno_between=0;

  for (int i=0; i<n; i++){
    int temp_c=c(i);
    for (int j=i+1; j<n; j++){
      if(r(j)<=r(i) || y(j)==y(i)) continue;
      double wij=w(i)*w(j);
      if (c(j)==c(i)){
        deno_by(temp_c)=deno_by(temp_c)+wij;
        nume_by(temp_c)=nume_by(temp_c)+sign(y(j)-y(i))*wij;
      }
      else{
        deno_between=deno_between+wij;
        nume_between=nume_between+sign(y(j)-y(i))*wij;
      }
    }
  }

  double nume_within=sum(nume_by);
  double deno_within=sum(deno_by);

  return List::create(
    _["deno"] = deno_within+deno_between,
    _["strat"] = (nume_within+nume_between)/(deno_within+deno_between),
    _["weight_by"] = deno_by/deno_within,
    _["strat_by"] = nume_by/deno_by,
    _["weight_w"] = deno_within/(deno_within+deno_between),
    _["strat_w"] = nume_within/deno_within,
    _["weight_b"] = deno_between/(deno_within+deno_between),
    _["strat_b"] = nume_between/deno_between
  );
}
