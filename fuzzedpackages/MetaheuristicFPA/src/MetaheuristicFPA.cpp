#include <stdio.h>
#include <RcppArmadilloExtensions/sample.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


NumericVector simpleBounds(NumericVector x, NumericVector LowerV, NumericVector UpperV){
  for(int i=0;i<x.size();i++){
    if(x[i]<LowerV[i]) x[i] = LowerV[i];
    if(x[i]>UpperV[i]) x[i] = UpperV[i];
  }
  return x;
}

NumericVector getLevy(NumericVector V, NumericVector U,double beta){
  NumericVector ret(V.size());
  V = Rcpp::abs(V);
  for(int i=0;i<V.size();i++){
    ret[i] = U[i]/(pow(V[i],(1/beta)));

  }

  return ret;
}


//[[Rcpp::export]]
List fpa_optim(double N=25,double p=0.8,double beta=1.5,double eta=0.1,int maxiter=2000,bool randEta=false,double gloMin=0,double objfit=1e-05, double D=2, double Lower=-10, double Upper=10, Function FUN=NULL){

  NumericVector LowerV(D);
  NumericVector UpperV(D);
  for(int i=0;i<D;i++){
    LowerV[i]=Lower;
    UpperV[i]=Upper;
  }
  NumericVector Gbest(D);
  NumericVector Fitness(N);
  NumericMatrix XX(N,D);
  NumericMatrix YY(N,D);

  //set to infinity
  double FitnessBest = 999999999999999999;
  double minIndex = 0;

  for(int i=0;i<N;i++){
    for(int j=0;j<D;j++){
      XX(i,j)=LowerV[j]+(UpperV[j]-LowerV[j])*R::runif(0,1);
    }
    Fitness[i] = Rcpp::as<double>(FUN(XX.row(i)));

  }

  for(int i=0;i<N;i++){
    if(Fitness[i]<FitnessBest){
      FitnessBest = Fitness[i];
      minIndex = i;
    }
  }

  Gbest = XX.row(minIndex);
  NumericVector dataIter;
  YY = clone(XX);

  //algorithm start
  int iter = 0;
  bool stop = false;

  double term1 = R::gammafn(1+beta)*sin(M_PI*beta/2);
  double term2 = R::gammafn((1+beta)/2)*beta*pow(2,(beta-1)/2);
  double sigma = pow(term1/term2,1/beta);



  NumericVector levy;
  NumericVector indexAcak(N);
  for(int i=0;i<N;i++) indexAcak[i] = i;

  while(iter<maxiter && !stop){
    for(int i=0;i<N;i++){
      double dumpP = R::runif(0,1);

      if(dumpP<p){
        //gloMin pollination
        NumericVector V = Rcpp::rnorm(D,0,1);
        NumericVector U = Rcpp::rnorm(D,0,sigma);
        levy = getLevy(V,U,beta);

        for(int j=0;j<D;j++){
          if(randEta) YY(i,j) = XX(i,j)+fabs(R::rnorm(0,eta))*levy[j]*(Gbest[j]-XX(i,j));
          else YY(i,j) = XX(i,j)+eta*levy[j]*(Gbest[j]-XX(i,j));
        }

      }
      else{
        //local pollination
        double eps = R::runif(0,1);

        NumericVector IndexJk;

        //select 2 random value
        do{
          IndexJk = Rcpp::RcppArmadillo::sample(indexAcak,2,true);

        } while (IndexJk[0]==i || IndexJk[1]==i);


        for(int j=0;j<D;j++){
          YY(i,j) = XX(i,j)+eps*(XX(IndexJk[0],j)-XX(IndexJk[1],j));
        }
      }

      YY.row(i) = simpleBounds(YY.row(i),LowerV,UpperV);

      double Fnew = Rcpp::as<double>(FUN(YY.row(i)));


      if(Fnew<Fitness[i]){
        Fitness[i] = Fnew;
        XX.row(i) = YY.row(i);
      }

      if(Fnew<FitnessBest){
        FitnessBest = Fnew;
        Gbest = YY.row(i);
      }


    }
    if(fabs(FitnessBest-gloMin)<objfit) stop = true;
    else iter++;

    dataIter.push_back(FitnessBest);

  }


  return Rcpp::List::create(Rcpp::Named("min_fitness") = FitnessBest,
                      Rcpp::Named("best_solution") = Gbest,
                      Rcpp::Named("Iteration") = iter);

}
