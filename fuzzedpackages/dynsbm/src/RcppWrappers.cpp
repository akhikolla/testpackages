/*
  This file is part of dynsbm.

  dysbm is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  dynsbm is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with dynsbm.  If not, see <http://www.gnu.org/licenses/>
*/
#include<Rcpp.h>
#include<DynSBMBinary.h>
#include<DynSBMDiscrete.h>
#include<DynSBMGaussian.h>
#include<EM.h>
#include<string>
#include<iostream>

using namespace dynsbm;
using namespace Rcpp;
using namespace std;
#ifdef _OPENMP
#include<omp.h>
#endif

// [[Rcpp::export]]
List dynsbmcore(int T, int N, int Q,
		NumericVector Yasvector, const Rcpp::IntegerMatrix & present,
		std::string edgetype, int K,
		IntegerVector clustering,
		int nbit = 20,
		int nbthreads = 1,
		bool isdirected = false,
		bool withselfloop = false,
		bool frozen = false) {
#ifdef _OPENMP
  omp_set_num_threads(nbthreads);
#endif
  ////////////////////////
  ////////////////////////
  if (edgetype=="binary"){
    EM<DynSBMBinary,int> em(T,N,Q,present,isdirected,withselfloop);
    int*** Y;
    allocate3D<int>(Y,T,N,N);
    int p=0;
    for(int j=0; j<N; j++){
      for(int i=0; i<N; i++){
	for(int t=0; t<T; t++){
	  Y[t][i][j] = int(Yasvector[p]);
	  p++;
	}
      }
    }
    em.initialize(as<vector<int> >(clustering),Y,frozen);
    int nbiteff = em.run(Y,nbit,10,frozen);
    NumericMatrix trans(Q,Q);
    for(int q=0;q<Q;q++) for(int l=0;l<Q;l++) trans[l+q*Q] = em.getModel().getTrans(q,l);
    IntegerMatrix membership(N,T);
    for(int t=0;t<T;t++){
      std::vector<int> groups = em.getModel().getGroupsByMAP(t);
      for(int i=0;i<N;i++) membership[i+t*N] = groups[i]+1;
    }
    Rcpp::NumericVector betadims(3);
    betadims[0] = T; betadims[1] = Q; betadims[2] = Q;
    Rcpp::Dimension d(betadims); // get the dim object
    Rcpp::NumericVector beta(d);  // create vec. with correct dims
    for(int t=0;t<T;t++){
      for(int q=0;q<Q;q++){
	for(int l=0;l<Q;l++){
	  beta[l*(Q*T)+q*T+t]= 1-(em.getModel().getBeta(t,q,l)); // cf. paper
	}}}
    double lkl = em.getModel().modelselectionLoglikelihood(Y);
    deallocate3D<int>(Y,T,N,Q);
    return List::create(Rcpp::Named("trans") = trans,
			Rcpp::Named("membership") = membership,
			Rcpp::Named("beta") = beta,
			Rcpp::Named("loglikelihood") = lkl,
			Rcpp::Named("iter") = nbiteff,
			Rcpp::Named("directed") = isdirected,
			Rcpp::Named("self.loop") = withselfloop);
  } else{
    ////////////////////////
    ////////////////////////
    if (edgetype=="discrete"){
      EM<DynSBMDiscrete,int> em(T,N,Q,present,isdirected,withselfloop);
      const_cast<DynSBMDiscrete&>(em.getModel()).setK(K);
      int*** Y;
      allocate3D<int>(Y,T,N,N);
      int p=0;
      for(int j=0; j<N; j++){
	for(int i=0; i<N; i++){
	  for(int t=0; t<T; t++){
	    Y[t][i][j] = int(Yasvector[p]);
	    p++;
	  }
	}
      }
      em.initialize(as<vector<int> >(clustering),Y,frozen);
      int nbiteff = em.run(Y,nbit,10,frozen);
      NumericMatrix trans(Q,Q);
      for(int q=0;q<Q;q++) for(int l=0;l<Q;l++) trans[l+q*Q] = em.getModel().getTrans(q,l);
      IntegerMatrix membership(N,T);
      for(int t=0;t<T;t++){
	std::vector<int> groups = em.getModel().getGroupsByMAP(t);
	for(int i=0;i<N;i++) membership[i+t*N] = groups[i]+1;
      }
      Rcpp::NumericVector betadims(3);
      betadims[0] = T; betadims[1] = Q; betadims[2] = Q;
      Rcpp::Dimension d(betadims); // get the dim object
      Rcpp::NumericVector beta(d);  // create vec. with correct dims
      for(int t=0;t<T;t++){
        for(int q=0;q<Q;q++){
          for(int l=0;l<Q;l++){
            beta[l*(Q*T)+q*T+t]= 1-em.getModel().getBeta(t,q,l); // cf. paper
	  }}}
      Rcpp::NumericVector gammadims(4);
      gammadims[0] = T; gammadims[1] = Q; gammadims[2] = Q; gammadims[3] = K;
      Rcpp::Dimension d2(gammadims);                // get the dim object
      Rcpp::NumericVector gamma(d2);             // create vec. with correct dims
      for(int t=0;t<T;t++){
        for(int q=0;q<Q;q++){
          for(int l=0;l<Q;l++){
            for(int k=0;k<K;k++){
              gamma[k*Q*Q*T+l*(Q*T)+q*T+t]= em.getModel().getMultinomproba(t,q,l,k);
	    }}}}
      double lkl = em.getModel().modelselectionLoglikelihood(Y);
      deallocate3D<int>(Y,T,N,Q);
      return List::create(Rcpp::Named("trans") = trans,
			  Rcpp::Named("membership") = membership,
			  Rcpp::Named("beta") = beta,
			  Rcpp::Named("gamma") = gamma,
			  Rcpp::Named("loglikelihood") = lkl,
			  Rcpp::Named("iter") = nbiteff,
			  Rcpp::Named("directed") = isdirected,
			  Rcpp::Named("self.loop") = withselfloop);
    } else{
      ////////////////////////
      ////////////////////////
      // (edgetype=="continuous"){
      EM<DynSBMGaussian,double> em(T,N,Q,present,isdirected,withselfloop);
      double*** Y;
      allocate3D<double>(Y,T,N,N);
      int p=0;
      for(int j=0; j<N; j++){
	for(int i=0; i<N; i++){
	  for(int t=0; t<T; t++){
	    Y[t][i][j] = double(Yasvector[p]);
	    p++;
	  }
	}
      }
      em.initialize(as<vector<int> >(clustering),Y,frozen);
      int nbiteff = em.run(Y,nbit,10,frozen);
      NumericMatrix trans(Q,Q);
      for(int q=0;q<Q;q++) for(int l=0;l<Q;l++) trans[l+q*Q] = em.getModel().getTrans(q,l);
      IntegerMatrix membership(N,T);
      for(int t=0;t<T;t++){
	std::vector<int> groups = em.getModel().getGroupsByMAP(t);
	for(int i=0;i<N;i++) membership[i+t*N] = groups[i]+1;
      }
      Rcpp::NumericVector betadims(3);
      betadims[0] = T; betadims[1] = Q; betadims[2] = Q;
      Rcpp::Dimension d(betadims); // get the dim object
      Rcpp::NumericVector beta(d);  // create vec. with correct dims
      for(int t=0;t<T;t++){
        for(int q=0;q<Q;q++){
          for(int l=0;l<Q;l++){
            beta[l*(Q*T)+q*T+t]= 1-(em.getModel().getBeta(t,q,l)); // cf. paper
	  }}}
      Rcpp::NumericVector mu(d);  // create vec. with correct dims
      for(int t=0;t<T;t++){
        for(int q=0;q<Q;q++){
          for(int l=0;l<Q;l++){
            mu[l*(Q*T)+q*T+t]= em.getModel().getMu(t,q,l);
	  }}}
      NumericVector sigma(T);
      for(int t=0;t<T;t++) sigma[t] = em.getModel().getSigma(t);
      double lkl = em.getModel().modelselectionLoglikelihood(Y);
      deallocate3D<double>(Y,T,N,Q);
      return List::create(Rcpp::Named("trans") = trans,
			  Rcpp::Named("membership") = membership,
			  Rcpp::Named("beta") = beta,
			  Rcpp::Named("mu") = mu,
			  Rcpp::Named("sigma") = sigma,
			  Rcpp::Named("loglikelihood") = lkl,
			  Rcpp::Named("iter") = nbiteff,
			  Rcpp::Named("directed") = isdirected,
			  Rcpp::Named("self.loop") = withselfloop);
    }
  }
}
  
