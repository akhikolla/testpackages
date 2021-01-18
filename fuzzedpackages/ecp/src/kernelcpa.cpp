// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP srcGetV( SEXP K_ ) ;
SEXP srcGetBandwidth( SEXP X_, SEXP rows_ );
SEXP srcKcpa( SEXP II_, SEXP V_, SEXP H_ );
}


SEXP srcGetV( SEXP K_ ){
BEGIN_RCPP	

using namespace Rcpp;
#define SUM(A) std::accumulate(A.begin() , A.end() , 0.0)

NumericMatrix K(K_);
int N = K.nrow();
NumericMatrix V(N,N);

for(int i=0;i<N;++i)
	for(int j=i;j<N;++j)
		V(i,j) = V(j,i) = SUM(diag(K(Range(i,j),Range(i,j))))-SUM(K(Range(i,j),Range(i,j)))/(j-i+1);
return wrap(V);

END_RCPP
}


SEXP srcGetBandwidth( SEXP X_, SEXP rows_ ){
BEGIN_RCPP	

NumericMatrix X(X_);
NumericVector rows(rows_);
int N = rows.size();
NumericVector u(N*N,(double)0);
for(int i=0;i<N;++i){
	for(int j=0;j<N;++j){
		u[i*N+j] = sum((X(rows[i],_)-X(rows[j],_))*(X(rows[i],_)-X(rows[j],_)));
	}
}
return wrap(u);

END_RCPP
}


SEXP srcKcpa( SEXP II_, SEXP V_, SEXP H_ ){
BEGIN_RCPP	

NumericMatrix II(II_);
NumericMatrix V(V_);
IntegerMatrix H(H_);
int N = V.nrow();
int L = H.nrow();

for(int k=1;k<L;++k){
	for(int i=k;i<N;++i){
		for(int j=k-1;j<i;++j){
			double tmp = II(k-1,j) + V(j+1,i);
			if(tmp < II(k,i)){
				II(k,i) = tmp;
				H(k,i) = j+1;//fix indexing differences between R and C++
			}
		}
	}
}
return List::create(II,H);

END_RCPP
}
