#include "gaussian_utils.h"


double _Soft(double b, double lambda) {
		double b1 = fabs(b) - lambda;
		int s = (b <= 0) ? -1 : 1;
		double b2 = (b1 <= 0) ? 0 : b1;
		return (s * b2);
}


double cycle_gaussian_lasso (const arma::vec & Y, const arma::mat & X,  const arma::mat & L,  const arma::vec & diagL,
		arma::vec & b_est, const arma::uvec & active_set, const arma::vec &multiplier, double lambda1,double lambda2,bool standardizeX,int initial) {
		
		int n=  X.n_rows;
		int p = X.n_cols;
		arma::vec rr(n);
		double bj, zz, xi, bb, aa, obj, p1, p2, p3;
		arma::vec E(n);
		
		if(standardizeX){
			if(initial==1) rr=Y;
			else rr = Y - X*b_est(arma::span(1,p)); 
			for (int j=0; j<p; j++) {
				if (active_set(j) == 0)  continue;
	        	bj = b_est(j+1);        	
				zz = arma::dot(X.col(j), rr) / n + bj;
				xi = lambda2 * arma::dot(b_est(arma::span(1,p)), L.col(j));
				bb = zz + xi;
				aa = lambda2 * diagL(j) +1;
				if(multiplier(j)!=0) b_est(j+1) = _Soft(bb,lambda1*multiplier(j))/aa;
				else b_est(j+1) = zz  ;
				rr = rr - (b_est(j+1) - bj) * X.col(j);

			}
		}else{
			if(initial==1){
				b_est(0)=arma::mean(Y);
				rr=Y;
			}else{
				rr = Y -X*b_est(arma::span(1,p)); 
				b_est(0)=arma::mean(rr);
				rr= Y -b_est(0)-X*b_est(arma::span(1,p)); 
			}
			for (int j=0; j<p; j++) {
				if (active_set(j) == 0)  continue;
	        	bj = b_est(j+1);
				zz = arma::dot(X.col(j), rr) / n + bj*arma::dot(X.col(j),X.col(j))/n;           //update
				xi = lambda2 * arma::dot(b_est(arma::span(1,p)), L.col(j));
				bb = zz + xi;
				aa = lambda2 * diagL(j) + arma::dot(X.col(j),X.col(j))/n ;  			//update
				if(multiplier(j)!=0) b_est(j+1) = _Soft(bb,lambda1*multiplier(j))/aa;
				else b_est(j+1) = zz  ;
				rr = rr - (b_est(j+1) - bj) * X.col(j);
			}
		}	
		
		E=b_est(0)+X*b_est(arma::span(1,p));
	    p1=0.5*arma::mean((Y-E)%(Y-E));
	    p2=arma::sum(arma::abs(b_est(arma::span(1,p))));
	    p3=arma::dot(  b_est(arma::span(1,p)), L*b_est(arma::span(1,p)));
	    obj=p1+p2+p3;
	    return obj;

}


double cycle_gaussian_MCP (const arma::vec & Y,const arma::mat & X, const arma::mat & L, const arma::vec & diagL, 
		arma::vec & b_est, const arma::uvec & active_set,const arma::vec &multiplier,double lambda1, double lambda2, double gamma,bool standardizeX,int initial) {

		int n=  X.n_rows;
		int p = X.n_cols;
		arma::vec rr(n);
		double bj, zz, xi, bb;
		double t1,t2;
		double obj,p1,p2,p3;
		arma::vec E(n);

		if(standardizeX){
			if(initial==1) rr=Y;
			else rr = Y - X*b_est(arma::span(1,p)); 
			for (int j=0; j<p; j++) {
				if (active_set(j) == 0)  continue;
	       		bj = b_est(j+1);
				zz = arma::dot(X.col(j), rr) / n + bj;
				xi = lambda2 * arma::dot(b_est(arma::span(1,p)), L.col(j));
				bb = zz + xi;

				if(multiplier(j)!=0) {
					t1 = 1 + lambda2 * diagL(j);
					t2 = gamma * lambda1 * t1;
					if (fabs(bb) <= t2) {
						b_est(j+1) = _Soft(bb, lambda1*multiplier(j)) * gamma / (gamma * t1 - 1);
					} else {
						b_est(j+1) = bb / t1; 
					}
				}
				else b_est(j+1) = zz  ;
				rr = rr - (b_est(j+1) - bj) * X.col(j);
			}
		}else{
			if(initial==1){
				b_est(0)=arma::mean(Y);
				rr=Y;
			}else{
				rr = Y -X*b_est(arma::span(1,p)); 
				b_est(0)=arma::mean(rr);
				rr= Y -b_est(0)-X*b_est(arma::span(1,p)); 
			}
			
			for (int j=0; j<p; j++) {
				if (active_set(j) == 0)  continue;
	       		bj = b_est(j+1);
				zz = arma::dot(X.col(j), rr) / n + bj*arma::dot(X.col(j),X.col(j))/n;           //update
				xi = lambda2 * arma::dot(b_est(arma::span(1,p)), L.col(j));
				bb = zz + xi;

				if(multiplier(j)!=0) {
					t1 = arma::dot(X.col(j),X.col(j))/n + lambda2 * diagL(j);
					t2 = gamma * lambda1 * t1;
					if (fabs(bb) <= t2) {
						b_est(j+1) = _Soft(bb, lambda1*multiplier(j)) * gamma / (gamma * t1 - 1);
					} else {
						b_est(j+1) = bb / t1; 
					}
				}
				else b_est(j+1) = zz  ;
				rr = rr - (b_est(j+1) - bj) * X.col(j);
			}
		}	
		
		E=b_est(0)+X*b_est(arma::span(1,p));
	    p1=0.5*arma::mean((Y-E)%(Y-E));
	    p2=arma::sum(arma::abs(b_est(arma::span(1,p))));
	    p3=arma::dot(  b_est(arma::span(1,p)), L*b_est(arma::span(1,p)));
	    obj=p1+p2+p3;
	    return obj;
}




