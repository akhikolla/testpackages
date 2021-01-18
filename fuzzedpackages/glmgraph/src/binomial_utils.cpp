#include "binomial_utils.h"


	double cycle_binomial_lasso (const arma::vec & Y, const arma::mat & X, const arma::mat & L, const arma::vec & diagL,  
		arma::vec & b_est,arma::vec & b_est0, const arma::uvec & active_set, const arma::vec &multiplier, double lambda1,double lambda2,int initial) {
		
			int n=  X.n_rows;
			int p = X.n_cols;
			arma::vec eta(n);
			arma::vec w(n);
			arma::vec pi(n);
			arma::vec r(n);
			arma::uvec temp1;
			arma::vec shift(n);
			
			arma::vec Z(n);
			arma::vec E(n);
			
			double bj,abj,dj,bb,aa,tau,delta,obj,p1,p2,p3;
			
			if(initial==1){
				double ybar=arma::mean(Y);
				bj=log(ybar/(1-ybar));
				b_est(0)=bj;
				eta=eta.ones()*bj;
			}else{
				bj=b_est(0);
				eta=bj+X*b_est(arma::span(1,p));				
			}
			
			
			pi=1/(1+exp(-eta));
			w=pi%(1-pi);
			
			
			temp1 = arma::find(eta > 10);
			if(temp1.n_elem>0) {
				pi.elem( temp1 )=arma::ones<arma::vec>(temp1.n_elem);
				w.elem( temp1 )=0.0001*arma::ones<arma::vec>(temp1.n_elem);
			}
			temp1 = arma::find(eta < -10);
			if(temp1.n_elem>0) {
				pi.elem( temp1 )=arma::zeros<arma::vec>(temp1.n_elem);
				w.elem( temp1 )=0.0001*arma::ones<arma::vec>(temp1.n_elem);
			}
			
			//update intercept
			delta=arma::mean(w);
			tau=delta*b_est(0)+arma::mean(Y-pi);
			b_est(0)=tau/delta;
			shift=(b_est(0)-bj)*arma::ones<arma::vec>(n);
			eta=eta+shift;
			
		
			for (int j=0; j<p; j++) {
				if (active_set(j) == 0)  continue;
				
				bj = b_est(j+1);
				pi=1/(1+exp(-eta));
				w=pi%(1-pi);
				
				temp1 = arma::find(eta > 10);
				if(temp1.n_elem>0) {
					pi.elem( temp1 )=arma::ones<arma::vec>(temp1.n_elem);
					w.elem( temp1 )=0.0001*arma::ones<arma::vec>(temp1.n_elem);
				}
				temp1 = arma::find(eta < -10);
				if(temp1.n_elem>0) {
					pi.elem( temp1 )=arma::zeros<arma::vec>(temp1.n_elem);
					w.elem( temp1 )=0.0001*arma::ones<arma::vec>(temp1.n_elem);
				}
				
				delta=arma::mean(X.col(j)%X.col(j)%w);
				tau=delta*bj+arma::mean(X.col(j)%(Y-pi));
				abj=lambda2 * arma::dot(b_est(arma::span(1,p)), L.col(j));
				//dj=arma::sum(A.col(j)); 
				dj=diagL(j);
				bb=tau+abj;
				aa=delta+lambda2*dj;
				if(multiplier(j)!=0){
					b_est(j+1)=_Soft(bb,lambda1*multiplier(j))/aa;
				}else{
					b_est(j+1)= tau/delta ;
				}
				

				shift=(b_est(j+1)-bj)*X.col(j);
				eta=eta+shift;
				

	        }
	        E=b_est0(0)+X*b_est0(arma::span(1,p));
	        Z=b_est(0)+X*b_est(arma::span(1,p))+(Y-pi)/w;
	        p1=0.5*arma::mean(w%(Z-E)%(Z-E));
	        p2=arma::sum(arma::abs(b_est(arma::span(1,p))));
	        p3=arma::dot(  b_est(arma::span(1,p)), L*b_est(arma::span(1,p)));
	        obj=p1+p2+p3;
	        return obj;
}	




double cycle_binomial_MCP (const arma::vec & Y,const arma::mat & X, const arma::mat & L, const arma::vec &diagL, 
		arma::vec & b_est,arma::vec & b_est0, const arma::uvec & active_set,const arma::vec &multiplier,double lambda1, double lambda2, double gamma,const std::string &mcpapproach,int initial) {

			int n=  X.n_rows;
			int p = X.n_cols;
			double bj,abj,dj,bb,tau,delta,obj,p1,p2,p3;
			double M=0.25;		
			double aa=0;
			
			arma::vec eta(n);
			arma::vec w(n);
			arma::vec pi(n);
			arma::vec shift(n);
			arma::uvec temp1;
			
			arma::vec Z(n);
			arma::vec E(n);
			
			if(initial==1){
				double ybar=arma::mean(Y);
				bj=log(ybar/(1-ybar));
				b_est(0)=bj;
				eta=eta.ones()*bj;
			}else{
				bj=b_est(0);
				eta=bj+X*b_est(arma::span(1,p));				
			}
			
			pi=1/(1+exp(-eta));
			w=pi%(1-pi);
			
			temp1 = arma::find(eta > 10);
			if(temp1.n_elem>0) {
				pi.elem( temp1 )=arma::ones<arma::vec>(temp1.n_elem);
				w.elem( temp1 )=0.0001*arma::ones<arma::vec>(temp1.n_elem);
			}
			temp1 = arma::find(eta < -10);
			if(temp1.n_elem>0) {
				pi.elem( temp1 )=arma::zeros<arma::vec>(temp1.n_elem);
				w.elem( temp1 )=0.0001*arma::ones<arma::vec>(temp1.n_elem);
			}
			
			
			//update intercept
			
			
			if(mcpapproach=="mmcd"){
				tau=M*b_est(0)+arma::mean(Y-pi);
				b_est(0)=tau/M;
			}else if(mcpapproach=="adaptive" || mcpapproach=="original"){
				delta=arma::mean(w);
				tau=delta*b_est(0)+arma::mean(Y-pi);
				b_est(0)=tau/delta;
			}
			
			delta=arma::mean(w);
			tau=delta*b_est(0)+arma::mean(Y-pi);
			b_est(0)=tau/delta;
				
				
			shift=(b_est(0)-bj)*arma::ones<arma::vec>(n);
			eta=eta+shift;
			

			for (int j=0; j<p; j++) {
				
				if (active_set(j) == 0)  continue;
								
				bj = b_est(j+1);
				pi=1/(1+exp(-eta));
				w=pi%(1-pi);
				
				temp1 = arma::find(eta > 10);
				if(temp1.n_elem>0) {
					pi.elem( temp1 )=arma::ones<arma::vec>(temp1.n_elem);
					w.elem( temp1 )=0.0001*arma::ones<arma::vec>(temp1.n_elem);
				}
				temp1 = arma::find(eta < -10);
				if(temp1.n_elem>0) {
					pi.elem( temp1 )=arma::zeros<arma::vec>(temp1.n_elem);
					w.elem( temp1 )=0.0001*arma::ones<arma::vec>(temp1.n_elem);
				}
								
								
				dj=diagL(j);
				abj=lambda2 * arma::dot(b_est(arma::span(1,p)), L.col(j));
				M=0.25*arma::mean(X.col(j)%X.col(j));
				delta=arma::mean(X.col(j)%X.col(j)%w);

				if(mcpapproach=="mmcd"){
					tau=M*bj+arma::mean(X.col(j)%(Y-pi));
					aa=M+lambda2*dj;
				}else if(mcpapproach=="adaptive" || mcpapproach=="original"){
					tau=delta*bj+arma::mean(X.col(j)%(Y-pi));
					aa=delta+lambda2*dj;
				}
				bb=tau+abj;
									
				
				if(multiplier(j)!=0){
					if(mcpapproach=="mmcd"){
						if( fabs(bb)<=(aa*gamma*lambda1) ){
							b_est(j+1)=_Soft(bb,lambda1)/(aa-1/gamma); 
						}else{
							b_est(j+1)=bb/aa;
						}
					}else if(mcpapproach=="adaptive"){
						if(fabs(bb)<=(gamma*lambda1)){
							b_est(j+1)=_Soft(bb,lambda1)/aa/(1-1/gamma);
						}else{
							b_est(j+1)=bb/aa;
						}
					}else if(mcpapproach=="original"){
						if(fabs(bb)<=(aa*gamma*lambda1)){
							b_est(j+1)=_Soft(bb,lambda1)/(aa-1/gamma);
						}else{
							b_est(j+1)=bb/aa;
						}
					}
				}else{
					b_est(j+1)= tau/delta;
				}
				
				shift=(b_est(j+1)-bj)*X.col(j);
				eta=eta+shift;
	     }	
	     
	      	E=b_est0(0)+X*b_est0(arma::span(1,p));
	        Z=b_est(0)+X*b_est(arma::span(1,p))+(Y-pi)/w;
	        p1=0.5*arma::mean(w%(Z-E)%(Z-E));
	        p2=arma::sum(arma::abs(b_est(arma::span(1,p))));
	        p3=arma::dot(  b_est(arma::span(1,p)), L*b_est(arma::span(1,p)));
	        obj=p1+p2+p3;
	        return obj;
	     	       
}




















