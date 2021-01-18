#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// soft thresholding
double ST1(double z,double gam){
	double sparse=0;
	if(z>0 && gam<fabs(z)) return(z-gam);
	if(z<0 && gam<fabs(z)) return(z+gam);
	if(gam>=fabs(z)) return(sparse);
	else return(0);
	}

// [[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]

List M_Update_Rcpp(arma::mat Mean, arma::mat D, arma::mat Uinv, arma::mat Vinv, 
	arma::vec nc, arma::mat weightmat, double rho, double lambda, int C, int r, double mtol){

	bool notconverged = true; 

	int p = Mean.n_cols; 
	double alpha = 1; 
	double alphatemp; 
	double rresid; 
	double sresid;
	double out; 
	int K = C*(C-1)/2; 
	int Cinner = C-1; 

	double N = sum(nc(span::all)); 
	vec pihat = nc/N; 

	int k = 1; 
	int d = 1; 
	mat lambdaweight(K*r,p);  

	mat Mtemp(C*r, p); 
	mat Q(K*r, p); 
	mat Dtemp(K*r, p); 
	mat Dnew = D;

	mat Mnew(C*r, p); 

	mat Gnew(K*r, p);
	mat Gtemp(K*r, p); 

	mat Dhatnew = D; 
	mat Dhattemp; 
	mat righttemp =  zeros<mat>(r, p);

	lambdaweight = lambda*weightmat; 
	
	do{		
	 	/* Mean updates; note that we use Kr x p mean matrix */
	 	/* Mean update for class 1; indexing is different and can be faster */
	 	righttemp  = zeros<mat>(r, p); 
	 	
	 	for(int j=0;  j < Cinner; ++j){
	 		righttemp = righttemp + Dhatnew(span(j*r, (j+1)*r - 1), span::all); 
	 	}
	 	
	 	Mtemp(span(0, r-1), span::all) = Mean(span(0,r-1), span::all) + (.5/pihat(0))*Uinv*righttemp*Vinv; 

	 	/* Mean update for class C; */
		righttemp = zeros<mat>(r, p); 
	 	
	 	for(int i=0;  i < Cinner; ++i){
	 		d = (Cinner)*i + Cinner - 1 - (i*(i+1))/2; 
	 		righttemp = righttemp - Dhatnew(span(d*r, (d+1)*r - 1), span::all); 
	 	}

	 	Mtemp(span(Cinner*r, (Cinner+1)*r-1), span::all) = Mean(span(Cinner*r, (Cinner+1)*r-1), span::all) 
	 	+ (.5/pihat(Cinner))*Uinv*righttemp*Vinv; 


	 	/* Mean for classes 2 to C-1; need to loop to add D+ and subtract D- (see manuscript) */
	 	if(C > 2){
	 		/* fix mean index betwen 1 and Cinner*/
			for(int j=1; j < Cinner; ++j){
				
				/* clear righttemp */
				righttemp = zeros<mat>(r, p); 

				// subtract columns where j > i; meaning use j as a column index //
				for(int i = 0; i < j; ++i){
					d = Cinner*i + j - 1 - (i*(i+1))/2; 
					righttemp = righttemp - Dhatnew(span(d*r, (d+1)*r - 1), span::all); 
				}


				// add rows where j < i; meaning use j as a row index  // 			
				for(int i=j+1; i < (Cinner+1); ++i){
					d = Cinner*j + i - 1 - (j*(j+1))/2; 
					righttemp = righttemp + Dhatnew(span(d*r, (d+1)*r - 1), span::all); 
				}

				Mtemp(span(j*r, (j+1)*r - 1), span::all) =  Mean(span(j*r, (j+1)*r - 1), span::all) 
				+ (.5/pihat(j))*Uinv*righttemp*Vinv;
		
			}
		}
		
	 	
	 	/* compute entire Q mat (Kr x p) for fast soft thresholding */
	 	for(int i=0; i < Cinner; ++i){
		 	for(int j=i+1; j < Cinner+1; ++j){

				d = Cinner*i + j - 1 - (i*(i+1))/2; 	 
				Q(span(d*r, (d+1)*r -1), span::all) = -(1/rho)*Dhatnew(span(d*r, (d+1)*r -1), span::all) 
					- Mtemp(span(j*r, (j+1)*r-1), span::all) + Mtemp(span(i*r, (i+1)*r-1), span::all); 

		 	}
		}
	



	 	/* perform soft-thresholding of Q matrix with scaled weights */
	 	for(int i=0; i < (K*r); ++i){
	 		for(int j=0; j<p; ++j){
	 			out = Q(i,j); 
	 			Gtemp(i,j) = ST1(out, lambdaweight(i,j)/rho); 
	 		}
	 	}

	 	/* update D; use same loop at for Q */
	 	for(int i=0; i < Cinner; ++i){
	 		for(int j=i+1; j < Cinner+1; ++j){

				d = Cinner*i + j - 1 - (i*(i+1))/2; 	 

				Dtemp(span(d*r, (d+1)*r -1), span::all) = Dhatnew(span(d*r, (d+1)*r -1), span::all) + 
				rho*(Gtemp(span(d*r, (d+1)*r -1), span::all) + Mtemp(span(j*r, (j+1)*r-1), span::all) 
					- Mtemp(span(i*r, (i+1)*r-1), span::all)); 

	 		}
	 	}	

	 	/* change this number to impose fixed iteration restarts */
	  	if(k%10000 == 0){ 
	 		alphatemp = 1; 
	 		Dhattemp = Dtemp; 

	 	}  else {
	 	
	 		alphatemp = (1 + sqrt(1 + 4*pow(alpha,2)))/2; 
	  		Dhattemp = Dtemp + ((alpha - 1)/(alphatemp))*(Dtemp - Dnew);  
	
		}

	 if(k > 10){
	 	 rresid = 0; 
	 	 	for(int i=0; i < Cinner ; ++i){
	 			for(int j=i+1; j < Cinner+1; ++j){
				d = Cinner*i + j - 1 - (i*(i+1))/2; 	
	 			rresid =  rresid + norm(Mtemp(span(i*r, (i+1)*r-1), span::all) 
	 				- Mtemp(span(j*r, (j+1)*r-1), span::all) - Gtemp(span(d*r, (d+1)*r -1), span::all), "fro")/(K*r*p); 
	 		}
	 	}	

	 	sresid = norm(Dtemp - Dnew, "fro")/norm(Dnew, "fro"); 

	 if(rresid < mtol && sresid < mtol){
		 notconverged = false;
	 } else {
	 	 notconverged = true;
	 }

	 }

	 Mnew = Mtemp;
	 Gnew = Gtemp; 
	 Dnew = Dtemp; 
	 alpha = alphatemp; 
	 Dhatnew = Dhattemp; 

	 k++; 

	 if(k > 10000){
	  	notconverged = false; 
	 }
	} while(notconverged); 
	
  return List::create(Named("M") = wrap(Mnew), 
  					  Named("G") = wrap(Gnew), 
  					  Named("D") = wrap(Dnew), 
  					  Named("k") = wrap(k), 
  					  Named("r_resid") = wrap(rresid), 
  					  Named("d_resid") = wrap(sresid),
  					  Named("rho") = wrap(rho), 
  					  Named("lambda") = wrap(lambda));
}
