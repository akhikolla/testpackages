#include "binomial_utils.h" 

using namespace Rcpp ;

RcppExport SEXP cdfit_binomial(SEXP X_, SEXP y_, SEXP L_, SEXP diagL_, SEXP penalty_, SEXP lambda1_, SEXP lambda2_,SEXP gamma_, SEXP multiplier_, SEXP eps_, SEXP max_iter_, SEXP dfmax_,SEXP mcpapproach_, SEXP warn_) {


	arma::mat X= Rcpp::as<arma::mat>(X_);
	arma::vec Y= Rcpp::as<arma::vec>(y_);
	arma::mat L= Rcpp::as<arma::mat>(L_);
	arma::vec diagL=Rcpp::as<arma::vec>(diagL_);

	std::string penalty = Rcpp::as<std::string>(penalty_);
	arma::vec lambda2=Rcpp::as<arma::vec>(lambda2_);	
	double eps=Rcpp::as<double>(eps_);

	int max_iter=Rcpp::as<int>(max_iter_);
	double gamma=Rcpp::as<double>(gamma_);
	
	arma::vec multiplier=Rcpp::as<arma::vec>(multiplier_);
	int dfmax=Rcpp::as<int>(dfmax_);
	std::string mcpapproach = Rcpp::as<std::string>(mcpapproach_);
	bool warn=Rcpp::as<bool>(warn_);

	int q = lambda2.n_elem;
	
	arma::vec lambda1= Rcpp::as<arma::vec>(lambda1_);
	arma::field<arma::vec> lambda1s(q);
	for(int i=0;i<q;i++) lambda1s(i)=lambda1;
	

	double _lambda1, _lambda2;
	int i, j, iter,breakpoint;
	int p = X.n_cols;
	
	arma::field<arma::mat> b_est_all(q);
	arma::field<arma::vec> b_iter_all(q);
	

	arma::vec b_est1(p+1);
	arma::vec b_est2(p+1);
	arma::vec b_est0(p+1);
	arma::uvec active_set1(p);
	arma::uvec active_set2(p);
	
	
	double obj1=0,obj2=0;
	
	
	for (j=0; j<q; j++) {

	    arma::mat b_estj(p+1,lambda1s(j).n_elem);
	    arma::vec b_iterj=arma::vec(lambda1s(j).n_elem);
	    b_iterj.fill(NA_REAL);
	    b_estj.fill(NA_REAL);

	    _lambda2 = lambda2(j);
	    
	    for (i=0; i<(int) lambda1s(j).n_elem; i++) {
	    
	
	        _lambda1 = lambda1s(j)(i);  
	        if (i == 0) b_est1.zeros();	
	        active_set1.ones();
	        iter=0;
	    	breakpoint=0;
	    	
	    	
	        if(penalty=="lasso")  obj1=cycle_binomial_lasso (Y, X, L, diagL, b_est1, b_est0, active_set1, multiplier, _lambda1, _lambda2,1);
	    	else if(penalty=="MCP")  obj1=cycle_binomial_MCP (Y, X, L, diagL, b_est1, b_est0, active_set1, multiplier, _lambda1, _lambda2, gamma,mcpapproach,1);
	        active_set1 = (arma::abs(b_est1(arma::span(1,p))) > arma::datum::eps   ); 
	        if(arma::all(active_set1==0)) continue;
			
			
	        while (TRUE) {															 
	        	active_set1 = (arma::abs(b_est1(arma::span(1,p))) > arma::datum::eps   ); 
	    	    b_est2=b_est1+1;	
	            while (TRUE) {	
	            	iter += 1;
	            	if( std::abs( (obj1-obj2)/obj1 ) < eps ) break;
	            	b_est2 = b_est1;
	            	obj2=obj1;
	            	if (iter >= max_iter) {
	            		if(warn){
	            			Rcpp::Rcout  << "Exceeds the maximum iteration for lambda1: " << _lambda1 <<", lambda2: " << _lambda2;
	            			Rcpp::Rcout  << ". Further lambda1 will not be considered." << std::endl;
	            		}
	            		breakpoint=1;
	            		break;
	            	}
	            	b_est0=b_est1;
	            	if(penalty=="lasso")  obj1=cycle_binomial_lasso (Y, X, L, diagL, b_est1, b_est0, active_set1, multiplier, _lambda1, _lambda2, 0);
	    			else if(penalty=="MCP")  obj1=cycle_binomial_MCP (Y, X, L, diagL, b_est1, b_est0, active_set1, multiplier, _lambda1, _lambda2, gamma,mcpapproach, 0);
	    			
	            }
	           	if(breakpoint==1) break;
	           	if(penalty=="lasso")  obj1=cycle_binomial_lasso (Y, X, L, diagL, b_est1, b_est0, arma::ones<arma::uvec>(p), multiplier, _lambda1, _lambda2, 0);
	    		else if(penalty=="MCP")  obj1=cycle_binomial_MCP (Y, X, L, diagL, b_est1, b_est0, arma::ones<arma::uvec>(p), multiplier, _lambda1, _lambda2, gamma, mcpapproach,0);
	           
	            active_set2 = (arma::abs(b_est1(arma::span(1,p))) > arma::datum::eps   ); 
	            if (arma::accu(active_set2 != active_set1) <1) break;
			}   
			if(breakpoint==1) break;
			
	    	if ( (int) arma::sum(   arma::abs( b_est1(arma::span(1,p)) ) > arma::datum::eps   ) > dfmax   )   {
	        	if(warn){
	        		Rcpp::Rcout  << "Exceeds dfmax "<<dfmax<<", for lambda1: " << _lambda1 <<", lambda2:" << _lambda2;
	            	Rcpp::Rcout  << ". Further lambda1 will not be considered." << std::endl;
	            }
	        	break;
	        }
	        
			b_iterj(i)=i;
			double intercept=b_est1(0);
	        b_est1(arma::find(arma::abs(b_est1) < eps)).zeros();
	        b_est1(0)=intercept;
	        b_estj.col(i) = b_est1;
	    }
	    b_est_all(j) = b_estj;
	    b_iter_all(j)= b_iterj;
	    
	}
	
	Rcpp::List res=List::create(b_est_all,b_iter_all);
	return Rcpp::wrap(res);
	
}




