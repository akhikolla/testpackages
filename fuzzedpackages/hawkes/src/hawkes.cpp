#include "hawkes_utils.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
 
 
int getDimension(SEXP lambda0){
  Rcpp::NumericVector lambda0_internal(lambda0);  
  return lambda0_internal.size();
}
 
int attribute(double alea, double I_star,const arma::vec& m_lambda)
{
  int index = 0;
	double cumul = m_lambda(0);
	while (alea > (cumul / I_star))
	{
		index = index + 1;
		cumul = cumul + m_lambda(index);
	}
	return (index);
}

// [[Rcpp::export]]
std::vector<std::vector<double> > simulateHawkes(SEXP lambda0,SEXP alpha,SEXP beta,SEXP horizon)
{
  RNGScope scope;
  Rcpp::NumericVector tempUnif(1);
  
  int dimension = getDimension(lambda0);
  double m_horizon = as<double>(horizon);
  
  std::vector<std::vector<double> > history;
  for (int i=0;i<dimension;i++)
  {
    std::vector<double> a;
    history.push_back(a);
  }
  if (dimension == 1)
  {
    double m_lambda0 = as<double>(lambda0);
    double m_alpha = as<double>(alpha);
    double m_beta = as<double>(beta);
    if(m_beta<m_alpha){
      stop("Unstable. You must have alpha < beta");
    }
    
    double lambda_star = m_lambda0;
	  double dlambda = 0.0,t=0;
	  //first event
    tempUnif = runif(1);
	  double U = tempUnif(0);
	  double s = -(1.0 / lambda_star) * log(U);
	  if (s <= m_horizon)
	  {
		  history[0].push_back(s);
		  dlambda = m_alpha;
		  t = s;
	  }
  	else
  	{
  		return (history);
  	}
  	//general routine
  	while (true)
  	{
  		lambda_star = m_lambda0+dlambda*exp(-m_beta*(s-t));
      tempUnif = runif(1);
  		U = tempUnif(0);
  		s = s - (1.0 / lambda_star) * log(U);
  		if (s > m_horizon)
  		{
  			return (history);
  		}
      tempUnif = runif(1);
  		double D = tempUnif(0);
  		if (D <= (m_lambda0+dlambda*exp(-m_beta*(s-t))) / lambda_star)
  		{
  			history[0].push_back(s);
  			dlambda = dlambda*exp(-m_beta*(s-t)) + m_alpha;
  			t=s;
  		}
  	}
  }
  else{
    arma::mat dlambda(dimension,dimension,arma::fill::zeros);
    Rcpp::NumericVector lambda0_internal(lambda0);
    Rcpp::NumericMatrix alpha_internal(alpha);
    Rcpp::NumericVector beta_internal(beta);
    
    arma::vec m_lambda0(lambda0_internal.begin(),dimension,false);
    arma::mat m_alpha(alpha_internal.begin(),dimension,dimension,false);
    arma::vec m_beta(beta_internal.begin(),dimension,false);
    arma::mat m_beta_matrix = diagmat(m_beta);
    arma::vec m_lambda(dimension);
    
    
    
    checkStability(m_beta_matrix,m_alpha);
    
  	double lambda_star = 0.0;
	  
  	double t=0;
  	for (int i = 0; i < dimension; i++)
  	{
  		lambda_star += m_lambda0(i);
  		m_lambda(i) = m_lambda0(i);
  	}
    
  	//first event
    tempUnif = runif(1);
  	double U = tempUnif(0);
  	double s = -(1.0 / lambda_star) * log(U);
    
    
  	if (s <= m_horizon)
  	{
        tempUnif = runif(1);
  		double D = tempUnif(0);
  		int n0 = attribute(D, lambda_star,m_lambda);
  		history[n0].push_back(s);
  
  		for (int i=0;i<dimension;i++)
  		{
  			dlambda(i,n0) = m_alpha(i,n0);
  			m_lambda(i) = m_lambda0(i)+m_alpha(i,n0);
  		}
  	}
  	else
  	{
  		return (history);
  	}
  	t=s;
    //general routine
  	lambda_star = 0;
  	for (int i = 0; i < dimension; i++)
  	{
  		lambda_star = lambda_star + m_lambda(i);
  	}
  	while (true)
  	{
        tempUnif = runif(1);
  		U = tempUnif(0);
  		s = s - (1.0 / lambda_star) * log(U);
  		if (s <= m_horizon)
  		{
            tempUnif = runif(1);
  			double D = tempUnif(0);
  			double I_M = 0.0;
  			for (int i = 0; i < dimension; i++)
  			{
  				double dl = 0.0;
  				for (int j = 0; j < dimension; j++)
  				{
  					dl += dlambda(i,j)*exp(-m_beta(i)*(s-t));
  				}
  				m_lambda[i] = m_lambda0(i)+dl;
  				I_M = I_M + m_lambda(i);
  			}
  			if (D <= (I_M / lambda_star))
  			{
  				int n0 = attribute(D, lambda_star,m_lambda);
  				history[n0].push_back(s);
  				lambda_star=0.0;
  				for (int i=0;i<dimension;i++)
  				{
  					double dl=0.0;
  					for (int j = 0; j < dimension; j++)
  					{
  						dlambda(i,j) = dlambda(i,j)*exp(-m_beta(i)*(s-t));
  						if (n0==j)
  						{
  							dlambda(i,n0) += m_alpha(i,n0);
  						}
  						dl +=dlambda(i,j);
  					}
  					lambda_star+= m_lambda0(i)+dl;
  				}
  				t=s;
  			}
  			else
  			{
  				lambda_star = I_M;
  			}
  
  		}
  		else
  		{
  			return (history);
  		}
  	}
  }
}

// [[Rcpp::export]]
std::vector<double>  jumpMean(SEXP lambda0,SEXP alpha,SEXP beta,SEXP tau)
{
  int dimension = getDimension(lambda0);
  double m_tau = as<double>(tau);
  
  std::vector<double> res;
  
  if (dimension == 1)
  {
    double m_lambda0 = as<double>(lambda0);
    double m_alpha = as<double>(alpha);
    double m_beta = as<double>(beta);
    if(m_beta<m_alpha){
      stop("Unstable. You must have alpha < beta");
    }
    double lambda = m_beta*m_lambda0/(m_beta-m_alpha);
    res.push_back(lambda*m_tau); 
    return res;
  }else{
    Rcpp::NumericVector lambda0_internal(lambda0);
    Rcpp::NumericMatrix alpha_internal(alpha);
    Rcpp::NumericVector beta_internal(beta);
    
    arma::vec m_lambda0(lambda0_internal.begin(),dimension,false);
    arma::mat m_alpha(alpha_internal.begin(),dimension,dimension,false);
    arma::vec m_beta(beta_internal.begin(),dimension,false);
    arma::mat m_beta_matrix = diagmat(m_beta);
    
    checkStability(m_beta_matrix,m_alpha);
    
    arma::mat beta_minus_alpha = m_beta_matrix- m_alpha;
    arma::mat matrixBetaMinusAlpha_inv = inv(beta_minus_alpha);
    arma::mat temp = matrixBetaMinusAlpha_inv*m_beta_matrix;
	  arma::vec means = temp*m_lambda0;
	  means = means*m_tau;
    for(int i=0;i<dimension;i++){
      res.push_back(means[i]);
    }
	  return res;
  }
}
// [[Rcpp::export]]
arma::mat jumpVariance(SEXP lambda0,SEXP alpha,SEXP beta,SEXP tau)
{
  int dimension = getDimension(lambda0);
  double m_tau = as<double>(tau);
  
  arma::mat res = arma::zeros(dimension,dimension);
  
  if (dimension == 1)
  {
    double m_lambda0 = as<double>(lambda0);
    double m_alpha = as<double>(alpha);
    double m_beta = as<double>(beta);
    if(m_beta<m_alpha){
      stop("Unstable. You must have alpha < beta");
    }
    double gamma = m_lambda0/(1-m_alpha/m_beta);
    double kappa = 1./(1.-m_alpha/m_beta);
    double v=gamma*(m_tau*kappa*kappa+(1-kappa*kappa)*(1-exp(-(m_beta-m_alpha)*m_tau))/(m_beta-m_alpha));
    res(0,0)= v; 
    return res;
  }else{
    Rcpp::NumericVector lambda0_internal(lambda0);
    Rcpp::NumericMatrix alpha_internal(alpha);
    Rcpp::NumericVector beta_internal(beta);
    
    arma::vec m_lambda0(lambda0_internal.begin(),dimension,false);
    arma::mat m_alpha(alpha_internal.begin(),dimension,dimension,false);
    arma::vec m_beta(beta_internal.begin(),dimension,false);
    arma::mat m_beta_matrix = diagmat(m_beta);
    
    checkStability(m_beta_matrix,m_alpha);
    
    arma::mat temp0 = computeC5(m_lambda0,m_alpha,m_beta_matrix,m_tau);
    arma::vec tempVec = expectedStationaryLambda(m_lambda0,m_alpha,m_beta_matrix,0);
    arma::mat temp = grandLambdaInfini(m_lambda0,m_alpha,m_beta_matrix,0)+(m_alpha*vectorToDiagonalMatrix(tempVec));
    arma::mat J1 = temp0*temp;
    arma::vec tempVec2 = expectedStationaryLambda(m_lambda0,m_alpha,m_beta_matrix,m_tau);
    res = J1+arma::trans(J1)+(vectorToDiagonalMatrix(tempVec2)*m_tau);
  
	  return res;
  }
}
// [[Rcpp::export]]
arma::mat jumpAutocorrelation(SEXP lambda0,SEXP alpha,SEXP beta,SEXP tau,SEXP lag)
{
  int dimension = getDimension(lambda0);
  double m_tau = as<double>(tau);
  double m_lag = as<double>(lag);
  
  arma::mat res = arma::zeros(dimension,dimension);
  
  if (dimension == 1)
  {
    
    double m_alpha = as<double>(alpha);
    double m_beta = as<double>(beta);
    if(m_beta<m_alpha){
      stop("Unstable. You must have alpha < beta");
    }
    double alpha2=m_alpha*m_alpha,beta2=m_beta*m_beta,beta3=beta2*m_beta;
	  double delta = m_lag;
	  double num = (exp(m_alpha* delta - m_beta* (delta + 
     2 *m_tau))* pow((exp(m_alpha *m_tau) - 
     exp(m_beta* m_tau)),2)* m_alpha *(m_alpha - 
     2 *m_beta));
   
   double den =(2* (-alpha2 + 
   exp((m_alpha - m_beta) *m_tau)* m_alpha *(m_alpha - 2* m_beta) + 
   2 *m_alpha *m_beta + m_alpha* beta2* m_tau - beta3* m_tau));
   
    res(0,0)= num/den; 
    return res;
  }else{
    Rcpp::NumericVector lambda0_internal(lambda0);
    Rcpp::NumericMatrix alpha_internal(alpha);
    Rcpp::NumericVector beta_internal(beta);
    
    arma::vec m_lambda0(lambda0_internal.begin(),dimension,false);
    arma::mat m_alpha(alpha_internal.begin(),dimension,dimension,false);
    arma::vec m_beta(beta_internal.begin(),dimension,false);
    arma::mat m_beta_matrix = diagmat(m_beta);
    
    checkStability(m_beta_matrix,m_alpha);
    
    arma::mat variance = jumpVariance(lambda0,alpha,beta,tau);
   
     
    arma::mat temp0 = computeC2(m_lambda0,m_alpha,m_beta_matrix,m_tau)*computeC0(m_lambda0,m_alpha,m_beta_matrix,m_lag);
    temp0 = temp0*computeC2(m_lambda0,m_alpha,m_beta_matrix,m_tau);
    arma::vec tempVec = expectedStationaryLambda(m_lambda0,m_alpha,m_beta_matrix,0);
    arma::mat temp = grandLambdaInfini(m_lambda0,m_alpha,m_beta_matrix,0)+
    (m_alpha*vectorToDiagonalMatrix(tempVec));
    res = temp0*temp;
    for (int i = 0; i < dimension; i++){
        for (int j = 0; j < dimension; j++){
            res(i,j)=res(i,j)/sqrt(variance(i,i)*variance(j,j));
        }
    }
  
  
    return res;
  }
}

// [[Rcpp::export]]
double  likelihoodHawkes(SEXP lambda0,SEXP alpha,SEXP beta,SEXP history)
{
  int dimension = getDimension(lambda0);
  
  
  double res = 0;
  
  if (dimension == 1)
  {
    double m_lambda0 = as<double>(lambda0);
    double m_alpha = as<double>(alpha);
    double m_beta = as<double>(beta);
    Rcpp::NumericVector m_history(history);
    double m_T = m_history[m_history.size()-1];
    
    if(m_beta<m_alpha){
      stop("Unstable. You must have alpha < beta");
    }
    double *A = new double[m_history.size()];
    A[0] = 0;
  	for (int i = 1; i < m_history.size(); i++)
  	{
  		A[i] = (1.0+A[i-1])*exp(-m_beta * (m_history[i] - m_history[i-1])) ;
  	}

  	double sum = 0.0;
  	for (int i = 0; i < m_history.size(); i++)
  	{
  		sum = sum+   (1 - exp(-m_beta * (m_T - m_history[i]))) ;
  	}
  	sum = (m_alpha / m_beta) * sum;
  	res = - m_lambda0 * m_T - sum;
  	
  	
  	for (int i = 0; i < m_history.size(); i++)
  	{
  		res = res + log(m_lambda0+m_alpha*A[i]);
  	}
  	delete [] A;
	
	  return (-res);
  }else{
    Rcpp::NumericVector lambda0_internal(lambda0);
    Rcpp::NumericMatrix alpha_internal(alpha);
    Rcpp::NumericVector beta_internal(beta);
    arma::vec m_lambda0(lambda0_internal.begin(),dimension,false);
    arma::mat m_alpha(alpha_internal.begin(),dimension,dimension,false);
    arma::vec m_beta(beta_internal.begin(),dimension,false);
    Rcpp::List m_history(history);
    
    double m_T = 0;
    for (int n = 0; n < dimension; n++)
    {
  	  m_T = std::max(as<Rcpp::NumericVector>(m_history[n])[as<Rcpp::NumericVector>(m_history[n]).size()-1],m_T);
  	}
    
    for (int m = 0; m < dimension; m++)
  	{
  		double sum = 0.0;
     
      double *Rdiag = new double[as<Rcpp::NumericVector>(m_history[m]).size()];
    	double *RNonDiag = new double[as<Rcpp::NumericVector>(m_history[m]).size()];
	    int * index=new int[dimension];
    	for (int n = 0; n < dimension; n++)
    	{
    		index[n] = 0;
    	}
    	Rdiag[0] = 0;
    	RNonDiag[0] = 0;
    	for (int i = 1; i <as<Rcpp::NumericVector>(m_history[m]).size(); i++)
    	{
    		Rdiag[i] = (1.0+Rdiag[i-1])*exp(-m_beta(m) * (as<Rcpp::NumericVector>(m_history[m])[i] -as<Rcpp::NumericVector>(m_history[m])[i-1])) ;
    	}
    	for (int i = 1; i < as<Rcpp::NumericVector>(m_history[m]).size(); i++)
    	{
    		RNonDiag[i] = (RNonDiag[i-1])*exp(-m_beta(m) * (as<Rcpp::NumericVector>(m_history[m])[i] - as<Rcpp::NumericVector>(m_history[m])[i-1])) ;
    		for (int n = 0; n < dimension; n++)
    		{
    			if (m==n)
    			{
    				continue;
    			}
    			for (int k = index[n]; k < as<Rcpp::NumericVector>(m_history[n]).size(); k++)
    			{
    				if (as<Rcpp::NumericVector>(m_history[n])[k] >= as<Rcpp::NumericVector>(m_history[m])[i-1])
    				{
    					if (as<Rcpp::NumericVector>(m_history[n])[k] < as<Rcpp::NumericVector>(m_history[m])[i])
    					{
    						RNonDiag[i] += exp(-m_beta(m) * (as<Rcpp::NumericVector>(m_history[m])[i] - as<Rcpp::NumericVector>(m_history[n])[k]));
    					}
    					else
    					{
    						index[n] = k;
    						break;
    					}
    				}
    			}
    		}
    	}
    	for (int n = 0; n < dimension; n++)
    	{
    		for (int k = 0; k < as<Rcpp::NumericVector>(m_history[n]).size(); k++)
    		{
    			sum = sum + (m_alpha(m,n) / m_beta(m)) *
    				(1-exp(-m_beta(m) * (m_T - as<Rcpp::NumericVector>(m_history[n])[k])));//Beta diagonal
    
    		}
    	}

  	  res = res - m_lambda0(m) * m_T - sum;
  
  
    	for (int i = 0; i < as<Rcpp::NumericVector>(m_history[m]).size(); i++)
    	{
    		sum = m_lambda0(m);
    		for (int n = 0; n < dimension; n++)		
    		{
    
    			if(m==n)
    			{
    				sum = sum+m_alpha(m,n)*Rdiag[i];
    			}
    			else
    			{
    				sum = sum +m_alpha(m,n)*RNonDiag[i];
    			}
    		}
    		res = res+log(sum);
    	}
    	delete[] Rdiag;
    	delete[] RNonDiag;
    	delete[] index;
    }
   
	  return (-res);
  }
}
