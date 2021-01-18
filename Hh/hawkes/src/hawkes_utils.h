#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat vectorToDiagonalMatrix(arma::vec& v)
{
  arma::mat res=arma::zeros(v.n_elem,v.n_elem);
	for (unsigned int i=0;i<v.n_elem;i++)
	{
		res(i,i) = v(i);
	}
	return(res);
}
arma::mat expm_pad(arma::mat& H, double t = 1.0, const int p = 6)
{
  //adapted from https://www.dbtsai.com/blog/2008-11-25-matrix-exponential/
	const unsigned int n = H.n_rows;
	arma::mat I = arma::eye(n,n);
	if(t==0.0)
		return I;
	arma::mat U(n,n),H2(n,n),P(n,n),Q(n,n);
	double norm = 0.0;
// Calcuate Pade coefficients
	arma::vec c(p+1);
	c(0)=1;  
	for(int i = 0; i < p; ++i) 
		c(i+1) = c(i) * ((p - i)/((i + 1.0) * (2.0 * p - i)));
// Calcuate the infinty norm of H, which is defined as the largest row sum of a matrix
	for(unsigned int i=0; i<n; ++i) 
	{
		double temp = 0.0;
		for(unsigned int j = 0; j < n; j++)
			temp += std::abs(H(i, j)); 
		norm = t * std::max<double>(norm, temp);
	}
// If norm = 0, and all H elements are not NaN or infinity but zero, 
// then U should be identity.
	if (norm == 0.0) {
		bool all_H_are_zero = true;
		for(unsigned int i = 0; i < n; i++)
			for(unsigned int j = 0; j < n; j++)
				if( H(i,j) != double(0.0) ) 
					all_H_are_zero = false; 
		if( all_H_are_zero == true ) return I;
// Some error happens, H has elements which are NaN or infinity. 
		stop("Null input error in the template expm_pad.");
		
	}
// Scaling, seek s such that || H*2^(-s) || < 1/2, and set scale = 2^(-s)
 	int s = 0;
	double scale = 1.0;
	if(norm > 0.5) {
		s = std::max<int>(0, static_cast<int>((log(norm) / log(2.0) + 2.0)));
		scale /= double(std::pow(2.0, s));
		U = (scale * t) * H; // Here U is used as temp value due to that H is const
	}
	else
		U = H;

// Horner evaluation of the irreducible fraction, see the following ref above.
// Initialise P (numerator) and Q (denominator) 
	H2 = U*U;
	Q= ( c(p)*I );
	P=( c(p-1)*I );
	unsigned int odd = 1;
	for( unsigned int k = p - 1; k > 0; --k) {
		( odd == 1 ) ?
			( Q = (  (Q* H2) + c(k-1) * I ) ) :
			( P = (  (P* H2) + c(k-1) * I ) ) ;
		odd = 1 - odd;
	}
	( odd == 1 ) ? ( Q =  (Q* U) ) : ( P =  (P* U) );
	Q -= P;
// In origine expokit package, they use lapack ZGESV to obtain inverse matrix,
// and in that ZGESV routine, it uses LU decomposition for obtaing inverse matrix.
// Since in ublas, there is no matrix inversion template, I simply use the build-in
// LU decompostion package in ublas, and back substitute by myself.

// Implement Matrix Inversion
	H2=inv(Q);
	(odd == 1) ? 
		( U = ( -(I + double(2.0) * (H2* P))) ):
		( U = (   I + double(2.0) * (H2* P) ) );
// Squaring 
	for(int i = 0; i < s; ++i)
		U = ( (U*U));
	return U;
}
arma::mat computeC0(arma::vec& lambda,arma::mat& alpha,arma::mat& beta,double tau)
{
  arma::mat temp = alpha-beta;
  temp = tau*temp;
  arma::mat res = expm_pad(temp);
  
	return (res);
}
arma::mat computeC1(arma::vec& lambda,arma::mat& alpha,arma::mat& beta,double tau)
{
  arma::mat I = arma::eye(lambda.n_elem,lambda.n_elem);
	
	arma::mat alphaMinusBeta = alpha-beta;
	arma::mat matrixAlphaMinusBeta_inv = inv(alphaMinusBeta);
	arma::mat temp1 = beta*lambda;
	arma::mat temp = computeC0(lambda,alpha,beta,tau)-I;
  temp = temp*temp1;
	arma::mat res = matrixAlphaMinusBeta_inv*temp;
	
	return (res);
}
arma::mat computeC2(arma::vec& lambda,arma::mat& alpha,arma::mat& beta,double tau)
{
	arma::mat I = arma::eye(lambda.n_elem,lambda.n_elem);
	arma::mat alphaMinusBeta = alpha-beta;
	arma::mat matrixAlphaMinusBeta_inv = inv(alphaMinusBeta);
	arma::mat res = matrixAlphaMinusBeta_inv*(computeC0(lambda,alpha,beta,tau)- I);
	return (res);
}
arma::mat computeC3(arma::vec& lambda,arma::mat& alpha,arma::mat& beta,double tau)
{
	arma::mat I = arma::eye(lambda.n_elem,lambda.n_elem);
	
	arma::mat alphaMinusBeta = alpha-beta;
	arma::mat matrixAlphaMinusBeta_inv = inv(alphaMinusBeta);
	
	arma::mat betaL = beta*lambda;
	arma::mat temp1 = matrixAlphaMinusBeta_inv* matrixAlphaMinusBeta_inv;
	arma::mat temp2 = temp1*(computeC0(lambda,alpha,beta,tau)- I);
	arma::mat res1 = temp2*betaL;
	
	arma::mat res2 = matrixAlphaMinusBeta_inv*betaL;
	res2 = res2* tau;
	
	arma::mat res = res1-res2;
	return (res);
}
arma::mat computeC4(arma::vec& lambda,arma::mat& alpha,arma::mat& beta,double tau)
{
	arma::mat I = arma::eye(lambda.n_elem,lambda.n_elem);
	arma::mat alphaMinusBeta = alpha-beta;
	arma::mat matrixAlphaMinusBeta_inv = inv(alphaMinusBeta);
	arma::mat inv2 = matrixAlphaMinusBeta_inv*matrixAlphaMinusBeta_inv;
	arma::mat inv3 = matrixAlphaMinusBeta_inv*inv2;
	arma::mat temp = inv3*(computeC0(lambda,alpha,beta,tau)-I);
	arma::mat res = (((matrixAlphaMinusBeta_inv*(-tau * tau / 2.0))-(inv2*tau))+temp);
	return (res);
}
arma::mat computeC5(arma::vec& lambda,arma::mat& alpha,arma::mat& beta,double tau)
{
	arma::mat I = arma::eye(lambda.n_elem,lambda.n_elem);
	
	arma::mat alphaMinusBeta = alpha-beta;
	arma::mat matrixAlphaMinusBeta_inv = inv(alphaMinusBeta);
	arma::mat inv2 = matrixAlphaMinusBeta_inv*matrixAlphaMinusBeta_inv;
	arma::mat temp = inv2*(computeC0(lambda,alpha,beta,tau)-I);
	arma::mat res = (matrixAlphaMinusBeta_inv*(-tau))+ temp;
	
	return (res);
}
arma::vec expectedStationaryLambda(arma::vec& lambda,arma::mat& alpha,arma::mat& beta,double tau)
{
	arma::mat betaMinusalpha = beta-alpha;
	arma::mat matrixBetaMinusAlpha_inv = inv(betaMinusalpha);
	arma::mat temp = matrixBetaMinusAlpha_inv*beta;
	arma::vec res = temp*lambda;
	return (res);
	
}
arma::mat solveLyapounov(arma::mat&A,arma::mat& Q)
{
  //AX+XA'+Q=0
  arma::mat I = arma::eye(A.n_rows,A.n_rows);
  arma::mat AA = arma::kron(I, A)+ arma::kron(A, I);
  arma::mat QQ = reshape(Q,Q.n_rows*Q.n_cols,1);
  arma::mat mQQ = QQ*(-1);
  arma::mat invAA = inv(AA);
  arma::mat resTemp = invAA*mQQ;
  arma::mat result = arma::zeros(A.n_rows,A.n_cols);
  
  for (unsigned int i=0;i<A.n_rows;i++)
  {
      for (unsigned int j =0;j<A.n_cols;j++)
      {
          result(i, j) = resTemp(j * A.n_rows + i, 0);
      }
  }
  
  return (result);
}
arma::mat grandLambdaInfini(arma::vec& lambda,arma::mat& alpha,arma::mat& beta,double tau)
{
  arma::vec tempVec = expectedStationaryLambda(lambda,alpha,beta,tau);
	arma::mat temp = vectorToDiagonalMatrix(tempVec);
	arma::mat temp2 = alpha*temp;
	arma::mat Q = temp2*arma::trans(alpha);
	arma::mat A = alpha- beta;
	arma::mat res = solveLyapounov(A, Q);
  
	return (res);
}
void checkStability(arma::mat& m_beta_matrix,arma::mat& m_alpha){
  arma::mat beta_minus_alpha = m_beta_matrix- m_alpha;
  arma::cx_vec eigval=arma::eig_gen(beta_minus_alpha);
  arma::vec eigval_real=real(eigval);
  
  if(eigval_real.min()<0){
    stop("Unstable. beta - alpha must have eigenvalues with strictly positive real part.");
  }
}
