#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
//' Returns a vector containing the upper triangle of a matrix, without the diagonal.
//'
//' @param mat a square matrix
//' @return Returns a vector containing the upper triangle of a matrix, without the diagonal.
//' @seealso unvectorize
//' @export
//'
//' @examples
//' vectorize(matrix(1:9,3,3))
// [[Rcpp::export]]
  NumericVector vectorize(const NumericMatrix& mat){

    int p = mat.nrow();
    NumericVector vect(p*(p-1)/2);
    int index = 0;
    for(int irow=0; irow<p-1; ++irow){
	    for(int icol=irow+1; icol<p; ++icol){
	    	vect(index) = mat(irow,icol);
		    index = index+1;
	    }
    }
    return(vect);
  }


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
//' Returns an upper-triangle matrix, without the diagonal, containing the elements of a given vector.
//'
//' @param vect A vector containing the upper triangle of a matrix, without the diagonal
//' @return Returns an upper-triangle matrix where each entry is given by the vector containing the upper triangle of a matrix, without the diagonal.
//' @seealso vectorize
//' @export
//'
//' @examples
//' unvectorize(1:10)
// [[Rcpp::export]]
  NumericMatrix unvectorize(const NumericVector& vect){

    int p = (int) floor((1+sqrt((double) 1+8*vect.length())/2));
    NumericMatrix mat(p,p);
    int index = 0;
    for(int irow=0; irow<p-1; ++irow){
	    for(int icol=irow+1; icol<p; ++icol){
	    	mat(irow,icol) = vect(index);
		    index = index+1;
	    }
    }
    return(mat);
  }



//' Returns the theoretical covariance of empirical correlations.
//'
//' @param r a correlation matrix
//' @return Returns the theoretical covariance of empirical correlations.
//' @export
//' @references Aitkin, M. A. (1969). Some tests for correlation matrices. Biometrika, 443-446.
//' @seealso covDcorNorm
//'
//' @examples
//' p <- 10
//' corr_theo <- diag(1,p)
//' corr_theo[2:p,] <- 0.3
//' corr_theo[,2:p] <- 0.3
//' covDcor(corr_theo)
// [[Rcpp::export]]
  NumericMatrix covDcor(const NumericMatrix& r) {

    int n2 = r.nrow()*(r.ncol()-1)/2;
    NumericMatrix w(n2,n2);
	
    int ligne = -1;
    int colonne = -1;
    for(int i=0; i<r.nrow()-1; ++i){
      for(int j=i+1; j<r.nrow(); ++j){
	      ligne = ligne+1;
	      colonne = -1;
	      for(int k=0; k<r.ncol()-1; ++k){
	        for(int l=k+1; l<r.ncol(); ++l){
	    	    colonne = colonne+1;
  		      if( ((i==k)&&(j==l)) || ((i==l)&&(j==k)) ){
  		         w(ligne,colonne) = (1.0-r(i,j)*r(i,j))*(1.0-r(i,j)*r(i,j));
		        }
              if( ((i==k)&&(j!=l)) ){
		           w(ligne,colonne) = -0.5*r(i,j)*r(i,l)*(1.0-r(i,j)*r(i,j)-r(i,l)*r(i,l)-r(j,l)*r(j,l))+r(j,l)*(1.0-r(i,j)*r(i,j)-r(i,l)*r(i,l));
    		    }
		        if( ((i==l)&&(j!=k)) ){
		           w(ligne,colonne) =  -0.5*r(i,j)*r(i,k)*(1.0-r(i,j)*r(i,j)-r(i,k)*r(i,k)-r(j,k)*r(j,k))+ r(j,k)*(1.0-r(i,j)*r(i,j)-r(i,k)*r(i,k));
    		    }
		        if( ((i!=k)&&(j==l)) ){
		           w(ligne,colonne) = -0.5*r(j,i)*r(j,k)*(1.0-r(j,i)*r(j,i)-r(j,k)*r(j,k)-r(i,k)*r(i,k)) + r(i,k)*(1.0-r(j,i)*r(j,i)-r(j,k)*r(j,k));
		        }
		        if( ((i!=l)&&(j==k)) ){
 		           w(ligne,colonne) = -0.5*r(j,i)*r(j,l)*(1.0-r(j,i)*r(j,i)-r(j,l)*r(j,l)-r(i,l)*r(i,l)) + r(i,l)*(1.0-r(j,i)*r(j,i)-r(j,l)*r(j,l));
		        }
    		    if( (i!=k)&&(i!=l)&&(j!=k)&&(j!=l) ){
		          w(ligne,colonne) = 0.5*r(i,j)*r(k,l)*(r(i,k)*r(i,k)+r(i,l)*r(i,l)+r(j,k)*r(j,k)+r(j,l)*r(j,l)) + r(i,k)*r(j,l) + r(i,l)*r(j,k) - r(i,k)*r(j,k)*r(k,l) - r(i,j)*r(i,k)*r(i,l) - r(i,j)*r(j,k)*r(j,l) - r(i,l)*r(j,l)*r(k,l);
  		      }
	        }
        }
      }
    }


    return w;
  }


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
//' Returns the theoretical covariance of empirical correlations.
//'
//' @param r a correlation matrix
//' @return Returns the theoretical covariance of 2nd order statistics, 
//'         \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}.
//' @export
//' @seealso covDcor
//'
//' @examples
//' p <- 10
//' corr_theo <- diag(1,p)
//' corr_theo[2:p,] <- 0.3
//' corr_theo[,2:p] <- 0.3
//' covD2nd(corr_theo)
// [[Rcpp::export]]
Rcpp::NumericMatrix covD2nd(const Rcpp::NumericMatrix& r) {
  
  int n2 = r.nrow()*(r.ncol()-1)/2;
  Rcpp::NumericMatrix w(n2,n2);
  
  int ligne = -1;
  int colonne = -1;
  for(int i=0; i<r.nrow()-1; ++i){
    for(int j=i+1; j<r.nrow(); ++j){
      ligne = ligne+1;
      colonne = -1;
      for(int k=0; k<r.ncol()-1; ++k){
        for(int l=k+1; l<r.ncol(); ++l){
          colonne = colonne+1;
          w(ligne,colonne) = (r(i,j)*r(k,l)+r(i,k)*r(j,l)+r(i,l)*r(j,k))/sqrt(1+r(i,j)*r(i,j))/sqrt(1+r(k,l)*r(k,l));
        }
      }
    }
  }
  
  return w;
}
