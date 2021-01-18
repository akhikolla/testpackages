
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
using namespace Rcpp;

#include <RcppEigen.h>

/*
// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}
*/


// [[Rcpp::export]]
Rcpp::List marshal_faa_di_bruno(const Eigen::MatrixXd& Jf,
				const Eigen::MatrixXd& Hf,
				const Eigen::MatrixXd& Jg,
				const Eigen::MatrixXd& Hg)
			      
			      
{

  // domain dimension is n
  // co domain dimension is m
  // co co domain dimension is l


  // determine dimensions from Jacobians
  unsigned int nf = Jf.cols();
  unsigned int mf = Jf.rows();
  unsigned int ng = Jg.cols();
  unsigned int mg = Jg.rows();  

  // are they compatible ?
  if(mg != nf) { /* problem !! */ }

  // create the resulting Jacobian
  Eigen::MatrixXd Jfog = Jf*Jg;

  // use faa di bruno's formula to obtain the (stacked) Hessian
  Eigen::MatrixXd Hfog(ng,mf*ng);
  Hfog.fill(0);

  // loop over all of the second order partial derivatives
  for (unsigned int m = 0; m < mf; m++)
    {
      for(unsigned int i = 0; i < ng; i++)
	{
	  for(unsigned int j = 0; j < ng; j++)
	    {
	      for(unsigned int k = 0; k < nf; k++)
		{
		  Hfog(i,j+m*ng) = Hfog(i,j+m*ng) + Jf(m,k)*Hg(i,j+k*ng);
		  for(unsigned int l = 0; l < nf; l++)
		    {
		      Hfog(i,j+m*ng) = Hfog(i,j+m*ng) + Hf(k,l+m*ng)*Jg(k,i)*Jg(l,j);		      
		    }
		}	   
	    }
	}
    }

  /*
    for (unsigned int m = 0; m < mf; m++)
    {
      for(unsigned int i = 0; i < ng; i++)
	{
	  for(unsigned int j = 0; j < ng; j++)
	    {
	      for(unsigned int k = 0; k < nf; k++)
		{
		  Hfog(i,j+m*ng) = Hfog(i,j+m*ng) + Jf(m,k)*Hg(i,j+k*ng);
		}	   
	    }
	}
    } 
  */
  return Rcpp::List::create(Rcpp::Named("Jfog") = Jfog,Rcpp::Named("Hfog") = Hfog);

  
}

