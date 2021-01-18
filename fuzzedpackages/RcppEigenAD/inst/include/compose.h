
#ifndef ___COMPOSE_H___
#define ___COMPOSE_H___


#include "cppad.eigen.h"
#include "adfunc.h"
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>


template <class T>
boost::tuple<MATRIX(T),MATRIX(T),TENSOR(T) > compose(const boost::tuple<MATRIX(T),MATRIX(T),TENSOR(T) >& f, 
						     const boost::tuple<MATRIX(T),MATRIX(T),TENSOR(T) >& g)

{

  // NB - THIS IS WORK IN PROGRESS - ONLY APPLIES TH THE CASE   g:R^n->R^m and f:R^m->R

  typedef  MATRIX(T) Matrix_Type_T;
  typedef  TENSOR(T) Tensor_Type_T;

  // get jacobians and hessians for f and g
  MATRIX(T) jac_g = boost::get<1>(g); // g.get<1>();
  MATRIX(T) jac_f = boost::get<1>(f); //f.get<1>(); 
  MATRIX(T) hess_f = boost::get<2>(f)(0,0); // f.get<2>()(0,0); // know that there is only one hessian for f - see note above)

  // work out the dimensions (using the jacobian of g)
  unsigned int n = jac_g.cols();
  unsigned int m = jac_g.rows();

  // chain rule for jacobians
  MATRIX(T) jacobian = jac_f*jac_g;

  // implement Faa di Bruno's formula for second derivatives (y \in R only for now - see note above)
  // configure result
  MATRIX(T) hessian(n,n);
  hessian.setZero();
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < n; j++)
	{
	  for(int k = 0; k < m; k++)
	    {
	      MATRIX(T) hess_g_k = boost::get<2>(g)(0,k); // g.get<2>()(0,k);
	      hessian(i,j) = hessian(i,j) + jac_f(0,k)*hess_g_k(i,j); 
	      for(int l = 0; l < m; l++)
		{
		  hessian(i,j) = hessian(i,j) + hess_f(k,l)*jac_g(k,i)*jac_g(l,j);
		}
	    }
	} 
    }

  TENSOR(T) hessians(1,1);
  hessians(0,0) = hessian;
  boost::tuple<Matrix_Type_T,Matrix_Type_T,Tensor_Type_T> result = boost::make_tuple(boost::get<0>(f),jacobian,hessians);
  // boost::tuple<Mat,Mat,Tensor> result = boost::make_tuple(f.get<0>(),jacobian,hessians);
  return result;
}

template <class T>
boost::tuple<MATRIX(T),MATRIX(T),TENSOR(T)> compose(const boost::tuple<MATRIX(T),MATRIX(T),TENSOR(T)>& f, 
						    const TRIPLE(T)& g)
{

  typedef  MATRIX(T) Matrix_Type_T;
  typedef  TENSOR(T) Tensor_Type_T;

  // find out how many ys there are
  unsigned int num_ys = boost::get<0>(g).size();
  // find out the dimension
  unsigned int dimensions = 0;
  for(unsigned int i = 0; i < num_ys; i++)
    {
      dimensions += boost::get<0>(g)[i].rows()*boost::get<0>(g)[i].cols();
    }
  MATRIX(T) y_star(1,dimensions);
  unsigned int cursor = 0;
  for(unsigned int i = 0; i < num_ys; i++)
    {
      unsigned int n_rows = boost::get<0>(g)[i].rows();
      unsigned int n_cols = boost::get<0>(g)[i].cols();
      for(unsigned int r = 0; r < n_rows; r++)
	{
	  for(unsigned int c = 0; c < n_cols; c++)
	    {
	      y_star(0,cursor++) = boost::get<0>(g)[i](r,c);
	    }
	}
    }
  boost::tuple<MATRIX(T),MATRIX(T),TENSOR(T)> g_star = boost::make_tuple(y_star,boost::get<1>(g),boost::get<2>(g));
  return compose(f,g_star);
}

template <class T>
TRIPLE(T) compose(const TRIPLE(T)& f, 
	       const TRIPLE(T)& g)
{
  typedef  MATRIX(T) Matrix_Type_T;
  typedef  TENSOR(T) Tensor_Type_T;
  typedef  MULTIARG(T) MultiArg_Type_T;
  boost::tuple<Matrix_Type_T,Matrix_Type_T,Tensor_Type_T> f_star;
  boost::get<1>(f_star) = boost::get<1>(f);
  boost::get<2>(f_star) = boost::get<2>(f);
  boost::tuple<Matrix_Type_T,Matrix_Type_T,Tensor_Type_T> composite = compose(f_star,g);
  return boost::make_tuple<MultiArg_Type_T,Matrix_Type_T,Tensor_Type_T>(boost::get<0>(f),boost::get<1>(composite),boost::get<2>(composite));
}

template <class T1, class T2>
TRIPLE(T1) compose(const FUNCTION(T2)& f, 
		  const FUNCTION(T2)& g,
		  const MULTIARG(T1)& xs)
{
  typedef TRIPLE(T1) Triple_Type_T1;
  Triple_Type_T1 ad_info_g = adfunc(g,xs);
  Triple_Type_T1 ad_info_f = adfunc(f,boost::get<0>(ad_info_g));
  return compose(ad_info_f,ad_info_g);
}


template <class T1, class T2>
TRIPLE(T1) compose(const FUNCTION(T2)& f, 
		   const TRIPLE(T1)& g)
{
  typedef TRIPLE(T1) Triple_Type_T1;
  Triple_Type_T1 ad_info_f = adfunc(f,boost::get<0>(g));
  return compose(ad_info_f,g);
}


#endif
