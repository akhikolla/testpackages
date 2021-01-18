

#ifndef ___ADFUNC_H___
#define ___ADFUNC_H___

#include "cppad.eigen.h"
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>



template <class T1,class T2>
TRIPLE(T1) adfunc
(const FUNCTION(T2)& f,const MULTIARG(T1)& xs)
{
  typedef MULTIARG(T1) MultiArg_Type_1;
  typedef MATRIX(T1) Matrix_Type_1; 
  typedef MATRIX(T2) Matrix_Type_2; 
  typedef TENSOR(T1) Tensor_Type_1;
  typedef TRIPLE(T1) Triple_Type_1;

  unsigned int num_xs = xs.size();
  typename MultiArg_Type_1::const_iterator it_xs = xs.begin();
  typename MultiArg_Type_1::const_iterator it_xs_end = xs.end();

  // determine size of AD vector for use with ADfun and INdependent
  unsigned int xs_AD_vec_size = 0;
  while(it_xs != it_xs_end)
    {
      xs_AD_vec_size += (*it_xs).rows()*(*it_xs).cols() ;
      it_xs++;
    }
  // create an AD vector for use with ADfun and Independent
  std::vector<T2> a_vec_x(xs_AD_vec_size);
  // create a vector representing the point at which the jacobian and hessian will be calculated
  std::vector<T1> eval_point(xs_AD_vec_size);
  // populate the vectors
  it_xs = xs.begin();
  unsigned int cursor = 0;
  while(it_xs != it_xs_end)
    {
      unsigned int x_nrows = (*it_xs).rows();
      unsigned int x_ncols = (*it_xs).cols();
      for(int row = 0; row < x_nrows; row++)
	{
	  for(int col = 0; col < x_ncols; col++)
	    { 
	      a_vec_x[cursor] = eval_point[cursor] = (*it_xs)(row,col);
	      cursor++;
	    }
	}
      it_xs++;
    }
 
  // start the tape recording
  CppAD::Independent(a_vec_x);  
  // wrap the vector into an aMat for dispatching to f
  std::vector<Matrix_Type_2 > a_xs(num_xs);
  it_xs = xs.begin();
  cursor = 0;
  unsigned int xs_entry = 0;
  while(it_xs != it_xs_end)
    {
      unsigned int x_nrows = (*it_xs).rows();
      unsigned int x_ncols = (*it_xs).cols();
      Matrix_Type_2 a_mat_x(x_nrows,x_ncols);
      for(int row = 0; row < x_nrows; row++)
	{
	  for(int col = 0; col < x_ncols; col++)
	    { 
	      a_mat_x(row,col) = a_vec_x[cursor++];
	    }
	}
      a_xs[xs_entry++] = a_mat_x;
      *it_xs++;
    }
  
  // dispatch to f
  std::vector<Matrix_Type_2 > a_ys = f(a_xs);

  // unwrap the result into a vector for the tape and a matrix for the result
  unsigned int num_a_ys = a_ys.size();
  typename std::vector<Matrix_Type_2 >::const_iterator it_a_ys = a_ys.begin();
  typename std::vector<Matrix_Type_2 >::const_iterator it_a_ys_end = a_ys.end();
  // determine size of AD vector for use with ADfun and INdependent
  unsigned int a_ys_AD_vec_size = 0;
  while(it_a_ys != it_a_ys_end)
    {
      a_ys_AD_vec_size += (*it_a_ys).rows()*(*it_a_ys).cols() ;
      it_a_ys++;
    }
  // create an AD vector for use with ADfun and Independent
  std::vector<T2> a_vec_y(a_ys_AD_vec_size);

 
  it_a_ys = a_ys.begin();
  cursor = 0;
  while(it_a_ys != it_a_ys_end)
    {
      unsigned int y_nrows = (*it_a_ys).rows();
      unsigned int y_ncols = (*it_a_ys).cols();
      for(int row = 0; row < y_nrows; row++)
	{
	  for(int col = 0; col < y_ncols; col++)
	    { 
	      a_vec_y[cursor++] = (*it_a_ys)(row,col);
	    }
	}
      it_a_ys++;
    }

  // process the tape
  // CppAD::ADFun<T1> f_tape(a_vec_x, a_vec_y);
  CppAD::ADFun<T1> f_tape; 
  f_tape.Dependent(a_vec_x, a_vec_y);

  // calculate the Jacobian
  std::vector<T1> v_jacobian = f_tape.Jacobian(eval_point);  // TODO -will need to do this across additional directions later

  // unwrap into the jacobian
  Matrix_Type_1 jacobian(a_ys_AD_vec_size,xs_AD_vec_size);
  cursor = 0;
  for(int row = 0; row < a_ys_AD_vec_size; row++)
    {
      for(int col = 0; col < xs_AD_vec_size; col++)
        { 
	  jacobian(row,col) = v_jacobian[cursor++];
    
        }
    }
  
  Tensor_Type_1 hessians(1,a_ys_AD_vec_size);
  for(int d = 0; d < a_ys_AD_vec_size; d++)
    {
      std::vector<T1> v_hessian = f_tape.Hessian(eval_point,d);  
      Matrix_Type_1 hessian(xs_AD_vec_size,xs_AD_vec_size);
      cursor = 0;
      for(int row = 0; row < xs_AD_vec_size; row++)
	{
	  for(int col = 0; col < xs_AD_vec_size; col++)
	    {
	      hessian(row,col) = v_hessian[cursor++]; 
	    } 
	  hessians(0,d) = hessian;
	}
    }
  
  // result matrices
  unsigned int num_ys = a_ys.size();
  MultiArg_Type_1  ys(num_ys);
  for(unsigned int i = 0; i < num_ys; i++)
    {
      ys[i] = Convert<Matrix_Type_1>(a_ys[i]);
    }
  Triple_Type_1 result = boost::make_tuple(ys,jacobian,hessians); 
  return result;

}




template <class T1,class T2>
  TRIPLE(T1) adfunc_prior(const FUNCTION(T2)& f,const MULTIARG(T1)& xs, std::vector<T1>& eval_point_copy)
{
  typedef MULTIARG(T1) MultiArg_Type_1;
  typedef MATRIX(T1) Matrix_Type_1; 
  typedef MATRIX(T2) Matrix_Type_2; 
  typedef TENSOR(T1) Tensor_Type_1;
  typedef TRIPLE(T1) Triple_Type_1;

  unsigned int num_xs = xs.size();
  typename MultiArg_Type_1::const_iterator it_xs = xs.begin();
  typename MultiArg_Type_1::const_iterator it_xs_end = xs.end();

  // determine size of AD vector for use with ADfun and INdependent
  unsigned int xs_AD_vec_size = 0;
  while(it_xs != it_xs_end)
    {
      xs_AD_vec_size += (*it_xs).rows()*(*it_xs).cols() ;
      it_xs++;
    }
  // create an AD vector for use with ADfun and Independent
  std::vector<T2> a_vec_x(xs_AD_vec_size);
  // create a vector representing the point at which the jacobian and hessian will be calculated
  std::vector<T1> eval_point(xs_AD_vec_size);
  // populate the vectors
  it_xs = xs.begin();
  unsigned int cursor = 0;
  while(it_xs != it_xs_end)
    {
      unsigned int x_nrows = (*it_xs).rows();
      unsigned int x_ncols = (*it_xs).cols();
      for(int row = 0; row < x_nrows; row++)
	{
	  for(int col = 0; col < x_ncols; col++)
	    { 
	      a_vec_x[cursor] = eval_point[cursor] = (*it_xs)(row,col);
	      cursor++;
	    }
	}
      it_xs++;
    }
 
  // start the tape recording
  CppAD::Independent(a_vec_x);
  
  // wrap the vector into an aMat for dispatching to f
  std::vector<Matrix_Type_2 > a_xs(num_xs);
  it_xs = xs.begin();
  cursor = 0;
  unsigned int xs_entry = 0;
  while(it_xs != it_xs_end)
    {
      unsigned int x_nrows = (*it_xs).rows();
      unsigned int x_ncols = (*it_xs).cols();
      Matrix_Type_2 a_mat_x(x_nrows,x_ncols);
      for(int row = 0; row < x_nrows; row++)
	{
	  for(int col = 0; col < x_ncols; col++)
	    { 
	      a_mat_x(row,col) = a_vec_x[cursor++];
	    }
	}
      a_xs[xs_entry++] = a_mat_x;
      *it_xs++;
    }
  
  // dispatch to f
  std::vector<Matrix_Type_2 > a_ys = f(a_xs);

  // unwrap the result into a vector for the tape and a matrix for the result
  unsigned int num_a_ys = a_ys.size();
  typename std::vector<Matrix_Type_2 >::const_iterator it_a_ys = a_ys.begin();
  typename std::vector<Matrix_Type_2 >::const_iterator it_a_ys_end = a_ys.end();
  // determine size of AD vector for use with ADfun and INdependent
  unsigned int a_ys_AD_vec_size = 0;
  while(it_a_ys != it_a_ys_end)
    {
      a_ys_AD_vec_size += (*it_a_ys).rows()*(*it_a_ys).cols() ;
      it_a_ys++;
    }
  // create an AD vector for use with ADfun and Independent
  std::vector<T2> a_vec_y(a_ys_AD_vec_size);

 
  it_a_ys = a_ys.begin();
  cursor = 0;
  while(it_a_ys != it_a_ys_end)
    {
      unsigned int y_nrows = (*it_a_ys).rows();
      unsigned int y_ncols = (*it_a_ys).cols();
      for(int row = 0; row < y_nrows; row++)
	{
	  for(int col = 0; col < y_ncols; col++)
	    { 
	      a_vec_y[cursor++] = (*it_a_ys)(row,col);
	    }
	}
      it_a_ys++;
    }

  // process the tape
  CppAD::ADFun<T1> f_tape; 
  f_tape.Dependent(a_vec_x, a_vec_y);

  // this is PRIOR mark posterior as independent
  Independent(eval_point);
  eval_point_copy = eval_point; // retain a copy for use by posterior

  // calculate the Jacobian
  std::vector<T1> v_jacobian = f_tape.Jacobian(eval_point);  // TODO -will need to do this across additional directions later

  // unwrap into the jacobian
  Matrix_Type_1 jacobian(a_ys_AD_vec_size,xs_AD_vec_size);
  cursor = 0;
  for(int row = 0; row < a_ys_AD_vec_size; row++)
    {
      for(int col = 0; col < xs_AD_vec_size; col++)
        { 
	  jacobian(row,col) = v_jacobian[cursor++];
    
        }
    }
  
  Tensor_Type_1 hessians(1,a_ys_AD_vec_size);  
  for(int d = 0; d < a_ys_AD_vec_size; d++)
    {
      std::vector<T1> v_hessian = f_tape.Hessian(eval_point,d);  
      Matrix_Type_1 hessian(xs_AD_vec_size,xs_AD_vec_size);
      cursor = 0;
      for(int row = 0; row < xs_AD_vec_size; row++)
	{
	  for(int col = 0; col < xs_AD_vec_size; col++)
	    {
	      hessian(row,col) = v_hessian[cursor++]; 
	    } 
	  hessians(0,d) = hessian;
	}
    }
 
  // result matrices
  unsigned int num_ys = a_ys.size();
  MultiArg_Type_1  ys(num_ys);
  for(unsigned int i = 0; i < num_ys; i++)
    {
      ys[i] = Convert<Matrix_Type_1>(a_ys[i]);
    }
  Triple_Type_1 result = boost::make_tuple(ys,jacobian,hessians); 
  return result;

}


template <class T1,class T2>
  TRIPLE(T1) adfunc_posterior(const FUNCTION(T2)& f,const MULTIARG(T2)& xs,std::vector<T2>& eval_point_copy)
{
  typedef MULTIARG(T1) MultiArg_Type_1;
  typedef MULTIARG(T2) MultiArg_Type_2;
  typedef MATRIX(T1) Matrix_Type_1; 
  typedef MATRIX(T2) Matrix_Type_2; 
  typedef TENSOR(T1) Tensor_Type_1;
  typedef TRIPLE(T1) Triple_Type_1;

  unsigned int num_xs = xs.size();

  typename MultiArg_Type_2::const_iterator it_xs = xs.begin();
  typename MultiArg_Type_2::const_iterator it_xs_end = xs.end();

  // determine size of AD vector for use with ADfun and INdependent
  unsigned int xs_AD_vec_size = 0;
  while(it_xs != it_xs_end)
    {
      xs_AD_vec_size += (*it_xs).rows()*(*it_xs).cols() ;
      it_xs++;
    }
  // create an AD vector for use with ADfun and Independent
  std::vector<T2> a_vec_x(xs_AD_vec_size);
  // create a vector representing the point at which the jacobian and hessian will be calculated 
  // populated from the PRIOR evaluation point
  std::vector<T1> eval_point(eval_point_copy.size());
  for(unsigned int i = 0; i < eval_point_copy.size(); i++)
    {
      eval_point[i] = Value(eval_point_copy[i]);
    }


  // populate the vectors
  it_xs = xs.begin();
  unsigned int cursor = 0;
  while(it_xs != it_xs_end)
    {
      unsigned int x_nrows = (*it_xs).rows();
      unsigned int x_ncols = (*it_xs).cols();
      for(int row = 0; row < x_nrows; row++)
	{
	  for(int col = 0; col < x_ncols; col++)
	    { 
	      a_vec_x[cursor] = (*it_xs)(row,col);
	      cursor++;
	    }
	}
      it_xs++;
    }
 
  // POSTERIOR the tape is already recording
  
  // wrap the vector into an aMat for dispatching to f
  std::vector<Matrix_Type_2 > a_xs(num_xs);
  it_xs = xs.begin();
  cursor = 0;
  unsigned int xs_entry = 0;
  while(it_xs != it_xs_end)
    {
      unsigned int x_nrows = (*it_xs).rows();
      unsigned int x_ncols = (*it_xs).cols();
      Matrix_Type_2 a_mat_x(x_nrows,x_ncols);
      for(int row = 0; row < x_nrows; row++)
	{
	  for(int col = 0; col < x_ncols; col++)
	    { 
	      a_mat_x(row,col) = a_vec_x[cursor++];
	    }
	}
      a_xs[xs_entry++] = a_mat_x;
      *it_xs++;
    }
  
  // dispatch to f
  std::vector<Matrix_Type_2 > a_ys = f(a_xs);

  // unwrap the result into a vector for the tape and a matrix for the result
  unsigned int num_a_ys = a_ys.size();
  typename std::vector<Matrix_Type_2 >::const_iterator it_a_ys = a_ys.begin();
  typename std::vector<Matrix_Type_2 >::const_iterator it_a_ys_end = a_ys.end();
  // determine size of AD vector for use with ADfun and INdependent
  unsigned int a_ys_AD_vec_size = 0;
  while(it_a_ys != it_a_ys_end)
    {
      a_ys_AD_vec_size += (*it_a_ys).rows()*(*it_a_ys).cols() ;
      it_a_ys++;
    }
  // create an AD vector for use with ADfun and Independent
  std::vector<T2> a_vec_y(a_ys_AD_vec_size);

 
  it_a_ys = a_ys.begin();
  cursor = 0;
  while(it_a_ys != it_a_ys_end)
    {
      unsigned int y_nrows = (*it_a_ys).rows();
      unsigned int y_ncols = (*it_a_ys).cols();
      for(int row = 0; row < y_nrows; row++)
	{
	  for(int col = 0; col < y_ncols; col++)
	    { 
	      a_vec_y[cursor++] = (*it_a_ys)(row,col);
	    }
	}
      it_a_ys++;
    }

  // process the tape
  CppAD::ADFun<T1> f_tape;
  f_tape.Dependent(eval_point_copy, a_vec_y);


  // this is POSTERIOR - Independent already implemented in PRIOR

  // calculate the Jacobian
  std::vector<T1> v_jacobian = f_tape.Jacobian(eval_point);  // TODO -will need to do this across additional directions later

  // using the eval point for size of the domain - adjust the size of jacobian and hessians accordingly
  xs_AD_vec_size = eval_point_copy.size();
  // unwrap into the jacobian
  Matrix_Type_1 jacobian(a_ys_AD_vec_size,xs_AD_vec_size);
  cursor = 0;
  for(int row = 0; row < a_ys_AD_vec_size; row++)
    {
      for(int col = 0; col < xs_AD_vec_size; col++)
        { 
	  jacobian(row,col) = v_jacobian[cursor++];
    
        }
    }
  
  // unwrp into the hessians
  Tensor_Type_1 hessians(1,a_ys_AD_vec_size);  
  for(int d = 0; d < a_ys_AD_vec_size; d++)
    {
      std::vector<T1> v_hessian = f_tape.Hessian(eval_point,d);  
      Matrix_Type_1 hessian(xs_AD_vec_size,xs_AD_vec_size);
      cursor = 0;
      for(int row = 0; row < xs_AD_vec_size; row++)
	{
	  for(int col = 0; col < xs_AD_vec_size; col++)
	    {
	      hessian(row,col) = v_hessian[cursor++]; 
	    } 
	  hessians(0,d) = hessian;
	}
    }
 
  // result matrices
  unsigned int num_ys = a_ys.size();
  MultiArg_Type_1  ys(num_ys);
  for(unsigned int i = 0; i < num_ys; i++)
    {
      ys[i] = Convert<Matrix_Type_1>(a_ys[i]);
    }
  Triple_Type_1 result = boost::make_tuple(ys,jacobian,hessians); 
  return result;

}






#endif
