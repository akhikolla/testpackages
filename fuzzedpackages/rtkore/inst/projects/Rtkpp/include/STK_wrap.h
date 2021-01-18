/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp
 * created on: 22 août 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_wrap.h
 *  @brief In this file we define a stk++ wrapper converter for expressions.
 **/


#ifndef STK_WRAP_H
#define STK_WRAP_H

namespace STK
{
namespace hidden
{
// forward declaration
template<class Derived, int Structure_ = Traits<Derived>::structure_, int Orient_ =  Traits<Derived>::orient_>
struct WrapHelper;


/** @ingroup hidden
 * Specialization of WrapHelpher for vector_
 */
template<class Derived, int Orient_>
struct WrapHelper<Derived, Arrays::vector_, Orient_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Vector<Rtype_> Result;

  static SEXP wrapImpl(Derived const& vec)
  {
    Result res(vec.size());
    for(int i=vec.begin(), iRes=0; i< vec.end(); i++, iRes++)
    { res(iRes) = vec.elt(i);}
    return Rcpp::wrap(res);
  }
};

// specialization for vector_
template<class Derived, int Orient_>
struct WrapHelper<Derived, Arrays::point_, Orient_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Vector<Rtype_> Result;

  static SEXP wrapImpl(Derived const& vec)
  {
    Result res(vec.size());
    for(int i=vec.begin(), iRes=0; i< vec.end(); i++, iRes++)
    { res(iRes) = vec.elt(i);}
    return Rcpp::wrap(res);
  }
};

// specialization for diagonal_
template<class Derived, int Orient_>
struct WrapHelper<Derived, Arrays::diagonal_, Orient_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Matrix<Rtype_> Result;

  static SEXP wrapImpl(Derived const& d)
  {
    Result res(d.size(), d.size());
    for(int i=d.begin(), iRes=0; i< d.end(); i++, iRes++)
    {
      for(int j=d.begin(), jRes=0; j< d.end(); j++, jRes++)
      { res(iRes, jRes) = Type(0);}
      res(iRes, iRes) = d.elt(i);
    }
    return Rcpp::wrap(res);
  }
};

// specialization for array2D_
template<class Derived>
struct WrapHelper<Derived, Arrays::array2D_, Arrays::by_col_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Matrix<Rtype_> Result;

  static SEXP wrapImpl(Derived const& matrix)
  {
    Result res(matrix.sizeRows(), matrix.sizeCols());
    for(int j=matrix.beginCols(), jRes=0; j< matrix.endCols(); j++, jRes++)
    {
      for(int i=matrix.beginRows(), iRes=0; i< matrix.endRows(); i++, iRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
};
template<class Derived>
struct WrapHelper<Derived, Arrays::array2D_, Arrays::by_row_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Matrix<Rtype_> Result;

  static SEXP wrapImpl(Derived const& matrix)
  {
    Result res(matrix.sizeRows(), matrix.sizeCols());
    for(int i=matrix.beginRows(), iRes=0; i< matrix.endRows(); i++, iRes++)
    {
      for(int j=matrix.beginCols(), jRes=0; j< matrix.endCols(); j++, jRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
};

// specialization for array2D_
template<class Derived>
struct WrapHelper<Derived, Arrays::square_, Arrays::by_col_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Matrix<Rtype_> Result;

  static SEXP wrapImpl(Derived const& matrix)
  {
    Result res(matrix.size(), matrix.size());
    for(int j=matrix.begin(), jRes=0; j< matrix.end(); j++, jRes++)
    {
      for(int i=matrix.begin(), iRes=0; i< matrix.end(); i++, iRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
};
template<class Derived>
struct WrapHelper<Derived, Arrays::square_, Arrays::by_row_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Matrix<Rtype_> Result;

  static SEXP wrapImpl(Derived const& matrix)
  {
    Result res(matrix.size(), matrix.size());
    for(int i=matrix.begin(), iRes=0; i< matrix.end(); i++, iRes++)
    {
      for(int j=matrix.begin(), jRes=0; j< matrix.end(); j++, jRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
};
// specialization for upper_triangular_
template<class Derived>
struct WrapHelper<Derived, Arrays::upper_triangular_, Arrays::by_col_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Matrix<Rtype_> Result;

  static SEXP wrapImpl(Derived const& matrix)
  {
    Result res(matrix.sizeRows(), matrix.sizeCols());
    for(int j=matrix.beginCols(), jRes=0; j< matrix.endCols(); j++, jRes++)
    {
      for(int i=matrix.beginRows(), iRes=0; i<= j; i++, iRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
};
template<class Derived>
struct WrapHelper<Derived, Arrays::upper_triangular_, Arrays::by_row_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Matrix<Rtype_> Result;

  static SEXP wrapImpl(Derived const& matrix)
  {
    Result res(matrix.sizeRows(), matrix.sizeCols());
    for(int i=matrix.beginRows(), iRes=0; i< matrix.endRows(); i++, iRes++)
    {
      for(int j=i, jRes=iRes; j< matrix.endCols(); j++, jRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
};

// specialization for array2D_
template<class Derived>
struct WrapHelper<Derived, Arrays::lower_triangular_, Arrays::by_col_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Matrix<Rtype_> Result;

  static SEXP wrapImpl(Derived const& matrix)
  {
    Result res(matrix.sizeRows(), matrix.sizeCols());
    for(int j=matrix.beginCols(), jRes=0; j< matrix.endCols(); j++, jRes++)
    {
      for(int i=j, iRes=jRes; i< matrix.endRows(); i++, iRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
};
template<class Derived>
struct WrapHelper<Derived, Arrays::lower_triangular_, Arrays::by_row_>
{
  typedef typename Traits<Derived>::Type Type;
  enum
  {
    Rtype_ = RcppTraits<Type>::Rtype_
  };
  typedef Rcpp::Matrix<Rtype_> Result;

  static SEXP wrapImpl(Derived const& matrix)
  {
    Result res(matrix.sizeRows(), matrix.sizeCols());
    for(int i=matrix.beginRows(), iRes=0; i< matrix.endRows(); i++, iRes++)
    {
      for(int j=matrix.beginCols(), jRes=0; j<=i; j++, jRes++)
      { res(iRes, jRes) = matrix.elt(i,j);}
    }
    return Rcpp::wrap(res);
  }
};


} // namespace hidden

/** template method allowing to wrap expressions in a RMatrix and a SEXP
 *  @param expr the expression to wrap
 *  @return a copy of expr in a SEXP
 **/
template<typename Derived>
SEXP wrap( ExprBase<Derived> const& expr)
{ return hidden::WrapHelper<Derived>::wrapImpl(expr.asDerived());}

} // namespace STK


#endif /* STK_WRAP_H */
