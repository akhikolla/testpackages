/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Arrays
 * created on: 17 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_BinaryImpl.h
 *  @brief In this file we implement the Binary Implementation helper classes.
 **/


#ifndef STK_BINARYIMPL_H
#define STK_BINARYIMPL_H

#define EGAL(arg1, arg2) ((arg1::structure_ == int(Arrays::arg2)))

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  @brief Helper class giving the Structure of a binary operator
 **/
template<typename FunctorOp, typename Lhs, typename Rhs, int LStructure_, int RStructure_>
struct BinaryEltImpl;

//----------------
// Lhs array2D_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_      = Arrays::array2D_
       , binary_op_Kind_ = Arrays::binary_op_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_ = Arrays::binary_op_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_2D_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_2D_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
      };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_2D_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_2D_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_2D_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};



//----------------
// Lhs square_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_ = Arrays::binary_op_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_ = Arrays::binary_op_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_LowSym_
       , sizeRows_  = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_  = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};


//----------------
// lhs diagonal_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_Diag_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_Diag_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::upper_triangular_
       , binary_op_Kind_= Arrays::binary_op_Diag_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f( LType(), r.elt(i,j)) : f( LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::lower_triangular_
       , binary_op_Kind_= Arrays::binary_op_Diag_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f( LType(), r.elt(i,j)) : f( LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_Diag_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::upper_symmetric_
       , binary_op_Kind_= Arrays::binary_op_Diag_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f( LType(), r.elt(i,j)) : f( LType(), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::lower_symmetric_
       , binary_op_Kind_= Arrays::binary_op_Diag_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f( LType(), r.elt(i,j)) : f( LType(), r.elt(j,i));}
};



//-------------------------
// lhs upper_triangular_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_UpTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
      };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::upper_triangular_
       , binary_op_Kind_= Arrays::binary_op_UpTri_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f(l.elt(i,j), RType()) : f(LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::upper_triangular_
       , binary_op_Kind_= Arrays::binary_op_UpTri_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_UpTri_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f(l.elt(i,j), RType()) : f(LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpTri_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpTri_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpTri_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(j,i)) : f( LType(), r.elt(i,j));}
};



// lhs is lower_triangular_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_LowTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::lower_triangular_
       , binary_op_Kind_= Arrays::binary_op_LowTri_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f(l.elt(i,j), RType()) : f( LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_LowTri_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f(l.elt(i,j), RType()) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::lower_triangular_
       , binary_op_Kind_= Arrays::binary_op_LowTri_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(j,i)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(j,i));}
};



// lhs is symmetric_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_Sym_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
      };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_Sym_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
      };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_Sym_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_Sym_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_Sym_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_Sym_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_Sym_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_Sym_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};


//---------------------
// lhs upper_symmetric_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpSym_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpSym_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::upper_symmetric_
       , binary_op_Kind_= Arrays::binary_op_UpSym_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
      };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f(l.elt(i,j), RType()) : f(l.elt(j,i), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::upper_triangular_
       , binary_op_Kind_= Arrays::binary_op_UpSym_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(j,i), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpSym_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f(l.elt(i,j), RType()) : f(l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_UpSym_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::upper_symmetric_
       , binary_op_Kind_= Arrays::binary_op_UpSym_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j,i), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_UpSym_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(j,i)) : f( l.elt(j,i), r.elt(i,j));}
};


//-----------------------
// lhs lower_symmetric_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j, i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j, i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::lower_symmetric_
       , binary_op_Kind_= Arrays::binary_op_LowTri_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f(l.elt(i,j), RType()) : f(l.elt(j,i), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowSym_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f(l.elt(i,j), RType()) : f(l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::lower_triangular_
       , binary_op_Kind_= Arrays::binary_op_LowSym_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(j,i), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_LowSym_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j, i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_LowSym_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(j,i)) : f(l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::lower_symmetric_
       , binary_op_Kind_= Arrays::binary_op_LowSym_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(j,i), r.elt(j,i));}
};


//--------------------------------------------------------------------------------------------
// 1D case

//-----------------------
// Lhs diagonal_, Rhs 1D
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::diagonal_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f(LType(), RType());}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::point_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::diagonal_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeCols_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsOtherSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::vector_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::diagonal_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeRows_)
       , useForRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsOtherSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};

//-----------------------
// Lhs vector_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::vector_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_     = Arrays::diagonal_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = Rhs::sizeRows_ != UnknownSize ? int(Rhs::sizeRows_) : int(Lhs::sizeCols_)
       , sizeCols_      = Rhs::sizeCols_ != UnknownSize ? int(Rhs::sizeCols_) : int(Lhs::sizeCols_)
       , useForRows_    = Rhs::sizeRows_ != UnknownSize ? Arrays::useRhsSize_ : Arrays::useLhsSize_
       , useForCols_    = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::vector_, Arrays::vector_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::vector_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_      = 1
       , useForRows_    = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useForCols_    = Arrays::useLhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::vector_, Arrays::point_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_     = Arrays::vector_ // vector + point is vector (somehow an arbitrary choice)
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeCols_)
       , sizeCols_      = 1
       , useForRows_    = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsOtherSize_
       , useForCols_    = Arrays::useLhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};

//-----------------------
// Lhs point_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::point_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_     = Arrays::diagonal_ // point_ + diagonal_ is diagonal_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = Rhs::sizeRows_ != UnknownSize ? int(Rhs::sizeRows_) : int(Lhs::sizeCols_)
       , sizeCols_      = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_    = Rhs::sizeRows_ != UnknownSize ? Arrays::useRhsSize_ : Arrays::useLhsOtherSize_
       , useForCols_    = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::point_, Arrays::vector_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_     = Arrays::point_ // point + vector_ is point_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = 1
       , sizeCols_      = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeRows_)
       , useForRows_    = Arrays::useLhsSize_
       , useForCols_    = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsOtherSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::point_, Arrays::point_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::point_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = 1
       , sizeCols_      = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useForRows_    = Arrays::useLhsSize_
       , useForCols_    = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};

//-----------------------
// Lhs number_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::number_, Arrays::number_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_     = Arrays::number_
       , binary_op_Kind_= Arrays::binary_op_0D_
       , sizeRows_      = 1
       , sizeCols_      = 1
       , useForRows_    = Arrays::useLhsSize_
       , useForCols_    = Arrays::useLhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
  inline static result_type elt0Impl(FunctorOp const& f, Lhs const& l, Rhs const& r)
  { return f(l.elt(), r.elt());}
};

/** @ingroup hidden
 *  @brief implement the access to the rows of the BinaryOperator
 *  Possible cases are:
 *  - use lhs.rows()
 *  - use rhs.rows()
 *  - use lhs.cols()
 *  - use rhs.cols()
 **/
template< typename Lhs, typename Rhs, int Size_, int useForRows_>
struct BinaryRowsImpl;
/** @ingroup hidden
 *  @brief implement the access to the columns of the BinaryOperator
 *  Possible cases are:
 *  - use lhs.cols()
 *  - use rhs.cols()
 *  - use lhs.rows()
 **/
template< typename Lhs, typename Rhs, int Size_, int useForCols_>
struct BinaryColsImpl;
/** @ingroup hidden
  * @brief specialization for the case useLhsSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryRowsImpl< Lhs, Rhs,  Size_, Arrays::useLhsSize_>
{
  /** Type of the Range for the rows */
  typedef TRange<Size_> RowRange;
  /**  @return the range of the rows */
  inline static RowRange const& rowsImpl(Lhs const& lhs, Rhs const& rhs) { return lhs.rows();}
};
/** @ingroup hidden
  * @brief specialization for the case useRhsSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryRowsImpl<Lhs, Rhs, Size_, Arrays::useRhsSize_>
{
  /** Type of the Range for the rows */
  typedef TRange<Size_> RowRange;
  /**  @return the range of the rows */
  inline static RowRange const& rowsImpl(Lhs const& lhs, Rhs const& rhs) { return rhs.rows();}
};
/** @ingroup hidden
  * @brief specialization for the case useLhsOtherSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryRowsImpl<Lhs, Rhs, Size_, Arrays::useLhsOtherSize_>
{
  /** Type of the Range for the rows */
  typedef TRange<Size_> RowRange;
  /**  @return the range of the rows */
  inline static RowRange const& rowsImpl(Lhs const& lhs, Rhs const& rhs) { return lhs.cols();}
};
/** @ingroup hidden
  * @brief specialization for the case useRhsOtherSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryRowsImpl<Lhs, Rhs, Size_, Arrays::useRhsOtherSize_>
{
  /** Type of the Range for the rows */
  typedef TRange<Size_> RowRange;
  /**  @return the range of the rows */
  inline static RowRange const& rowsImpl(Lhs const& lhs, Rhs const& rhs) { return rhs.cols();}
};
/** @ingroup hidden
  * @brief specialization for the case useLhsSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryColsImpl<Lhs, Rhs, Size_, Arrays::useLhsSize_>
{
  /** Type of the Range for the columns */
  typedef TRange<Size_> ColRange;
  /**  @return the range of the columns */
  inline static ColRange const& colsImpl(Lhs const& lhs, Rhs const& rhs) { return lhs.cols();}
};
/** @ingroup hidden
  * @brief specialization for the case useRhsSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryColsImpl< Lhs, Rhs, Size_, Arrays::useRhsSize_>
{
  /** Type of the Range for the columns */
  typedef TRange<Size_> ColRange;
  /**  @return the range of the columns */
  inline static ColRange const& colsImpl(Lhs const& lhs, Rhs const& rhs) { return rhs.cols();}
};
/** @ingroup hidden
  * @brief specialization for the case useLhsOtherSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryColsImpl< Lhs, Rhs, Size_, Arrays::useLhsOtherSize_>
{
  /** Type of the Range for the columns */
  typedef TRange<Size_> ColRange;
  /**  @return the range of the columns */
  inline static ColRange const& colsImpl(Lhs const& lhs, Rhs const& rhs) { return lhs.rows();}
};
/** @ingroup hidden
  * @brief specialization for the case useRhsOtherSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryColsImpl< Lhs, Rhs,Size_, Arrays::useRhsOtherSize_>
{
  /** Type of the Range for the columns */
  typedef TRange<Size_> ColRange;
  /**  @return the range of the columns */
  inline static ColRange const& colsImpl(Lhs const& lhs, Rhs const& rhs) { return rhs.rows();}
};


} // end namespace hidden


} // namespace STK

#undef EGAL

#endif /* STK_BINARYOPERATORS_H */
