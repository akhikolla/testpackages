// Copyright (c) 1999-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_LAZY_EXACT_NT_H
#define CGAL_LAZY_EXACT_NT_H

#define CGAL_int(T)    typename First_if_different<int,    T>::Type
#define CGAL_double(T) typename First_if_different<double, T>::Type
#define CGAL_To_interval(T) To_interval<T>


#include <CGAL/number_type_basic.h>
#include <CGAL/assertions.h>

#include <CGAL/boost/iterator/transform_iterator.hpp> // for Root_of functor
#include <boost/operators.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>

#include <boost/mpl/if.hpp>
#include <boost/mpl/logical.hpp>

#include <CGAL/Interval_nt.h>
#include <CGAL/Handle.h>
#include <CGAL/NT_converter.h>

#include <CGAL/Profile_counter.h>
#include <CGAL/Lazy.h>

#include <CGAL/Sqrt_extension_fwd.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/tss.h>
#include <CGAL/type_traits.h>

#include <CGAL/IO/io.h>


/*
 * This file contains the definition of the number type Lazy_exact_nt<ET>,
 * where ET is an exact number type (must provide the exact operations needed).
 *
 * Lazy_exact_nt<ET> provides a DAG-based lazy evaluation, like LEDA's real,
 * Core's Expr, LEA's lazy rationals...
 *
 * The values are first approximated using Interval_base.
 * The exactness is provided when needed by ET.
 *
 * Lazy_exact_nt<ET> is just a handle to the abstract base class
 * Lazy_exact_nt_rep which has pure virtual methods .approx() and .exact().
 * From this class derives one class per operation, with one constructor.
 *
 * The DAG is managed by :
 * - Handle and Rep.
 * - virtual functions to denote the various operators (instead of an enum).
 *
 * Other packages with vaguely similar design : APU, MetaCGAL, LOOK.
 */

/*
 * TODO :
 * - Generalize it for constructions at the kernel level.
 * - Add mixed operations with ET too ?
 * - Interval refinement functionnality ?
 * - Separate the handle and the representation(s) in 2 files (?)
 *   maybe not a good idea, better if everything related to one operation is
 *   close together.
 * - Add a CT template parameter ?
 * - Add a string constant to provide an expression string (a la MetaCGAL) ?
 *   // virtual ostream operator<<() const = 0; // or string, like Core ?
 * - Have a template-expression (?) thing that evaluates a temporary element,
 *   and allocates stuff in memory only when really needs to convert to the
 *   NT.  (cf gmp++, and maybe other things, Blitz++, Synaps...)
 */

/*
 * Interface of the rep classes:
 * - .approx()      returns Interval_nt<> (assumes rounding=nearest).
 *                  [ only called from the handle, and declared in the base ]
 * - .exact()       returns ET, if not already done, computes recursively
 *
 * - .rafine_approx()   ??
 */

namespace CGAL {

template <class NT> class Lazy_exact_nt;


#ifdef CGAL_LAZY_KERNEL_DEBUG
template <typename ET>
inline
void
print_dag(const Lazy_exact_nt<ET>& l, std::ostream& os, int level=0)
{
  l.print_dag(os, level);
}
#endif

// Abstract base representation class for lazy numbers
template <typename ET>
struct Lazy_exact_nt_rep : public Lazy_exact_nt<ET>::Self_rep
{
  typedef typename Lazy_exact_nt<ET>::Self_rep  Base;

  Lazy_exact_nt_rep (const Interval_nt<false> & i)
      : Base(i) {}

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
  }
#endif
};

// int constant
template <typename ET>
struct Lazy_exact_Int_Cst : public Lazy_exact_nt_rep<ET>
{
  Lazy_exact_Int_Cst (int i)
      : Lazy_exact_nt_rep<ET>(double(i)) {}

  void update_exact() const { this->et = new ET((int)this->approx().inf()); }
};

// double constant
template <typename ET, typename X>
struct Lazy_exact_Cst : public Lazy_exact_nt_rep<ET>
{
  Lazy_exact_Cst (X x)
      : Lazy_exact_nt_rep<ET>(x), cste(x) {}

  void update_exact() const { this->et = new ET(cste); }

  private:
  X cste;
};

// Exact constant
template <typename ET>
struct Lazy_exact_Ex_Cst : public Lazy_exact_nt_rep<ET>
{
  Lazy_exact_Ex_Cst (const ET & e)
      : Lazy_exact_nt_rep<ET>(CGAL_NTS to_interval(e))
  {
    this->et = new ET(e);
  }

  void update_exact() const { CGAL_error(); }
};

// Construction from a Lazy_exact_nt<ET1> (which keeps the lazyness).
template <typename ET, typename ET1>
class Lazy_lazy_exact_Cst : public Lazy_exact_nt_rep<ET>
{
  mutable Lazy_exact_nt<ET1> l;

public:

  Lazy_lazy_exact_Cst (const Lazy_exact_nt<ET1> &x)
      : Lazy_exact_nt_rep<ET>(x.approx()), l(x)
  {
    this->set_depth(l.depth() + 1);
  }

  void update_exact() const
  {
    this->et = new ET(l.exact());
    this->at = l.approx();
    prune_dag();
  }

  void prune_dag() const { l = Lazy_exact_nt<ET1>(); }
};


// Unary  operations: abs, sqrt, square.
// Binary operations: +, -, *, /, min, max.

// Base unary operation
template <typename ET>
struct Lazy_exact_unary : public Lazy_exact_nt_rep<ET>
{
  mutable Lazy_exact_nt<ET> op1;

  Lazy_exact_unary (const Interval_nt<false> &i, const Lazy_exact_nt<ET> &a)
      : Lazy_exact_nt_rep<ET>(i), op1(a)
  {
    this->set_depth(op1.depth() + 1);
  }

  void prune_dag() const { op1 = Lazy_exact_nt<ET>(); }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      msg(os, level, "Unary number operator:");
      CGAL::print_dag(op1, os, level+1);
    }
  }
#endif
};

// Base binary operation
template <typename ET, typename ET1 = ET, typename ET2 = ET>
struct Lazy_exact_binary : public Lazy_exact_nt_rep<ET>
{
  mutable Lazy_exact_nt<ET1> op1;
  mutable Lazy_exact_nt<ET2> op2;

  Lazy_exact_binary (const Interval_nt<false> &i,
		     const Lazy_exact_nt<ET1> &a, const Lazy_exact_nt<ET2> &b)
      : Lazy_exact_nt_rep<ET>(i), op1(a), op2(b)
  {
    this->set_depth((std::max)(op1.depth(), op2.depth()) + 1);
  }

  void prune_dag() const
  {
    op1 = Lazy_exact_nt<ET1>();
    op2 = Lazy_exact_nt<ET2>();
  }

#ifdef CGAL_LAZY_KERNEL_DEBUG
  void
  print_dag(std::ostream& os, int level) const
  {
    this->print_at_et(os, level);
    if(this->is_lazy()){
      msg(os, level, "Binary number operator:");
      CGAL::print_dag(op1, os, level+1);
      CGAL::print_dag(op2, os, level+1);
    }
  }
#endif
};

// Here we could use a template class for all operations (STL provides
// function objects plus, minus, multiplies, divides...).  But it would require
// a template template parameter, and GCC 2.95 seems to crash easily with them.

// Macro for unary operations
#define CGAL_LAZY_UNARY_OP(OP, NAME)                                     \
template <typename ET>                                                   \
struct NAME : public Lazy_exact_unary<ET>                                \
{                                                                        \
  typedef typename Lazy_exact_unary<ET>::AT::Protector P;                \
  NAME (const Lazy_exact_nt<ET> &a)                                      \
      : Lazy_exact_unary<ET>((P(), OP(a.approx())), a) {}                \
                                                                         \
  void update_exact() const                                              \
  {                                                                      \
    this->et = new ET(OP(this->op1.exact()));                            \
    if (!this->approx().is_point())                                      \
      this->at = CGAL_NTS to_interval(*(