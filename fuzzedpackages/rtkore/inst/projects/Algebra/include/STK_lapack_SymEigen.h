/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff

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
 * Project:  stkpp::Algebra
 * created on: 20 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_lapack_SymEigen.h
 *  @brief In this file we define the enclosing class of the syevr lapck routine.
 **/


#ifndef STK_LAPACK_SYMEIGEN_H
#define STK_LAPACK_SYMEIGEN_H

#include "STK_ISymEigen.h"
#include "STK_lapack_Util.h"

namespace STK
{

// forward declaration
namespace lapack
{
template<class SquareArray> class SymEigen;
}

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the lapack::SymEigen class.
 **
 **/
template<class SquareArray_>
struct AlgebraTraits< lapack::SymEigen<SquareArray_> >
{
  typedef SquareArray_ SquareArray;
};

} // namespace hidden


namespace lapack
{
/** @ingroup Algebra
 *  {
 *    @class SymEigen
 *    @brief  SymEigen computes the eigenvalues and optionally the
 *    eigenvectors of a symmetric real matrix using the syevr Lapack routine.
 *
 *  The decomposition of a symmetric matrix require
 *  - Input:  A symmetric matrix A of size (n,n)
 *  - Output:
 *     -# P Array of size (n,n).
 *     -# D Vector of dimension n
 *     -# \f$ A = PDP' \f$
 *  The matrix A can be copied or overwritten by the class.
 *
 *  The 2-norm (operator norm) of the matrix is given. if the 2-norm is less
 *  than the arithmetic precision of the type @c Real, the rank is not
 *  full.
 *  Thus the user can be faced with a deficient rank matrix and with a norm and
 *  a determinant very small (but not exactly 0.0).
 *
 *  @note Default values for syevr lapack parameters are
 *  @code
 *  jobz_ = 'V'; RANGE_ = 'A'; UPLO_ = 'U';
 *   VL_ = 0.0; VU_ = 0.0; IL_ = 0; IU_ = 0;
 *  @endcode
 *
 *  @sa STK::ISymEigen, STK::SymEigen
 */
template<class SquareArray>
class SymEigen: public ISymEigen< SymEigen<SquareArray> >
{
  public:
    typedef ISymEigen< SymEigen<SquareArray> > Base;
    using Base::eigenValues_;
    using Base::eigenVectors_;
    using Base::range_;
    using Base::norm_;
    using Base::rank_;
    using Base::det_;

    /** @brief Default Constructor */
    SymEigen();
    /** @brief Constructor
     *  @param data reference on a symmetric square matrix
     *  @param ref @c true if we overwrite the data set, @c false otherwise
     *  @note data can be a reference if and only if it is a CSquareX
     *  @note Default values for syevr lapack parameters are
     *  @code
     *  jobz_ = 'V'; RANGE_ = 'A'; UPLO_ = 'U';
     *   VL_ = 0.0; VU_ = 0.0; IL_ = 0; IU_ = 0;
     *  @endcode
     */
    SymEigen( SquareArray const& data, bool ref =false);
    /** @brief Constructor
     *  @param data reference on a symmetric square expression
     */
    template<class Derived>
    SymEigen( ArrayBase<Derived> const& data);
    /** @brief copy constructor
     *  @param eigen the SymEigen to copy
     */
    SymEigen( SymEigen const& eigen);
    /** Destructor. */
    virtual ~SymEigen() {}

    // setters
    /** @param jobz If jobz ='N': Compute eigenvalues only; If jobz = 'V': Compute
     * eigenvalues and eigenvectors.
     **/
    inline void setJobz(char jobz) { jobz_ = jobz;}
    /** @param range range of the eigenvalues to be found.
     *  If range = 'A': all eigenvalues will be found.
     *  If range = 'V': all eigenvalues in the half-open interval  (VL_,VU_] will be found.
     *  If range = 'I': the IL_-th through IU_-th eigenvalues will be found.
     **/
    inline void setRange(char range) { RANGE_ = range;}
    /** @param uplo values to used in A.
     * If uplo = 'U':  Upper triangle of A is stored;
     * If uplo = 'L':  Lower triangle of A is stored.
     **/
    inline void setUplo(char uplo) { UPLO_ = uplo;}
    /** @param[in] vl,vu lower and upper bounds of the interval to be searched
     * for eigenvalues.
     * Not referenced if RANGE_ = 'A' or 'I'.
    **/
    inline void setVlAndVu(Real const& vl, Real const& vu) { VL_ = vl; VU_ = vu;}
    /** @param il, iu
     *  If RANGE_='I', the indices (in ascending order)
     *  of the smallest and largest eigenvalues to be returned.
     *  1 <= IL_ <= IU_ <= NL, if NL > 0; IL_ = 1 and IU_ = 0 if NL = 0. Not
     *  referenced if RANGE_ = 'A' or 'V'.
    **/
    inline void setIlAndIu(int il, int iu) { IL_ = il; IU_ = iu;}

    /** @brief clone pattern */
    inline virtual SymEigen* clone() const { return new SymEigen(*this);}
    /** @brief Run eigenvalues decomposition
     *  Launch SYEVR LAPACK routine to perform the eigenvalues decomposition.
     *  @return @c true if no error occur, @c false otherwise
     */
    bool runImpl();

  private:
    /** Lapack pptions */
    char jobz_, RANGE_, UPLO_;
    Real VL_, VU_;
    int IL_, IU_;
};
/** @} */



/* @brief Default Constructor */
template<class SquareArray>
SymEigen<SquareArray>::SymEigen(): Base()
                                 , jobz_('V'), RANGE_('A'), UPLO_('U')
                                 , VL_(0.0), VU_(0.0), IL_(0), IU_(0)
{}
/* @brief Constructor
 *  @param data reference on a symmetric square matrix
 *  @param ref @c true if we overwrite the data set, @c false otherwise
 */
template<class SquareArray>
SymEigen<SquareArray>::SymEigen( SquareArray const& data, bool ref)
                              : Base(data)
                               , jobz_('V'), RANGE_('A'), UPLO_('U')
                               , VL_(0.0), VU_(0.0), IL_(0), IU_(0)
{}
/* @brief Constructor
 *  @param data reference on a symmetric square matrix
 *  @param ref @c true if we overwrite the data set, @c false otherwise
 */
template<>
inline SymEigen<CSquareX>::SymEigen( CSquareX const& data, bool ref)
                                  : Base(data, ref)
                                   , jobz_('V'), RANGE_('A'), UPLO_('U')
                                   , VL_(0.0), VU_(0.0), IL_(0), IU_(0)
{}

template<class SquareArray>
template<class Derived>
SymEigen<SquareArray>::SymEigen( ArrayBase<Derived> const& data)
                              : Base(data)
                               , jobz_('V'), RANGE_('A'), UPLO_('U')
                               , VL_(0.0), VU_(0.0), IL_(0), IU_(0)
{}
/** @brief copy constructor
 *  @param eigen the SymEigen to copy
 */
template<class SquareArray>
SymEigen<SquareArray>::SymEigen( SymEigen const& eigen)
                              : Base(eigen)
                               , jobz_(eigen.jobz_), RANGE_(eigen.RANGE_), UPLO_(eigen.UPLO_)
                               , VL_(eigen.VL_), VU_(eigen.VU_), IL_(eigen.IL_), IU_(eigen.IU_)
 {}

/* @brief Run eigen decomposition
 *  Launch SYEVR LAPACK routine to perform the eigenvalues decomposition.
 *  @return @c true if no error occur, @c false otherwise
 */
template<class SquareArray>
bool SymEigen<SquareArray>::runImpl()
{
#ifdef STK_ALGEBRA_VERY_VERBOSE
  stk_cout << _T("Enter in SymEigen::run\n");
#endif
  /* copy square matrix with the original data set. */
  CSquareX data_ = eigenVectors_;
  // shift data sets
  data_.shift(0);
  eigenVectors_.shift(0);
  eigenValues_.shift(0);
  this->SupportEigenVectors_.shift(0);
  /* set default behavior */
  Real absTol = 0.0; // let Lapack chose the correct tolerance
  // get optimal size necessary for work
  Real work; // work is just one place, get the optimal size for work
  // iwork is just on place, get the optimal size for iWork
  int iwork, lWork =-1, liwork =-1; // workspace variable

  int info = 1;
#ifdef STK_ALGEBRA_DEBUG
  stk_cout << _T("Data dimensions: ") << data_.rows() << " " << data_.cols() << "\n";
  stk_cout << _T("eigenValues_ dimensions: ") << eigenValues_.rows() << " " << eigenValues_.cols() << "\n";
  stk_cout << _T("eigenVectors_ dimensions: ") << eigenVectors_.rows() << " " << eigenVectors_.cols() << "\n";
  stk_cout << _T("Options: ") << jobz_ << " " << RANGE_ << " " << UPLO_ << "\n";
#endif
  info = syevr( jobz_, RANGE_, UPLO_
              , range_.size(), data_.p_data(), range_.size()
              , VL_, VU_, IL_, IU_
              , absTol, &rank_,  eigenValues_.p_data()
              , eigenVectors_.p_data(), range_.size(), this->SupportEigenVectors_.p_data()
              , &work, lWork, &iwork, liwork);
  // check any error
  if (info!=0)
  {
    if (info>0)
    { this->msg_error_ = STKERROR_NO_ARG(SymEigen::run,internal error);
      return false;
    }
    this->msg_error_= STKERROR_1ARG(SymEigen::run,-info,error parameter);
    return false;
  }
#ifdef STK_ALGEBRA_DEBUG
  stk_cout << _T("Size needed:") << (int)work << " " << iwork << " " << "\n";
#endif
  // get results and allocate space
  lWork = (int)work;
  liwork = iwork;
  Real* p_work = new Real[lWork];
  int* p_iwork = new int[liwork];

  // Call SYEVR with the optimal block size
  info = syevr( jobz_, RANGE_, UPLO_
              , range_.size(), data_.p_data(), range_.size()
              , VL_, VU_, IL_, IU_
              , absTol, &rank_, eigenValues_.p_data()
              , eigenVectors_.p_data(), range_.size(), this->SupportEigenVectors_.p_data()
              , p_work, lWork, p_iwork, liwork);
  // recover memory
  delete[] p_work;
  delete[] p_iwork;

  // finalize
  data_.shift(range_.begin());
  eigenVectors_.shift(range_.begin());
  eigenValues_.shift(range_.begin());
  this->SupportEigenVectors_.shift(range_.begin());
  this->finalizeStep();
  this->hasRun_ = true;
  // return the result of the computation
  if (!info) return true;
  if (info>0)
  { this->msg_error_ = STKERROR_NO_ARG(SymEigen::run,internal error);
    return false;
  }
  this->msg_error_= STKERROR_1ARG(SymEigen::run,-info,error parameter);
  return false;
}

} // namespace lapack

} // namespace STK

#endif /* STK_LAPACK_SYMEIGEN_H */
