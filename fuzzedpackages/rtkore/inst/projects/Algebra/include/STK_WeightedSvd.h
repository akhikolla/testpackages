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
 * Project:  stkpp::Algebra
 * created on: 10 août 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_WeightedSvd.h
 *  @brief In this file we define the WeightedSvd class.
 **/

#ifndef STK_WEIGHTEDSVD_H
#define STK_WEIGHTEDSVD_H

namespace STK
{
// forward declaration
template<class Array, class WRows, class WCols> class WeightedSvd;

namespace hidden
{
/** @ingroup hidden
 *  Specialization for the WeightedSvd class.
 **/
template<class Array_, class WRows, class WCols>
struct AlgebraTraits< WeightedSvd<Array_, WRows, WCols> >
{
  typedef CArrayXX ArrayU;
  typedef CVectorX ArrayD;
  typedef CArrayXX ArrayV;
};
} // namespace hidden

/** @brief This class perform a weighted svd decomposition */
template<class Array, class WRows, class WCols>
class WeightedSvd: public ISvd< WeightedSvd<Array, WRows, WCols> >
{
  public:
    typedef ISvd< WeightedSvd<Array, WRows, WCols> > Base;
    using Base::U_;
    using Base::D_;
    using Base::V_;
    /** default constructor.
     *  The matrix U_ will be weigthed by this constructor on the fly.
     *  The decomposition will be achieved using lapack (if present) or STK++
     *  implementations (otherwise) when the methode @© run is called.
     *  @param a the matrix to decompose
     *  @param wrows, wcols weights of the rows and columns
     *  @param dim the number of left and right eigenVectors required
     **/
    WeightedSvd( Array const& a, WRows const& wrows, WCols const& wcols, int dim)
              : Base(a, false, (dim>0) ? true:false, (dim>0) ? true:false)
               , wrows_(wrows), wcols_(wcols), dim_(dim)

    {
      if (wrows.range() != U_.rows()) { wrows_.resize(U_.rows()).setValue(1./U_.sizeRows());}
      else                            { wrows_ /= wrows.sum();}
      if (wcols.range() != U_.cols()) { wcols_.resize(U_.cols()).setOne();}
      dim_ = std::min(dim_, U_.sizeRows());
      dim_ = std::min(dim_, U_.sizeCols());
    }
    /** destructor */
    virtual ~WeightedSvd() {}
    /** run the weighted svd */
    virtual bool run()
    {
      U_ = wrows_.diagonalize() * U_ * wcols_.diagonalize();
#ifdef STKUSELAPACK
      lapack::Svd solver(U_, false, this->withU_, this->withV_);
      // if there is no cv, fall back to STK++ svd
      if (!solver.run())
      {
        Svd<CArrayXX> dec(U_, true, this->withU_, this->withV_);
        if (!dec.run()) return false;
        else
        {
          U_.move(dec.U_);
          D_ = dec.D_;
          V_ = dec.V_;
        }
      }
      else
      {
        U_.move(solver.U_);
        D_.move(solver.D_);
        V_.move(solver.V_);
      }
#else
      Svd<CArrayXX> solver(U_, true, this->withU_, this->withV_);
      if (!solver.run()) return false;
      else
      {
        U_.move(dec.U_);
        D_ = dec.D_;
        V_ = dec.V_;
      }

#endif
      // weight back
      return true;
    }
  private:
    /** rows weights*/
    CVectorX wrows_;
    /** columns weights*/
    CPointX wcols_;
    /** number of eigenvectors (left and right)*/
    int dim_;
};

} // namespace STK

#endif /* STK_WEIGHTEDSVD_H */
