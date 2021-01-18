/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016 Serge Iovleff

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: Oct 24, 2013
 * Author:   Serge Iovleff
 **/

/** @file STK_HDMatrixGaussianModel.h
 *  @brief In this file we define the HDMatrixGaussianModel class
 **/

#ifndef STK_HDMATRIXGAUSSIANMODEL_H
#define STK_HDMATRIXGAUSSIANMODEL_H

#include "../STK_IMixtureDensity.h"
#include "STK_HDMatrixGaussianParameters.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<int IdRow_, int IdCol_, class Array>
class HDMatrixGaussianModel;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the HDMatrixGaussianModel traits policy. */
template<int IdRow_, int IdCol_, class Array_>
struct MixtureTraits< HDMatrixGaussianModel<IdRow_, IdCol_, Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a MixturGaussian_sjk model*/
  typedef HDMatrixModelParameters<Array_> Parameters;
};

} // namespace hidden


/** @ingroup Clustering
 *  Main class for HD matrix valued matrix models
 **/
template<int IdRow_, int IdCol_, class Array_>
class HDMatrixGaussianModel: public IMixtureDensity<HDMatrixGaussianModel<IdRow_, IdCol_, Array_> >
{
  public:

    /** Base class */
    typedef IMixtureDensity<HDMatrixGaussianModel<IdRow_, IdCol_, Array_> > Base;
    /** Type of the structure storing the parameters of a MixturGaussian_sjk model*/
    typedef HDMatrixModelParameters<Array_> Parameters;

    using Base::nbCluster;
    using Base::nbSample;
    using Base::param_;
    using Base::p_data;

    /** constructor
     *  @param nbCluster number of clusters
     **/
    HDMatrixGaussianModel( int nbCluster)
                         : Base(nbCluster)
                         , isRowAj_(hidden::HDCovarianceChooser<IdRow_>::isAj_)
                         , isRowAk_(hidden::HDCovarianceChooser<IdRow_>::isAk_)
                         , isRowBk_(hidden::HDCovarianceChooser<IdRow_>::isBk_)
                         , isRowQk_(hidden::HDCovarianceChooser<IdRow_>::isQk_)
                         , isRowDk_(hidden::HDCovarianceChooser<IdRow_>::isDk_)
                         , isColAj_(hidden::HDCovarianceChooser<IdCol_>::isAj_)
                         , isColAk_(hidden::HDCovarianceChooser<IdCol_>::isAk_)
                         , isColBk_(hidden::HDCovarianceChooser<IdCol_>::isBk_)
                         , isColQk_(hidden::HDCovarianceChooser<IdCol_>::isQk_)
                         , isColDk_(hidden::HDCovarianceChooser<IdCol_>::isDk_)
    {}
    /** constructor
     *  @param nbCluster number of clusters
     **/
    HDMatrixGaussianModel( HDMatrixGaussianModel const& model)
                         : Base(model)
                         , nbRow_(model.nbRow_)
                         , nbCol_(model.nbCol_)
                         , isRowAj_(model.isRowAj_), isRowAk_(model.isRowAk_)
                         , isRowBk_(model.isRowBk_)
                         , isRowQk_(model.isRowQk_)
                         , isRowDk_(model.isRowDk_)
                         , isColAj_(model.isColAj_), isColAk_(model.isColAk_)
                         , isColBk_(model.isColBk_)
                         , isColQk_(model.isColQk_)
                         , isColDk_(model.isColDk_)
    {}
    /** destructor */
    ~HDMatrixGaussianModel(){}

    /** Initialize randomly the parameters of the Gaussian mixture. The centers
     *  will be selected randomly among the data set and the standard-deviation
     *  will be set to 1.
     */
    void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk);
    /** Compute the weighted mean and the common standard deviation. */
    bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk);
    /** @return the number of free parameters of the model */
    int computeNbFreeParameters() const;

  protected:
    /** number of rows and columns of the data */
    int nbRow_, nbCol_;
    /** Structure of the covariance matrices between clusters for the rows */
    bool isRowAj_, isRowAk_, isRowBk_, isRowQk_, isRowDk_;
    /** Structure of the covariance matrices between clusters for the columns */
    bool isColAj_, isColAk_, isColBk_, isColQk_, isColDk_;

  private:
    /** Update parameters for models with free orientation by rows */
    bool runFreeOrientationByRow( CArrayXX const* const& p_tik, CPointX const* const& p_tk);
    /** Update parameters for models with common orientation by rows*/
    bool runCommonOrientationByRow( CArrayXX const* const& p_tik, CPointX const* const& p_tk);
    /** Update parameters for models with free orientation by columns */
    bool runFreeOrientationByCol( CArrayXX const* const& p_tik, CPointX const* const& p_tk);
    /** Update parameters for models with common orientation by columns */
    bool runCommonOrientationByCol( CArrayXX const* const& p_tik, CPointX const* const& p_tk);
};

template<int IdRow_, int IdCol_, class Array_>
int HDMatrixGaussianModel<IdRow_, IdCol_, Array_>::computeNbFreeParameters() const
{
  int sum = nbRow_*nbCol_;
  int maxRowDk = param_.rowDk_.maxElt();
  int maxColDk = param_.colDk_.maxElt();
  // for Bk
  sum += (isRowBk_)? nbCluster :1;
  sum += (isColBk_)? nbCluster :1;
  // for Ajk
  sum += (isRowAk_) ? ((isRowAj_) ? param_.rowDk_.sum() : nbCluster_) : maxRowDk;
  sum += (isColAk_) ? ((isColAj_) ? param_.colDk_.sum() : nbCluster_) : maxColDk;
  // for Qk
  if (isRowQk_)
  {
    int dk = param_.rowDk_[k];
    sum += dk*nbRow_ - (dk*(dk+1))/2;
  }
  else
  { sum += maxRowDk*nbRow - (maxRowDk*(maxRowDk+1))/2;}
  if (isColQk_)
  {
    int dk = param_.colDk_[k];
    sum += dk*nbCol_ - (dk*(dk+1))/2;
  }
  else
  { sum += maxColDk*nbRow - (maxColDk*(maxColDk+1))/2;}
  //
  return sum;
}
/* Compute the weighted mean and the common standard deviation. */
template<int IdRow_, int IdCol_, class Array_>
bool HDMatrixGaussianModel<IdRow_, IdCol_, Array_>::run( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{

  return true;
}

template<int IdRow_, int IdCol_, class Array_>
void HDMatrixGaussianModel<IdRow_, IdCol_, Array_>::randomInit( CArrayXX const* const&  p_tik, CPointX const* const& p_tk)
{
  VectorXi indexes(p_data()->rows());
  for(int i=indexes.begin(); i< indexes.end(); ++i) { indexes[i] = i;}
  Range rind(p_data()->rows());
  // sample mean between individuals without repetition
  for (int k= p_tk->begin(); k < p_tk->end(); ++k)
  {
    // random number in [0, end-k[
    int i = Law::UniformDiscrete::rand(rind.begin(), rind.end()-1);
    // get ith individuals
    param_.meank_[k] = p_data()->row(indexes[i]);
    // exchange it with nth
    indexes.swap(i, rind.lastIdx());
    // decrease
    rind.decLast(1);
  }
  //
  Real rowB = Law::Exponential::rand(1);
  Real colB = Law::Exponential::rand(1);;
  Real rowA = rowB + Law::Exponential::rand(1);
  Real colA = colB + Law::Exponential::rand(1);
  //
  for (int k= p_tk->begin(); k < p_tk->end(); ++k)
  {
    Real v = Law::Normal::rand(0,1), w = Law::Normal::rand(0,1);
    param_.rowQk_[k]  = v * Const::Identity(nbRow_);
    param_.colQk_[k]  = w * Const::Identity(nbCol_);
    param_.rowAjk_[k] = rowA;
    param_.rowBk_[k]  = rowB;
    param_.colAjk_[k] = colA;
    param_.colBk_[k]  = colB;
  }
}

/* Update parameters for models with free orientation */
template<int IdRow_, int IdCol_, class Array_>
bool HDMatrixGaussianModel<IdRow_, IdCol_, Array_>::runFreeOrientationByRow( CArrayXX const* const& p_tik, CPointX const* const& p_tk)
{
  CSquareX w;
  STK::lapack::SymEigen<CSquareX> decomp(w);
  for (int k=p_tk->begin(); k< p_tk->end(); ++k)
  {

    // compute the unbiased covariance function
    Stat::Covariance(*p_data(), p_tik()->col(k), false);
  }
}

} // namespace STK

#endif /* STK_HDMATRIXGAUSSIANMODEL_H */
