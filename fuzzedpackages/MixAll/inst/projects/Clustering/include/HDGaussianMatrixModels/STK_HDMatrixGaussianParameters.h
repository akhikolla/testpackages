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

/** @file  STK_HDMatrixGaussianParameters.h
 *  @brief In this file we define the Parameters classes for High Dimensional
 *  Matrix Gaussian mixture models
 **/

#ifndef STK_HDMATRIXGAUSSIANPARAMETERS_H
#define STK_HDMATRIXGAUSSIANPARAMETERS_H

#include <Arrays/include/STK_Array1D.h>
#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_CArraySquare.h>
#include <Arrays/include/STK_CArrayPoint.h>
#include <Arrays/include/STK_CArrayVector.h>

#include <STatistiK/include/STK_Stat_Online.h>

namespace STK
{

/** @ingroup Clustering
 *  @brief structure storing the parameters of the HD matrix valued mixture
 **/
template<class Array_>
struct HDMatrixModelParameters
{
    /** array of size nbCluster with for each cluster the mean matrix */
    Array1D<Array_> meank_;

    /** array of size nbCluster with the rotation matrix for rows */
    Array1D<CSquareX> rowQk_;
    /** array of size nbCluster with inertia in low dimensional space for rows */
    Array1D<CVectorX> rowAjk_;
    /** array of size nbCluster with remaining variance noise for rows */
    Array1D<Real> rowBk_;

    /** array of size nbCluster with dimension of low dimensional space for rows */
    CVectorXi rowDk_;

    /** array of size nbCluster with the rotation matrix for columns */
    Array1D<CSquareX> colQk_;
    /** array of size nbCluster with inertia in low dimensional space for columns */
    Array1D<CVectorX> colAjk_;
    /** array of size nbCluster with variance noise for columns */
    Array1D<Real> colBk_;

    /** array of size nbCluster with dimension of low dimensional space for columns */
    CVectorXi colDk_;

    /** Array of size nbCluster of the mean statistics */
    Array1D< Stat::Online<Array_, Real> > statMeank_;

    /** Array  of size nbCluster of the rotation matrix statistics for rows */
    Array1D< Stat::Online<CSquareX, Real> > statRowQk_;
    /** Array of size nbCluster of the rotation matrix statistics for rows */
    Array1D< Stat::Online<CVectorX, Real> > statRowAjk_;
    /** Array of size nbCluster with the variance noise statistics for rows */
    Array1D< Stat::Online<Real, Real> > statRowBk_;

    /** Array of the rotation matrix statistics for rows */
    Array1D< Stat::Online<CSquareX, Real> > statColQk_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<CVectorX, Real> > statColAjk_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > statColBk_;

    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    HDMatrixModelParameters( int nbCluster)
                           : meank_(nbCluster)
                           , rowQk_(nbCluster), rowAjk_(nbCluster)
                           , rowBk_(nbCluster), rowDk_(nbCluster)
                           , colQk_(nbCluster), colAjk_(nbCluster)
                           , colBk_(nbCluster), colDk_(nbCluster)
                           , statMeank_(nbCluster)
                           , statRowQk_(nbCluster), statRowAjk_(nbCluster)
                           , statRowBk_(nbCluster)
                           , statColQk_(nbCluster), statColAjk_(nbCluster)
                           , statColBk_(nbCluster)
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    HDMatrixModelParameters( HDMatrixModelParameters const& param)
                           : meank_(param.meank_)
                           , rowQk_(param.rowQk_), rowAjk_(param.rowAjk_)
                           , rowBk_(param.rowBk_), rowDk_(param.rowDk_)
                           , colQk_(param.colQk_), colAjk_(param.colAjk_)
                           , colBk_(param.colBk_), colDk_(param.colDk_)
                           , statMeank_(param.statMeank_)
                           , statRowQk_(param.statRowQk_), statRowAjk_(param.statRowAjk_)
                           , statRowBk_(param.statRowBk_)
                           , statColQk_(param.statColQk_), statColAjk_(param.statColAjk_)
                           , statColBk_(param.statColBk_)
    {}
    /** destructor */
    ~HDMatrixModelParameters() {}
    /** copy operator.
     *  @param param the parameters to copy.
     **/
    HDMatrixModelParameters& operator=( HDMatrixModelParameters const& param)
    {
      meank_ = param.meank_;
      rowQk_ = param.rowQk_; rowAjk_ = param.rowAjk_;
      rowBk_ = param.rowBk_; rowDk_  = param.rowDk_;
      colQk_ = param.colQk_; colAjk_ = param.colAjk_;
      colBk_ = param.colBk_; colDk_  = param.colDk_;
      statMeank_ = param.statMeank_;
      statRowQk_ = param.statRowQk_; statRowAjk_ = param.statRowAjk_;
      statRowBk_ = param.statRowBk_;
      statColQk_ = param.statColQk_; statColAjk_ = param.statColAjk_;
      statColBk_ = param.statColBk_;
    }

    /** resize and initialize the set of parameter.
     *  @param range range of the variables in the data set
     **/
    void resize(Range const& rangeRows, Range const& rangeCols);

    /** update statistics of the parameters. */
    void updateStatistics()
    {
      for(int k=statMeank_.begin(); k<statMeank_.end(); ++k)
      {
        statMeank_[k].update(meank_[k]);
        statRowQk_[k].update(rowQk_[k]);
        statRowAjk_[k].update(rowAjk_[k]);
        statRowBk_[k].update(rowBk_[k]);
        statColQk_[k].update(colQk_[k]);
        statColAjk_[k].update(colAjk_[k]);
        statColBk_[k].update(colBk_[k]);
      }
    }
    /** Set the computed statistics */
    void setStatistics()
    {
      for(int k=statMeank_.begin(); k<statMeank_.end(); ++k)
      {
        meank_[k]  = statMeank_[k].mean();
        rowQk_[k]  = statRowQk_[k].mean();
        rowAjk_[k] = statRowAjk_[k].mean();
        rowBk_[k]  = statRowBk_[k].mean();
        colQk_[k]  = statColQk_[k].mean();
        colAjk_[k] = statColAjk_[k].mean();
        colBk_[k]  = statColBk_[k].mean();
      }
      releaseStatistics();
    }
    /** Release the computed statistics */
    void releaseStatistics()
    {
      for(int k=statMeank_.begin(); k<statMeank_.end(); ++k)
      {
        statMeank_[k].release();
        statRowQk_[k].release();
        statRowAjk_[k].release();
        statRowBk_[k].release();
        statColQk_[k].release();
        statColAjk_[k].release();
        statColBk_[k].release();
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the mean and
     *  sd parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters( ExprBase<Array> const& params)
    {}
//    {
//      for(int k=mean_.begin(), kp= params.beginRows(); k<mean_.end(); ++k, kp+=2)
//      {
//        mean_[k] = params.row(kp);
//        sigma_[k] = params.row(kp+1);
//      }
//    }
};

} // namespace STK

#endif /* STK_HDGAUSSIANPARAMETERS_H */
