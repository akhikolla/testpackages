/*--------------------------------------------------------------------*/
/*     Copyright (C) 2013-2013  Serge Iovleff, Quentin Grimonprez

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

    Contact : quentin.grimonprez@inria.fr
*/

/*
 * Project:  MPAGenomics::
 * created on: 4 juin 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Cvlars.h
 *  @brief In this file, we define the @c Cvlars class.
 **/


#ifndef CVLARS_H_
#define CVLARS_H_

namespace HD
{
/**
 * Cross-validation in order to choose the best result of the path of lars algorithm.
 */
  class Cvlars
  {
    public:
      /**
       * Constructor with no index ( it will be a sequence from 0 to 1 by 0.01)
       * @param X matrix of data, a row=a individual
       * @param y response
       * @param nbFolds number of folds
       * @param maxSteps number of maximum step to do
       * @param intercept if true, there is an intercept in the model
       * @param eps epsilon (for 0)
       */
      Cvlars(STK::CArrayXX const& X, STK::CVectorX const& y, int nbFolds, int maxSteps, bool intercept = true, STK::Real eps = STK::Arithmetic<STK::Real>::epsilon());

      /**
       * Constructor
       * @param X matrix of data, a row=a individual
       * @param y response
       * @param nbFolds number of folds
       * @param index elements to test for cross validation (l1norm fraction or lambda)
       * @param lambdaMode if true index contains lambda values, else it contains real between 0 and 1 (ratio (norm coefficient)/max(norm coefficient) for which we compute the prediction error)
       * @param maxSteps number of maximum step to do
       * @param intercept if true, there is an intercept in the model
       * @param eps epsilon (for 0)
       */
      Cvlars(STK::CArrayXX const& X, STK::CVectorX const& y, int nbFolds, std::vector<double> const& index, bool lambdaMode, int maxSteps, bool intercept = true, STK::Real eps = STK::Arithmetic<STK::Real>::epsilon());

      /**
       * run a k-fold cross validation
       */
      void run();
#ifdef _OPENMP
      /**
       * run a k-fold cross validation (parallelized version)
       */
      void run2();
#endif
      //getter
      /** @return return the prediction error for each index*/
      inline STK::CVectorX const& cv() const {return cv_;}
      /** @return return the standard deviation of prediction error for each index*/
      inline STK::CVectorX const& cvError() const {return cvError_;}
      /** @return return the index*/
      inline std::vector<double> const& index() const {return index_;}
      /**@param partition of the */
      void setPartition(std::vector<int> const& partition);

    private:
      /**
       * create a random partition in k folds
       */
      void partition();
      /**
       * run cross validation for folds from idxStartFold to idxEndFold
       * @param idxStartFold index of the first fold
       * @param idxEndFold index of the last fold
       */
      void subrun(int idxStartFold,int idxEndFold);

    private:
      ///pointer on the data
      STK::CArrayXX const* p_X_;
      ///pointer on the response
      STK::CVectorX const* p_y_;
      ///repartition of the sample into k-folds
      std::vector<int> partition_;
      ///size of each fold
      std::vector<int> sizePartition_;
      //std::vector<STK::Range> foldRange_;
      ///vector with real between 0 and 1 (ratio (norm coefficient)/max(norm coefficient) for which we compute the prediction error)
      std::vector<double> index_;
      ///if true, index is fraction, else it's lambda
      bool lambdaMode_;
      ///residuals
      STK::CArrayXX residuals_;
      ///criterion
      STK::CVectorX cv_;
      ///criterion error
      STK::CVectorX cvError_;
      ///number of folds
      int nbFolds_;
      ///number of sample
      int n_;
      ///number of variables
      int p_;
      ///maximum number of steps for the lars algorithm
      int maxSteps_;
      ///numerical zero
      STK::Real eps_;
      /// if true, there is an intercept in the model
      bool intercept_;
  };
}

#endif /* CVLARS_H_ */
