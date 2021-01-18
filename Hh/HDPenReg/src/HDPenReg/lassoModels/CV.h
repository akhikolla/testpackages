/*--------------------------------------------------------------------*/
/*     Copyright (C) 2013-2013  Serge Iovleff, Quentin Grimonprez

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public
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
 * created on: 8 oct. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file CV.h
 *  @brief In this file, definition of the abstract class CV.
 **/


#ifndef CV_H_
#define CV_H_

#include "IAlgo.h"

namespace HD
{
/**
 * Cross-validation class
 * Abstract class. Derived class must implement the runModel method for computing the cverror for all
 * values of index on a specific data test
 */
  class CV
  {

    public:

      /**default constructor*/
      CV();
      /**
       * Constructor
       * @param X matrix of data, a row=a individual
       * @param y response
       * @param nbFolds number of folds
       * @param index vector with value to test
       */
      CV(STK::ArrayXX const& X, STK::VectorX const& y, int nbFolds, std::vector<double> const& index);

      /**destructor*/
      virtual ~CV() {}

      /** run a k-fold cross validation*/
      void run();
      /** run a k-fold cross validation (parallelized version)*/
      void run2();
      void run3();

      //setter
      /**set the data*/
      inline void setX(STK::ArrayXX const& X) {p_X_ = &X;}
      /**set the response*/
      inline void setY(STK::VectorX const& y) {p_y_ = &y;}
      /**set the number of folds for the cross validation*/
      inline void setNbFolds(int const& nbFolds) {nbFolds_ = nbFolds;}
      /**set the index to test*/
      inline void setIndex(std::vector<double> const& index) {index_ = index;}

      //getter
      /** @return return the prediction error for each index*/
      inline STK::VectorX const& cv() const {return cv_;}
      /** @return return the standard deviation of prediction error for each index*/
      inline STK::VectorX const& cvError() const {return cvError_;}
      /** @return return the index*/
      inline std::vector<double> const& index() const {return index_;}

    protected:
      /**initialize the class after setting all parameters*/
      void initializeCV();

      /** create a random partition in k folds */
      void partition();
      /**
       * run cross validation for folds from idxStartFold to idxEndFold (non parallelized version)
       * @param idxStartFold index of the first fold
       * @param idxEndFold index of the last fold
       */
      void subrun(int idxStartFold,int idxEndFold);

      /**
       * run model on all values of index for the given control data
       * @param i index of the actual delete folder
       * @param XTest data test
       * @param yTest response test
       * @param p_XControl pointer to the data control
       * @param p_yControl pointer to the response control
       */
      virtual void runModel(int i, STK::ArrayXX const& XTest, STK::VectorX const& yTest, STK::ArrayXX const* p_XControl, STK::VectorX const* p_yControl) = 0;

    protected:
      ///pointer on the data
      STK::ArrayXX const* p_X_;
      ///pointer on the response
      STK::VectorX const* p_y_;
      ///repartition of the sample into k-folds
      std::vector<int> partition_;
      ///size of each fold
      std::vector<int> sizePartition_;
      ///vector with real between 0 and 1 (ratio (norm coefficient)/max(norm coefficient) for which we compute the prediction error)
      std::vector<double> index_;
      ///measure of good fitness : mse, auc
      STK::ArrayXX measure_;
      ///criterion
      STK::VectorX cv_;
      ///criterion error
      STK::VectorX cvError_;
      ///number of folds
      int nbFolds_;
      ///number of sample
      int n_;
      ///number of variables
      int p_;
  };
}


#endif /* CV_H_ */
