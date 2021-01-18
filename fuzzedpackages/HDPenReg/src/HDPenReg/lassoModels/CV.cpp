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

/** @file CV.cpp
 *  @brief In this file, implementation of the method of the class @c CV .
 **/


#include "CV.h"
#include <cstdlib>
#include <algorithm>

namespace HD
{

  /*default cosntructor*/
  CV::CV()
        : p_X_(0), p_y_(0), partition_()
        , sizePartition_(), index_()
        , measure_(), cv_()
        , cvError_(), nbFolds_(0)
        , n_(0), p_(0)
  {
  }
  /*
   * Constructor
   * @param X matrix of data, a row=a individual
   * @param y response
   * @param nbFolds number of folds
   * @param index vector with real between 0 and 1 (ratio (norm coefficient)/max(norm coefficient) for which we compute the prediction error)
   */
  CV::CV(STK::ArrayXX const& X, STK::VectorX const& y, int nbFolds, std::vector<double> const& index)
                : p_X_(&X)
                , p_y_(&y)
                , partition_(X.sizeRows())
                , sizePartition_(nbFolds,0)
                , index_(index)
                , measure_(index.size(),nbFolds)
                , cv_(index.size(),0)
                , cvError_(index.size())
                , nbFolds_(nbFolds)
                , n_(X.sizeRows())
                , p_(X.sizeCols())
  {
    //create the partition
    partition();
  }

  /*
   * create a random partition in k folds
   */
  void CV::partition()
  {
    //fill the container with the index of folds
    for(int i = 0 ; i< n_ ;i++)
    {
      partition_[i] = i%nbFolds_;
      sizePartition_[i%nbFolds_]++;
    }
    //make a random rearrangement
    srand(time(NULL));
    random_shuffle(partition_.begin(),partition_.end());
  }

  /*initialize containers and create partition*/
  void CV::initializeCV()
  {
    n_ = p_X_->sizeRows();
    p_ = p_X_->sizeCols();
    partition_.resize(p_X_->sizeRows());
    sizePartition_.resize(nbFolds_);
    measure_.resize(index_.size(),nbFolds_);
    cv_.resize(index_.size());
    cvError_.resize(index_.size());
    partition();
  }

  /*
   * run the cross validation
   */
  void CV::run()
  {
    //search the first and last fold with the same size
    std::vector<int> startIndex(1,0),endIndex(1,nbFolds_-1);
    int k = 0;
    for(int i = 1; i < nbFolds_; i++)
    {
      if(sizePartition_[i]!= sizePartition_[startIndex[k]])
      {
        startIndex.push_back(i);
        endIndex[k] = i-1;
        endIndex.push_back(nbFolds_-1);
        k++;
      }
    }
    //run for each size of fold
    for(int i = 0; i < (int) startIndex.size(); i++)
      subrun(startIndex[i],endIndex[i]);
    // compute mean prediction error for each index
    STK::VectorX one(nbFolds_,1);
    cv_ = (measure_ * one) / nbFolds_;
    // compute mean standard deviation of cv_ for each index
    for(int i = 0; i < (int) index_.size(); i++)
      for(int j = 0; j < nbFolds_; j++)
        measure_(i,j) -= cv_[i];
      //residuals_.row(i) -= cv_[i];
    measure_ = measure_.square();
    cvError_ = (measure_ * one)/(nbFolds_-1)/nbFolds_;
    cvError_ = cvError_.sqrt();

  }

/*
 * run cross validation for folds from idxStartFold to idxEndFold
 * @param idxStartFold index of the first fold
 * @param idxEndFold index of the last fold
 */
void CV::subrun(int idxStartFold,int idxEndFold)
{
  //create test and control container
  STK::ArrayXX XControl( n_ - sizePartition_[idxStartFold], p_);
  STK::VectorX yControl( n_ - sizePartition_[idxStartFold] );
  STK::ArrayXX XTest(sizePartition_[idxStartFold], p_);
  STK::VectorX yTest(sizePartition_[idxStartFold] );
  STK::VectorX yPred(sizePartition_[idxStartFold] );

  for(int i = idxStartFold ; i <= idxEndFold ; i++)
  {
    //fill the container
    int index = 1;
    int index2 = 1;
    for(int j = 0; j < n_; j++)
    {
      if(partition_[j-1] != i)
      {
        yControl[index] = (*p_y_)[j];
        XControl.row(index)=p_X_->row(j);

        index++;
      }
      else
      {
        yTest[index2] = (*p_y_)[j];
        XTest.row(index2)=p_X_->row(j);
        index2++;
      }
    }
    runModel(i,XTest,yTest,&XControl,&yControl);
  }
}

/*parallelized version of run*/
void CV::run2()
 {
#ifdef HD_CVDEBUG
       std::cout << "Entering CV::run2\n";
       for (int i=0; i<sizePartition_.size(); i++)
       { std::cout << sizePartition_[i] << " ";}
       std::cout << "\n";
#endif
   //run for each size of fold
   //create test and control container
//#pragma omp parallel
   {
//     #pragma omp for schedule(dynamic,1)
     for(int i = 0; i < nbFolds_ ; i++)
     {
#ifdef HD_CVDEBUG
       std::cout << "---------->Creating Fold " << i << "\n";
#endif
       STK::ArrayXX XControl( n_ - sizePartition_[i], p_);
       STK::VectorX yControl( n_ - sizePartition_[i] );
       STK::ArrayXX XTest(sizePartition_[i], p_);
       STK::VectorX yTest(sizePartition_[i] );
       //fill the container
       int index = 0;
       int index2 = 0;
       for(int j = 0; j < n_; j++)
       {
         if(partition_[j] != i)
         {
           yControl[index] = (*p_y_)[j];
           XControl.row(index)=p_X_->row(j);
           index++;
         }
         else
         {
           yTest[index2] = (*p_y_)[j];
           XTest.row(index2)=p_X_->row(j);
           index2++;
         }
       }
#ifdef HD_CVDEBUG
       std::cout << "Fold build. Call runModel\n";
#endif
       runModel(i,XTest,yTest,&XControl,&yControl);
     }
   }//end parallel
#ifdef HD_CVDEBUG
       std::cout << "measure_ =\n" << measure_ << "\n";
#endif
   // compute mean prediction error for each index
   cv_ = STK::sumByRow(measure_) / nbFolds_;
#ifdef HD_CVDEBUG
       std::cout << "cv_ =" << cv_ << "\n";
#endif
   // compute mean standard deviation of cv_ for each index
   for(int i = 0; i < (int) index_.size(); i++)
     measure_.row(i) -= cv_[i];
   measure_ = measure_.square();
   cvError_ = STK::sumByRow(measure_)/(nbFolds_-1)/nbFolds_;
   cvError_ = cvError_.sqrt();
 }

void CV::run3()
{
   //search the first and last fold with the same size
   std::vector<int> startIndex(1,0),endIndex(1,nbFolds_-1);
   int k = 0;
   for(int i = 1; i < nbFolds_; i++)
   {
     if(sizePartition_[i]!= sizePartition_[startIndex[k]])
     {
       startIndex.push_back(i);
       endIndex[k] = i-1;
       endIndex.push_back(nbFolds_-1);
       k++;
     }
   }
   //run for each size of fold
   //create test and control container
     for(int i = 0; i < nbFolds_ ; i++)
     {
       STK::ArrayXX XControl( n_ - sizePartition_[i], p_);
       STK::VectorX yControl( n_ - sizePartition_[i] );
       STK::ArrayXX XTest(sizePartition_[i], p_);
       STK::VectorX yTest(sizePartition_[i] );
       //fill the container
       int index = 1;
       int index2 = 1;
       for(int j = 0; j < n_; j++)
       {
         if(partition_[j] != i)
         {
           yControl[index] = (*p_y_)[j];
           XControl.row(index)=p_X_->row(j);
           index++;
         }
         else
         {
           yTest[index2] = (*p_y_)[j];
           XTest.row(index2)=p_X_->row(j);
           index2++;
         }
       }

       runModel(i,XTest,yTest,&XControl,&yControl);
     }
   // compute mean prediction error for each index
   STK::VectorX one(nbFolds_,1);
   cv_ = (measure_ * one) / nbFolds_;

   // compute mean standard deviation of cv_ for each index
   for(int i = 0; i < (int) index_.size(); i++)
     measure_.row(i) -= cv_[i];

   measure_ = measure_.square();
   cvError_ = (measure_ * one)/(nbFolds_-1)/nbFolds_;
   cvError_ = cvError_.sqrt();
}
}//end namespace HD
