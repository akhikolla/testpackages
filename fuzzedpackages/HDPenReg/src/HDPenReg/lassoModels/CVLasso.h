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
 * created on: 9 oct. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file CVLasso.h
 *  @brief In this file, definition of the class @c CVLasso derived from class @c CV .
 **/


#ifndef CVLASSO_H_
#define CVLASSO_H_

#include "IMeasure.h"
#include "CV.h"
#include "EM.h"
#include "Lasso.h"
#include "LogisticLasso.h"

namespace HD
{
/**
 * Class derived from @c CV, implementing the cross validation for @c Lasso with @c EM algorithm.
 * This class contains setters and implementation of the pure virtual runModel method from @c CV.
 */
  template<class LassoModel>
  class CVLasso : public CV
  {
    public:
      /**default constructor*/
      CVLasso() : CV(), eps_(1e-5), threshold_(1e-8), epsCG_(1e-8), maxStep_(1000), burn_(30), p_typeMeasure_(0){};

      /**initialize containers and class*/
      void initialize() {initializeCV();};

      /**set the epsilon for the convergence of the @c EM algorithm*/
      inline void setEps(STK::Real const& eps) {eps_ = eps;}
      /**set the maximum number of step of the @c EM algorithm*/
      inline void setMaxStep(int const& maxStep) {maxStep_ = maxStep;}
      /**set the number of burn steps of the @c EM algorithm*/
      inline void setBurn(int const& burn) {maxStep_ = burn;}
      /**set the epsilon for the convergene of the conjugate gradient (@c CG)*/
      inline void setEpsCG(STK::Real const& epsCG) {epsCG_ = epsCG;}
      /**set the threshold of the @c LassoSolver*/
      inline void setThreshold(STK::Real const& threshold) {threshold_ = threshold;}
      /**set the type of measure for evaluate the model*/
      inline void setTypeMeasure(IMeasure* p_typeMeasure) {p_typeMeasure_ = p_typeMeasure;}


    protected:
      /**
       * run lasso on all values of index for the given control data
       * @param i index of the actual delete folder
       * @param XTest data test
       * @param yTest response test
       * @param p_XControl pointer to the data control
       * @param p_yControl pointer to the response control
       */
      void runModel(int i, STK::ArrayXX const& XTest, STK::VectorX const& yTest, STK::ArrayXX const* p_XControl, STK::VectorX const* p_yControl)
      {
#ifdef HD_CVDEBUG
         std::cout << "Entering CVlasso::runModel with i=" << i << "\n";
#endif
        STK::VectorX yPred(sizePartition_[i]);
        //create em algorithm
        EM algo(maxStep_, burn_, eps_);
        //create model
        LassoModel lasso(p_XControl, p_yControl, index_[0], threshold_, epsCG_);
        //run the lasso on all value of index
        for(int s = 0 ; s < (int) index_.size(); s++)
        {
#ifdef HD_CVDEBUG
        	std::cout << "In CVlasso::runModel lauching algo with lambda= " << index_[s] <<"\n";
#endif
          //set the new value of lambda to test
          lasso.setLambda(index_[s]);
          //initialize the model
          lasso.initializeBeta();
          //run algorithm
          if (!algo.run(&lasso))
          {

      #ifdef HD_CVDEBUG
            std::cout << "\nIn CVlasso::runModel. An error occur in algo.run(&lasso).\nWhat: " << algo.error() << "\n";
      #endif
          }
          //we compute the prediction of the y associated to XTest
          yPred = XTest * lasso.beta();
          //compute the residuals
          measure_(s,i) = p_typeMeasure_->measure(yTest,yPred);
#ifdef HD_CVDEBUG
      std::cout << "yTest = " << yTest.transpose() << "\n";
      std::cout << "yPred = " << yPred.transpose() << "\n";
      std::cout << "measure_(" << s  <<"," << i << ") = "<< measure_(s,i) << "\n";
#endif
        }
#ifdef HD_CVDEBUG
         std::cout << "Terminating CVlasso::runModel\n";
#endif
      }

    private:
      /// eps for EM algorithm convergence
      STK::Real eps_;
      /// threshold for set value to 0 in LassoPenalty
      STK::Real threshold_;
      /// eps for conjugate gradient convergence
      STK::Real epsCG_;
      /// maximum number of step for EM algorithm
      int maxStep_;
      /// burn for EM algorithm
      int burn_;
      ///type of measure
      IMeasure* p_typeMeasure_;
  };


}


#endif /* CVLASSO_H_ */
