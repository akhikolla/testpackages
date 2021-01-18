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

/** @file CVFusedLasso.h
 *  @brief In this file, definition of class @c FusedLasso1D and @c FusedLasso2D derived from @c CV .
 **/


#ifndef CVFUSEDLASSO_H_
#define CVFUSEDLASSO_H_


#include "IMeasure.h"
#include "CV.h"
#include "EM.h"
#include "FusedLasso.h"
#include "LogisticFusedLasso.h"

namespace HD
{
/**
 * Cross validation for lasso with EM algorithm.
 * Optimization of either lambda1 or lambda2
 */
  template<class FusedLassoModel>
  class CVFusedLasso1D : public CV
  {
    public:
      /**default constructor*/
      CVFusedLasso1D(): CV(), lambda_(1.), optimL1_(false), eps_(1e-6),threshold_(1e-10), epsCG_(1e-8), maxStep_(1000), burn_(30), p_typeMeasure_(0) {};

      /**set optimL1 parameter
       * @param optimL1 If true, optimization of lambda1
       */
      inline void setOptimL1(bool optimL1) {optimL1_ = optimL1;}
      /** set lambda_
       * @param lambda Value of the fixed parameter (lambda1 if optimL1=false)
       */
      inline void setLambda(STK::Real lambda) {lambda_ = lambda;}
      /**set the eps for th EM
       * @param eps eps for the convergence of the em
       */
      inline void setEps(STK::Real const& eps) {eps_ = eps;}
      /**set the eps for the CG
       * @param eps eps for the convergence of the CG
       */
      inline void setEpsCG(STK::Real const& epsCG) {epsCG_ = epsCG;}
      /**set the threshold
       * @param threshold threshold for 0
       */
      inline void setThreshold(STK::Real const& threshold) {threshold_ = threshold;}
      /**set maxStep
       * @param maxStep maximum number of step of the em
       */
      inline void setMaxStep(int const& maxStep) {maxStep_ = maxStep;}
      /**set the burn
       * @param burn number of burn step for the em
       */
      inline void setBurn(int const& burn) {maxStep_ = burn;}
      /**set the type of measure for evaluate the model*/
      inline void setTypeMeasure(IMeasure* p_typeMeasure) {p_typeMeasure_ = p_typeMeasure;}

      /**initialize CV*/
      void initialize() {initializeCV();};

    protected:
      void runModel(int i, STK::ArrayXX const& XTest, STK::VectorX const& yTest, STK::ArrayXX const* p_XControl, STK::VectorX const* p_yControl)
      {
        STK::VectorX yPred(sizePartition_[i] );

        //create model
        FusedLassoModel fusedlasso(p_XControl, p_yControl,index_[0],lambda_,threshold_,epsCG_);
        fusedlasso.setEps(threshold_);
        EM algo(maxStep_,burn_,eps_);
        //run the model for all values of lambda
        for(int s = 0 ; s < (int) index_.size(); s++)
        {
          //set the value of lambda to test
          if(optimL1_)
            fusedlasso.setLambda1(index_[s]);
          else
            fusedlasso.setLambda2(index_[s]);

          //initialize the model
          fusedlasso.initializeModel();
          //run the algo
          algo.run(&fusedlasso);

          //we compute the prediction of the y associated to XTest
          STK::VectorX yPred = XTest * fusedlasso.beta();
          //compute the residuals
//          stk_cout<<(fusedlasso.p_penalty())->lambda1()<<"   "<<(fusedlasso.p_penalty())->lambda2()<<"  "<<index_[s-1]<<"  "<<i+1<<"   "<<(yPred-yTest).square().sum()/sizePartition_[i]<<std::endl;
          measure_(s,i) = p_typeMeasure_->measure(yTest,yPred);

        }
      }

//      void runModel(int i, STK::ArrayXX const& XTest, STK::VectorX const& yTest, STK::ArrayXX const* p_XControl, STK::VectorX const* p_yControl)
//      {
//        #pragma omp parallel
//        {
//          #pragma omp for schedule(dynamic,1)
//          for(int s = 1 ; s <= (int) index_.size(); s++)
//          {
//            STK::VectorX yPred(sizePartition_[i] );
//            FusedLasso fusedlasso;
//            fusedlasso.setP_x(p_XControl);
//            fusedlasso.setP_y(p_yControl);
//            fusedlasso.setLambda1(lambda_);
//            fusedlasso.setLambda2(lambda_);
//            fusedlasso.setEps(1e-12);
//            fusedlasso.setThreshold(1e-8);
//
//            EM p_algo2_(200,10,1e-5);
//            if(optimL1_)
//              fusedlasso.setLambda1(index_[s-1]);
//            else
//              fusedlasso.setLambda2(index_[s-1]);
//
//            fusedlasso.initializeModel();
//            p_algo2_.run(&fusedlasso);
//
//            //we compute the prediction of the y associated to XTest
//            yPred = XTest * fusedlasso.beta();
//            //compute the residuals
//            residuals_(s,i+1) = (yPred-yTest).square().sum()/sizePartition_[i];
//  //          stk_cout<<residuals_(s,i+1)<<"     "<<p_fusedLasso_->p_penalty()->lambda1()<<"  sol:  "<<p_fusedLasso_->beta()<<std::endl;
//
//          }
//        }
//      }


    private:
      /// value of the lambda to not optimize
      STK::Real lambda_;
      /// true if we optimize lambda1
      bool optimL1_;
      /// eps for EM algorithm convergence
      STK::Real eps_;
      /// threshold for set value to 0 in FusedLassoSolver
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

  /**
   * Cross validation for lasso with EM algorithm.
   * Optimization on a grid of values for lambda1 and lambda2
   */
  template<class FusedLassoModel>
  class CVFusedLasso2D : public CV
  {
    public:
      /**default constructor*/
      CVFusedLasso2D() : CV(), indexL2_(), eps_(1e-6),threshold_(1e-10), epsCG_(1e-8), maxStep_(1000), burn_(30), p_typeMeasure_(0){};

      /**set the eps for th EM
       * @param eps eps for the convergence of the em
       */
      inline void setEps(STK::Real const& eps) {eps_ = eps;}
      /**set the eps for the CG
       * @param eps eps for the convergence of the CG
       */
      inline void setEpsCG(STK::Real const& epsCG) {epsCG_ = epsCG;}
      /**set the threshold
       * @param threshold threshold for 0
       */
      inline void setThreshold(STK::Real const& threshold) {threshold_ = threshold;}
      /**set maxStep
       * @param maxStep maximum number of step of the em
       */
      inline void setMaxStep(int const& maxStep) {maxStep_ = maxStep;}
      /**set the burn
       * @param burn number of burn step for the em
       */
      inline void setBurn(int const& burn) {maxStep_ = burn;}
      /** set indexL2
       * @param indexL2 values of lambda2 to test
       */
      inline void setIndexL2(std::vector<STK::Real> const& indexL2) {indexL2_ = indexL2;}
      /**set the type of measure for evaluate the model*/
      inline void setTypeMeasure(IMeasure* p_typeMeasure) {p_typeMeasure_ = p_typeMeasure;}

      /**initialize class and containers*/
      void initialize()
      {
        initializeCV();
        measure_.resize(index_.size() * indexL2_.size(), nbFolds_);
        cv_.resize(index_.size() * indexL2_.size());
        cvError_.resize(index_.size() * indexL2_.size());
      };

    protected:

      void runModel( int i, STK::ArrayXX const& XTest, STK::VectorX const& yTest
                   , STK::ArrayXX const* p_XControl, STK::VectorX const* p_yControl)
      {
#ifdef HD_CVDEBUG
         std::cerr << "Entering CVFusedlasso::runModel with i=" << i << "\n";
#endif
        //create model
        FusedLassoModel fusedlasso(p_XControl, p_yControl, index_[0], indexL2_[0], threshold_, epsCG_);
        fusedlasso.setEps(threshold_);
        EM algo(maxStep_,burn_,eps_);
        //run on all the grid of lambda1 and lambda2
        for(int s = 0 ; s < (int) index_.size(); s++)
        {
          //set lambda1
          fusedlasso.setLambda1(index_[s]);
          //run for all lambda2
          for(int j = 0; j < (int) indexL2_.size(); j++)
          {
#ifdef HD_CVDEBUG
            std::cerr << "---------------------------------------------\n";
            std::cerr << "In CVlasso::runModel lauching algo with lambda1= "
                      << index_[s] << ", lambda2=" << indexL2_[j] <<"\n";
#endif
            //set the values of lambda2
            fusedlasso.setLambda2(indexL2_[j]);
            //initialize the model
            fusedlasso.initializeModel();
            //run algorithm
            if (!algo.run(&fusedlasso))
            {

        #ifdef HD_CVDEBUG
              std::cout << "\nIn CVFusedlasso::runModel. An error occur in algo.run(&lasso).\nWhat: " << algo.error() << "\n";
        #endif
            }
#ifdef HD_CVDEBUG
            std::cerr << "Algo terminated.\n";
#endif
           //we compute the prediction of the y associated to XTest
            STK::VectorX yPred = XTest * fusedlasso.beta();
            //compute the residuals
            measure_(s*indexL2_.size() + j, i) = p_typeMeasure_->measure(yTest,yPred);
#ifdef HD_CVDEBUG
      std::cout << "measure_(" << s*indexL2_.size() + j  <<"," << i << ") = "<< measure_(s*indexL2_.size() + j, i) << "\n";
#endif
          }
        }
      }

    private:
      ///values to test for lambda2
      std::vector<STK::Real> indexL2_;
      /// eps for EM algorithm convergence
      STK::Real eps_;
      /// threshold for set value to 0 in FusedLassoSolver
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


#endif /* CVFUSEDLASSO_H_ */
