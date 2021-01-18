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
 * created on: 30 mai 2013
 * Author:   Quentin Grimonprez
 **/

/** @file EM.h
 *  @brief In this file, definition of the EM class .
 **/


#ifndef EM_H_
#define EM_H_

#include "IAlgo.h"

namespace HD
{
  /**
   * This class runs an EM algorithm on a @c PenalizedModels object.
   * The stopping criterion is the convergence of the completed
   * loglikelihood or the number of iterations.
   */
  class EM : public IAlgo
  {
    public:
      /**
       * Constructor
       * @param maxStep maximal number of steps of the algorithm
       * @param eps threshold for convergence of the completed loglikelihood
       */
      EM( int maxStep = 1, int burn = 1, STK::Real eps=1e-5)
        : IAlgo(maxStep, eps), step_(0), burn_(burn), llc_(-std::numeric_limits<STK::Real>::max())
      {}
      /**@return the number of step of the algorithm*/
      inline int step() const { return step_;};
      /** @param burn new burn period */
      inline void setBurn(int burn) { burn_ = burn;}
      /** run the EM algorithm on a PenalizedModels object
       *  @param model pointer to a PenalizedModels object
       */
      template<class Model>
      bool run(PenalizedModels<Model>* model)
      {
#ifdef HD_DEBUG
            std::cout << "\n\nentering EM::run(model)"
                      << ", burn_= " << burn_
                      << ", eps_= "  << eps_
                      << ", maxStep_ = "<< maxStep_ << std::endl;
#endif
        try
        {
          //initialization
          step_ = 0;
          llc_ = -std::numeric_limits<STK::Real>::max();
          // we stop the burning after convergence of the completed log-likelihood
          // or after reaching the maximal number of burning step
          runBurn(model, burn_);
          int div = 0;
          while (step_ < maxStep_ && div < 2)
          {
            step_++;
            model->eStep(true);
            model->mStep(true);
            //difference between the log-likelihood of 2 successive steps.
            STK::Real llcOld=llc_;
            llc_ = model -> lnLikelihood();
            STK::Real diff = std::abs((llc_-llcOld)/llcOld);
#ifdef HD_DEBUG
            std::cout << "Selection step= " << step_
                      << ", llcOld= " << llcOld << ", llc_= "    << llc_
                      << ", diff ="   << diff <<std::endl;
#endif
            // convergence
            if (diff<eps_) break;
            // divergence
            if (llc_ < llcOld)
            {
              div++;
#ifdef HD_DEBUG
                std::cout << "Divergence in selection step. Trying Burning step..." << std::endl;
#endif
              if (!runBurn(model, 3)) break;
            }
          }
          step_++;
          model->eStep(true);
          model->mStep(true);
#ifdef HD_DEBUG
          //difference between the log-likelihood of 2 successive steps.
          STK::Real llcOld=llc_;
          llc_ = model -> lnLikelihood();
          STK::Real diff = std::abs((llc_-llcOld)/llcOld);
          std::cout << "Selection (final) step= " << step_
                    << ", llcOld= " << llcOld << ", llc_= "    << llc_
                    << ", diff ="   << diff <<std::endl;
#endif
        }
        catch(const STK::Exception& e)
        {
          msg_error_ = e.error();
#ifdef HD_DEBUG
          std::cout << "An error occur in EM::run(model). What: " << msg_error_ << std::endl;
#endif
          return false;
        }
        return true;
      }
      /** @return the last error message**/
      inline STK::String const& error() const { return msg_error_;}

    private:
      /** run the EM algorithm on a PenalizedModels object
       *  @param model pointer to a PenalizedModels object
       */
      template<class Model>
      bool runBurn(PenalizedModels<Model>* model, int burn)
      {
#ifdef HD_DEBUG
        std::cout << "\n----------------------------"<< std::endl;
        std::cout <<   "Entering runBurn with burn = " << burn << std::endl;
#endif
        // we stop the burning after convergence of the completed log-likelihood
        // or after reaching the maximal number of burning step
        // or in case of divergence
        int step = 0, div = 0;
        STK::Real diff;
        while (step < burn)
        {
          step_++; step++;
          model -> eStep(false);
          model -> mStep(false);
          //difference between the log-likelihood of 2 successive steps.
          STK::Real llcOld=llc_;
          llc_ = model -> lnLikelihood();
          // update divergence
          if (llc_ < llcOld) {div++;}
          diff = std::abs((llc_-llcOld)/(llcOld));
#ifdef HD_DEBUG
          std::cout << "runBurn step= " << step_
                    << ", llcOld= "  << llcOld
                    << ", llc_= "    << llc_
                    << ", diff ="    << diff
                    << ", div_= "    << div <<std::endl;
#endif
          if (diff<eps_ || div > 2) break;
        }
#ifdef HD_DEBUG
          std::cout << "runBurn terminated at step = " << step
                    << ", with llc_ = " << llc_ << std::endl;
          std::cout << "----------------------------"<< std::endl<< std::endl;
#endif
        // if we don't converge and get divergence, then return false
        return (div>2 && diff>eps_) ? false : true;
      }
      ///step of the algorithm
      int step_;
      ///burn period
      int burn_;
      /// current llc
      STK::Real llc_;
      ///last error message
      STK::String msg_error_;
  };

}

#endif /* EM_H_ */
