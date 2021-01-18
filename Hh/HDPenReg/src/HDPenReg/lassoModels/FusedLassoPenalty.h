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
 * created on: 21 juin 2013
 * Author:   Quentin Grimonprez
 **/

/** @file FusedLassoPenalty.h
 *  @brief In this file, definition of the @c FusedLassoPenalty class .
 **/


#ifndef FUSEDLASSOPENALTY_H_
#define FUSEDLASSOPENALTY_H_

#include "PenalizedModels.h"
#include "IPenalty.h"

namespace HD
{
/**
 * FusedLassoMultiplicator associated to the FusedLassoPenalty for the @c CG.
 */
  struct FusedLassoMultiplicator
  {
      /** constructor
       *  @param p_x pointer to the current data
       *  @param p_tridiag pointer to the penalty matrix
       *  @param p_sigma2 pointer to the sigma2 value
       */
      FusedLassoMultiplicator( STK::ArrayXX const* p_x = 0
                             , STK::VectorX const* p_mainDiagonal = 0
                             , STK::VectorX const* p_offDiagonal = 0
                             , STK::Real const* p_sigma2 = 0)
                             : p_x_(p_x)
                             , p_mainDiagonal_(p_mainDiagonal), p_offDiagonal_(p_offDiagonal)
                             , p_sigma2_(p_sigma2)
      {}
      /**
       * multiplicator for the CG in the mstep
       * @param x a vector
       * @return the product M*x
       */
      STK::VectorX operator()(STK::VectorX &x) const
      {       //a = (sigma2*B+tX*X)*x
#ifdef HD_DEBUG
//   std::cout << "x.range() =" << x.range() << "\n";
//   std::cout << "tridiagMult(x).range() =" << tridiagMult(x).range() << "\n";
//   std::cout << "(*p_x_).rows() =" << (*p_x_).rows() << "\n";
#endif
        STK::VectorX a = tridiagMult(x) + (p_x_->transpose() * ((*p_x_) * x));
        return   a ;
      }
      /**
       * Product of a the tridiagonal matrix penalty with a vector x
       * @param x vector
       * @return the product penalty matrix * x
       */
      STK::VectorX tridiagMult(STK::VectorX &x) const
      {
        STK::VectorX a(x.range());
#ifdef HD_DEBUG
        //check the size
        if(x.range() != p_mainDiagonal_->range())
        {
          std::cerr << "In FusedLassoMultiplicator::tridiagMult\n"
                    << "x.range() != p_mainDiagonal_->range()\n";
          throw(STK::out_of_range("size mismatch."));
        }
#endif
        if(x.size() > 0)
        {
          //matrix of size 1, no offDiagonal
          if(x.size() == 1)
          {  a.front() = p_mainDiagonal_->front() * x.front();}
          else
          {
            //first element of the product
            a.front() = (*p_mainDiagonal_).front() * x.front() +  (*p_offDiagonal_).front() * x[x.begin()+1];
            // intermediate elements
            for(int i = x.begin()+1; i < x.lastIdx(); i++)
              a[i] =(*p_offDiagonal_)[i-1] * x[i-1] + (*p_mainDiagonal_)[i] * x[i] +  (*p_offDiagonal_)[i] * x[i+1] ;
            //last element of the product
            a.back() = x.back() * (*p_mainDiagonal_).back() + x[x.lastIdx()-1] * (*p_offDiagonal_)[x.lastIdx()-1] ;
          }
        }
        return a;
      }
      ///pointer to the current data
      STK::ArrayXX const* p_x_;
      ///pointer to the current data
      STK::VectorX const* p_mainDiagonal_;
      ///pointer to the current data
      STK::VectorX const* p_offDiagonal_;
      ///pointer to the sigma2 value
      STK::Real const* p_sigma2_;
  };


  /**
   * Class FusedLassoPenalty derived from @c IPenalty
   * This class contains the penalty term of a fused lasso in a EM algorithm and the way to update the penalty
   *
   */
  class FusedLassoPenalty : public IPenalty
  {
    public:
      typedef STK::Array2DVector<STK::Range> Segment;
      /** Constructor
       *  @param lambda1 penalization parameter for the l1-norm of the estimates
       *  @param lambda2 penalization parameter for the l1-norm of the difference between successive estimates
       *  @param eps epsilon to add to denominator of fraction to avoid zeros.
       */
      FusedLassoPenalty(STK::Real lambda1, STK::Real lambda2, STK::Real eps = 1e-8);
      /** Copy constructor
       *  @param penalty LassoPenalty object to copy
       */
      FusedLassoPenalty(FusedLassoPenalty const& penalty);
      /** destructor */
      virtual ~FusedLassoPenalty() {};
      /**clone*/
      FusedLassoPenalty* clone() const ;

      //getter
      /**@return lambda1_ parameter of the lasso penalty */
      inline STK::Real const& lambda1() const {return lambda1_;}
      /**@return lambda2_ parameter of the fusion penalty */
      inline STK::Real const& lambda2() const {return lambda2_;}
      /**@return sigma2 variance of the response*/
      inline STK::Real const& sigma2() const { return sigma2_;}
      /**@return a pointer to sigma2*/
      inline STK::Real const*  p_sigma2() const { return &sigma2_;}
      /**@return */
      inline STK::Real const& eps() const {return eps_;}
      /**@return the main diagonal of the penalty matrix.*/
      inline STK::VectorX const& mainDiagonal() const {return mainDiagonal_;}
      /**@return a pointer to the main diagonal of the penalty matrix.*/
      inline STK::VectorX const* p_mainDiagonal() const {return &mainDiagonal_;}
      /**@return the off diagonal of the penalty matrix.*/
      inline STK::VectorX const& offDiagonal() const {return offDiagonal_;}
      /**@return a pointer to the off diagonal of the penalty matrix.*/
      inline STK::VectorX const* p_offDiagonal() const {return &offDiagonal_;}

      //setter
      /**
       * Set the parameter of the penalty value.
       * @param lambda1 new value of lambda1
       */
      inline void setLambda1(STK::Real const& lambda1) {lambda1_ = lambda1;}
      /**
       * Set the parameter of the penalty value.
       * @param lambda2 new value of lambda2
       */
      inline void setLambda2(STK::Real const& lambda2) {lambda2_ = lambda2;}
      /**
       * Set eps
       * @param eps epsilon to add to the denominator
       */
      inline void setEps(STK::Real const& eps) {eps_ = eps;}

      /** update sigma2 and the fused lasso penalty
       *  @param beta current estimates
       *  @param normResidual ||y-X*beta||_2^2
       */
      //void update(STK::VectorX const& beta, STK::Real const& normResidual){};
      virtual void update(STK::VectorX const& beta);
      virtual void update(STK::VectorX const& beta, Segment const& segment);
      /** @param beta current estimates
       *  @return t(beta) * matrixB * beta
       */
      STK::Real penaltyTerm(STK::VectorX const& beta) const;

    protected:
      /** update sigma2_
       *  @param beta current estimates
       *  @param normResidual ||y-X*beta||_2^2
       */
      //void updateSigma2(STK::VectorX const& beta, STK::Real const& normResidual);
      /** initialization in the constructor */
      //void initialization();

    private:
      /// parameter associated with the l1 norm of estimates
      STK::Real lambda1_;
      /// parameter associated with the l1 norm of the difference of successive estimates
      STK::Real lambda2_;
      /// main diagonal of the penalty matrix
      STK::VectorX mainDiagonal_;
      /// off diagonal of the penalty matrix
      STK::VectorX  offDiagonal_;
      /// variance
      STK::Real sigma2_;
      /// eps to add to the denominator in order to avoid 0
      STK::Real eps_;
  };
}

#endif /* FUSEDLASSOPENALTY_H_ */
