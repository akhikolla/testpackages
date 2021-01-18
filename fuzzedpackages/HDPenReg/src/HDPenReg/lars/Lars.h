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
 * created on: 12 févr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Lars.h
 *  @brief In this file, we define the class @c Lars.
 **/


#ifndef LARS_H_
#define LARS_H_

/** @namespace HD
 *  The namespace HD contains all the lars class and methods.
 **/
namespace HD
{

/**
 * Class for running the LARS algorithm.
 *
 * Let \f$ X\f$ a matrix of size \f$ n\times p\f$, \f$ y\f$ a vector of length \f$ n\f$ and \f$\lambda\f$ a positive real.
 *
 * The lasso problem is to find \f$ \hat{\beta}\f$ such that
 *
 * \f$ \hat{\beta}=argmin_{\beta} \|\mathbf{y}-\mathbf{X}\mathbf{\beta}\|_2^2 + \lambda\|\beta\|_1\f$
 *
 *
 * The LARS algorithm solves the lasso problem for all values of lambda.
 *
 */
  class Lars
  {
    public:
      //constructors
      /**
       * Constructor
       * @param X matrix of data, a row=a individual
       * @param y response
       * @param intercept if true there is an intercept in the model
       */
      Lars(STK::CArrayXX const& X, STK::CVectorX const& y, bool intercept = true);
      /**
       * Constructor
       * @param X matrix of data, a row=a individual
       * @param y response
       * @param maxSteps number of maximum step to do
       * @param intercept if true there is an intercept in the model
       * @param eps epsilon (for 0)
       */
      Lars( STK::CArrayXX const& X
          , STK::CVectorX const& y
          , int maxSteps
          , bool intercept=true
          , STK::Real eps =STK::Arithmetic<STK::Real>::epsilon());

      //getters
      /**@return path of the coefficients*/
      inline Path const& path() const {return path_;}
      /**
       * @param i step
       * @return the Pathstate object : the state of the path at the step i
       */
      inline PathState const& path(int i) const {return path_.states(i);}
      /**
       * @param i index of the step
       * @param j index of the coefficients
       * @return the value of the j-th coefficient at the step i
       */
      inline STK::Real coefficient(int i,int j) const {return path_.varCoeff(i,j);}
      /**
       * @param i index of the step
       * @param j index of the coefficients
       * @return the value of the j-th coefficient at the step i
       */
      inline int varIdx(int i,int j) const {return path_.varIdx(i,j);}
      /**
       * @param i index of the step
       * @return the value of l1norm at the i-th step
       */
      inline STK::Real l1norm(int i) const {return path_.l1norm(i);}
      /** @return the vector of l1norm*/
      inline STK::VectorX const l1norm() const {return path_.l1norm();}
      /** @return the vector of lambda*/
      inline std::vector<STK::Real> const lambda() const {return path_.lambda();}
      /** @return the historic of add and drop variable*/
      inline std::vector< std::pair<std::vector<int> ,std::vector<int> > > evolution() const {return path_.evolution();}
      /**@return Number of step of the algorithm*/
      inline int step() const {return step_;}
      /** @return the intercept of the solution*/
      inline STK::Real mu() const {return mu_;}
      /** @return the ignored variable*/
      inline STK::CArrayVector<bool> toIgnore() const {return toIgnore_;}
      /** @return msg_error_*/
      inline std::string msg_error() const {return msg_error_;}
      /**@return muX_*/
      inline STK::CVectorX muX() const {return muX_;}
      /**@return muX_[i]*/
      inline STK::Real muX(int i) const {return muX_[i];}

      //methods
      /** run lars algorithm*/
      void run();

      /**
       * predict the path for a ratio fraction = l1norm/l1normmax
       * @param X new data for predict the response
       * @param index index (lambda or fraction) where the response is estimated.
       * @param lambdaMode if TRUE, index corresponds to a value of lambda, if FALSE, index is a real between 0 and 1
       * corresponding to ratio between the l1 norm of estimates to calculate and l1 norm max of solution
       * @param yPred container for the predicted response (will be modified)
       */
      void predict(STK::CArrayXX const& X, STK::Real index, bool lambdaMode, STK::CVectorX &yPred);

   protected:
      /**
       * initialization of algorithm
       */
      void initialization();

      /**
       * search non active variable with the greatest correlation
       * @param Cmax correlation max
       * @param newId a vector containing the index of variable to potentially add
       */
      void computeAddSet(STK::Real Cmax, std::vector<int>& newId) const;

      /**
       * update the QR decomposition of Xi
       * @param idxVar index of active variable to add
       * @param signC sign of correlation of active variable
       * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
       */
      void updateR(int idxVar, STK::VectorXi &signC, std::pair<bool,std::vector<int> > &action);

      /**
       * compute inv(Xi'*Xi)*1 from qr decomposition
       * @param Gi1 for stock inv(Xi'*Xi)*1
       * @param signC sign of correlation of active variable
       */
      void computeGi1(STK::CVectorX &Gi1, STK::VectorXi const& signC) const;

      /**
       * compute Cmax
       * @return Cmax the correlation max
       */
      STK::Real computeCmax();

      /**
       * add Cmax to the vector of lambda
       * @return Cmax the correlation max of the step
       */
      inline void addCmax(STK::Real const& Cmax) {path_.addLambda(Cmax);};

      /** dropStep
       * downdate qr decomposition,  X and signC
       * @param idxVar index of active variable to drop
       * @param signC sign of correlation of active variable
       */
      void dropStep(std::vector<int> const& idxVar, STK::VectorXi &signC);

      /**
       * Compute gammahat for the update of coefficient in add case
       * @param Aa norm of the inverse of G
       * @param a X' * equiangular vector
       * @param Cmax correlation max
       * @return gammaHat a real
       */
      STK::Real computeGamHat(STK::Real const& Aa, STK::CVectorX const& a, STK::Real Cmax) const;

      /**
       * Compute gammaTilde for the update of coefficient in drop case
       * @param w Aa*Gi1 @see computeGi1
       * @param idxMin we stock the index (in the activeVariable vector) of the variable with the min value
       * @return gammatilde a real
       */
      STK::Real computeGamTilde(STK::CVectorX const& w, std::vector<int> &idxMin) const;

      /**
       * Update the coefficient of the path
       * @param gamma gammaHat or gammaTilde
       * @param w Aa*Gi1 @see computeGi1
       * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
       * @param isAddCase true if we add a variable
       * @param dropId id top potentially drop
       */
      void updateBeta(STK::Real gamma, STK::CVectorX const& w, std::pair<bool,std::vector<int> > action, bool isAddCase, std::vector<int> dropId);

      /**
       * first step
       * @param Cmax correlation max
       * @param newId vector of index of active variable to add
       * @param signC sign of correlation of active variable
       * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
       * @param Aa norm of the inverse of G
       * @param Gi1 for stock inv(Xi'*Xi)*1
       * @param w Aa*Gi1
       * @param u unit vector making equal angles with the column of Xi
       * @param a X' * equiangular vector       * @param Gi1
       * @param gam the step for update coefficients
       * @return
       */
      bool firstStep( STK::Real &Cmax
                    , std::vector<int> &newId
                    , STK::VectorXi &signC
                    , std::pair<bool,std::vector<int> > &action
                    , STK::Real &Aa
                    , STK::CVectorX &Gi1, STK::CVectorX &w, STK::CVectorX &u, STK::CVectorX &a
                    , STK::Real &gam);

      /**
       * updateR only for the first step
       * @see updateR
       */
      void firstUpdateR(int idxVar, STK::VectorXi &signC, std::pair<bool,std::vector<int> > &action);
      /**
       * Compute the coefficients for a given value of lambda
       * @param state1 state of a lars step
       * @param state2 state of the next lars step
       * @param evolution difference between the 2 lars step
       * @param lambda abscissa to compute ordinates
       * @param coef value of coefficients for lambda
       */
      void computeCoefficients( PathState const& state1
                              , PathState const& state2
                              , std::pair<std::vector<int>
                              , std::vector<int> > const& evolution
                              , STK::Real const& lambda
                              , STK::Array2DVector< std::pair<int,STK::Real> > &coeff);
    private:
      ///number of individuals
      int n_;
      ///number of variables
      int p_;
      /// maximal number of steps
      int maxSteps_;
      ///covariate size n*p
      STK::CArrayXX X_;
      ///response size p*1
      STK::CVectorX y_;
      ///mean of each covariate of X
      STK::CVectorX muX_;
      /// path solution
      Path path_;
      /// current active variable (non zero coefficient)
      STK::CArrayVector<bool> isActive_;
      ///index of variables to ignore because it causes singularity
      STK::CArrayVector<bool> toIgnore_;
      ///number of active variables
      int nbActiveVariable_;
      ///number of ignored variables
      int nbIgnoreVariable_;
      ///index of active variables
      STK::VectorXi activeVariables_;
      ///current step
      int step_;
      /// Beta_0
      STK::Real mu_;
      ///eps for zero approximation
      STK::Real eps_;
      ///X reduced to covariates of active Set
      STK::ArrayXX Xi_;
      ///qr decomposition of Xi
      STK::lapack::Qr qrX_;
      /// vector of correlation (size p*1)
      STK::CVectorX c_;
      /// if true, there is an intercept in the model
      bool intercept_;
      ///last error message
      std::string msg_error_;
  };

}//end namespace

#endif /* LARS_H_ */
