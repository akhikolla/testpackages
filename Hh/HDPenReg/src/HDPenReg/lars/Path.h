/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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
 * Project:  testC::
 * created on: 1 févr. 2013
 * Author:   grimonprez, serge.iovleff@stkpp.org
 **/

/** @file Path.h
 *  @brief In this file, we define the class @c Path.
 **/


#ifndef PATH_H_
#define PATH_H_

namespace HD
{
/**
 * Path solution of a lars algorithm. It contains the @c PathState at each step of the algorithm and the evolution
 * of the active variable (added or deleted variable) at each step, and the lambda parameter of the lars.
 */
  class Path
  {
    public:
      /** default constructor
       * @param maxSizePath maximal size of the path
       */
      Path(int maxSizePath);

      //getters
      /**@return lambda_*/
      inline std::vector< STK::Real > const& lambda() const {return lambda_;}
      /**@return lambda_[i]*/
      inline STK::Real const lambda(int i) const {return lambda_[i];}
      /**@return states_*/
      inline std::vector< PathState > const& states() const {return states_;}
      /**@return states_[i]*/
      inline PathState const& states(int i) const {return states_[i];}
      /**@return last state*/
      inline STK::Array1D< std::pair<int,STK::Real> > const& lastState() const
      { return states_.back().coefficients();}
      /**@return l1norm of state i*/
      inline STK::Real const l1norm(int i) const {return states_[i].l1norm();}
      /**@return l1norm*/
      inline STK::VectorX const l1norm() const
      {
        STK::VectorX l1norm(states_.size());
        for(int i = 0, il1norm=l1norm.begin(); i < l1norm.end(); i++, il1norm++)
        { l1norm[il1norm] = states_[i].l1norm();}
        return l1norm;
      }
      /**@return coefficient j of state i*/
      inline STK::Real varCoeff(int i,int j) const {return states_[i].varCoeff(j);}
      /**@return index variable j of state i*/
      inline int varIdx(int i, int j) const {return states_[i].varIdx(j);}
      /**@return last coefficients*/
      inline STK::Real lastVarCoeff(int i) const {return states_.back().varCoeff(i);}
      /**@return last index variable*/
      inline int lastVarIdx(int i) const {return states_.back().varIdx(i);}
      /**@return last step change*/
      inline std::pair<std::vector<int> ,std::vector<int> > const& lastStep() const {return evolution_.back();}
      /**@return size of path*/
      inline int size() const {return states_.size();}
      /**@return size of path*/
      inline int size(int i) const {return states_[i].size();}
      /** print states_[i]*/
      inline void print(int i) const {states_[i].printCoeff();}
      /** print evolution_*/
      //inline void printEvolution() const {for(std::vector< std::pair<int,std::vector<int> > >::const_iterator it=evolution_.begin();it!=evolution_.end();it++){std::cout<<(*it).first<<"   "<<(*it).second<<std::endl;}}
      /** @return evolution_*/
      inline std::vector< std::pair<std::vector<int> ,std::vector<int> > > evolution() const {return evolution_;}
      /** @return evolution_[i]*/
      inline std::pair<std::vector<int> ,std::vector<int> > evolution(int i) const {return evolution_[i];}

      //methods
      /** @brief get coefficient associates to a lambda value.
       *  @param lambda is the norm value for which we want values of coefficient
       *  @return a vector containing pair<int,double>=(index of non zero coefficient,coefficient)
       */
      STK::Array1D< std::pair<int,STK::Real> >  coeff(STK::Real lambda) const;

      /** @brief Add coefficients of a LARS step to the actual path
       * @param indexVariables Array2DVector containing the index of active variables
       * @param coefficients Array2DVector containing the value of estimates for the active variables
       * @param idxVarAdd index of the new variable (0 if no new variable)
       * @param idxVarDrop index of the delete variable variable (0 if no delete variable)
       */
      void addCoeff(STK::VectorXi const& indexVariables,STK::VectorX const& coefficients,int idxVarAdd,int idxVarDrop);

      /**
       * update of the coefficients of the previous state with a new variable
       * @param w direction of the update
       * @param gamma step of the update
       * @param addIdxVar index of the variable to add
       */
      void addCaseUpdate(STK::Real gamma, STK::CVectorX const& w, std::vector<int> const& addIdxVar);

      /**
       * update of the coefficients of the previous step
       * @param w direction of the update
       * @param gamma step of the update
       */
      void update(STK::Real gamma, STK::CVectorX const& w);

      /**
       * update of the coefficients of the previous state with a variable to drop and a variable to add
       * @param w direction of the update
       * @param gamma step of the update
       * @param addIdxVar index of the variable to add
       * @param dropIdxVar index of the delete variable
       * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
       */
      void addWithDropCaseUpdate(STK::Real gamma, STK::CVectorX const& w, std::vector<int> const& addIdxVar, std::vector<int> const& dropIdxVar, std::vector<int> const& dropIdx);

      /**
       * update of the coefficients of the previous state with a variable to drop
       * @param w direction of the update
       * @param gamma step of the update
       * @param dropIdxVar index of the delete variable
       * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
       */
      void dropAfterDropCaseUpdate(STK::Real gamma, STK::CVectorX const& w, std::vector<int> const& dropIdxVar, std::vector<int> const& dropIdx);

      /**
       * add an element at the end of the vector of the correlation max
       * @param lambda correlation max to add at the end of the vector
       */
      void addLambda(STK::Real const& lambda);

      /**
       * Compute the coefficients for a given value of lambda
       * @param state1 state of a lars step
       * @param state2 state of the next lars step
       * @param evolution difference between the 2 lars step
       * @param lambda abscissa to compute ordinates
       * @return value of coefficients for lambda
       */
      STK::Array1D< std::pair<int,STK::Real> > computeCoefficients( PathState const& state1
                                                                  , PathState const& state2
                                                                  , std::pair<std::vector<int>
                                                                  , std::vector<int> > const& evolution
                                                                  , STK::Real const& lambda) const;

    private:
      ///path solution
      std::vector< PathState > states_;
      ///the first index is the index of the add variable, the second the index of the drop variable
      std::vector< std::pair<std::vector<int> ,std::vector<int> > > evolution_;
      ///vector containing the lambda (correlation max at each step)
      std::vector< STK::Real > lambda_;
  };
}

#endif /* PATH_H_ */
