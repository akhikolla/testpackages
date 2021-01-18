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
 * Project:  MPAGenomics::
 * created on: 4 févr. 2013
 * Author:   grimonprez, serge.iovleff@stkpp.org
 **/

/** @file Path.cpp
 *  @brief In this file, implementation of the methods of @c Path .
 **/

#include "../larsRmain.h"

using namespace STK;
using namespace std;

namespace HD
{
//Constructors
Path::Path(int maxSizePath)
{
  states_.reserve(maxSizePath);
  evolution_.reserve(maxSizePath);
  states_.push_back(PathState());
  lambda_.reserve(maxSizePath);
}

  //Methods
/* @brief Add coefficients of a LARS step to the actual path
 * @param indexVariables Array2DVector containing the index of active variables
 * @param coefficients Array2DVector containing the value of estimates for the active variables
 * @param idxVarAdd index of the new variable (0 if no new variable)
 * @param idxVarDrop index of the delete variable variable (0 if no delete variable)
 */
void Path::addCoeff(VectorXi const& indexVariables,VectorX const& coefficients,int idxVarAdd,int idxVarDrop)
{
  states_.push_back(PathState());
  states_.back().update(indexVariables,coefficients);
  vector<int> idxVarAddVect(1,idxVarAdd), idxVarDropVect(1,idxVarDrop);
  evolution_.push_back(make_pair(idxVarAddVect,idxVarDropVect));
}

/* @brief get coefficient associates to a l1norm value.
 *  @param l1norm is the norm value for which we want values of coefficient
 *  @return a vector containing pair<int,double>=(index of non zero coefficient,coefficient)
 */
Array1D< pair<int,Real> > Path::coeff(Real l1norm) const
{
  Array1D< pair<int,Real> > coeff;
  if(l1norm!=0)
  {
    //search off the interval of m_l1norm containing l1norm
    int indexl1norm=0;
    while( ( states_[indexl1norm].l1norm() < l1norm ) && ( indexl1norm < (int)states_.size()-1 ) )
    { indexl1norm++;}
    if(l1norm==states_[indexl1norm].l1norm())//if l1norm is equal to an actual l1norm
    {  coeff=states_[indexl1norm].coefficients();}
    else
    {
      if( indexl1norm == (int) states_.size()-1 )//last m_l1norm, we return the last coefficient
        coeff=states_[indexl1norm].coefficients();
      else
        coeff=computeCoefficients( states_[indexl1norm-1]
                                 , states_[indexl1norm]
                                 , evolution_[indexl1norm-1]
                                 , l1norm);
    }
  }
  return coeff;
}

/*
 * update of the coefficients of the previous state with a new variable
 * @param w direction of the update
 * @param gamma step of the update
 * @param addIdxVar index of the variable to add
 */
void Path::addCaseUpdate(Real gamma, CVectorX const& w, std::vector<int> const& addIdxVar)
{
  //copy the previous states
  states_.push_back(states_.back());
  //update of evolution with the new index
  vector<int> vide;
  evolution_.push_back(make_pair(addIdxVar,vide));
  //update of the coefficients
  states_.back().addUpdate(w,gamma,addIdxVar);
}

/*
 * update of the coefficients of the previous step
 * @param w direction of the update
 * @param gamma step of the update
 */
void Path::update(Real gamma, CVectorX const& w)
{
  //copy the previous states
  states_.push_back(states_.back());
  //update of evolution with the new index
  vector<int> vide;
  evolution_.push_back(make_pair(vide,vide));
  //update of the coefficients
  states_.back().update(w,gamma);
}

/*
 * update of the coefficients of the previous state with a variable to drop and a variable to add
 * @param w direction of the update
 * @param gamma step of the update
 * @param addIdxVar index of the variable to add
 * @param dropIdxVar index of the delete variable
 * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
 */
void Path::addWithDropCaseUpdate(Real gamma, CVectorX const& w, std::vector<int> const& addIdxVar, std::vector<int> const& dropIdxVar, std::vector<int> const& dropIdx)
{
  //copy the previous states
  states_.push_back(states_.back());
  //update of evolution with the new index
  evolution_.push_back(make_pair(addIdxVar,dropIdxVar));
  //update of the coefficients
  states_.back().addWithDropUpdate(w, gamma, addIdxVar, dropIdx);
}

/*
 * update of the coefficients of the previous state with a variable to drop
 * @param w direction of the update
 * @param gamma step of the update
 * @param dropIdxVar index of the delete variable
 * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
 */
void Path::dropAfterDropCaseUpdate(Real gamma, CVectorX const& w, std::vector<int> const& dropIdxVar, std::vector<int> const& dropIdx)
{
  //copy the previous states
  states_.push_back(states_.back());
  //update of evolution with the new index
  vector<int> vide;
  evolution_.push_back(make_pair(vide,dropIdxVar));
  //update of the coefficients
  states_.back().dropAfterDropUpdate(w,gamma,dropIdx);
}


/*
 * add an element at the end of the vector of the correlation max
 * @param lambda correlation max to add at the end of the vector
 */
void Path::addLambda(Real const& lambda)
{ lambda_.push_back(lambda);}

Array1D< pair<int,Real> > Path::computeCoefficients( PathState const& state1
                                                   , PathState const& state2
                                                   , pair<std::vector<int>, std::vector<int> > const& evolution
                                                   , Real const& l1norm) const
{
  int maxSize = state1.size() + evolution.first.size();
  Array1D< pair<int,Real> > coeff(Range(1,maxSize));
  if(evolution.second.size() == 0)
  {//no drop variable
    int j;
    for(j=coeff.begin(); j <= state1.size(); j++)
      coeff[j]=make_pair(state1.varIdx(j),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(j), state2.varCoeff(j)));
    //add variable case
    if(evolution.first.size() != 0)
    {
      for(int i = coeff.begin(); i < (int) evolution.first.size(); i++)
        coeff[j]=make_pair( evolution.first[i], computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, 0., state2.varCoeff(j)));
    }
  }
  else
  {
    //delete variable case
    int i = coeff.begin();
    for( int j = 0; j < (int) evolution.second.size(); j++)
    {
      //while we don't meet the delete variable, variable has the same index in the two sets
      while(evolution.second[j]!=state1.varIdx(i))
      {
        coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i), state2.varCoeff(i)));
        i++;
      }
      //compute coefficient for the delete variable
      coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i),0.));
      i++;
    }
    //compute coefficient for the other variable
    while(i < state1.size()+1)
    {
      coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i), state2.varCoeff(i-1)));
      i++;
    }
    //drop with an add variable
    if(evolution.first.size()!=0)
    {
      for(int j = 0; j < (int) evolution.first.size(); j++)
        coeff[i]=make_pair(evolution.first[j], computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, 0., state2.varCoeff(i-1)));
    }
  }
  return coeff;
}


}//end namespace




