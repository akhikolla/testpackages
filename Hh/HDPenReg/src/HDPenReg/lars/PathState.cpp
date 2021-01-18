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
 * created on: 5 févr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file PathState.cpp
 *  @brief Code of methods associates to @c PathState.
 **/

#include "../larsRmain.h"

using namespace STK;
using namespace std;

namespace HD
{
  //constructors
/* default constructor*/
PathState::PathState(): coefficients_(Range(1,0)), l1norm_(0.) {}

/*
 * constructor with reserve size
 * @param nbMaxVariable maximal number of variable to potentially stock
 */
PathState::PathState(int nbMaxVariable): coefficients_(Range(1,0)), l1norm_(0.)
{
  coefficients_.reserve(nbMaxVariable);
  coefficients_.shift(1); // in case
}

/* update coefficients with values specified in parameters
 * @param indexVariables vector containing the index of active variables
 * @param coefficients vector containing the value of estimates for active variables
 */
void PathState::update(VectorXi const& indexVariables, VectorX const& coefficients)
{
  //resize the container
  coefficients_.resize(coefficients.range());
  l1norm_=0;
  //fill the container and compute l1norm
  for(int i=coefficients_.begin(); i<coefficients_.end(); i++)
  {
    coefficients_[i]=make_pair(indexVariables[i], coefficients[i]);
    l1norm_ += std::abs(coefficients[i]);
  }
}

/*print coefficients*/
void PathState::printCoeff() const
{
  for(int i=coefficients_.begin(); i<coefficients_.end(); i++)
    std::cout << coefficients_[i].first<<"        ";
  std::cout<<std::endl;
  for(int i=coefficients_.begin(); i<coefficients_.end(); i++)
    std::cout << coefficients_[i].second<<" ";
  std::cout<<std::endl;
}

/*
 * update of the coefficients of the previous step
 * @param w direction of the update
 * @param gamma step of the update
 */
void PathState::update(CVectorX const& w, Real gamma)
{
  l1norm_=0;
  for(int i=coefficients_.begin(); i<coefficients_.end(); i++)
  {
    coefficients_[i].second += gamma*w[i];
    l1norm_ += std::abs(coefficients_[i].second);
  }
}

/*
 * update of the coefficients of the previous state with a variable to drop and a variable to add
 * @param w direction of the update
 * @param gamma step of the update
 * @param addIdxVar index of the variable to add
 * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
 */
void PathState::addWithDropUpdate(CVectorX const& w, Real gamma, vector<int> const& addIdxVar, vector<int> const& dropIdx)
{
  //update the other variable
  l1norm_=0;
//  if(dropIdx[0]!=1)
//  {
    for(int i=coefficients_.begin(); i < dropIdx[0]; i++)
    {
      coefficients_[i].second += gamma*w[i];
      l1norm_ += std::abs(coefficients_[i].second);
    }
//  }
  if(dropIdx.size() > 1)
  {
    for(int j = 0; j < (int) dropIdx.size(); j++)
    {
      for(int i = dropIdx[j]+1; i < dropIdx[j]; i++)
      {
        coefficients_[i].second += gamma*w[i];
        l1norm_ += std::abs(coefficients_[i].second);
      }
    }
  }
  if(dropIdx.back()!=coefficients_.lastIdx())
  {
    for(int i=dropIdx.back()+1; i < coefficients_.end(); i++)
    {
      coefficients_[i].second += gamma*w[i];
      l1norm_ += std::abs(coefficients_[i].second);
    }
  }
  //add the new variable
  for(int i = 0; i < (int) addIdxVar.size(); i++)
  {
    coefficients_.pushBack(1);
    coefficients_.back()=make_pair(addIdxVar[i],gamma*w[coefficients_.lastIdx()]);
    l1norm_ += std::abs(coefficients_.back().second);
  }
  //delete the variable to delete
  for(int i = dropIdx.size()-1; i >= 0; i--)
    coefficients_.erase(dropIdx[i]);
}

/*
 * update of the coefficients of the previous state with a variable to drop
 * @param w direction of the update
 * @param gamma step of the update
 * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
 */
void PathState::dropAfterDropUpdate(CVectorX const& w, Real gamma, vector<int> const& dropIdx)
{
  l1norm_=0;
  //update the other coefficient
  if(dropIdx[0]!=1)
  {
    for(int i = coefficients_.begin(); i < dropIdx[0]; i++)
    {
      coefficients_[i].second += gamma*w[i];
      l1norm_ += std::abs(coefficients_[i].second);
    }
  }
  if(dropIdx.size() > 1)
  {
    for(int j = 0; j < (int) dropIdx.size(); j++)
    {
      for(int i = dropIdx[j]+1; i < dropIdx[j]; i++)
      {
        coefficients_[i].second += gamma*w[i];
        l1norm_ += std::abs(coefficients_[i].second);
      }
    }
  }
  if(dropIdx.back()!=coefficients_.lastIdx())
  {
    for(int i = dropIdx.back() + 1; i < coefficients_.end(); i++)
    {
      coefficients_[i].second += gamma*w[i];
      l1norm_ += std::abs(coefficients_[i].second);
    }
  }
  //delete the variables to delete
  for(int i = dropIdx.size()-1; i >= 0; i--)
    coefficients_.erase(dropIdx[i]);
}

/*
 * update of the coefficients of the previous state with a new variable
 * @param w direction of the update
 * @param gamma step of the update
 * @param addIdxVar index of the variable to add
 */
void PathState::addUpdate(CVectorX const& w, Real gamma, vector<int> const& addIdxVar)
{
  //update previous coefficients
  l1norm_=0;
  for(int i=coefficients_.begin(); i<coefficients_.end(); i++)
  {
    coefficients_[i].second += gamma * w[i];
    l1norm_ += std::abs(coefficients_[i].second);
  }
  for(int i = 0; i < (int) addIdxVar.size(); i++)
  {
    coefficients_.pushBack(1);
    coefficients_.back() = make_pair(addIdxVar[i], gamma * w[coefficients_.lastIdx()]);
    l1norm_ += std::abs(coefficients_.back().second);
  }
}

}//end namespace



