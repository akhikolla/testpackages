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
 * created on: 6 févr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file functions.cpp
 *  @brief Contains utilities functions for lars algorithm.
 **/

#include "../larsRmain.h"

using namespace std;
using namespace HD;
using namespace STK;

/*
 * (x1,y1) and (x2,y2) are points of an affine function, we wants to calculate the image of x3.
 * @param x1 abscissa of the first point
 * @param x2 abscissa of the second point
 * @param x3 abscissa of the point we want ordinate
 * @param y1 ordinate of the first point
 * @param y2 ordinate of the second point
 * @return ordinate of x3
 */
Real computeOrdinate(Real x1,Real x2,Real x3,Real y1,Real y2)
{
  return y1 + (y2-y1) * ((x3-x1)/(x2-x1)) ;
}

/*
 * Compute the coefficients for a given value of l1norm
 * Use with move()
 * @param state1 state of a lars step
 * @param state2 state of the next lars step
 * @param evolution difference between the 2 lars step
 * @param l1norm abscissa to compute ordinates
 * @return value of coefficients for l1norm
 */
//STK::Array1D< pair<int,Real> > computeCoefficients(PathState const& state1,PathState const& state2,pair<std::vector<int> ,std::vector<int> > const& evolution, Real const& l1norm)
//{
//  int maxSize = state1.size() + evolution.first.size();
//
//  STK::Array1D< pair<int,Real> > coeff(maxSize);
//
//  if(evolution.second.size() == 0)
//  {//no drop variable
//    int j(1);
//    for(j=1; j <= state1.size(); j++)
//      coeff[j]=make_pair(state1.varIdx(j),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(j), state2.varCoeff(j)));
//
//    //add variable case
//    if(evolution.first.size() != 0)
//    {
//      for(int i = 0; i < (int) evolution.first.size(); i++)
//        coeff[j]=make_pair( evolution.first[i], computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, 0., state2.varCoeff(j)));
//    }
//  }
//  else
//  {
//    //delete variable case
//    int i = 1;
//    for( int j = 0; j < (int) evolution.second.size(); j++)
//    {
//      //while we don't meet the delete variable, variable has the same index in the two sets
//      while(evolution.second[j]!=state1.varIdx(i))
//      {
//        coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i), state2.varCoeff(i)));
//        i++;
//      }
//      //compute coefficient for the delete variable
//      coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i),0.));
//      i++;
//    }
//    //compute coefficient for the other variable
//    while(i < state1.size()+1)
//    {
//      coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i), state2.varCoeff(i-1)));
//      i++;
//    }
//
//    //drop with an add variable
//    if(evolution.first.size()!=0)
//    {
//      for(int j = 0; j < (int) evolution.first.size(); j++)
//        coeff[i]=make_pair(evolution.first[j], computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, 0., state2.varCoeff(i-1)));
//    }
//
//  }
//
//  return coeff;
//}

//void computeCoefficients(PathState const& state1,PathState const& state2,pair<std::vector<int> ,std::vector<int> > const& evolution, Real const& l1norm, STK::Array2DVector< pair<int,Real> > &coeff)
//{
//  //STK::Array2DVector< pair<int,Real> > coeff(std::max(state1.size(),state2.size()));
//  int maxSize = state1.size() + evolution.first.size();
//
//  coeff.resize1D(Range(1,maxSize));
//
//  if(evolution.second.size()==0)
//  {//no drop variable
//    int j(1);
//    for(j = 1; j <= state1.size(); j++)
//      coeff[j]=make_pair(state1.varIdx(j),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(j), state2.varCoeff(j)));
//
//    //add variable case
//    if(evolution.first.size()!=0)
//    {
//      for(int i = 0; i < (int) evolution.first.size(); i++)
//        coeff[j]=make_pair( evolution.first[i], computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, 0., state2.varCoeff(j)));
//    }
//  }
//  else
//  {
//    //delete variable case
//    int i = 1;
//
//    for( int j = 0; j < (int) evolution.second.size(); j++)
//    {
//      //while we don't meet the delete variable, variable has the same index in the two sets
//      while(evolution.second[j]!=state1.varIdx(i))
//      {
//        coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i), state2.varCoeff(i)));
//        i++;
//      }
//
//      //compute coefficient for the delete variable
//      coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i),0.));
//      i++;
//    }
//
//
//    //compute coefficient for the other variable
//    while(i < state1.size()+1)
//    {
//      coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, state1.varCoeff(i), state2.varCoeff(i-1)));
//      i++;
//    }
//
//    //drop with an add variable
//    if(evolution.first.size()!=0)
//    {
//      for(int j = 0; j < (int) evolution.first.size(); j++)
//        coeff[i]=make_pair(evolution.first[j], computeOrdinate(state1.l1norm(), state2.l1norm(), l1norm, 0., state2.varCoeff(i-1)));
//
//    }
//
//  }
//}



void print(STK::Array2DVector< pair<int,Real> > const& state)
{
  for(int j=state.begin();j<state.end();j++)
    cout<<setw(8)<<state[j].first<<"  ";
  cout<<endl;
  for(int j=state.begin();j<state.end();j++)
    cout<<setw(8)<<state[j].second<<"  ";
  cout<<endl;
}


bool import(std::string adressFichier,int n,int p,STK::CArrayXX &data)
{
  std::ifstream flux(adressFichier.c_str());

  if (flux.is_open())//si le fichier est ouvert
  {
    Real real;
    int i=data.beginRows(),j=data.beginCols();
    while(flux>>real)//on lit le fichier données par données
    {
      data(i,j)=real;
      j++;
      if(j==data.endCols())//bout de ligne, on passe à la ligne
      {
        j=data.beginCols();
        i++;
      }
    }
    return true;
  }
  else
    return false;
}

bool import(std::string adressFichier,int n,STK::CVectorX &data)
{
  std::ifstream flux(adressFichier.c_str());

  if (flux.is_open())//si le fichier est ouvert
  {
    Real real;
    int i=data.begin();
    while(flux>>real)//on lit le fichier entier par entier
    {
      data[i]=real;
      i++;
    }
    return true;
  }
  else
    return false;
}
