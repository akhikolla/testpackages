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

/** @file functions.h
 *  @brief Contains utilities functions for lars algorithm.
 **/

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

/**
 * (x1,y1) and (x2,y2) are points of an affine function, we wants to calculate the image of x3.
 * @param x1 abscissa of the first point
 * @param x2 abscissa of the second point
 * @param x3 abscissa of the point we want ordinate
 * @param y1 ordinate of the first point
 * @param y2 ordinate of the second point
 * @return ordinate of x3
 */
STK::Real computeOrdinate(STK::Real x1,STK::Real x2,STK::Real x3,STK::Real y1,STK::Real y2);

///**
// * Compute the coefficients for a given value of lambda
// * @param state1 state of a lars step
// * @param state2 state of the next lars step
// * @param evolution difference between the 2 lars step
// * @param lambda abscissa to compute ordinates
// * @return value of coefficients for lambda
// */
//STK::Array1D< std::pair<int,STK::Real> > computeCoefficients(HD::PathState const& state1,HD::PathState const& state2,std::pair<std::vector<int> ,std::vector<int> > const& evolution, STK::Real const& lambda);
//void computeCoefficients(HD::PathState const& state1,HD::PathState const& state2,std::pair<std::vector<int> ,std::vector<int> > const& evolution, STK::Real const& lambda, STK::Array2DVector< std::pair<int,STK::Real> > &coeff);

void print(STK::Array2DVector< std::pair<int,STK::Real> > const& state);
bool import(std::string adressFichier,int n,int p,STK::CArrayXX &data);
bool import(std::string adressFichier,int n,STK::CVectorX &data);

template <typename T>
std::ostream& operator<<( std::ostream &flux, std::vector<T> vect)
{
  for (unsigned int i(0);i<vect.size();i++)
      flux << vect[i]<< " ";
  return flux;
}

#endif /* FUNCTIONS_H_ */
