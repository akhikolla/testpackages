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
 * created on: 23 mai 2013
 * Author:   Quentin Grimonprez
 **/

/** @file IPenalty.h
 *  @brief In this file, we define the interface class @c IPenalty.
 **/


#ifndef PENALTY_H_
#define PENALTY_H_

#include <RTKpp.h>

namespace HD
{
  /** Interface base class for penalty that can be applied to  a penalized
   *  regression model. This is essentially a functor that given a current
   *  estimation beta will compute the value of the penalty.
   **/
  class IPenalty
  {
    public:
      /**default constructor*/
      inline IPenalty() {}
      /**copy constructor*/
      inline IPenalty(IPenalty const& penalty) {}
      /**destructor*/
      inline virtual ~IPenalty() {}
      /**clone*/
      virtual IPenalty* clone() const = 0;
      /** update penalty
       *  @param beta current estimates
       */
      virtual void update(STK::VectorX const& beta) = 0;
      /** penalty term
       *  @param beta current estimates
       *  @return t(beta) * penalty * beta
       */
      virtual STK::Real penaltyTerm(STK::VectorX const& beta) const = 0;
      /** @return sigma2 */
      virtual STK::Real const& sigma2() const = 0;
  };
}


#endif /* PENALTY_H_ */
