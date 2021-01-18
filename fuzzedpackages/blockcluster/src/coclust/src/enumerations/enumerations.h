
/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2015  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>

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

    Contact : parmeet.bhatia@inria.fr , bhatia.parmeet@gmail.com
*/


/** @file enumerations.h
 *  @brief Defines all the enumerations used in CoClust Project.
 **/


#ifndef ENUMERATIONS_H_
#define ENUMERATIONS_H_

/** @brief Enumeration for Data-type. */
enum DataType
{
  binary_      = 1,
  contingency_ = 2,
  continuous_  = 3,
  categorical_ = 4
};

/**@ brief Enumeration for Algorithms.*/
enum Algorithm
{
  bem_    = 1,
  bcem_   = 2,
  bsem_   = 3,
  bgibbs_ = 4
};

/**@brief Enumeration for Stopping Criteria. */
enum StopCriteria
{
  parameter_  = 1,
  likelihood_ = 2
};

/** @brief Enumeration for Model Initialization. */
enum Initialization
{
  e_CEMInit_    = 1,
  e_EMInit_     = 2,
  e_RandomInit_ = 3
};

/** @brief Enumeration for all Models. */
enum Model
{
  pi_rho_epsilon_     = 1,
  pik_rhol_epsilon_   = 2,
  pi_rho_epsilonkl_   = 3,
  pik_rhol_epsilonkl_ = 4,
  pi_rho_unknown_     = 5,
  pik_rhol_unknown_   = 6,
  pi_rho_known_       = 7,
  pik_rhol_known_     = 8,
  pi_rho_sigma2_      = 9,
  pik_rhol_sigma2_    = 10,
  pi_rho_sigma2kl_    = 11,
  pik_rhol_sigma2kl_  = 12,
  pi_rho_multi_       = 13,
  pik_rhol_multi_     = 14
};

#endif /* ENUMERATIONS_H_ */
