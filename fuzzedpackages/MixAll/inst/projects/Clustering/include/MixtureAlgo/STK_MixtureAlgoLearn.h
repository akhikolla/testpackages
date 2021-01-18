/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_MixtureAlgoLearn.h
 *  @brief In this file we define learning mixture algorithms
 **/

#ifndef STK_MIXTUREALGOLEARN_H
#define STK_MIXTUREALGOLEARN_H

#include "STK_IMixtureAlgoLearn.h"

namespace STK
{

/** @ingroup Clustering
 *  @brief Implementation of the ImputeAlgo learning algorithm.
 *  The Impute algorithm call alternatively the steps:
 *  - imputeStep()
 *  - paramUpdateStep()
 *  until the maximum number of iterations is reached or the variation of the
 *  log-likelihood is less than the tolerance.
 **/
class ImputeAlgo: public IMixtureAlgoLearn
{
  public:
    /** default constructor */
    inline ImputeAlgo(): IMixtureAlgoLearn() {}
    /** Copy constructor.
     *  @param algo the algorithm to copy */
    inline ImputeAlgo( ImputeAlgo const& algo): IMixtureAlgoLearn(algo) {}
    /** destructor */
    inline virtual ~ImputeAlgo(){}
    /** clone pattern */
    inline virtual ImputeAlgo* clone() const { return new ImputeAlgo(*this);}
    /** run the algorithm on the model calling the eStep and mStep of the model
     *  until the maximal number of iteration is reached or the variation
     *  of the lnLikelihood is less than epsilon.
     * @return @c true if no error occur, @c false otherwise
     **/
    virtual bool run();
};

/** @ingroup Clustering
 *  @brief Implementation of the SimulAlgo learning algorithm.
 *  The SimulAlgo algorithm calls alternatively steps:
 *  - samplingStep()
 *  - paramUpdateStep()
 *  - storeIntermediateResults(iter)
 *  until the maximum number of iterations is reached.
 **/
class SimulAlgo: public IMixtureAlgoLearn
{
  public:
    /** default constructor */
    inline SimulAlgo(): IMixtureAlgoLearn() {}
    /** Copy constructor.
     *  @param algo the algorithm to copy */
    inline SimulAlgo( SimulAlgo const& algo): IMixtureAlgoLearn(algo) {}
    /** destructor */
    inline virtual ~SimulAlgo(){}
    /** clone pattern */
    inline virtual SimulAlgo* clone() const { return new SimulAlgo(*this);}
    /** run the algorithm on the model calling sStep, mStep and eStep of the
     *  model until the maximal number of iteration is reached.
     *  @return @c true if no error occur, @c false otherwise.
     **/
    virtual bool run();
};
} // namespace STK

#endif /* STK_MIXTUREALGOLEARN_H */
