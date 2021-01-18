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

/** @file STK_MixtureAlgoPredict.h
 *  @brief In this file we define algorithms for predicting in a mixture model.
 **/

#ifndef STK_MIXTUREALGOPREDICT_H
#define STK_MIXTUREALGOPREDICT_H

#include "STK_IMixtureAlgoPredict.h"

namespace STK
{

/** @ingroup Clustering
 *  @brief Implementation of the EMPredict algorithm.
 *  EMPredict algorithm start calling an initializationStep and then
 *  calls until convergence steps:
 *  - imputationStep()
 *  - eStep()
 *  until the maximum number of iterations or the threshold is reached.
 *  This is the counterpart of EMSAlgo for prediction.
 **/
class EMPredict: public IMixtureAlgoPredict
{
  public:
    /** default constructor */
    inline EMPredict(): IMixtureAlgoPredict() {}
    /** Copy constructor.
     *  @param algo the algorithm to copy */
    inline EMPredict( EMPredict const& algo): IMixtureAlgoPredict(algo) {}
    /** destructor */
    inline virtual ~EMPredict(){}
    /** clone pattern */
    inline virtual EMPredict* clone() const { return new EMPredict(*this);}
    /** run the algorithm on the model until the maximal number of iteration
     * or the threshold is reached.
     *  @return @c true if no error occur, @c false otherwise.
     **/
    virtual bool run();
};

/** @ingroup Clustering
 *  @brief Implementation of the SemiSEMPredict algorithm.
 *  SemiSEMPredict algorithm starts calling an initializationStep and then
 *  calls until the number of iteration is reached steps:
 *  - samplingStep()
 *  - eStep()
 *  - storeIntermediateResults(iter).
 *  It ends by a call to setParameterStep().
 *  This is the counterpart of SemiSEMSAlgo for prediction.
 **/
class SemiSEMPredict: public IMixtureAlgoPredict
{
  public:
    /** default constructor */
    inline SemiSEMPredict(): IMixtureAlgoPredict() {}
    /** Copy constructor.
     *  @param algo the algorithm to copy */
    inline SemiSEMPredict( SemiSEMPredict const& algo): IMixtureAlgoPredict(algo) {}
    /** destructor */
    inline virtual ~SemiSEMPredict(){}
    /** clone pattern */
    inline virtual SemiSEMPredict* clone() const { return new SemiSEMPredict(*this);}
    /** run the algorithm on the model until the maximal number of iteration is
     *  reached.
     *  @return @c true if no error occur, @c false otherwise.
     **/
    virtual bool run();
};

} // namespace STK

#endif /* STK_MIXTUREALGOPREDICT_H */
