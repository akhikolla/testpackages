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

/** @file STK_MixtureAlgo.h
 *  @brief In this file we define mixture algorithms
 **/

#ifndef STK_MIXTUREALGO_H
#define STK_MIXTUREALGO_H

#include "STK_IMixtureAlgo.h"

namespace STK
{

/** @ingroup Clustering
 *  @brief Implementation of the EM algorithm.
 *  The EM algorithm call alternatively the steps:
 *  - paramUpdateStep()
 *  - eStep()
 *  until the maximum number of iterations is reached or the variation of the
 *  ln-likelihood is less than the tolerance.
 **/
class EMAlgo: public IMixtureAlgo
{
  public:
    /** default constructor */
    inline EMAlgo(): IMixtureAlgo() {}
    /** Copy constructor.
     *  @param algo the algorithm to copy */
    inline EMAlgo( EMAlgo const& algo): IMixtureAlgo(algo) {}
    /** destructor */
    inline virtual ~EMAlgo(){}
    /** clone pattern */
    inline virtual EMAlgo* clone() const { return new EMAlgo(*this);}
    /** run the algorithm on the model calling the eStep and mStep of the model
     *  until the maximal number of iteration is reached or the variation
     *  of the lnLikelihood is less than epsilon.
     * @return @c true if no error occur, @c false otherwise
     **/
    virtual bool run();
};

/** @ingroup Clustering
 *  @brief Implementation of the SEM algorithm.
 *  The CEM algorithm calls alternatively the steps:
 *  - cStep()
 *  - paramUpdateStep()
 *  - eStep()
 *  until the maximum number of iterations is reached or the variation of the
 *  ln-likelihood is less than the tolerance.
 **/
class CEMAlgo: public IMixtureAlgo
{
  public:
    /** default constructor */
    inline CEMAlgo(): IMixtureAlgo() {}
    /** Copy constructor.
     *  @param algo the algorithm to copy */
    inline CEMAlgo( CEMAlgo const& algo): IMixtureAlgo(algo) {}
    /** destructor */
    inline virtual ~CEMAlgo(){}
    /** clone pattern */
    inline virtual CEMAlgo* clone() const { return new CEMAlgo(*this);}
    /** run the algorithm on the model calling cStep, mStep and eStep of the
     *  model until the maximal number of iteration is reached or the variation
     *  of the lnLikelihood is less than epsilon.
     *  @return @c true if no error occur, @c false otherwise
     **/
    virtual bool run();
};

/** @ingroup Clustering
 *  @brief Implementation of the SEM algorithm.
 *  The SEM algorithm calls alternatively the steps:
 *  - sStep()
 *  - samplingStep()
 *  - pStep()
 *  - paramUpdateStep()
 *  - eStep()
 *  - storeIntermediateResults(iter)
 *  until the maximum number of iterations is reached.
 **/
class SEMAlgo: public IMixtureAlgo
{
  public:
    /** default constructor */
    inline SEMAlgo(): IMixtureAlgo() {}
    /** Copy constructor.
     *  @param algo the algorithm to copy */
    inline SEMAlgo( SEMAlgo const& algo): IMixtureAlgo(algo) {}
    /** destructor */
    inline virtual ~SEMAlgo(){}
    /** clone pattern */
    inline virtual SEMAlgo* clone() const { return new SEMAlgo(*this);}
    /** run the algorithm on the model calling sStep, mStep and eStep of the
     *  model until the maximal number of iteration is reached.
     *  @return @c true if no error occur, @c false otherwise.
     **/
    virtual bool run();
};

/** @ingroup Clustering
 *  @brief Implementation of the SemiSEM algorithm.
 *  The SemiSEM algorithm calls alternatively the steps:
 *  - samplingStep()
 *  - pStep()
 *  - paramUpdateStep()
 *  - eStep()
 *  - storeIntermediateResults(iter)
 *  until the maximum number of iterations is reached.
 **/
class SemiSEMAlgo: public IMixtureAlgo
{
  public:
    /** default constructor */
    inline SemiSEMAlgo(): IMixtureAlgo() {}
    /** Copy constructor.
     *  @param algo the algorithm to copy */
    inline SemiSEMAlgo( SemiSEMAlgo const& algo): IMixtureAlgo(algo) {}
    /** destructor */
    inline virtual ~SemiSEMAlgo(){}
    /** clone pattern */
    inline virtual SemiSEMAlgo* clone() const { return new SemiSEMAlgo(*this);}
    /** run the algorithm on the model calling sStep, mStep and eStep of the
     *  model until the maximal number of iteration is reached.
     *  @return @c true if no error occur, @c false otherwise.
     **/
    virtual bool run();
};

} // namespace STK

#endif /* STK_MIXTUREALGO_H */
