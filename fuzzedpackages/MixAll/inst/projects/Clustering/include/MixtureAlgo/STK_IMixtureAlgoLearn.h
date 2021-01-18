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

/** @file STK_IMixtureAlgoLearn.h
 *  @brief In this file we define the interface base class for mixture algorithms.
 **/

#ifndef STK_IMIXTURELEARNALGO_H
#define STK_IMIXTURELEARNALGO_H

#include <Sdk.h>
#include "../STK_Clust_Util.h"

namespace STK
{
// forward declaration
class IMixtureLearner;

/** @ingroup Clustering
 * Interface base class for the learning algorithms.
 * All algorithms are runners applying on a model instance given by pointer
 * and have to implement the run method.
 **/
class IMixtureAlgoLearn: public IRunnerBase
{
  protected:
    /** default constructor */
    IMixtureAlgoLearn();
    /** Copy constructor.
     *  @param algo the algorithm to copy */
    IMixtureAlgoLearn( IMixtureAlgoLearn const& algo);

  public:
    /** destructor */
    virtual ~IMixtureAlgoLearn();

    // getters
    /** @return the maximal number of iteration of the algorithm */
    inline int nbIterMax() const { return nbIterMax_; }
    /** @return the epsilon of the algorithm */
    inline int epsilon() const { return epsilon_;}

    // setters
    /** set model */
    void setModel(IMixtureLearner* p_model);
    /** set maximal number of iterations */
    inline void setNbIterMax(int nbIterMax) { nbIterMax_ = nbIterMax; }
    /** set tolerance value */
    inline void setEpsilon(Real epsilon) { epsilon_ = epsilon; }

  protected:
    /** pointer on the mixture model */
    IMixtureLearner* p_model_;
    /** maximal number of iterations of the algorithm */
    int nbIterMax_;
    /** tolerance of the algorithm. */
    Real epsilon_;
};


} // namespace STK

#endif /* STK_IMIXTURELEARNALGO_H */
