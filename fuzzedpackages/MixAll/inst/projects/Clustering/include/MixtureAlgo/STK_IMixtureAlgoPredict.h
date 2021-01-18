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

/** @file STK_IMixtureAlgoPredict.h
 *  @brief In this file we define the interface base class for mixture predicting algorithms.
 **/

#ifndef STK_IMIXTUREALGOPREDICT_H
#define STK_IMIXTUREALGOPREDICT_H

#include <Sdk.h>
#include "../STK_Clust_Util.h"

namespace STK
{
// forward declaration
class IMixtureComposer;

/** @ingroup Clustering
 * @brief Interface base class for predicting algorithms
 *
 * Predicting algorithms are decomposed in two stages :
 * - A burning stage using SEM algorithm
 * - An estimation stage using EM or SemiSEM algorithm
 * If there is no missing values, the algorithm simplify in a posterior allocation stage.
 *
 * All algorithms are runners applying on a IMixtureComposer model instance given by pointer
 * and have to implement the run method.
 **/
class IMixtureAlgoPredict: public IRunnerBase
{
  protected:
    /** default constructor */
    IMixtureAlgoPredict();
    /** Copy constructor.
     *  @param algo the algorithm to copy */
    IMixtureAlgoPredict( IMixtureAlgoPredict const& algo);

  public:
    /** destructor */
    virtual ~IMixtureAlgoPredict();

    // getters
    /** @return the maximal number of iteration of the algorithm */
    inline int nbIterBurn() const { return nbIterBurn_; }
    /** @return the maximal number of iteration of the algorithm */
    inline int setNbIterLong() const { return nbIterLong_; }
    /** @return the epsilon of the algorithm */
    inline int epsilon() const { return epsilon_;}

    // setters
    /** set model */
    void setModel(IMixtureComposer* p_model);
    /** set number of burning iterations */
    inline void setNbIterBurn(int nbIterBurn) { nbIterBurn_ = nbIterBurn; }
    /** set number of long iterations */
    inline void setNbIterLong(int nbIterLong) { nbIterLong_ = nbIterLong; }
    /** set tolerance value */
    inline void setEpsilon(Real epsilon) { epsilon_ = epsilon; }

  protected:
    /** pointer on the mixture model */
    IMixtureComposer* p_model_;
    /** Number of burning iterations of the algorithm */
    int nbIterBurn_;
    /** maximal number of iterations of the algorithm */
    int nbIterLong_;
    /** tolerance of the algorithm. */
    Real epsilon_;

    /** predict class labels when there is no missing values.
     * In this case, there is no algorithm to do. Just call:
     * @code
     *   initializeStep();
     *   eStep(); // will compute tik and zi values
     *   finalizeStep();
     * @endcode
     * @return @c true if no error occur, @c false otherwise */
    bool predictBayesClassifier();
    /** Perform burn step using SEM algorithm
     * @code
     *   initializeStep();
     *   eStep(); // will compute tik and zi values
     *   finalizeStep();
     * @endcode
     * @return @c true if no error occur, @c false otherwise */
    bool burnStep();

};

} // namespace STK

#endif /* STK_IMIXTUREALGO_H */
