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

/** @file STK_IMixtureAlgo.h
 *  @brief In this file we define the interface base class for mixture algorithms.
 **/

#ifndef STK_IMIXTUREALGO_H
#define STK_IMIXTUREALGO_H

#include <Sdk.h>
#include "../STK_Clust_Util.h"

namespace STK
{
// forward declaration
class IMixtureComposer;

/** @ingroup Clustering
 * Interface base class for the algorithms.
 * All algorithms are runners applying on a model instance given by pointer
 * and have to implement the run method.
 *
 * All algorithms start with an paramUpdateStep(), so user have to provide an instance of
 * the model with initial parameters values.
 **/
class IMixtureAlgo: public IRunnerBase
{
  protected:
    /** default constructor */
    IMixtureAlgo();
    /** Copy constructor.
     *  @param algo algorithm to copy */
    IMixtureAlgo( IMixtureAlgo const& algo);

  public:
    /** destructor */
    virtual ~IMixtureAlgo();

    // getters
    /** @return the maximal number of iteration of the algorithm */
    inline int nbIterMax() const { return nbIterMax_;}
    /** @return the epsilon of the algorithm */
    inline int epsilon() const { return epsilon_;}
    /** @return the threshold of the algorithm */
    inline Real threshold() const { return threshold_;}

    // setters
    /** set model */
    void setModel(IMixtureComposer* p_model);
    /** set maximal number of iterations */
    inline void setNbIterMax(int nbIterMax) { nbIterMax_ = nbIterMax;}
    /** set tolerance value */
    inline void setEpsilon(Real epsilon) { epsilon_ = epsilon;}
    /** set threshold value */
    inline void setThreshold(Real threshold) { threshold_ = threshold;}

  protected:
    /** pointer on the mixture model */
    IMixtureComposer* p_model_;
    /** number of iterations of the algorithm */
    int nbIterMax_;
    /** tolerance of the algorithm. */
    Real epsilon_;
    /** Minimal number of individuals. If the expected number of individuals
     *  is under this number, the algorithm will stop and return false.
     **/
    Real threshold_;
};


} // namespace STK

#endif /* STK_IMIXTUREALGO_H */
