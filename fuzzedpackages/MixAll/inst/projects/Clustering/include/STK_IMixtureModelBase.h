/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff

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
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_IMixtureModelBase.h
 *  @brief In this file we define the interface base class for mixture models.
 **/

#ifndef STK_IMIXTUREMODELBASE_H
#define STK_IMIXTUREMODELBASE_H

#include <cmath>
#include "STK_Clust_Util.h"

namespace STK
{
/** @ingroup Clustering
 *  @brief Base class for all Mixture model.
 *  Let @e X be an arbitrary measurable space and let
 *  \f$ \mathbf{x} =\{{\mathbf x}_1,...,{\mathbf x}_n\}\f$
 *  be @e n independent vectors in @e X. If each \f${\mathbf x}_i\f$
 *  arises from a probability distribution with density
 *  \f[
 *    f({\mathbf x}_i|\theta) = \sum_{k=1}^K p_k h({\mathbf x}_{i}| \lambda_{k},\alpha)
 *  \f]
 *  where the \f$p_k\f$'s are the mixing proportions, \f$h(\cdot| \lambda_{k},\alpha)\f$
 *  denotes a @e d-dimensional distribution parameterized by \f$\lambda_k\f$ and
 *  \f$\alpha\f$, it is said that we observe a mixture model..
 *  @sa IMixture, IMixtureModel, MixtureComposer
 */
class IMixtureModelBase
{
  public:
    /** default constructor
     * @param nbCluster number of cluster */
    IMixtureModelBase( int nbCluster);
    /** copy constructor.
     *  @note the pointer are initialized to 0.
     *  @param model the model to clone
     **/
    IMixtureModelBase( IMixtureModelBase const& model);
    /** destructor */
    ~IMixtureModelBase();
    /** @return the number of cluster */
    inline int nbCluster() const { return nbCluster_;}
    /** @return the total available observations */
    inline int nbSample() const { return nbSample_;}
    /** @return the Log of the total available observations */
    inline Real lnNbSample() const
    { return (nbSample_ <= 0) ? -Arithmetic<Real>::infinity() : std::log((Real)nbSample_);}

  protected:
    /** Set the number of sample of the model
     *  @param nbSample number of sample of the model
     * */
    inline void setNbSample( int nbSample) { nbSample_ = nbSample;}

  private:
    /** number of cluster. */
    int nbCluster_;
    /** total available samples */
    int nbSample_;
};

} // namespace STK

#endif /* STK_IMIXTUREMODELBASE_H */
