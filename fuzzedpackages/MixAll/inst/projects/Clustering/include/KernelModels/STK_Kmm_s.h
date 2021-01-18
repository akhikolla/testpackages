/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016 Serge Iovleff

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: Oct 24, 2014
 * Author:   Serge Iovleff
 **/

/** @file STK_Kmm_s.h
 *  @brief In this file we define the Kmm_s class
 **/

#ifndef STK_KMM_S_H
#define STK_KMM_S_H

#include "STK_KmmBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
class Kmm_s;

namespace hidden
{
/** @ingroup Clustering
 *  MixtureTraits class for the Kmm_s traits policy. */
template<>
struct MixtureTraits< Kmm_s >
{
  typedef CArrayXX Array;
  typedef Array::Type       Type;
  /** Type of the structure storing the parameters of a Kmm_s model*/
  typedef ModelParameters<Clust::Kmm_s_> Parameters;
};

} // namespace Clust

/** @ingroup Clustering
 *  The Gaussian mixture model @c Kmm_s is an isotrope Gaussian
 *  mixture model on a kernel space. It has a density function of the form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k
 *    \sum_{k=1}^K p_k \left(\frac{1}{\sqrt{2\pi}\sigma}\right)^{d_k}
 *    \exp\left\{ -\frac{\|\phi(x)-m_k\|^2}{2\sigma^2}  \right\}
 * \f]
 * where \f$ \phi \f$ denote a feature mapping from the original space to an RKHS.
 *
 * In a Kmm_s model, the data set refer to the Gram's matrix.
 **/
class Kmm_s: public KmmBase<Kmm_s >
{
  public:
    typedef KmmBase<Kmm_s> Base;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Kmm_s( int nbCluster);
    /** copy constructor
     *  @param model The model to copy
     **/
    Kmm_s( Kmm_s const& model);
    /** destructor */
    ~Kmm_s();

    /** @return the number of free parameters of the model */
    int computeNbFreeParameters() const;

    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    Real lnComponentProbability(int i, int k) const;

    /** Initialize randomly the variances of the Gaussian kernel mixture. */
    void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
    /** update the variances. */
    bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk);
};

} // namespace STK

#endif /* STK_KMM_S_H */
