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
 * created on: 23 oct. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_KmmBridge.h
 *  @brief In this file we define the bridge class between the kernel mixtures
 *  and the IMixture interface.
 **/

#ifndef STK_KMMBRIDGE_H
#define STK_KMMBRIDGE_H

#include "STK_Kmm_s.h"
#include "STK_Kmm_sk.h"
#include "../STK_IMixtureBridge.h"

namespace STK
{
// forward declaration
template<int Id, class Data> class KmmBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Kmm_sk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< KmmBridge< Clust::Kmm_sk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Kmm_sk Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Kmm_sk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  // class of mixture
  enum
  {
    idMixtureClass_ = Clust::Kmm_
  };
};

/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Kmm_s_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< KmmBridge< Clust::Kmm_s_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Kmm_s Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Kmm_s_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  // class of mixture
  enum
  {
    idMixtureClass_ = Clust::Kmm_
  };
};

} // namespace hidden

} // namespace STK

namespace STK
{
/** @ingroup Clustering
 *  @brief template implementation of the IMixture interface allowing
 *  to bridge a STK++ kernel mixture with the composer.
 *
 *  This class inherit from the interface IMixtureBridge.
 *
 *  @note The p_data_ DataManager is wrapping the Gram matrix.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureDensity class.
 * @tparam Data is any container storing the data.
 */
template<int Id, class Data>
class KmmBridge: public IMixtureBridge< KmmBridge<Id,Data> >
{
  public:
    typedef IMixtureBridge< KmmBridge<Id,Data> > Base;
    typedef typename hidden::MixtureBridgeTraits< KmmBridge<Id,Data> >::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits< KmmBridge<Id,Data> >::Parameters Parameters;
    typedef typename Data::Type Type;

    // parameters type to get
    using Base::mixture_;
    using Base::p_tik;
    using Base::p_tk;

    /** default constructor.
     *  @param p_dataij pointer on the data set (should be zero, never used)
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    KmmBridge( Data* p_dataij, String const& idData, int nbCluster)
             : Base( p_dataij, idData, nbCluster)
    { initializeBridge();}
    /** copy constructor */
    KmmBridge( KmmBridge const& bridge): Base(bridge) { initializeBridge();}
    /** destructor */
    virtual ~KmmBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual KmmBridge* clone() const { return new KmmBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual KmmBridge* create() const
    {
      KmmBridge* p_bridge = new KmmBridge( mixture_, this->idData(), this->nbCluster());
      p_bridge->setDim(mixture_.param().dim_);
      p_bridge->setKernel(mixture_.p_kernel());
      return p_bridge;
    }
    /** set the dimension of the kernel mixture model */
    inline void setDim(Real const& dim) { mixture_.setDim(dim);}
    /** set the dimension of the kernel mixture model using row vector */
    template<class Vector>
    inline void setDim(ExprBase<Vector> const& dim)
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Vector);
      mixture_.setDim(dim.asDerived());
    }
    /** set the kernel */
    inline void setKernel(Kernel::IKernel const* p_kernel) { mixture_.setKernel(p_kernel);}

    /** @brief do nothing for kernel mixture models */
    virtual void imputationStep() {}
    /** @brief do nothing for kernel mixture models*/
    virtual void samplingStep() {}

    /** @return a safe value for the jth variable
     *  @param j index of the column with the safe value needed */
    Type safeValue( int j) const { return Type();}

  private:
    /** This function will be used in order to initialize the mixture model
     *  using informations stored by the DataBridge.
     *
     *  In the kernel bridge the mixture use the intermediary array @c dik_ computed
     *  by the bridge. The initialization step consist in resizing the array
     *  and to set the pointer to @c mixture_.
     **/
    void initializeBridge() {}
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    KmmBridge( Mixture const& mixture, String const& idData, int nbCluster)
             : Base( mixture, idData, nbCluster)
    {}
};

} // namespace STK

#endif /* STK_KMMBRIDGE_H */
