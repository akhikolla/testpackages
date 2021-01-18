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
 * created on: 15 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_CategoricalBridge.h
 *  @brief In this file we define the bridge classes between the categorical
 *  mixtures and the composer.
 **/

#ifndef STK_CATEGORICALBRIDGE_H
#define STK_CATEGORICALBRIDGE_H

#include "STK_Categorical_pjk.h"
#include "STK_Categorical_pk.h"
#include "../STK_IMixtureBridge.h"

namespace STK
{

// forward declaration
template<int Id, class Data> class CategoricalBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Categorical_pjk model
 **/
template<class Data_>
struct MixtureBridgeTraits< CategoricalBridge<Clust::Categorical_pjk_, Data_> >
{
  /** Structure storing the Data */
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the Mixture model */
  typedef Categorical_pjk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Categorical_pjk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Categorical_
  };
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the Categorical_pk model
 **/
template<class Data_>
struct MixtureBridgeTraits< CategoricalBridge< Clust::Categorical_pk_, Data_> >
{
  /** Structure storing the Data */
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Categorical_pk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Categorical_pk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Categorical_
  };
};

} // namespace hidden

/** @ingroup Clustering
 *  @brief template implementation of the IMixtureBridge interface allowing
 *  to bridge a STK++ mixture with the composer.
 *
 *  This class inherit from the IMixtureBridge.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureDensity class.
 */
template<int Id, class Data>
class CategoricalBridge: public IMixtureBridge< CategoricalBridge<Id,Data> >
{
  public:
    // Base class
    typedef IMixtureBridge< CategoricalBridge<Id,Data> > Base;
    // type of Mixture
    typedef typename hidden::MixtureBridgeTraits< CategoricalBridge<Id,Data> >::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits< CategoricalBridge<Id,Data> >::Parameters Parameters;
    typedef typename hidden::MixtureBridgeTraits< CategoricalBridge<Id,Data> >::MissingIndexes MissingIndexes;
    // type of data
    typedef typename Data::Type Type;
    // class of mixture
    enum
    {
      idMixtureClass_ = Clust::Categorical_
    };

    typedef typename MissingIndexes::const_iterator ConstIterator;
    using Base::mixture_;
    using Base::p_dataij_;
    using Base::p_tik;
    using Base::v_missing_;
    using Base::removeMissing;

    /** default constructor.
     *  - Remove the missing values from the data set,
     *  - initialize the mixture by setting the data set
     *  @param p_dataij pointer on the data set that will be used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    CategoricalBridge( Data* p_dataij, String const& idData, int nbCluster);
    /** copy constructor */
    CategoricalBridge( CategoricalBridge const& bridge);
    /** destructor */
    virtual ~CategoricalBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual CategoricalBridge* clone() const { return new CategoricalBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual CategoricalBridge* create() const
    {
      CategoricalBridge* p_bridge = new CategoricalBridge( mixture_, this->idData(), this->nbCluster());
      p_bridge->p_dataij_ = p_dataij_;
      p_bridge->mixture_.setData(*p_dataij_);
      p_bridge->v_missing_ = v_missing_;
      return p_bridge;
    }
    /** Compute a safe value for the jth variable by counting the most present
     *  modality.
     *  @return a safe value for the jth variable
     *  @param j index of the column with the safe value needed
     **/
    Type safeValue( int j) const;

  private:
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    CategoricalBridge( Mixture const& mixture, String const& idData, int nbCluster)
                     : Base(mixture, idData, nbCluster)
    {}
};

/* default constructor.
 *  - Remove the missing values from the data set,
 *  - initialize the mixture by setting the data set
 *  @param p_dataij pointer on the Data set that will be used by the bridge.
 *  @param idData id name of the mixture model
 *  @param nbCluster number of cluster
 **/
template<int Id, class Data>
CategoricalBridge<Id,Data>::CategoricalBridge( Data* p_dataij, String const& idData, int nbCluster)
                                             : Base( p_dataij, idData, nbCluster)
{
  removeMissing();
  mixture_.setData(*p_dataij_);
}
/* copy constructor */
template<int Id, class Data>
CategoricalBridge<Id,Data>::CategoricalBridge( CategoricalBridge const& bridge): Base(bridge)
{/*mixture_.setData(*p_dataij_);*/}

/* Compute a safe value for the jth variable by counting the most present
 *  modality.
 *  @return a safe value for the jth variable
 *  @param j index of the column with the safe value needed
 **/
template<int Id, class Data>
typename Data::Type CategoricalBridge<Id,Data>::safeValue( int j) const
{
  int lmin = p_dataij_->col(j).safe().minElt(), lmax = p_dataij_->col(j).safe().maxElt();
  Array2DVector<int> count(Range(lmin, lmax, 0), 0);
  for (int i= p_dataij_->beginRows(); i < p_dataij_->endRows(); ++i)
  {
    if (!Arithmetic<int>::isNA(p_dataij_->elt(i,j)))
      count[p_dataij_->elt(i,j)]++;
  }
  int l; count.maxElt(l);
  return l;
}

} // namespace STK

#endif /* STK_CATEGORICALBRIDGE_H */
