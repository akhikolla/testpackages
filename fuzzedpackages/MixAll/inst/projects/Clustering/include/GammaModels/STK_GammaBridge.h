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

/** @file STK_GammaBridge.h
 *  @brief In this file we define the bridge classes between the diagonal
 *  Gaussian mixtures and the composer.
 **/

#ifndef STK_GAMMABRIDGE_H
#define STK_GAMMABRIDGE_H

#include "STK_Gamma_a_bjk.h"
#include "STK_Gamma_a_bk.h"
#include "STK_Gamma_aj_bjk.h"
#include "STK_Gamma_aj_bk.h"
#include "STK_Gamma_ajk_b.h"
#include "STK_Gamma_ajk_bj.h"
#include "STK_Gamma_ajk_bjk.h"
#include "STK_Gamma_ajk_bk.h"
#include "STK_Gamma_ak_b.h"
#include "STK_Gamma_ak_bj.h"
#include "STK_Gamma_ak_bjk.h"
#include "STK_Gamma_ak_bk.h"
#include "../STK_IMixtureBridge.h"

namespace STK
{

// forward declaration
template<int Id, class Data> class GammaBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ajk_bjk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ajk_bjk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_ajk_bjk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ajk_bjk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ajk_bk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ajk_bk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_ajk_bk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ajk_bk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ajk_bj_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ajk_bj_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_ajk_bj<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ajk_bj_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ajk_b_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ajk_b_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_ajk_b<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ajk_b_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ak_bjk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ak_bjk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_ak_bjk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ak_bjk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ak_bk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ak_bk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_ak_bk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ak_bk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ak_bj_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ak_bj_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_ak_bj<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ak_bj_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_ak_b_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ak_b_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_ak_b<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ak_b_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_aj_bjk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_aj_bjk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_aj_bjk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_aj_bjk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_aj_bk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_aj_bk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_aj_bk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_aj_bk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_a_bjk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_a_bjk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_a_bjk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_a_bjk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the Gamma_a_bk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_a_bk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Gamma_a_bk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_a_bk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};

} // namespace hidden

/** @ingroup Clustering
 *  @brief template implementation of the IMixtureBridge interface allowing
 *  to bridge a STK++ mixture with the composer.
 *
 *  This class inherit from the interface IMixture and delegate almost
 *  all the treatments to the wrapped class.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureDensity class.
 */
template<int Id, class Data>
class GammaBridge: public IMixtureBridge< GammaBridge<Id,Data> >
{
  public:
    // Base class
    typedef IMixtureBridge< GammaBridge<Id,Data> > Base;
    // type of Mixture
    typedef typename hidden::MixtureBridgeTraits< GammaBridge<Id,Data> >::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits< GammaBridge<Id,Data> >::Parameters Parameters;
    // type of data
    typedef typename Data::Type Type;
    // class of mixture
    enum
    {
      idMixtureClass_ = Clust::Gamma_
    };
    typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
    using Base::mixture_;
    using Base::p_dataij_;
    using Base::p_tik;
    using Base::v_missing_;
    using Base::removeMissing;

    /** default constructor. Remove the missing values from the data set and
     *  initialize the mixture by setting the data set.
     *  @param p_dataij pointer on the data set that will be used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    GammaBridge( Data* p_dataij, String const& idData, int nbCluster)
               : Base( p_dataij, idData, nbCluster)
    {
      removeMissing(); // remove missing from data only once at creation
      mixture_.setData(*p_dataij_);
    }
    /** copy constructor */
    GammaBridge( GammaBridge const& bridge): Base(bridge)
    {}
//    { mixture_.setData(*p_dataij_);}
    /** destructor */
    virtual ~GammaBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual GammaBridge* clone() const { return new GammaBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual GammaBridge* create() const
    {
      GammaBridge* p_bridge = new GammaBridge( mixture_, this->idData(), this->nbCluster());
      p_bridge->p_dataij_ = p_dataij_;
      p_bridge->mixture_.setData(*p_dataij_);
      p_bridge->v_missing_ = v_missing_;
      return p_bridge;
    }
    /** @return a safe value for the jth variable
     *  @param j index of the column with the safe value needed */
    Type safeValue( int j) const
    { return p_dataij_->col(j).safe(1).mean();}

  private:
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    GammaBridge( Mixture const& mixture, String const& idData, int nbCluster)
              : Base(mixture, idData, nbCluster)
    {}
};

} // namespace STK

#endif /* STK_GAMMABRIDGE_H */
