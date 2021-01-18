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

/** @file STK_PoissonBridge.h
 *  @brief In this file we define the bridge class between the Poisson mixture
 *  models and the composer.
 **/

#ifndef STK_POISSONBRIDGE_H
#define STK_POISSONBRIDGE_H

#include "STK_Poisson_ljk.h"
#include "STK_Poisson_ljlk.h"
#include "STK_Poisson_lk.h"
#include "../STK_IMixtureBridge.h"

namespace STK
{

// forward declaration
template<int Id, class Data> class PoissonBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the STK::Poisson_ljk
 *  mixture model
 **/
template<class Data_>
struct MixtureBridgeTraits< PoissonBridge<Clust::Poisson_ljk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the Mixture model */
  typedef Poisson_ljk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Poisson_ljk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  // class of mixture
  enum
  {
    idMixtureClass_ = Clust::Poisson_
  };
};

/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the STK::Poisson_lk
 *  mixture model
 **/
template<class Data_>
struct MixtureBridgeTraits< PoissonBridge<Clust::Poisson_lk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the Mixture model */
  typedef Poisson_lk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Poisson_lk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  // class of mixture
  enum
  {
    idMixtureClass_ = Clust::Poisson_
  };
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the STK::Poisson_ljlk
 *  mixture model
 **/
template<class Data_>
struct MixtureBridgeTraits< PoissonBridge< Clust::Poisson_ljlk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef Poisson_ljlk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Poisson_ljlk_> Parameters;
  /** Type of the array storing missing values indexes */
  typedef std::vector<std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;
  // class of mixture
  enum
  {
    idMixtureClass_ = Clust::Poisson_
  };
};

} // namespace hidden

/** @ingroup Clustering
 *  @brief template implementation of the IMixtureBridge interface allowing
 *  to bridge a STK++ Poisson mixture with the composer.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureDensity class.
 */
template<int Id, class Data>
class PoissonBridge: public IMixtureBridge< PoissonBridge<Id,Data> >
{
  public:
    // Base class
    typedef IMixtureBridge< PoissonBridge<Id,Data> > Base;
    typedef typename hidden::MixtureBridgeTraits< PoissonBridge<Id,Data> >::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits< PoissonBridge<Id,Data> >::Parameters Parameters;
    typedef typename Data::Type Type;
    // class of mixture
    enum
    {
      idMixtureClass_ = Clust::Poisson_
    };
    typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
    using Base::mixture_;
    using Base::p_dataij_;
    using Base::v_missing_;
    using Base::removeMissing;

    /** default constructor. Remove the missing values from the data set and
     *  initialize the mixture by setting the data set.
     *  @param p_dataij pointer on the data set used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    PoissonBridge( Data* p_dataij, String const& idData, int nbCluster)
                 : Base(p_dataij, idData, nbCluster)
    {
      removeMissing(); // remove missing from data only once at creation
      mixture_.setData(*p_dataij_);
    }
    /** copy constructor */
    PoissonBridge( PoissonBridge const& bridge): Base(bridge) {}
    /** destructor */
    virtual ~PoissonBridge()
    {
#ifdef STK_MIXTURE_DEBUG_CREATE
        stk_cout << _T("Entering PoissonBridge::~PoissonBridge()\n");
        stk_cout << _T("this =") << this << _T("\n");
#endif

    }
  private:
    /** private constructor used in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    PoissonBridge( Mixture const& mixture, String const& idData, int nbCluster)
                 : Base(mixture, idData, nbCluster)
    {}
  public:
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual PoissonBridge* clone() const { return new PoissonBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of PoissonBridge as that of calling object.
     */
    virtual PoissonBridge* create() const
    {
#ifdef STK_MIXTURE_DEBUG_CREATE
        stk_cout << _T("Entering PoissonBridge::create()\n");
        stk_cout << _T("this =") << this << _T("\n");
#endif

      PoissonBridge* p_bridge = new PoissonBridge( mixture_, this->idData(), this->nbCluster());
#ifdef STK_MIXTURE_DEBUG_CREATE
        stk_cout << _T("p_bridge created\n");
        stk_cout << _T("p_bridge =") << p_bridge << _T("\n");
#endif
      p_bridge->p_dataij_ = p_dataij_;
      p_bridge->mixture_.setData(*p_dataij_);
      p_bridge->v_missing_ = v_missing_;
      return p_bridge;
    }

    /** @return a safe value for the jth variable
     *  @param j index of the column with the safe value needed */
    Type safeValue( int j) const
    {
      int lmin = p_dataij_->col(j).safe().minElt(), lmax = p_dataij_->col(j).safe().maxElt();
      if (lmax -lmin > 10)
      { return Real(p_dataij_->col(j).safe().sum())/p_dataij_->sizeRows();}
      Array2DVector<int> count(Range(lmin, lmax, 0), 0);
      for (int i= p_dataij_->beginRows(); i < p_dataij_->endRows(); ++i)
      {
        if (!Arithmetic<int>::isNA(p_dataij_->elt(i,j)))
          count[p_dataij_->elt(i,j)]++;
      }
      int l; count.maxElt(l);
      return l;
    }
};

} // namespace STK

#endif /* STK_POISSONBRIDGE_H */
