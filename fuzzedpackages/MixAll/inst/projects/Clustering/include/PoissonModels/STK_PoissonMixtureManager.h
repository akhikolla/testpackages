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
 * created on: 15 mars 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_PoissonMixtureManager.h
 *  @brief In this file we define the PoissonMixtureManager class.
 **/


#ifndef STK_POISSONMIXTUREMANAGER_H
#define STK_POISSONMIXTUREMANAGER_H

#include "../STK_Clust_Util.h"
#include "../STK_IMixtureManager.h"

#include "STK_PoissonBridge.h"


#define STK_CREATE_MIXTURE(Data, Bridge) \
          Data* p_data = new Data(idData); \
          p_handler()->getData(idData, p_data->dataij() ); \
          registerDataBridge(p_data); \
          return new Bridge( &(p_data->dataij()), idData, nbCluster);

namespace STK
{

// forward declaration
template<class DataHandler> class PoissonMixtureManager;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization for PoissonMixtureManager
 **/
template <class DataHandler_>
struct MixtureManagerTraits <PoissonMixtureManager<DataHandler_> >
{
  /** type of data */
  typedef DataHandler_ DataHandler;
  /** type of data */
  typedef Integer Type;
  /** Type of the array storing missing values indexes */
  typedef std::vector< std::pair<int,int> > MissingIndexes;
  /** Type of the array storing missing values */
  typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;

  // All data handlers will store and return a specific container for
  // the data they handle. The DataHandlerTraits class allow us to know the
  // type of these containers when data is of type Type.
  /** */
  /** type of the data set */
  typedef typename DataHandlerTraits<DataHandler, Type>::Data Data;
  // Classes wrapping the Real and Integer containers
  /** class wrapping the data set */
  typedef DataBridge<Data>  DataBridgeType;
};

} // namespace hidden

/** @ingroup Clustering
 *  @brief A mixture manager is a factory class for injection dependency in the
 *  STK++ derived class of the MixtureComposer class.
 *
 *  It allows to handle all the creation and initialization stuff needed by the
 *  (bridged) mixture models of the stkpp library.
 *
 *  @tparam DataHandler is any concrete class from the interface DataHandlerBase
 */
template<class DataHandler>
class PoissonMixtureManager: public IMixtureManager< PoissonMixtureManager<DataHandler> >
{
  public:
    typedef typename hidden::MixtureManagerTraits< PoissonMixtureManager >::Type Type;
    typedef typename hidden::MixtureManagerTraits< PoissonMixtureManager >::MissingIndexes MissingIndexes;
    typedef typename hidden::MixtureManagerTraits< PoissonMixtureManager >::MissingValues MissingValues;
    typedef typename hidden::MixtureManagerTraits< PoissonMixtureManager >::Data Data;
    typedef typename hidden::MixtureManagerTraits< PoissonMixtureManager >::DataBridgeType DataBridgeType;

    typedef IMixtureManager< PoissonMixtureManager > Base;
    using Base::getDataBridge;
    using Base::getIdModel;
    using Base::registerDataBridge;
    using Base::p_handler;

    // All Poisson bridges
    typedef PoissonBridge<Clust::Poisson_ljk_,  Data> MixtureBridge_ljk;
    typedef PoissonBridge<Clust::Poisson_lk_,   Data> MixtureBridge_lk;
    typedef PoissonBridge<Clust::Poisson_ljlk_, Data> MixtureBridge_ljlk;

    /** Default constructor, need an instance of a DataHandler.  */
    PoissonMixtureManager(DataHandler const& handler): Base(&handler) {}
    /** destructor */
    ~PoissonMixtureManager() {}

    /** get the missing values from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param missing structure to return with the missing values
     **/
    void getMissingValuesImpl(IMixture* p_mixture, MissingValues& missing) const
    {
      Clust::Mixture idModel = getIdModel(p_mixture->idData());
      if (idModel == Clust::unknown_mixture_) return;
      // up-cast... (Yes it's bad....;)...)
      switch (idModel)
      {
        case Clust::Poisson_ljk_:
        { static_cast<MixtureBridge_ljk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Poisson_lk_:
        { static_cast<MixtureBridge_lk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Poisson_ljlk_:
        { static_cast<MixtureBridge_ljlk*>(p_mixture)->getMissingValues(missing);}
        break;
        default: // idModel is not implemented
        break;
      }
    }
    /** get the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param param the array to return with the parameters
     **/
    void getParametersImpl(IMixture* p_mixture, ArrayXX& param) const
    {
      Clust::Mixture idModel = getIdModel(p_mixture->idData());
      if (idModel == Clust::unknown_mixture_) return;
      // up-cast... (Yes it's bad....;)...)
      switch (idModel)
      {
        case Clust::Poisson_ljk_:
        { static_cast<MixtureBridge_ljk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Poisson_lk_:
        { static_cast<MixtureBridge_lk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Poisson_ljlk_:
        { static_cast<MixtureBridge_ljlk*>(p_mixture)->getParameters(param);}
        break;
        default: // idModel is not implemented
        break;
      }
    }
    /** set the parameters to an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param param the array with the parameters to set
     **/
    void setParametersImpl(IMixture* p_mixture, ArrayXX const& param) const
    {
      Clust::Mixture idModel = getIdModel(p_mixture->idData());
      if (idModel == Clust::unknown_mixture_) return;
      // up-cast... (Yes it's bad....;)...)
      switch (idModel)
      {
        case Clust::Poisson_ljk_:
        { static_cast<MixtureBridge_ljk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Poisson_lk_:
        { static_cast<MixtureBridge_lk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Poisson_ljlk_:
        { static_cast<MixtureBridge_ljlk*>(p_mixture)->setParameters(param);}
        break;
        default: // idModel is not implemented
        break;
      }
    }
    /** create a concrete mixture and initialize it.
     *  @param modelName a valid model name
     *  @param idData Id of the data
     *  @param nbCluster number of cluster of the model
     **/
    IMixture* createMixtureImpl(String const& modelName, String const& idData, int nbCluster)
    {
      Clust::Mixture idModel = Clust::stringToMixture(modelName);
      return createMixtureImpl(idModel, idData, nbCluster);
    }

  private:
    /** create a concrete mixture and initialize it.
     *  @param idModel Id number of the model
     *  @param idData Id name of the data
     *  @param nbCluster number of cluster of the model
     **/
    IMixture* createMixtureImpl(Clust::Mixture idModel, String const& idData, int nbCluster)
    {
      DataBridgeType* p_dataBridge = new DataBridgeType(idData);
      p_handler()->getData(idData, p_dataBridge->dataij() );
      registerDataBridge(p_dataBridge);
      switch (idModel)
      {
        case Clust::Poisson_ljk_:
        {
          MixtureBridge_ljk* p_mixt = new MixtureBridge_ljk( &(p_dataBridge->dataij()), idData, nbCluster);
          return p_mixt;
        }
        break;
        case Clust::Poisson_lk_:
        { return new MixtureBridge_lk( &(p_dataBridge->dataij()), idData, nbCluster);}
        break;
        case Clust::Poisson_ljlk_:
        { return new MixtureBridge_ljlk( &(p_dataBridge->dataij()), idData, nbCluster);}
        break;
        default:
          return 0; // 0 if idModel is not implemented
          break;
      }
      return 0; // 0 if idModel is not a STK++ model
    }
};

} // namespace STK

#undef STK_CREATE_MIXTURE

#endif /* STK_POISSONMIXTUREMANAGER_H */
