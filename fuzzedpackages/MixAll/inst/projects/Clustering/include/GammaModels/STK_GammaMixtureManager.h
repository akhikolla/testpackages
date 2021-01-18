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

/** @file STK_GammaMixtureManager.h
 *  @brief In this file we define the GammaMixtureManager class.
 **/


#ifndef STK_GAMMAMIXTUREMANAGER_H
#define STK_GAMMAMIXTUREMANAGER_H

#include "../STK_IMixtureManager.h"
#include "../STK_Clust_Util.h"

#include "STK_GammaBridge.h"

#define STK_CREATE_MIXTURE(Data, Bridge) \
          Data* p_data = new Data(idData); \
          registerDataBridge(p_data); \
          p_handler()->getData(idData, p_data->dataij() ); \
          return new Bridge( &(p_data->dataij()), idData, nbCluster);

namespace STK
{

// forward declaration
template<class DataHandler> class GammaMixtureManager;

namespace hidden
{
/** @ingroup hidden
 *  Specialization for GammaMixtureManager */
template <class DataHandler_>
struct MixtureManagerTraits< GammaMixtureManager<DataHandler_> >
{
  /** type of data */
  typedef DataHandler_ DataHandler;
  /** type of data */
  typedef Real Type;
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

}
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
class GammaMixtureManager: public IMixtureManager<GammaMixtureManager<DataHandler> >
{
  public:
    typedef typename hidden::MixtureManagerTraits< GammaMixtureManager >::Type Type;
    typedef typename hidden::MixtureManagerTraits< GammaMixtureManager >::MissingIndexes MissingIndexes;
    typedef typename hidden::MixtureManagerTraits< GammaMixtureManager >::MissingValues MissingValues;
    typedef typename hidden::MixtureManagerTraits< GammaMixtureManager >::Data Data;
    typedef typename hidden::MixtureManagerTraits< GammaMixtureManager >::DataBridgeType DataBridgeType;

    typedef IMixtureManager<GammaMixtureManager> Base;
    using Base::registerDataBridge;
    using Base::getDataBridge;
    using Base::getIdModel;
    using Base::p_handler;

    // All Gamma bridges
    typedef GammaBridge<Clust::Gamma_ajk_bjk_, Data> MixtureBridge_ajk_bjk;
    typedef GammaBridge<Clust::Gamma_ajk_bk_,  Data> MixtureBridge_ajk_bk;
    typedef GammaBridge<Clust::Gamma_ajk_bj_,  Data> MixtureBridge_ajk_bj;
    typedef GammaBridge<Clust::Gamma_ajk_b_,   Data> MixtureBridge_ajk_b;
    typedef GammaBridge<Clust::Gamma_ak_bjk_,  Data> MixtureBridge_ak_bjk;
    typedef GammaBridge<Clust::Gamma_ak_bk_,   Data> MixtureBridge_ak_bk;
    typedef GammaBridge<Clust::Gamma_ak_bj_,   Data> MixtureBridge_ak_bj;
    typedef GammaBridge<Clust::Gamma_ak_b_,    Data> MixtureBridge_ak_b;
    typedef GammaBridge<Clust::Gamma_aj_bjk_,  Data> MixtureBridge_aj_bjk;
    typedef GammaBridge<Clust::Gamma_aj_bk_,   Data> MixtureBridge_aj_bk;
    typedef GammaBridge<Clust::Gamma_a_bjk_,   Data> MixtureBridge_a_bjk;
    typedef GammaBridge<Clust::Gamma_a_bk_,    Data> MixtureBridge_a_bk;

    /** Default constructor, need an instance of a DataHandler.  */
    GammaMixtureManager(DataHandler const& handler): Base(&handler) {}
    /** destructor */
    ~GammaMixtureManager() {}

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
        case Clust::Gamma_ajk_bjk_:
        { static_cast<MixtureBridge_ajk_bjk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gamma_ajk_bk_:
        { static_cast<MixtureBridge_ajk_bk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gamma_ajk_bj_:
        { static_cast<MixtureBridge_ajk_bj*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gamma_ajk_b_:
        { static_cast<MixtureBridge_ajk_b*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gamma_ak_bjk_:
        { static_cast<MixtureBridge_ak_bjk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gamma_ak_bk_:
        { static_cast<MixtureBridge_ak_bk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gamma_ak_bj_:
        { static_cast<MixtureBridge_ak_bj*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gamma_ak_b_:
        { static_cast<MixtureBridge_ak_b*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gamma_aj_bjk_:
        { static_cast<MixtureBridge_aj_bjk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gamma_aj_bk_:
        { static_cast<MixtureBridge_aj_bk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gamma_a_bjk_:
        { static_cast<MixtureBridge_a_bjk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gamma_a_bk_:
        { static_cast<MixtureBridge_a_bk*>(p_mixture)->getMissingValues(missing);}
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
        case Clust::Gamma_ajk_bjk_:
        { static_cast<MixtureBridge_ajk_bjk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ajk_bk_:
        { static_cast<MixtureBridge_ajk_bk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ajk_bj_:
        { static_cast<MixtureBridge_ajk_bj*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ajk_b_:
        { static_cast<MixtureBridge_ajk_b*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ak_bjk_:
        { static_cast<MixtureBridge_ak_bjk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ak_bk_:
        { static_cast<MixtureBridge_ak_bk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ak_bj_:
        { static_cast<MixtureBridge_ak_bj*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_ak_b_:
        { static_cast<MixtureBridge_ak_b*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_aj_bjk_:
        { static_cast<MixtureBridge_aj_bjk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_aj_bk_:
        { static_cast<MixtureBridge_aj_bk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_a_bjk_:
        { static_cast<MixtureBridge_a_bjk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gamma_a_bk_:
        { static_cast<MixtureBridge_a_bk*>(p_mixture)->getParameters(param);}
        break;
        default: // idModel is not implemented
        break;
      }
    }
    /** set the parameters from an IMixture.
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
        case Clust::Gamma_ajk_bjk_:
        { static_cast<MixtureBridge_ajk_bjk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ajk_bk_:
        { static_cast<MixtureBridge_ajk_bk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ajk_bj_:
        { static_cast<MixtureBridge_ajk_bj*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ajk_b_:
        { static_cast<MixtureBridge_ajk_b*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ak_bjk_:
        { static_cast<MixtureBridge_ak_bjk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ak_bk_:
        { static_cast<MixtureBridge_ak_bk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ak_bj_:
        { static_cast<MixtureBridge_ak_bj*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_ak_b_:
        { static_cast<MixtureBridge_ak_b*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_aj_bjk_:
        { static_cast<MixtureBridge_aj_bjk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_aj_bk_:
        { static_cast<MixtureBridge_aj_bk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_a_bjk_:
        { static_cast<MixtureBridge_a_bjk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gamma_a_bk_:
        { static_cast<MixtureBridge_a_bk*>(p_mixture)->setParameters(param);}
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
    IMixture* createMixtureImpl(String const&  modelName, String const& idData, int nbCluster)
    {
      Clust::Mixture idModel = Clust::stringToMixture(modelName);
      return createMixtureImpl(idModel, idData, nbCluster);
    }

  private:
    /** create a concrete mixture and initialize it.
     *  @param idModel Id name of the model
     *  @param idData Id name of the data
     *  @param nbCluster number of cluster of the model
     **/
    IMixture* createMixtureImpl(Clust::Mixture idModel, String const& idData, int nbCluster)
    {
      switch (idModel)
      {
        // gamma models
        case Clust::Gamma_ajk_bjk_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_ajk_bjk)}
        break;
        case Clust::Gamma_ajk_bk_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_ajk_bk)}
        break;
        case Clust::Gamma_ajk_bj_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_ajk_bj)}
        break;
        case Clust::Gamma_ajk_b_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_ajk_b)}
        break;
        case Clust::Gamma_ak_bjk_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_ak_bjk)}
        break;
        case Clust::Gamma_ak_bk_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_ak_bk)}
        break;
        case Clust::Gamma_ak_bj_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_ak_bj)}
        break;
        case Clust::Gamma_ak_b_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_ak_b)}
        break;
        case Clust::Gamma_aj_bjk_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_aj_bjk)}
        break;
        case Clust::Gamma_aj_bk_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_aj_bk)}
        case Clust::Gamma_a_bjk_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_a_bjk)}
        break;
        case Clust::Gamma_a_bk_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_a_bk)}
        default:
          return 0; // 0 if idModel is not implemented
          break;
      }
      return 0; // 0 if idModel is not a STK++ model
    }
};

} // namespace STK

#undef STK_CREATE_MIXTURE

#endif /* STK_GAMMAMIXTUREMANAGER_H */
