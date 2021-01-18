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

/** @file STK_DiagGaussianMixtureManager.h
 *  @brief In this file we define the DiagGaussianMixtureManager class.
 **/


#ifndef STK_DIAGGAUSSIANMIXTUREMANAGER_H
#define STK_DIAGGAUSSIANMIXTUREMANAGER_H

#include "../STK_Clust_Util.h"
#include "../STK_IMixtureManager.h"

#include "STK_DiagGaussianBridge.h"

#define STK_CREATE_MIXTURE(Data, Bridge) \
          Data* p_data = new Data(idData); \
          p_handler()->getData(idData, p_data->dataij() ); \
          registerDataBridge(p_data); \
          return new Bridge( &(p_data->dataij()), idData, nbCluster);

namespace STK
{

// forward declaration
template<class DataHandler> class DiagGaussianMixtureManager;

namespace hidden
{
/** @ingroup hidden
 *  Specialization for Clust::Categorical_ mixtures */
template <class DataHandler_>
struct MixtureManagerTraits<DiagGaussianMixtureManager<DataHandler_> >
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

} // namespace hidden

/** @ingroup Clustering
 *  @brief A mixture manager is a factory class for injection dependency in the
 *  STK++ derived class of the IMixtureComposer class.
 *
 *  It allows to handle all the creation and initialization stuff needed by the
 *  (bridged) mixture models of the stkpp library.
 *
 *  @tparam DataHandler is any concrete class from the interface DataHandlerBase
 */
template<class DataHandler>
class DiagGaussianMixtureManager: public IMixtureManager<DiagGaussianMixtureManager<DataHandler> >
{
  public:
    typedef typename hidden::MixtureManagerTraits< DiagGaussianMixtureManager >::Type Type;
    typedef typename hidden::MixtureManagerTraits< DiagGaussianMixtureManager >::MissingIndexes MissingIndexes;
    typedef typename hidden::MixtureManagerTraits< DiagGaussianMixtureManager >::MissingValues MissingValues;
    typedef typename hidden::MixtureManagerTraits< DiagGaussianMixtureManager >::Data Data;
    typedef typename hidden::MixtureManagerTraits< DiagGaussianMixtureManager >::DataBridgeType DataBridgeType;

    typedef IMixtureManager< DiagGaussianMixtureManager > Base;
    using Base::registerDataBridge;
    using Base::getDataBridge;
    using Base::getIdModel;
    using Base::p_handler;

    // eiagonal Gaussian bridges
    typedef DiagGaussianBridge<Clust::Gaussian_sjk_,  Data> MixtureBridge_sjk;
    typedef DiagGaussianBridge<Clust::Gaussian_sk_,   Data> MixtureBridge_sk;
    typedef DiagGaussianBridge<Clust::Gaussian_sj_,   Data> MixtureBridge_sj;
    typedef DiagGaussianBridge<Clust::Gaussian_sjsk_, Data> MixtureBridge_sjsk;
    typedef DiagGaussianBridge<Clust::Gaussian_s_,    Data> MixtureBridge_s;

    /** Default constructor, need an instance of a DataHandler.  */
    DiagGaussianMixtureManager(DataHandler const& handler): Base(&handler) {}
    /** destructor */
    ~DiagGaussianMixtureManager() {}

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
        case Clust::Gaussian_sjk_:
        { static_cast<MixtureBridge_sjk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gaussian_sk_:
        { static_cast<MixtureBridge_sk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gaussian_sj_:
        { static_cast<MixtureBridge_sj*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gaussian_sjsk_:
        { static_cast<MixtureBridge_sjsk*>(p_mixture)->getMissingValues(missing);}
        break;
        case Clust::Gaussian_s_:
        { static_cast<MixtureBridge_s*>(p_mixture)->getMissingValues(missing);}
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
        // Gaussian models
        case Clust::Gaussian_sjk_:
        { static_cast<MixtureBridge_sjk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gaussian_sk_:
        { static_cast<MixtureBridge_sk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gaussian_sj_:
        { static_cast<MixtureBridge_sj*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gaussian_sjsk_:
        { static_cast<MixtureBridge_sjsk*>(p_mixture)->getParameters(param);}
        break;
        case Clust::Gaussian_s_:
        { static_cast<MixtureBridge_s*>(p_mixture)->getParameters(param);}
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
        // Gaussian models
        case Clust::Gaussian_sjk_:
        { static_cast<MixtureBridge_sjk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gaussian_sk_:
        { static_cast<MixtureBridge_sk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gaussian_sj_:
        { static_cast<MixtureBridge_sj*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gaussian_sjsk_:
        { static_cast<MixtureBridge_sjsk*>(p_mixture)->setParameters(param);}
        break;
        case Clust::Gaussian_s_:
        { static_cast<MixtureBridge_s*>(p_mixture)->setParameters(param);}
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
    { return createMixtureImpl(Clust::stringToMixture(modelName), idData, nbCluster);}


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
        // Gaussian models
        case Clust::Gaussian_sjk_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_sjk)}
        break;
        case Clust::Gaussian_sk_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_sk)}
        break;
        case Clust::Gaussian_sj_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_sj)}
        break;
        case Clust::Gaussian_sjsk_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_sjsk)}
        break;
        case Clust::Gaussian_s_:
        { STK_CREATE_MIXTURE(DataBridgeType, MixtureBridge_s)}
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

#endif /* STK_DIAGGAUSSIANMIXTUREMANAGER_H */
