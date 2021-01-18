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

/** @file STK_IMixtureManager.h
 *  @brief In this file we define the Interface class IMixtureManager
 **/


#ifndef STK_IMIXTUREMANAGER_H
#define STK_IMIXTUREMANAGER_H

#include <DManager/include/STK_DataHandlerBase.h>
#include <DManager/include/STK_DataBridge.h>

#include <Arrays/include/STK_Array2D.h> // for get and set parameters

#include "STK_IMixture.h"

namespace STK
{
/** @ingroup Clustering
 *  @brief Interface base class for mixture managers.
 *
 *  A mixture manager is a factory class for injection dependency in the
 *  STK++ derived class of the IMixtureComposer:
 *  - It handles all the creation and initialization stuff needed by mixture models,
 *  - It allows to get parameters and imputed missing values from specific mixtures,
 *  - It allows also to set parameters to a specific mixture,
 *  - all data set are enclosed in a DataBridge structure and stored in vector v_data_
 *
 *  The pseudo pure virtual method to implement in derived classes are
 *  @code
 *    void getMissingValuesImpl(
 *    void getParametersImpl(IMixture* p_mixture, ArrayXX& data) const;
 *    void setParametersImpl(IMixture* p_mixture, ArrayXX const& data) const;
 *    IMixture* createMixtureImpl(String const& modelName, String const& idData, int nbCluster);
 *  @endcode
 *
 *  @tparam DataHandler any concrete class from the interface STK::DataHandlerBase
 */
template<class Derived>
class IMixtureManager: public IRecursiveTemplate<Derived>
{
  public:
    typedef typename hidden::MixtureManagerTraits< Derived >::DataHandler DataHandler;
    typedef typename hidden::MixtureManagerTraits< Derived >::Type Type;
    typedef typename hidden::MixtureManagerTraits< Derived >::MissingIndexes MissingIndexes;
    typedef typename hidden::MixtureManagerTraits< Derived >::MissingValues MissingValues;
    typedef typename hidden::MixtureManagerTraits< Derived >::Data Data;
    typedef typename hidden::MixtureManagerTraits< Derived >::DataBridgeType DataBridgeType;

    /** Default constructor, need an instance of a DataHandler. */
    IMixtureManager(DataHandler const* const p_handler);
    /** destructor */
    ~IMixtureManager();

    /** @return constant pointer on the data handler */
    DataHandler const* const p_handler() const { return p_handler_;}

    /** Utility function allowing to find the idModel from the idData
     *  @param idData the id name of the data we want the idModel
     *  @return the idModel
     **/
    Clust::Mixture getIdModel( String const& idData) const;
    /** Utility function allowing to find the idModel name from the idData
     *  @param idData the id name of the data we want the idModel
     *  @return the idModel name
     **/
    String getIdModelName( String const& idData) const;

    /** @brief create a mixture and initialize it.
     *  This method get the modelName from the DataHandler and then delegate
     *  the concrete creation to derived class using the pseudo pure virtual method
     *   @c createMixture( modelName, idData, nbCluster).
     *  @param idData name of the model
     *  @param nbCluster number of cluster of the model
     *  @return 0 if the idData is not find, the result of
     *  @c createMixture( modelName, idData, nbCluster) otherwise.
     **/
    IMixture* createMixture(String const& idData, int nbCluster);
    /** @brief register a data bridge to the IMixtureManager.
     *  For each mixture created and registered, a data bridge is created
     *  and registered so that it will be deleted when the mixture itself is
     *  deleted.
     *  @param p_data a pointer on the data manager
     **/
    void registerDataBridge(IDataBridge* p_data);
    /** release a data bridge from v_data_.
     *  @param idData name of the data set to release
     **/
    void releaseDataBridge(String const& idData);

    /** get the wrapper for any kind of data set using its Id
     *  @param idData Id name of the data set attached to the mixture
     *  @return a constant reference on the array with the data set
     **/
    Data const& getData( String const& idData) const;

    // pure virtual methods
    /** get the missing values
     *  @param p_mixture pointer on the mixture
     *  @param missing array with the indexes and the missing values
     **/
    inline void getMissingValues(IMixture* p_mixture, MissingValues& missing) const
    {this->asDerived().getMissingValuesImpl(p_mixture, missing);}
    /** get the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param data the array to return with the parameters
     **/
    inline void getParameters(IMixture* p_mixture, ArrayXX& data) const
    { this->asDerived().getParametersImpl(p_mixture, data);}
    /** set the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param data the array with the parameters to set
     **/
    inline void setParameters(IMixture* p_mixture, ArrayXX const& data) const
    { this->asDerived().setParametersImpl(p_mixture, data);}

  protected:
    /** Utility lookup function allowing to find a DataBridge from its idData
     *  @param idData the id name of the mixture we want to get
     *  @return a pointer on the DataBridge
     **/
    IDataBridge* getDataBridge( String const& idData) const;

  private:
    /** A pointer on the concrete instance of the data handler */
    DataHandler const* const p_handler_;
    /** vector of pointers to the data components */
    std::vector<IDataBridge*> v_data_;

    /** create a concrete mixture and initialize it.
     *  @param modelName, idData strings with the Id name of the model and of the data
     *  @param nbCluster number of cluster of the model
     **/
    inline IMixture* createMixture(String const& modelName, String const& idData, int nbCluster)
    { return this->asDerived().createMixtureImpl(modelName, idData, nbCluster);}
};

/* Utility function allowing to find the idModel from the idData
 *  @param idData the id name of the data we want the idModel
 *  @return the idModel
 **/
template<class Derived>
Clust::Mixture IMixtureManager<Derived>::getIdModel( String const& idData) const
{
  String modelName;
  if (!p_handler()->getIdModelName( idData, modelName))
  {
#ifdef STK_MIXTURE_VERY_VERBOSE
    stk_cout << _T("In IMixtureManager::getIdModel, fail to get idData = ") << idData << _T("\n");
#endif
    return Clust::unknown_mixture_;
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("In IMixtureManager::getIdModel, success to get idData = ") << idData << _T("\n");
  stk_cout << _T("In IMixtureManager::getIdModel, modelName = ") << modelName << _T("\n");
#endif
  return Clust::stringToMixture(modelName);
}

/* Default constructor, need an instance of a DataHandler.  */
template<class Derived>
IMixtureManager<Derived>::IMixtureManager( DataHandler const* const p_handler)
                                         : p_handler_(p_handler)
{}
/* destructor */
template<class Derived>
IMixtureManager<Derived>::~IMixtureManager()
{
  typedef std::vector<IDataBridge*>::iterator DataIterator;
  for (DataIterator it = v_data_.begin() ; it != v_data_.end(); ++it)
  { delete (*it);}
}
/* Utility function allowing to find the idModel name from the idData
 *  @param idData the id name of the data we want the idModel
 *  @return the idModel name
 **/
template<class Derived>
String IMixtureManager<Derived>::getIdModelName( String const& idData) const
{
  String modelName;
  if (!p_handler_->getIdModelName( idData, modelName))
  {
#ifdef STK_MIXTURE_VERY_VERBOSE
    stk_cout << _T("In IMixtureManager::getIdModelName, fail to get idData = ") << idData << _T("\n");
#endif
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("In IMixtureManager::getIdModeName, success to get idData = ") << idData << _T("\n");
  stk_cout << _T("In IMixtureManager::getIdModel, modelName = ") << modelName << _T("\n");
#endif
  return modelName;
}

/* create a mixture and initialize it.
 *  @param idData name of the model
 *  @param nbCluster number of cluster of the model
 *  @return 0 if the idData is not find, the result of @c createMixture( idModel, idData, nbCluster)
 *  otherwise.
 **/
template<class Derived>
IMixture* IMixtureManager<Derived>::createMixture(String const& idData, int nbCluster)
{
  String modelName;
  if (!p_handler_->getIdModelName( idData, modelName)) { return 0;};
  return createMixture( modelName, idData, nbCluster);
}
/* @brief register a data manager to the IMixtureManager.
 *  For each mixture created and registered, a data manager is created
 *  and registered so that it will be deleted when the mixture itself is
 *  deleted.
 *  @param p_data a pointer on the data manager
 **/
template<class Derived>
void IMixtureManager<Derived>::registerDataBridge(IDataBridge* p_data)
{ v_data_.push_back(p_data);}
/* release a data set from v_data_.
 *  @param idData name of the data set to release
 **/
template<class Derived>
void IMixtureManager<Derived>::releaseDataBridge(String const& idData)
{
  typedef std::vector<IDataBridge*>::iterator DataIterator;
  for (DataIterator it = v_data_.begin(); it != v_data_.end(); ++it)
  { if ((*it)->idData() == idData) {delete (*it); v_data_.erase(it); break;}}
}

/* Utility lookup function allowing to find a DataBridge from its idData
 *  @param idData the id name of the mixture we want to get
 *  @return a pointer on the DataBridge
 **/
template<class Derived>
IDataBridge* IMixtureManager<Derived>::getDataBridge( String const& idData) const
{
  typedef std::vector<IDataBridge*>::const_iterator ConstDataIterator;
  for (ConstDataIterator it = v_data_.begin(); it != v_data_.end(); ++it)
  {  if ((*it)->idData() == idData) return (*it);}
  return 0;
}

/* get the wrapper for any kind of data set using its Id
 *  @param idData Id name of the data set attached to the mixture
 *  @return a constant reference on the array with the data set
 **/
template<class Derived>
typename hidden::MixtureManagerTraits< Derived >::Data const& IMixtureManager<Derived>::getData( String const& idData) const
{
  IDataBridge* p_bridge = getDataBridge(idData);
  if (p_bridge)
    return static_cast<DataBridgeType*>(p_bridge)->dataij();
  else
    STKRUNTIME_ERROR_1ARG(IMixtureManager::getData,idData,data bridge does not exist);
}

} // namespace STK


#endif /* STK_IMIXTUREMANAGER_H */
