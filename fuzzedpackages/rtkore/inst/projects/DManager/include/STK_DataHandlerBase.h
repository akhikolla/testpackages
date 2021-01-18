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
 * Project:  stkpp::DManager
 * created on: 15 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_DataHandler.h
 *  @brief In this file we define the interface base class DataHandlerBase.
 **/


#ifndef STK_DATAHANDLERBASE_H
#define STK_DATAHANDLERBASE_H

#include <map>
#include <string>

#include <STKernel/include/STK_Stream.h>
#include <Sdk/include/STK_IRecursiveTemplate.h>

namespace STK
{
namespace hidden
{
/** @ingroup hidden
 *  The DataHandlerTraits will give the type of container furnished by the
 *  concrete implementations of the DataHandlerBase class.
 *  @note In incoming version of STK++, DataHandlerBase could be renamed as
 *  IDataHandler
 **/
template<class DataHandler, typename Type> struct DataHandlerTraits;

} // namespace hidden

/** @ingroup DManager
 *  A class derived from a DataHandlerBase allows to store various data sets
 *  identified by an idData and an idModel.
 *
 *  A typical usage of a data handler is to store and handle various multiple
 *  sets of data. Each data set possess an id (idData) and can be retrieved by
 *  using it. To each data set is also associated an idModel allowing to know
 *  which kind of (statistical) model is applied in order to model the data set.
 *
 *  @note In incoming version of STK++ DataHandlerBase could be renamed as
 *  IDataHandler
 */
template<class Derived>
class DataHandlerBase: public IRecursiveTemplate<Derived>
{
  protected:
    /** default constructor */
    DataHandlerBase() {}

  public:
    typedef std::map<std::string, std::string> InfoMap;
    /** destructor */
    inline ~DataHandlerBase() {}
    /** @return the map with the idDatas and idModel of the models */
    inline InfoMap const& info() const { return info_;}
    /** @brief Add an info descriptor to the data handler.
     *  An info descriptor is a pair that allow to say that all columns of
     *  the data set(s) handled by the data handler and having the name "idData"
     *  are modeled by the model with model "idModel".
     *  @param idData can be any string given by the user for identifying data.
     *  @param idModel represent the idModel of a given model (can be defined
     *  inside or outside STK++).
     *
     *  @note If the pair (idData, idModel) already exists then addInfo will do nothing.
     *  @return @c false if there exists already an idData matched with an other
     *  idModel, @c true otherwise.
     **/
    bool addInfo(std::string const& idData, std::string const& idModel);
    /** @brief Giving the Id of a data set, find the Id of the model.
     *  @param idData can be any string given by the user for identifying data.
     *  @param idModel The Id of the model associated with the data
     *  (not modified if idData is not present in the map).
     *  @return @c true if there exists an idData in the InfoMap, @c false
     *  otherwise.
     **/
    bool getIdModelName(std::string const& idData, std::string& idModel) const;
    /** write infoMap on os */
    void writeInfo(ostream& os) const;

  protected:
    /** Store the informations  of the mixtures in the form (idData, idModel) with
     * - idData: an arbitrary idData for data.
     * - idModel: a string which represent a (statistical) model.
     * @sa stringToMixture */
    InfoMap info_;
};


/* @brief Add an info descriptor to the data handler.
 *  An info descriptor is a pair that allow to say that all columns of
 *  the data set(s) handled by the data handler and having the name "idData"
 *  are modeled by the model with model "idModel".
 *  @param idData can be any string given by the user for identifying data.
 *  @param idModel represent the idModel of a given model (can be defined
 *  inside or outside STK++).
 *
 *  @note If the pair (idData, idModel) already exists then addInfo will do nothing.
 *  @return @c false if there exists already an idData matched with an other
 *  idModel, @c true otherwise.
 **/
template<class Derived>
bool DataHandlerBase<Derived>::addInfo(std::string const& idData, std::string const& idModel)
{
  // parse descriptor file
  std::pair<InfoMap::iterator,bool> ret;
  // check if identifer is already present
  ret = info_.insert(std::pair<std::string,std::string>(idData, idModel));
  // if name already exists, check if there is incoherence
  if (ret.second==false)
  {
     if (ret.first->second != idModel)
     {
#ifdef STK_DMANAGER_DEBUG
       stk_cerr << _T("In DataHandlerBase::addInfo, There exists an idData with a different idModel.\n");
#endif
       return false;
     }
  }
  return true;
}

/* @brief Giving a the Id of a dataset, find the Id of the model.
 *  @param idData can be any string given by the user.
 *  @param idModel The Id of the model associated with the data
 *  (not modified if idData is not present in the map).
 *  @return @c false if there exists already an idData matched with an other
 *  idModel, @c true otherwise.
 **/
template<class Derived>
bool DataHandlerBase<Derived>::getIdModelName(std::string const& idData, std::string& idModel) const
{
  bool res = false;
  // find idData
  InfoMap::const_iterator it = info_.find(idData);
  if (it != info_.end()) { idModel = it->second; res = true;}
  return res;
}

/* write  info on os */
template<class Derived>
void DataHandlerBase<Derived>::writeInfo(ostream& os) const
{
  // show content
  for (InfoMap::const_iterator it=info_.begin(); it!=info_.end(); ++it)
  os << _T("IdData: ") << it->first << _T(", IdModel: ") << it->second << _T('\n');
}


} // namespace STK

#endif /* STK_DATAHANDLERBASE_H */
