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
 **/

/** @file STK_DataBridge.h
 *  @brief In this file we define the data wrapper class
 **/

#ifndef STK_DATABRIDGE_H
#define STK_DATABRIDGE_H

#include "STK_IDataBridge.h"

#include <Arrays/include/STK_ITContainer2D.h>

namespace STK
{

/** @ingroup Clustering
 *  @brief bridge a data set in order to handle its missing values.
 *
 * @tparam Data The data bridged by the DataBridge class
 */
template<class Data>
class DataBridge: public IDataBridge
{
  public:
    /** data set */
    Data dataij_;

    typedef typename hidden::Traits<Data>::Type Type;
    typedef std::vector<std::pair<int,int> > MissingIndexes;
    typedef std::vector< std::pair<std::pair<int,int>, Type > > MissingValues;

    using IDataBridge::v_missing_;

    /** default constructor. */
    DataBridge( std::string const& idData): IDataBridge(idData), dataij_() {}
    /** constructor with data. */
    DataBridge( std::string const& idData, Data const& dataij): IDataBridge(idData), dataij_(dataij) {}
    /** copy constructor (Warning: will copy the data set)
     *  @param bridge the DataBridge to copy
     **/
    DataBridge( DataBridge const& bridge): IDataBridge(bridge), dataij_(bridge.dataij_) {}
    /** destructor */
    virtual ~DataBridge() {}
    /** @return data set */
    Data const& dataij() const { return dataij_;}
    /** @return data set */
    Data& dataij() { return dataij_;}

    /** get (estimated) missing values of the data set */
    void getMissingValues( MissingValues& data) const;

 protected:
    /** utility function for lookup the data set and find missing values coordinates.
     *  @return the number of missing values
     **/
    virtual size_t findMissing();
};

template<class Data>
void DataBridge<Data>::getMissingValues( MissingValues& data) const
{
  data.resize(v_missing_.size());
  for(size_t i = 0; i< v_missing_.size(); ++i)
  {
    data[i].first  = v_missing_[i];
    data[i].second = dataij_(v_missing_[i].first, v_missing_[i].second);
  }
}

template<class Data>
std::vector< std::pair<int,int> >::size_type DataBridge<Data>::findMissing()
{
#ifdef STK_DMANAGER_VERBOSE
  stk_cout << _T("IDataBridge::Entering findMissing()\n");
#endif
  for (int j=dataij_.beginCols(); j< dataij_.endCols(); ++j)
  {
    for (int i=dataij_.beginRows(); i< dataij_.endRows(); ++i)
    {
      if (Arithmetic<Type>::isNA(dataij_(i,j)))
      { v_missing_.push_back(std::pair<int,int>(i,j));}
    }
  }
 return v_missing_.size();
#ifdef STK_DMANAGER_VERBOSE
  stk_cout << _T("IDataBridge::findMissing() terminated, nbMiss= ") << v_missing_.size() << _T("\n");
#endif
}

} // namespace STK

#endif /* STK_DATABRIDGE_H */
