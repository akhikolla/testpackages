/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 8 août 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_RDataHandler.cpp
 *  @brief In this file we implement the RDataHandler class.
 **/

#include "../inst/projects/MixAll/RDataHandler.h"

namespace STK
{

bool RDataHandler::addType(String const& idData, int Rtype)
{
  // parse descriptor file
  std::pair<InfoType::iterator,bool> ret;
  // check if identifier is already present
  ret = infoType_.insert(std::pair<String,int>(idData, Rtype));
  // if name already exists, check if there is incoherence
  if (ret.second==false)
  {
     if (ret.first->second != Rtype)
     {
#ifdef RTK_DEBUG
       stk_cerr << _T("In RDataHandler::addType, existing idData with a different Rtype.\n");
#endif
       return false;
     }
  }
  return true;
}


} // namespace STK
