/**
 * Project:  Rmixmod
 * created on: 4 avr. 2011
 * Purpose:  Create the main for the mixmod call.
 * Author:   
 *
 **/

/***************************************************************************
                             ClusteringMain.cpp  description
                             ------------------------
    copyright            : (C) MIXMOD Team - 2001-2003-2004-2005-2006-2007-2008-2009
    email                : mixmod@univ-fcomte.fr
 ***************************************************************************/

/***************************************************************************
    This file is part of MIXMOD

    MIXMOD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MIXMOD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MIXMOD.  If not, see <http://www.gnu.org/licenses/>.

 ***************************************************************************/
/** @file ClusteringMain.cpp
 *  @brief In this file we implement the wrapper .
 **/


#include <vector>
#include <string>

#include "mixmod/Utilities/Util.h"
#include "mixmod/Utilities/exceptions/OtherException.h"
#ifdef RMIXMOD_XML
#include "mixmod_iostream/IOStreamUtil.h"
#endif
#include <Rcpp.h>



/** This is the main method doing the interface between R and Mixmod for checking
 * the content of an XML file (i.e. clustering, learn or predict)
 *  The method will create a matrix in the format of mixmod and copy the data
 *  inside
 *
 * @param xem R S4 object 
 */
RcppExport SEXP xMain( SEXP xem )
{
  BEGIN_RCPP
 // wrap S4 object
  Rcpp::S4 mixmodObj(xem);
#ifdef RMIXMOD_XML
  Rcpp::StringVector XmlFile(mixmodObj.slot("xmlFile"));
  std::string xmlFile = "";
  if(XmlFile.size()>0)  xmlFile = Rcpp::as< std::vector<std::string> >(XmlFile)[0];

  switch(XEM::getProjectType(xmlFile)) {
  case XEM::ProjectType::Clustering:
    mixmodObj.slot("xmlType") = Rcpp::StringVector("clustering");
    break;
  case XEM::ProjectType::Learn:
    mixmodObj.slot("xmlType") = Rcpp::StringVector("learn");
    break;
  case XEM::ProjectType::Predict:
    mixmodObj.slot("xmlType") = Rcpp::StringVector("predict");    
  }
      
#else
    THROW(XEM::OtherException, XEM::xmlFeaturesNotAvailable);
#endif    

  // return final output
  return mixmodObj;
  END_RCPP  

}
