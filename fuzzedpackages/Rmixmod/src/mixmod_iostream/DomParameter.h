/***************************************************************************
                             SRC/MIXMOD_IOSTREAM/XEMDomParameter.h  description
    copyright            : (C) MIXMOD Team - 2001-2011
    email                : contact@mixmod.org
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

    All informations available on : http://www.mixmod.org                                                                                               
***************************************************************************/

#ifndef XEM_DOMPARAMETER_H
#define XEM_DOMPARAMETER_H


#include "mixmod_iostream/IOStreamUtil.h"
#include "mixmod/Kernel/Parameter/BinaryParameter.h"
#include "mixmod/Kernel/Parameter/GaussianParameter.h"
#include "mixmod/Clustering/ClusteringInput.h"
#include "mixmod/Clustering/ClusteringStrategy.h"
#include "mixmod/Clustering/ClusteringStrategyInit.h"
#include "mixmod_iostream/NodeOpOutput.h"

namespace XEM {

class DomParameter : public xmlpp::Document {

public:

	//constructor by default
	DomParameter();

	//destructor
	~DomParameter();

	DomParameter(string & sFilename);    
	DomParameter(ParameterDescription* parameterDescription, string sFilename);
	DomParameter(ClusteringInput* cInput, string sFilename);    
	DomParameter(PredictInput* cInput, string sFilename);    
    DomParameter(xmlpp::Element *element);

	///read Parameter
	ParameterDescription * readParameter(string sFilename);

private:
	
	//parameter
    xmlpp::Element *_root;
};

} //end namespace

#endif // XEM_DOMPARAMETER_H
