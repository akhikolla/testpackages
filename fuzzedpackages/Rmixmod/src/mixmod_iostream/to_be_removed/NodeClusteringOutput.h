/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMNodeClusteringOutput.h  description
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

#ifndef XEM_NODECLUSTERINGOUTPUT_H
#define XEM_NODECLUSTERINGOUTPUT_H

#include "mixmod_iostream/NodeOutput.h"

namespace XEM {

class LabelDescription;
class ProbaDescription;
class ParameterDescription;

///Output node in .mixmod file in case clustering
class NodeClusteringOutput : public NodeOutput {

public:

	///Constructor & Destructor
	NodeClusteringOutput();
	~NodeClusteringOutput();
	NodeClusteringOutput(ClusteringOutput* output, string& s);
	NodeClusteringOutput(xmlpp::Element * rootOutput);

	///writer node

	///read Node one <Output>
	ClusteringModelOutput * readClustering();

private:

	///read Parameter node
	ParameterDescription * readParameter(string sFilename);

	///read Label node
	unique_ptr<LabelDescription> readLabel(string sFilename);

	///read Proba node
	unique_ptr<ProbaDescription> readProba(string sFilename);

	///write output node 
	void writeOutput(ClusteringModelOutput* output, 
			const std::vector<CriterionName> & criterionName, string str, int64_t numOutput);
};

} //end namespace

#endif // XEM_NODECLUSTERINGOUTPUT_H
