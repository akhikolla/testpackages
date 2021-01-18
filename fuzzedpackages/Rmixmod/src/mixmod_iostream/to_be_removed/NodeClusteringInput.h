/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMNodeClusteringInput.h  description
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

#ifndef XEM_NODECLUSTERINGINPUT_H
#define XEM_NODECLUSTERINGINPUT_H

#include "mixmod_iostream/NodeInput.h"
#include "mixmod/Clustering/ClusteringStrategy.h"
#include "mixmod/Kernel/Algo/Algo.h"
#include "mixmod/Kernel/IO/Partition.h"

namespace XEM {

///Input node in .mixmod file in case clustering
class NodeClusteringInput : public NodeInput {

public:

	///Constructor & Destructor
	NodeClusteringInput();
	~NodeClusteringInput();
	NodeClusteringInput(ClusteringInput * input, string & s);
	NodeClusteringInput(xmlpp::Element * rootInput);

	///writer node
	void writeListModel(ClusteringInput * input);
	void writeNbClusterNode(ClusteringInput * input);
	void writeStrategyNode(ClusteringInput * input, string & s);
	void writeAlgoNode(xmlpp::Element *listAlgo, const Algo * algo);
	void writeInitNode(xmlpp::Element *strategyElement, ClusteringInput * input, string & s);
	void writeCriterionNode(ClusteringInput * input);
	void writePartitionNode(ClusteringInput * input, string & s);
	void writeWeightsNode(ClusteringInput * input, string & s);

	///read Node
	void readClustering(ClusteringInput & input);

	void readModelNode(ClusteringInput & input);
	void readStrategyNode(ClusteringInput & input);
	void readAlgoNode(ClusteringStrategy * strat, xmlpp::Element * n);
	void readInitNode(ClusteringStrategy * strat, xmlpp::Element * n);
	void readCriterionNode(ClusteringInput & input);
	void readPartitionNode(ClusteringInput & input);
	void readWeightsNode(ClusteringInput & input);
    static void setInitPartition(string sFilename, ClusteringStrategy * strat);
};

} //end namespace

#endif // XEM_DOMCLUSTERING_H
