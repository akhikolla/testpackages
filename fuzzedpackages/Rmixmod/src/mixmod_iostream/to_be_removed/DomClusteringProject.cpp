/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMDomClusteringProject.cpp  description
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

#include "mixmod_iostream/DomClusteringProject.h"
#include "mixmod_iostream/NodeClusteringInput.h"
#include "mixmod_iostream/NodeClusteringOutput.h"
#include "mixmod/Clustering/ClusteringMain.h"

namespace XEM {

DomClusteringProject::DomClusteringProject() : DomProject() {
}

DomClusteringProject::~DomClusteringProject() {
}

  DomClusteringProject::DomClusteringProject(xmlpp::Element *root) : DomProject(root) {
}

void DomClusteringProject::writeClustering(string& s, ClusteringMain * cMain) {
  _root->set_namespace_declaration("http://www.w3.org/2001/XMLSchema-instance", "xsi");
  _root->set_attribute("type", "Clustering", "xsi");
	//insertion of input
	NodeInput * inputNode = new NodeClusteringInput(cMain->getInput(), s);
	_root->import_node(inputNode->getRoot());

	//insertion of output
	if (cMain->getOutput() != NULL) {
		NodeOutput * outputNode = new NodeClusteringOutput(cMain->getOutput(), s);
		_root->import_node(outputNode->getRoot());
	}
}

void DomClusteringProject::readClustering(ClusteringInput* cInput) {
	//read the XML file which the case is clustering then fill the ClusteringInput

	//root in input
  xmlpp::Element* rootInput = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Input")); 
  NodeClusteringInput nodeInput(rootInput);

	//fill cInput from nodeInput (input node)
	nodeInput.readClustering(*cInput);
}

void DomClusteringProject::readClustering( ClusteringOutput * cOutput ) {
	//read the XML file which the case is clustering then fill the ClusteringOutput

	//root in input
  xmlpp::Element* rootOutput = dynamic_cast<xmlpp::Element*>(_root->get_first_child("ListOutput"));
  xmlpp::Element* nodeElementOutput;
  vector<ClusteringModelOutput *> vClusteringModelOutput;
  if (rootOutput) {
    //xmlpp::Node::NodeList
    auto children = rootOutput->get_children();
    for (auto it=children.begin(); it != children.end(); ++it){
      xmlpp::Element* nodeElementOutput = dynamic_cast<xmlpp::Element*>(*it);
      if(!nodeElementOutput) continue;
      NodeClusteringOutput nodeOutput(nodeElementOutput);
      ClusteringModelOutput * modelOutput = nodeOutput.readClustering();
      vClusteringModelOutput.push_back(modelOutput);
    }
    cOutput->setClusteringModelOutput(vClusteringModelOutput);
  } else {
    delete cOutput;
    cOutput = NULL;
  }
}

} //end namespace
