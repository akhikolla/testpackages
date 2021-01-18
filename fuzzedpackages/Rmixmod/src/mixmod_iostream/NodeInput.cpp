/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMNodeInput.cpp  description
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

#include "mixmod_iostream/NodeInput.h"
#include "mixmod_iostream/DomData.h"
#include <stdint.h>
#include <algorithm>
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/GaussianData.h"

namespace XEM {

  NodeInput::NodeInput() : xmlpp::Document() {
    _rootInput = create_root_node("Input");
  }
  
  NodeInput::~NodeInput() {
  }

  NodeInput::NodeInput(string & s) : xmlpp::Document() {
    set_internal_subset(s, "", "");
  }
  
  NodeInput::NodeInput( xmlpp::Element *rootInput) {
    _rootInput = create_root_node_by_import(rootInput);
  }

//common part between Clustering & AD
  NodeInput::NodeInput(Input * input, string & s) : xmlpp::Document() {

	//_rootInput = createElement( "Input" );
    _rootInput = create_root_node( "Input");
	writeDataNode(input, s);
	writeSelectVariableNode(input);
	writeSelectIndividualNode(input);
}

  void NodeInput::writeDataNode(Input * input, string & s) {

	//Data
    xmlpp::Element *data = _rootInput->add_child("Data");
	//dataFilename
    string dataFilename = s + ".mxd";
    data->add_child_text(dataFilename);
	DomData doc(input->getDataDescription(), s);
}

  void NodeInput::writeSelectVariableNode(Input * input) {

	//SelectVariable
    xmlpp::Element *select = _rootInput->add_child("SelectVariable");
	for (int64_t i = 0; i < input->getPbDimension(); ++i) {
      xmlpp::Element *variable = select->add_child("Variable");
      variable->add_child_text(std::to_string(i+1));
      //By default, all variables are selectionned 
	}

}

  void NodeInput::writeSelectIndividualNode(Input * input) {
  }

  DataDescription & NodeInput::readDataNode() {

	if (!_rootInput) throw;
    xmlpp::Element *elementData = dynamic_cast<xmlpp::Element*>(_rootInput->get_first_child("Data"));
    string filename;
    if (!elementData) throw;
      //absolute filename is the name of .mxd file
    filename = elementData->get_child_text()->get_content();
    ValidateSchema(filename, IOStreamXMLFile::Data);
    // throw IOStreamErrorType::badLoadXML;

    DomData data;

    //creation of dataDescription to construct the input and return it
    return (*data.readDataFile(filename));
  }

  vector<int64_t> NodeInput::readNbClusterNode() {

	//declaration of vector of nbCLuster
	vector<int64_t> vNbCluster(0);

	if (!_rootInput) throw;

    //node ListNbCluster
    xmlpp::Element *elementListNbCluster = dynamic_cast<xmlpp::Element*>(_rootInput->get_first_child("ListNbCluster"));
    //node Nbcluster
    int64_t compt = 0;
    auto children = elementListNbCluster->get_children();
    //cross the listNbCLuster node to find the different nbCluster
    for (auto it=children.begin(); it != children.end(); ++it){
      xmlpp::Element *elementNbCluster = dynamic_cast<xmlpp::Element*>(*it);
      if(!elementNbCluster) continue;
      //fill the vector of nbCluster
      int64_t nbCluster = std::stoll(elementNbCluster->get_child_text()->get_content());
      if (find(vNbCluster.begin(), vNbCluster.end(), nbCluster) == vNbCluster.end()) {
        vNbCluster.push_back(nbCluster);
      }
      compt++;
    }

    return vNbCluster;
  }

} //end namespace
