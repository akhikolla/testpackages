/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMNodeInput.h  description
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

#ifndef XEM_NODEINPUT_H
#define XEM_NODEINPUT_H

#include "mixmod_iostream/IOStreamUtil.h"

using namespace std;

namespace XEM {

///Main input node in .mixmod file
  class NodeInput : public xmlpp::Document {

public:

	///constructor by default : create an element Project in this xmlpp::Document
	NodeInput();

	///destructor
	~NodeInput();

	///constructor
	NodeInput(string & s);
	NodeInput(Input * input, string & s);
	NodeInput(xmlpp::Element *rootInput);

	///get
	const xmlpp::Element * getRoot()const;

protected:

	//write necessary node
	void writeDataNode(Input * input, string & s);
	void writeSelectVariableNode(Input * input);
	void writeSelectIndividualNode(Input * input);

	//read necessary node
	DataDescription & readDataNode();
	vector<int64_t> readNbClusterNode();

	//parameter
    xmlpp::Element *_rootInput;
};

  inline const xmlpp::Element * NodeInput::getRoot()const {
	return _rootInput;
}

} //end namespace

#endif // XEM_NODEINPUT_H
