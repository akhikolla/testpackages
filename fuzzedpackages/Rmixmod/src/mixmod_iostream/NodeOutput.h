/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMNodeOutput.h  description
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

#ifndef XEM_NODEOUTPUT_H
#define XEM_NODEOUTPUT_H

#include "mixmod_iostream/IOStreamUtil.h"

using namespace std;

namespace XEM {

///Main output node in .mixmod file
  class NodeOutput : public xmlpp::Document {

public:

	///constructor by default : create an element Project in this xmlpp::Document
	NodeOutput();

	///destructor
	~NodeOutput();

	///constructor
	NodeOutput(string & s);
	//NodeOutput(ClusteringOutput * output, string & s);
	NodeOutput(xmlpp::Element *rootOutput);

	///get
	const xmlpp::Element * getRoot()const;

protected:

	//write necessary node

	//read necessary node

	//parameter
    xmlpp::Element *_rootOutput;

};

  inline const xmlpp::Element * NodeOutput::getRoot()const {
	return _rootOutput;
}

} //end namespace

#endif // XEM_NODEOUTPUT_H
