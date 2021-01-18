/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMNodeOutput.cpp  description
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

#include "mixmod_iostream/NodeOutput.h"
#include "mixmod_iostream/DomData.h"
#include <stdint.h>
#include <algorithm>

namespace XEM {

  NodeOutput::NodeOutput() : xmlpp::Document() {
	_rootOutput = create_root_node( "ListOutput" );    
  }

  NodeOutput::~NodeOutput() {
  }
 
  NodeOutput::NodeOutput(string & s) : xmlpp::Document() {
    set_internal_subset(s, "", "");

  }
  
  NodeOutput::NodeOutput( xmlpp::Element *rootOutput) {
	_rootOutput = create_root_node_by_import(rootOutput);
  }

//common part between Clustering & DA
/*
  NodeOutput::NodeOutput(ClusteringOutput* output, string& s) : xmlpp::Document() {
	_rootOutput = create_root_node("ListOutput");        
  }
*/
} //end namespace
