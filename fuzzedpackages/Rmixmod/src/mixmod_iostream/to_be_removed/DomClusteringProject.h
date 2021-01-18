/***************************************************************************
                             SRC/MIXMOD_IOSTREAM/XEMDomClusteringProject.h  description
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

#ifndef XEM_DOMCLUSTERINGPROJECT_H
#define XEM_DOMCLUSTERINGPROJECT_H

#include "mixmod_iostream/DomProject.h"

namespace XEM {

///use to create .mixmod file in Clustering case
class DomClusteringProject : public DomProject {

public : 
    
	///constructor by default
	DomClusteringProject();    

	///destructor
	virtual ~DomClusteringProject();

	///constructor by initialization
	DomClusteringProject(xmlpp::Element *root);

	///fill the xmlpp::Document to create the .mixmod file from a ClusteringInput and ClusteringOutput
	void writeClustering(string & s, ClusteringMain * cMain);

	///read a XML file and fill ClusteringInput
	void readClustering(ClusteringInput * cInput);

	///read a XML file and fill ClusteringOutput
	void readClustering(ClusteringOutput * cOutput);
};

} //end namespace

#endif // XEM_DOMCLUSTERINGPROJECT_H
