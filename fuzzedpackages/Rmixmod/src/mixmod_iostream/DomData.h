/***************************************************************************
                             SRC/MIXMOD_IOSTREAM/XEMDomData.h  description
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

#ifndef XEM_DOMDATA_H
#define XEM_DOMDATA_H

#include "mixmod_iostream/IOStreamUtil.h"
#include "mixmod/Kernel/IO/DataDescription.h"
#include "mixmod/Clustering/ClusteringInput.h"

namespace XEM {

///create a .mxd file
class DomData : public xmlpp::Document {

public:
     
	//constructor by default
	DomData();

	//destructor
	~DomData();

	DomData(string & sFilename);
	DomData(const DataDescription & dataFile, string & sFilename);
    DomData(xmlpp::Element *root);

	DataDescription * readDataFile(const string & sFilename);

private:

	//write node in .mxd
    void writeListColumnNode(const vector<ColumnDescription *> & vColumnDescription);

    //parameter
    xmlpp::Element *_root;
};

} //end namespace

#endif // XEM_DOMDATA_H
