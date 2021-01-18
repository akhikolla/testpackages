/***************************************************************************
                             SRC/MIXMOD_IOSTREAM/XEMDomLabel.h  description
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

#ifndef XEM_DOMLABEL_H
#define XEM_DOMLABEL_H
#include "mixmod_iostream/IOStreamUtil.h"
#include "mixmod/Clustering/ClusteringInput.h"
#include "mixmod/Kernel/IO/LabelDescription.h"

namespace XEM {

  class DomLabel : public xmlpp::Document {
	
public:

	//constructor by default
	DomLabel();

	//destructor
	~DomLabel();

	DomLabel(string & sFilename);
	DomLabel(LabelDescription * labelDescription, string str);
	DomLabel(Partition * partition, string & sFilename);
	DomLabel(xmlpp::Element *root);

	unique_ptr<LabelDescription> readLabel(string sFilename); 
	unique_ptr<LabelDescription> readLabelAsData(string sFilename);   
	void writeListColumnNode(const vector<ColumnDescription *> & vColumnDescription);

private:

	//label
    xmlpp::Element *_root;
};

} //end namespace

#endif // XEM_DOMLABEL_H
