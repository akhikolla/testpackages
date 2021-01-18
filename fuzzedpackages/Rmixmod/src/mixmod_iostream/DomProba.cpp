/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMDomProba.cpp  description
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

#include "mixmod_iostream/DomProba.h"
#include "mixmod_iostream/IOStreamUtil.h"
#include "mixmod/Kernel/IO/ColumnDescription.h"

namespace XEM {

  DomProba::DomProba() {
  }

  DomProba::~DomProba() {
  }

  DomProba::DomProba(ProbaDescription* probaDescription, string sFilename) : xmlpp::Document() {
    _root = create_root_node("Data");
    _root->set_namespace_declaration("http://www.w3.org/2001/XMLSchema-instance", "xsi");

	//QDomText text;
    xmlpp::Element* new_elt = NULL;
	//Name
	if ( !probaDescription->getInfoName().empty() ) {
      new_elt = _root->add_child("Name");
      new_elt->add_child_text(probaDescription->getInfoName());
	}

	//NbSample
    new_elt = _root->add_child("NbSample");
    new_elt->add_child_text(std::to_string(probaDescription->getNbSample()));    

	//NbColumn
	//   QDomElement nbCol = createElement("NbColumn");  
	//   text = createTextNode(QString::number(probaDescription->getNbColumn()));
	//   nbCol.appendChild(text);        
	//   _root.appendChild(nbCol); 

	//format
    new_elt = _root->add_child("Format");
    new_elt->add_child_text(FormatNumericFileToString(probaDescription->getFormat()));    
    

	//datafilename
	probaDescription->saveNumericValues(getAbsolutePath(sFilename + ".txt"));
    new_elt = _root->add_child("DataFilename");
    new_elt->add_child_text(sFilename + ".txt");        

	writeListColumnNode(probaDescription->getAllColumnDescription());
	//write new file .mxd to describe data
    Glib::ustring filename = getAbsolutePath(sFilename + ".mxd");
    removeIfExists(filename);
    write_to_file(filename);    
}

  DomProba::DomProba(string& sFilename) : xmlpp::Document() {
    set_internal_subset(sFilename, "", "");
  }

  DomProba::DomProba(xmlpp::Element *root) {
    _root = create_root_node_by_import(root);
  }

  unique_ptr<ProbaDescription> DomProba::readProba(string sFilename) {
	//-------
	//load file in variable "doc"
	//-------
    xmlpp::DomParser parser;
    parser.parse_file(sFilename);
    xmlpp::Document *doc = parser.get_document();    
    xmlpp::Element* _root = doc->get_root_node();
	if ( _root->get_name() != "Data" ) return 0;
    
    //------------------------
    //Declaration of variables
    //------------------------
    xmlpp::Element *elementName, *elementNbSample, *elementNbColumn, 
      *elementFormat, *elementProbaFilename, *elementListColumn;

    //name
    elementName = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Name"));
    string sName = "";
    if (elementName) {
      sName = elementName->get_child_text()->get_content();
    }

    //nbSample
    elementNbSample = dynamic_cast<xmlpp::Element*>(_root->get_first_child("NbSample"));
    int64_t nbSample = std::stoll(elementNbSample->get_child_text()->get_content());
    //nbColumn
    
    //Format
    elementFormat = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Format"));
    FormatNumeric::FormatNumericFile format = 
      StringToFormatNumericFile(elementFormat->get_child_text()->get_content());

    //Parameter Filename
    elementProbaFilename = dynamic_cast<xmlpp::Element*>(_root->get_first_child("DataFilename"));
    string dataFilename = elementProbaFilename->get_child_text()->get_content();
    //ListColumn
    elementListColumn = dynamic_cast<xmlpp::Element*>(_root->get_first_child("ListColumn"));
    int64_t nbCluster = 0;
    
    //Each Column
    if (elementListColumn) {
      auto children = elementListColumn->get_children("Column");
      for (auto it=children.begin(); it != children.end(); ++it){
        xmlpp::Element* elementColumn = dynamic_cast<xmlpp::Element*>(*it);
        if(!elementColumn) continue;
        nbCluster++;
        if (elementColumn->get_attribute_value("type", "xsi") != "Quantitative")
          throw IOStreamErrorType::badElementInDataXML;
        
        //num of proba column ( = cluster)
        //	int64_t num = elementColumn.attribute("Num").toLongLong();

        //name [TODO: unused ?!]
        sName = elementColumn->get_attribute_value("Name");
      }
    }
    unique_ptr<ProbaDescription> newProbaDescription(new ProbaDescription(nbSample, nbCluster, format, dataFilename, sName));
    return newProbaDescription;
    
//  return new ProbaDescription(nbSample, nbCluster, format, dataFilename, sName);
    //}//

	//else return 0;
  }

void DomProba::writeListColumnNode(const vector< ColumnDescription* >& vColumnDescription) {
	//ListColumn
	int64_t nbColumn = vColumnDescription.size();
    if(!nbColumn) return;
    xmlpp::Element* listColumn = _root->add_child("ListColumn");
	for (int64_t i = 0; i < nbColumn; ++i) {
      xmlpp::Element *column = listColumn->add_child("Column");
      //position of variable 
      column->set_attribute("Num", std::to_string(vColumnDescription[i]->getIndex() + 1));
      //name
      if (!vColumnDescription[i]->getName().empty()) {
        column->set_attribute("Name", vColumnDescription[i]->getName());
      }

      //different cases
      IOStreamColumnType columnType = StringToColumnType(vColumnDescription[i]->editType());
      switch (columnType) {
      case (IOStreamColumnType::Quantitative):
        column->set_attribute("type", "Quantitative", "xsi");
        break;
      case (IOStreamColumnType::Qualitative):
        THROW(InputException, ColumnTypeNotValid);
        break;
      case (IOStreamColumnType::Individual):
        THROW(InputException, ColumnTypeNotValid);
        break;
      case (IOStreamColumnType::Weight):
        THROW(InputException, ColumnTypeNotValid);
        break;
      }
	}

}

} //end namespace
