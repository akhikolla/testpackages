/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMDomData.cpp  description
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

//#include <QTextStream>
#include "mixmod_iostream/DomData.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/QualitativeColumnDescription.h"
#include "mixmod/Kernel/IO/QuantitativeColumnDescription.h"
#include "mixmod/Kernel/IO/ColumnDescription.h"
#include "mixmod/Kernel/IO/IndividualColumnDescription.h"
#include "mixmod/Kernel/IO/WeightColumnDescription.h"

namespace XEM {

  DomData::DomData() : xmlpp::Document() {
  }

  DomData::~DomData() {
  }

  DomData::DomData(string & sFilename) : xmlpp::Document() {
    set_internal_subset(sFilename, "", "");
  }
 
  DomData::DomData(xmlpp::Element * root) {
	//_root = root;
    _root = create_root_node_by_import(root);
  }

  DomData::DomData(const DataDescription & dataFile, string & sFilename) {
    _root = create_root_node("Data");
    _root->set_namespace_declaration("http://www.w3.org/2001/XMLSchema-instance", "xsi");


	//QDomText text;
    xmlpp::Element* new_elt = NULL;
	//Name
	if ( !dataFile.getInfoName().empty() ) {
      new_elt = _root->add_child("Name");
      //new_elt->add_child_text(dataFile.getInfoName());
      new_elt->add_child_text(normalizeFilename(dataFile.getInfoName()));
      //new_elt->add_child_text(dataFile.getInfoName());            
        
	}

	//NbSample
    new_elt = _root->add_child("NbSample");
    new_elt->add_child_text(std::to_string(dataFile.getNbSample()));    
	//NbColumn
	//   QDomElement nbCol = createElement("NbColumn");  
	//   text = createTextNode(QString::number(dataFile.getNbColumn()));
	//   nbCol.appendChild(text);        
	//   _root.appendChild(nbCol); 

	//format
    new_elt = _root->add_child("Format");
    new_elt->add_child_text(FormatNumericFileToString(dataFile.getFormat()));        
	//datafilename
    new_elt = _root->add_child("DataFilename");
    new_elt->add_child_text(dataFile.getFileName());            
	writeListColumnNode(dataFile.getAllColumnDescription());

	//appendChild( _root );

	//write new file .mxd to describe data
    Glib::ustring filename = getAbsolutePath(sFilename + ".mxd");
    removeIfExists(filename);
    write_to_file(filename);
}

//with the filename (string), DomData writes the data 
//from the .mxd file to a struture Data of NvInput
DataDescription * DomData::readDataFile(const string & sFilename)
{
	// Data * res = NULL;//return the new Data for NvInput from .mxd file
  xmlpp::DomParser parser;
  parser.parse_file(sFilename);
  xmlpp::Document *doc = parser.get_document();
  //-------
  //Variable declaration
  //-------
  string iName = "";
  int64_t iNbSample;
  int64_t iNbColumn = 0;
  FormatNumeric::FormatNumericFile iFormat;
  //  int64_t * iTabModality;
  string iDataFilename = "";

  xmlpp::Element *listColumnNode = NULL;
  //initialization of (optional) ColumnDescription vector
  vector<ColumnDescription *> vColumnDescription;
  xmlpp::Element *root = doc->get_root_node();
  if ( root->get_name() != "Data" )
    // Since the file has been validated at this point, we should never return NULL
    return NULL;

  //bool if qualitative or not
  listColumnNode = dynamic_cast<xmlpp::Element*>(root->get_first_child("ListColumn"));

	//--------
	//assigment of main node value in main node "data" to declared variables
	//--------
	// TODO: since elements have to appear sequentially in a specific order, why do we use a loop ?
    xmlpp::Element *n = dynamic_cast<xmlpp::Element*>(root->get_first_child("Name"));
    if(n){
      iName = n->get_child_text()->get_content();
    }
    n = dynamic_cast<xmlpp::Element*>(root->get_first_child("NbSample"));
    if(n){
      iNbSample = std::stoll(n->get_child_text()->get_content());
    }
    
    n = dynamic_cast<xmlpp::Element*>(root->get_first_child("Format"));
    if(n){
      iFormat = StringToFormatNumericFile(n->get_child_text()->get_content());
    }    
    n = dynamic_cast<xmlpp::Element*>(root->get_first_child("DataFilename"));    
    if(n){ //TODO: manage RELATIVE_PATHS
      iDataFilename = n->get_child_text()->get_content();
    }  
	//--------
	//filling of Column in dataDescription (optional)
	//--------
	//if (!listColumnNode.isNull()) {
    if(listColumnNode){

      //n = listColumnNode.firstChildElement(); //1st ColumnDescription

      //cross each column
      auto children = listColumnNode->get_children();
      for (auto it=children.begin(); it != children.end(); ++it){
        xmlpp::Element* n = dynamic_cast<xmlpp::Element*>(*it);
        if(!n) continue;
        //Column Type
        IOStreamColumnType columnType = 
          StringToColumnType(n->get_attribute_value("type", "xsi"));

        //index of Column [0..nbColumn-1]
        int64_t index = std::stoll(n->get_attribute_value("Num")) - 1;
        
        string nameColumn = n->get_attribute_value("Name");

        switch (columnType) {
        case IOStreamColumnType::Qualitative:
          {
            //nbFactor
            int64_t nbFactor = std::stoll(n->get_attribute_value("NbFactor"));

            //creation of Column
            auto column = new QualitativeColumnDescription(index, nbFactor);
            column->setName(nameColumn);

            //assigment of column to vColumnDescription
            vColumnDescription.push_back(column);
            break;
          }
        case IOStreamColumnType::Quantitative:
          {
            //creation and assigment of column to dataDescription(index)
            ColumnDescription * column = new QuantitativeColumnDescription(index);
            column->setName(nameColumn);

            //assigment of column to vColumnDescription
            vColumnDescription.push_back(column);
            break;
          }
        case IOStreamColumnType::Individual:
          {
            //creation of Column
            IndividualColumnDescription* column = new IndividualColumnDescription(index);
            column->setName(nameColumn);
            //assigment of column to vColumnDescription
            vColumnDescription.push_back(column);
            xmlpp::Element* nn = dynamic_cast<xmlpp::Element*>(n->get_first_child("IndividualName"));
            if(!nn) break;
            auto individuals = nn->get_children("IndividualName");
            for (auto it2=individuals.begin(); it2 != individuals.end(); ++it2){
              xmlpp::Element* individualNode = dynamic_cast<xmlpp::Element*>(*it2);
              if(!individualNode) continue;

              //cross Column to find all individualdescription
              //creation of the individual description and filling
              IndividualDescription individual;
              individual.name = individualNode->get_child_text()->get_content();
              individual.num = std::stoll(individualNode->get_attribute_value("Num"));
              column->insertIndividualDescription(individual, index);

              //next individual
            }

            break;
          }
			case IOStreamColumnType::Weight:

              //creation and assigment of column to dataDescription(index)
              ColumnDescription * column = new WeightColumnDescription(index);
              column->setName(nameColumn);

              //assigment of column to vColumnDescription
              vColumnDescription.push_back(column);
              break;
        }

        //compute the number of column
        iNbColumn++;
        
        //next Column in ListColumn 
      }
	}

	//initialization of DataDescription
	return new DataDescription(iNbSample, iNbColumn, 
                               vColumnDescription, iFormat, iDataFilename, iName);
}

void DomData::writeListColumnNode(const vector<ColumnDescription *> & vColumnDescription) {

  //ListColumn
  int64_t nbColumn = 0;
  nbColumn = vColumnDescription.size();
  if(!nbColumn) return;
  xmlpp::Element* listColumn = _root->add_child("ListColumn");
                    
	for (int64_t i = 0; i < nbColumn; ++i) {
      xmlpp::Element *column = listColumn->add_child("Column");

		//position of variable 
      column->set_attribute("Num", std::to_string(vColumnDescription[i]->getIndex() + 1));

      //name
      if (!vColumnDescription[i]->getName().empty()) {
        //xmlpp::Element *name = column->add_child("Name");
        //name->add_child_text(vColumnDescription[i]->getName());
        column->set_attribute("Name", vColumnDescription[i]->getName());
      }

      //different cases
      IOStreamColumnType columnType = StringToColumnType(vColumnDescription[i]->editType());
      switch (columnType) {
      case (IOStreamColumnType::Quantitative):
        column->set_attribute("type", "Quantitative", "xsi");
        break;
      case (IOStreamColumnType::Qualitative):
		{
          column->set_attribute("type", "Qualitative", "xsi");

          QualitativeColumnDescription * QCD = 
            dynamic_cast<QualitativeColumnDescription *> (vColumnDescription[i]);
			
          //nbFactor
          column->set_attribute("NbFactor", std::to_string(QCD->getNbFactor()));

          //ListFactor [TODO: <Factor Name=... Num=.../> ...]
          xmlpp::Element *listFactorName = NULL; 
          for (int64_t j = 0 ; j < QCD->getVariableDescription().size(); ++j) {
            //factor j of variable i
            if (!QCD->getVariableDescription()[j].name.empty()) {
              if(!listFactorName) listFactorName = column->add_child("ListFactorName");
              xmlpp::Element *factorName = listFactorName->add_child("FactorName");
              //index of factor
              factorName->set_attribute("Num", 
                                        std::to_string(QCD->getVariableDescription()[j].num));
              factorName->add_child_text(QCD->getVariableDescription()[j].name);
            }
          }

          break;
		}
      case (IOStreamColumnType::Individual):
		{
          column->set_attribute("type", "Individual", "xsi");
          IndividualColumnDescription * ICD = 
            dynamic_cast<IndividualColumnDescription *> (vColumnDescription[i]);
          for (int64_t j = 0; j < ICD->getIndividualDescription().size(); ++j) {

            //creates individual nodes in .mxd
            xmlpp::Element *individualName = column->add_child("IndividualName");
				//index of individual
            individualName->set_attribute("Num", 
                                          std::to_string(ICD->getIndividualDescription()[j].num));
            individualName->add_child_text(ICD->getIndividualDescription()[j].name);

			}
			break;
		}
      case (IOStreamColumnType::Weight):
        column->set_attribute("type", "Weight", "xsi");
        break;
      }
	}

}

} //end namespace
