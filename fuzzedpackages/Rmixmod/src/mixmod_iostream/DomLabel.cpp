/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMDomLabel.cpp  description
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
#include "mixmod/Kernel/IO/Partition.h"
#include "mixmod_iostream/DomLabel.h"
#include "mixmod_iostream/IOStreamUtil.h"
#include "mixmod/Kernel/IO/QualitativeColumnDescription.h"
#include "mixmod/Kernel/IO/IndividualColumnDescription.h"

namespace XEM {

//constructor by default
  DomLabel::DomLabel() : xmlpp::Document() {
  }

//destructor
  DomLabel::~DomLabel() {
  }

  DomLabel::DomLabel(string & sFilename) : xmlpp::Document() {
    set_internal_subset(sFilename, "", "");
  }


  DomLabel::DomLabel(LabelDescription* labelDescription, string str) : xmlpp::Document() {
    _root = create_root_node("Label");
    _root->set_namespace_declaration("http://www.w3.org/2001/XMLSchema-instance", "xsi");
	//QDomText text;
    xmlpp::Element* new_elt = NULL;
	//Name
	if ( !labelDescription->getInfoName().empty() ) {
      new_elt = _root->add_child("Name");
      new_elt->add_child_text(labelDescription->getInfoName());      
	}

	//NbSample
    new_elt = _root->add_child("NbSample");
    new_elt->add_child_text(std::to_string(labelDescription->getNbSample()));        
    new_elt = _root->add_child("NbCluster");
    new_elt->add_child_text(std::to_string(labelDescription->getNbCluster()));        
    new_elt = _root->add_child("Format");
    new_elt->add_child_text(FormatNumericFileToString(labelDescription->getFormat()));           
    //new_elt = _root->add_child("Type");
    //new_elt->add_child_text("label");        
	//datafilename
	labelDescription->saveNumericValues(getAbsolutePath(str + ".txt"));
    new_elt = _root->add_child("Filename");
    new_elt->add_child_text(str + ".txt");               
	//writeListColumnNode(labelDescription->getAllColumnDescription());
	//write new file .mxd to describe data
    Glib::ustring filename = getAbsolutePath(str + ".mxl");
    removeIfExists(filename);
    write_to_file(filename);    
  }
  /*
  //Old style constructor, to be deleted ASAP
  //NB:  Labels as Data representation is not suitable because zeroes are not allowed as
  //qualitative data values but zeroes are valid label values in semi-supervized classification
  //TO BE REMOVED ASAP
  DomLabel::DomLabel(LabelDescription* labelDescription, string str) : xmlpp::Document() {
    _root = create_root_node("Data");
    _root->set_namespace_declaration("http://www.w3.org/2001/XMLSchema-instance", "xsi");
	//QDomText text;
    xmlpp::Element* new_elt = NULL;
	//Name
	if ( !labelDescription->getInfoName().empty() ) {
      new_elt = _root->add_child("Name");
      new_elt->add_child_text(labelDescription->getInfoName());      
	}

	//NbSample
    new_elt = _root->add_child("NbSample");
    new_elt->add_child_text(std::to_string(labelDescription->getNbSample()));        
	//NbColumn
	//   QDomElement nbCol = createElement("NbColumn");  
	//   text = createTextNode(QString::number(labelDescription->getNbColumn()));
	//   nbCol.appendChild(text);        
	//   _root.appendChild(nbCol); 

	//format
    new_elt = _root->add_child("Format");
    new_elt->add_child_text(FormatNumericFileToString(labelDescription->getFormat()));           
	//datafilename
	labelDescription->saveNumericValues(str + ".txt");
    new_elt = _root->add_child("DataFilename");
    new_elt->add_child_text(str + ".txt");               
	writeListColumnNode(labelDescription->getAllColumnDescription());

	//write new file .mxd to describe data
    Glib::ustring filename = str + ".mxl";
    removeIfExists(filename);
    write_to_file(filename);    
  }
  */
  DomLabel::DomLabel(Partition * partition, string & sFilename) {
    //_root = createElement( "Partition" );
    string tag = partition->getPartitionFile()._type==TypePartition::label ? "Label" : "Partition";
    string extn = partition->getPartitionFile()._type==TypePartition::label ? ".mxl" : ".mxd";
    _root = create_root_node(tag);
    xmlpp::Element* new_elt = NULL;
	//text

	//name TODO

	//NbSample
    new_elt = _root->add_child("NbSample");
    new_elt->add_child_text(std::to_string(partition->getNbSample()));
	//nbCluster
    new_elt = _root->add_child("NbCluster");
    new_elt->add_child_text(std::to_string(partition->getNbCluster()));    
	//Format
    new_elt = _root->add_child("Format");
    new_elt->add_child_text(FormatNumericFileToString(partition->getPartitionFile()._format));            
	//type
    //new_elt = _root->add_child("Type");
    //new_elt->add_child_text(TypePartitionToString(partition->getPartitionFile()._type));            
	//Filename
    new_elt = _root->add_child("Filename");
    new_elt->add_child_text(partition->getPartitionFile()._fileName);            

	//appendChild(_root);

	//write in new file .mxl to describe the partition
    Glib::ustring filename = sFilename + extn;
    removeIfExists(filename);
    write_to_file(filename);
    
  }

  DomLabel::DomLabel(xmlpp::Element *root) {
    _root = create_root_node_by_import(root);
  }
  
  unique_ptr<LabelDescription>  DomLabel::readLabelAsData(string sFilename) {
  //-------
  //load file in variable "doc"
  //-------
    xmlpp::DomParser parser;
    parser.parse_file(sFilename);
    xmlpp::Document *doc = parser.get_document();
    xmlpp::Element *_root = doc->get_root_node();
  
    if(_root->get_name() != "Data") throw IOStreamErrorType::badElementInDataXML;
    
    //------------------------
    //Declaration of variables
    //------------------------
    xmlpp::Element *elementName, *elementNbSample, *elementNbColumn, *elementFormat, 
				*elementLabelFilename, *elementListColumn ;

    //name
    elementName = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Name"));
    string sName = elementName ? elementName->get_child_text()->get_content() : "";
    //nbSample
    elementNbSample = dynamic_cast<xmlpp::Element*>(_root->get_first_child("NbSample"));
    int64_t nbSample = std::stoll(elementNbSample->get_child_text()->get_content());
    //nbColumn
    //     elementNbColumn = _root.namedItem("NbColumn").toElement();
    //     int64_t nbColumn = elementNbColumn.text().toLongLong();
  
    //Format
    elementFormat = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Format"));
    FormatNumeric::FormatNumericFile format =
      StringToFormatNumericFile(elementFormat->get_child_text()->get_content());
    //Parameter Filename
    elementLabelFilename = dynamic_cast<xmlpp::Element*>(_root->get_first_child("DataFilename"));
    string dataFilename = elementLabelFilename->get_child_text()->get_content();
    //ListColumn
    //elementListColumn = _root.namedItem("ListColumn").toElement();
    elementListColumn = dynamic_cast<xmlpp::Element*>(_root->get_first_child("ListColumn"));
    
    //initialization of ColumnDescription vector
    vector<ColumnDescription *> vColumnDescription;

    //NbColumn
    int64_t iNbColumn = 0;

    //cross each column
    if(elementListColumn){
      auto children = elementListColumn->get_children();
      for (auto it=children.begin(); it != children.end(); ++it){
        xmlpp::Element* n = dynamic_cast<xmlpp::Element*>(*it);
        if(!n) continue;
        //Column Type
        IOStreamColumnType columnType = 
          StringToColumnType(n->get_attribute_value("type","xsi"));

        //index of Column [0..nbColumn-1]
        //int64_t index = n.attribute("Num").toLongLong() - 1;
        int64_t index = std::stoll(n->get_attribute_value("Num")) - 1;
        string nameColumn = n->get_attribute_value("Name");
        switch (columnType) {
        case IOStreamColumnType::Qualitative:
          {
            //nbFactor
            int64_t nbFactor = std::stoll(n->get_attribute_value("NbFactor"));
            //creation of Column
            QualitativeColumnDescription * column = 
              new QualitativeColumnDescription(index, nbFactor);
            column->setName(nameColumn);

            //assigment of column to vColumnDescription
            vColumnDescription.push_back(column);
            break;
          }
        case IOStreamColumnType::Quantitative:
          {
            throw IOStreamErrorType::badElementInDataXML;
          }
        case IOStreamColumnType::Individual:
          {
            //creation of Column
            IndividualColumnDescription * column = new IndividualColumnDescription(index);
            column->setName(nameColumn);
            auto individuals = n->get_children("IndividualName");
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

            //assigment of column to vColumnDescription
            vColumnDescription.push_back(column);
            break;
          }
        case IOStreamColumnType::Weight:
          throw IOStreamErrorType::badElementInDataXML;
        }

        //compute the number of column
        iNbColumn++;

        //next Column in ListColumn 
      }//for...
    }//if...
    unique_ptr<LabelDescription> newLabelDescription(new LabelDescription(nbSample, iNbColumn, vColumnDescription, format, dataFilename, sName));
    return newLabelDescription;

//  return new LabelDescription(nbSample, iNbColumn, vColumnDescription, 
//                            format, dataFilename, sName);
  }

  
  
  unique_ptr<LabelDescription>  DomLabel::readLabel(string sFilename) {
  //-------
  //load file in variable "doc"
  //-------
    xmlpp::DomParser parser;
    parser.parse_file(sFilename);
    xmlpp::Document *doc = parser.get_document();
    xmlpp::Element *_root = doc->get_root_node();
  
    if(_root->get_name() != "Label"/*&&_root->get_name() != "Partition"*/) throw IOStreamErrorType::badElementInDataXML;
    
    //------------------------
    //Declaration of variables
    //------------------------
    xmlpp::Element *elementName, *elementNbSample, *elementNbCluster, *elementFormat, 
				*elementLabelFilename, *elementListColumn ;

    //name
    elementName = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Name"));
    string sName = elementName ? elementName->get_child_text()->get_content() : "";
    //nbSample
    elementNbSample = dynamic_cast<xmlpp::Element*>(_root->get_first_child("NbSample"));
    int64_t nbSample = std::stoll(elementNbSample->get_child_text()->get_content());
    //nbSample
    elementNbCluster = dynamic_cast<xmlpp::Element*>(_root->get_first_child("NbCluster"));
    int64_t nbCluster = std::stoll(elementNbCluster->get_child_text()->get_content());
  
    //Format
    elementFormat = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Format"));
    FormatNumeric::FormatNumericFile format =
      StringToFormatNumericFile(elementFormat->get_child_text()->get_content());
    //Parameter Filename
    elementLabelFilename = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Filename"));
    string dataFilename = elementLabelFilename->get_child_text()->get_content();
    //ListColumn
    //elementListColumn = _root.namedItem("ListColumn").toElement();
    //elementListColumn = dynamic_cast<xmlpp::Element*>(_root->get_first_child("ListColumn"));
    
    //initialization of ColumnDescription vector
    vector<ColumnDescription *> vColumnDescription(1);
    QualitativeColumnDescription * column = new QualitativeColumnDescription(0, nbCluster);
    string lab = "Label";
    column->setName(lab);
    vColumnDescription[0] = column;
    unique_ptr<LabelDescription> newLabelDescription(new LabelDescription(nbSample, 1, vColumnDescription, format, dataFilename, sName));
    return newLabelDescription;

//  return new LabelDescription(nbSample, iNbColumn, vColumnDescription, 
//                            format, dataFilename, sName);
  }
  
  void DomLabel::writeListColumnNode(const vector<ColumnDescription *> & vColumnDescription) {

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
        THROW(InputException, ColumnTypeNotValid);
        break;
      case (IOStreamColumnType::Qualitative):
		{
          column->set_attribute("type", "Qualitative", "xsi");
          //nbFactor
          QualitativeColumnDescription* QCD =
            dynamic_cast<QualitativeColumnDescription *> (vColumnDescription[i]);
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
        THROW(InputException, ColumnTypeNotValid);
        break;
      case (IOStreamColumnType::Weight):
        THROW(InputException, ColumnTypeNotValid);
        break;
      }
      
	}

  }

} //end namespace
