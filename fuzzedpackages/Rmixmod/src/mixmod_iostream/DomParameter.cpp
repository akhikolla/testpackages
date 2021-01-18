/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMDomParameter.cpp  description
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
#include "mixmod_iostream/DomParameter.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod_iostream/IOStreamUtil.h"
#include "mixmod/Kernel/IO/ParameterDescription.h"

namespace XEM {

  //constructor by default
  DomParameter::DomParameter() : xmlpp::Document() {
  }

  //destructor
  DomParameter::~DomParameter() {
  }

  DomParameter::DomParameter(string & sFilename) : xmlpp::Document() {
    set_internal_subset(sFilename, "", "");
  }

  DomParameter::DomParameter(ParameterDescription* parameterDescription, string sFilename) {
    _root = create_root_node("Parameter");
    _root->set_namespace_declaration("http://www.w3.org/2001/XMLSchema-instance", "xsi");
	//text
    xmlpp::Element* new_elt = NULL;
	//name 
	if ( !parameterDescription->getInfoName().empty() ) {
      new_elt = _root->add_child("Name");
      new_elt->add_child_text(parameterDescription->getInfoName());      
	}

	//NbVariable
    new_elt = _root->add_child("NbVariable");
    new_elt->add_child_text(std::to_string(parameterDescription->getNbVariable()));        
	//NbCluster
    new_elt = _root->add_child("NbCluster");
    new_elt->add_child_text(std::to_string(parameterDescription->getNbCluster()));       
	//Format
    new_elt = _root->add_child("Format");
    new_elt->add_child_text(FormatNumericFileToString(parameterDescription->getFormat()));        
	//Filename
	parameterDescription->saveNumericValues(getAbsolutePath(sFilename + ".txt"));
    new_elt = _root->add_child("ParameterFilename");
    new_elt->add_child_text(sFilename + ".txt");        
	//model
    new_elt = _root->add_child("Model");
    new_elt->add_child_text(ModelNameToString(parameterDescription->getModelType()->getModelName()));           
	if (isBinary(parameterDescription->getModelType()->getModelName())) {
      _root->set_attribute("type", "Qualitative", "xsi");
	}

	//TODO: HeterogeneousParameter handling is still highly experimental
	else if (isHeterogeneous(parameterDescription->getModelType()->getModelName())) {
      _root->set_attribute("type", "Composite", "xsi");
	}

	else {
      _root->set_attribute("type", "Quantitative", "xsi");
	}

    string file = getAbsolutePath(sFilename + ".mxp");
    removeIfExists(file);
    write_to_file(file);    
  }
  static void domParameterImpl(Input* cInput, Parameter * parameter, xmlpp::Element *root, string sFilename){
    root->set_namespace_declaration("http://www.w3.org/2001/XMLSchema-instance", "xsi");
	//text
    xmlpp::Element* new_elt = NULL;
	//name 
	/*if ( !parameterDescription->getInfoName().empty() ) {
      new_elt = root->add_child("Name");
      new_elt->add_child_text(parameterDescription->getInfoName());      
      }*/

	//NbVariable
    new_elt = root->add_child("NbVariable");
    new_elt->add_child_text(std::to_string(cInput->getPbDimension()));        
	//NbCluster
    new_elt = root->add_child("NbCluster");
    new_elt->add_child_text(std::to_string(cInput->getNbCluster(0)));       
	//Format
    new_elt = root->add_child("Format");
    new_elt->add_child_text(FormatNumericFileToString(parameter->getFormat()));        
	//Filename
	//parameterDescription->saveNumericValues(sFilename + ".txt");
    new_elt = root->add_child("ParameterFilename");
    std::string parFilename = parameter->getFilename();
    if(parFilename != ""){
      new_elt->add_child_text(normalizeFilename(parFilename));
    } else {
      std::string txtFile = sFilename + ".txt";
      string absPath = getAbsolutePath(txtFile);
      std::ofstream fo(absPath.c_str(), ios::out);
      parameter->edit(fo);
      new_elt->add_child_text(txtFile);
    }
	//model
    new_elt = root->add_child("Model");
    new_elt->add_child_text(ModelNameToString(parameter->getModelType()->getModelName()));           
	if (isBinary(parameter->getModelType()->getModelName())) {
      root->set_attribute("type", "Qualitative", "xsi");
	}

	//TODO: HeterogeneousParameter handling is still highly experimental
	else if (isHeterogeneous(parameter->getModelType()->getModelName())) {
      root->set_attribute("type", "Composite", "xsi");
	}

	else {
      root->set_attribute("type", "Quantitative", "xsi");
	}

  }
  
  DomParameter::DomParameter(xmlpp::Element *root) {
    _root = create_root_node_by_import(root);
  }

  DomParameter::DomParameter(ClusteringInput* cInput, string sFilename) {
    //ParameterDescription* parameterDescription = NULL;
    Parameter * parameter = cInput->getStrategy()->getStrategyInit()->getInitParameter(0);
    _root = create_root_node("Parameter");
    domParameterImpl(cInput, parameter, _root, sFilename);    
    string file = getAbsolutePath(sFilename + ".mxp");
    removeIfExists(file);
    write_to_file(file);    
  }
  DomParameter::DomParameter(PredictInput* pInput, string sFilename) {
    //ParameterDescription* parameterDescription = NULL;
    Parameter * parameter = pInput->getClassificationRule();
    _root = create_root_node("Parameter");
    domParameterImpl(pInput, parameter, _root, sFilename);    
    string file = getAbsolutePath(sFilename + ".mxp");
    removeIfExists(file);
    write_to_file(file);    
  }

  
  ParameterDescription * DomParameter::readParameter(string sFilename) {

	//-------
	//load file in this
	//-------
    xmlpp::DomParser parser;
    parser.parse_file(sFilename);
    xmlpp::Document *doc = parser.get_document();    
	_root = doc->get_root_node();//documentElement();

	if ( _root->get_name() != "Parameter" ) return 0;

    //------------------------
    //Declaration of variables
    //------------------------
    xmlpp::Element *elementName, *elementNbVariable, *elementNbCluster, *elementFormat, 
      *elementParameterFilename, *elementModel, *elementListNbFactor ;

    //if qualitative case, listNbFactor exists      
    if( _root->get_attribute_value("type", "xsi") == "Qualitative"){
      elementListNbFactor = dynamic_cast<xmlpp::Element*>(_root->get_first_child("ListNbFactor"));
    }
    
    //name
    elementName = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Name"));
    string sName = "";
    if (elementName) {
      sName = elementName->get_child_text()->get_content() ; 
    }
    
    //nbVariable
    elementNbVariable = dynamic_cast<xmlpp::Element*>(_root->get_first_child("NbVariable"));
    int64_t nbVariable = std::stoll(elementNbVariable->get_child_text()->get_content());
    
    //nbCluster
    elementNbCluster = dynamic_cast<xmlpp::Element*>(_root->get_first_child("NbCluster"));
    int64_t nbCluster = std::stoll(elementNbCluster->get_child_text()->get_content());
    
    //Format
    elementFormat = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Format"));
    FormatNumeric::FormatNumericFile format = 
      StringToFormatNumericFile(elementFormat->get_child_text()->get_content());
    
    //Parameter Filename
    elementParameterFilename = dynamic_cast<xmlpp::Element*>(_root->get_first_child("ParameterFilename"));
    string parameterFilename = elementParameterFilename->get_child_text()->get_content();
    
    //Model
    elementModel = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Model"));
    ModelName modelName = StringToModelName(elementModel->get_child_text()->get_content());
    
    //listNbFactor (if Qualitative case)
    if (isBinary(modelName)) {
      modelName = Binary_pk_Ekjh; // because ParameterDescription creates a BinaryEkjhParameter anyway
      return new ParameterDescription(nbCluster, nbVariable, Global::vNbFactor,
                                      format, parameterFilename, sName, modelName);
    }
    else if (isHeterogeneous(modelName)) {
      return new ParameterDescription(nbCluster, Global::nbVariables_binary,
                                      Global::nbVariables_gaussian, Global::vNbFactor, format,
                                      parameterFilename, sName, modelName);
    }
    else { //gaussian
      //if(!isGeneral(modelName)) modelName = Gaussian_pk_Lk_Ck; // because ParameterDescription creates a GaussianGeneralParameter object anyway
      return new ParameterDescription(nbCluster, nbVariable, format,
                                      parameterFilename, sName, modelName);
    }
    //}//


  }

} //end namespace
