/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMNodeClusteringOutput.cpp  description
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

#include "mixmod_iostream/NodeClusteringOutput.h"
#include "mixmod_iostream/DomParameter.h"
#include "mixmod_iostream/DomLabel.h"
#include "mixmod_iostream/DomProba.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/Criterion/CriterionOutput.h"
#include "mixmod/Clustering/ClusteringModelOutput.h"

namespace XEM {

  NodeClusteringOutput::NodeClusteringOutput() : NodeOutput() {
  }

  NodeClusteringOutput::~NodeClusteringOutput() {
  }

  NodeClusteringOutput::NodeClusteringOutput(ClusteringOutput * output, string & s) : NodeOutput() {
	for (int64_t i = 0; i < output->getNbClusteringModelOutput(); ++i) {
      writeOutput(output->getClusteringModelOutput(i), output->getCriterionName(), s, i + 1);
	}
  }

  NodeClusteringOutput::NodeClusteringOutput( xmlpp::Element *rootOutput ) : NodeOutput(rootOutput) {
  }

  //read the clustering to fill XEMOutput
  ClusteringModelOutput * NodeClusteringOutput::readClustering() {

	//model
    xmlpp::Element *elementModel = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("Model"));
	ModelName modelName = StringToModelName(elementModel->get_child_text()->get_content());
	ModelType modelType(modelName);

	//nbCluster
    xmlpp::Element *elementNbCluster = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("NbCluster"));
    int64_t nbCluster = std::stoll(elementNbCluster->get_child_text()->get_content());

	//error
    xmlpp::Element *elementError = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("Error"));
	if (elementError) {
      
      Exception error = Exception(elementError->get_child_text()->get_content());

      return new ClusteringModelOutput(modelType, nbCluster, error);
	}
	else {
      //Criterion
      xmlpp::Element *elementListCriterion =  dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("ListOutputCriterion"));
      auto children = elementListCriterion->get_children();
      vector<CriterionOutput*> vCriterionOutput;
      for (auto it=children.begin(); it != children.end(); ++it){
        xmlpp::Element * elementOutputCriterion = dynamic_cast<xmlpp::Element*>(*it);
        if(!elementOutputCriterion) continue;
        
        //criterion Name
        xmlpp::Element *cnElt = dynamic_cast<xmlpp::Element*>(elementOutputCriterion->get_first_child("CriterionName"));
        CriterionName criterionName = StringtoCriterionName(cnElt->get_child_text()->get_content());
        double criterionValue = 0;
        Exception * error = NOERROR.clone();
        xmlpp::Element *elementError = dynamic_cast<xmlpp::Element*>(elementOutputCriterion->get_first_child("Error"));
        if(elementError){
          //error
          delete error;
          error = new Exception(elementError->get_child_text()->get_content());
        }
        else {
          //criterion Value
          xmlpp::Element* cvElt = dynamic_cast<xmlpp::Element*>(elementOutputCriterion->get_first_child("CriterionValue"));
          if (IOMODE == IoMode::BINARY) {
            stringstream stream;
            uint64_t tmp;
            stream << hex << cvElt->get_child_text()->get_content(); 
            stream >> tmp;
            memcpy(&criterionValue, &tmp, sizeof(tmp));
          }
          else {
            criterionValue = std::stod(cvElt->get_child_text()->get_content()); 
          }
          //.toElement().text().toDouble();
        }

        vCriterionOutput.push_back(new CriterionOutput(criterionName, criterionValue, *error));
        delete error;
      }

      //likelihood
      xmlpp::Element *elementLikelihood = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("Likelihood"));
      double likelihood;
      if (IOMODE == IoMode::BINARY) {
        stringstream stream;
        uint64_t tmp;
        stream << hex << elementLikelihood->get_child_text()->get_content(); 
        stream >> tmp;
        memcpy(&likelihood, &tmp, sizeof(tmp));
      }
      else {
        likelihood = std::stod(elementLikelihood->get_child_text()->get_content());
      }
      //Parameter
      xmlpp::Element *elementParameter = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("Parameter"));

      string parameterFilename = elementParameter->get_child_text()->get_content();
      //TODO : take the XEMParameterDescription


      //Label
      xmlpp::Element *elementLabel = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("Label"));
      string labelFilename = elementLabel->get_child_text()->get_content();
      //TODO take labeldescription


      //Proba
      xmlpp::Element *elementProba = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("Proba"));
      string probaFilename = elementProba->get_child_text()->get_content();
      //TODO tke probadescription

      //ModelOutput(ModelType & modelType, int64_t nbCluster, 
      //vector<CriterionOutput> & criterionOutput, double likelihood, 
      //ParameterDescription & parameterDescription, LabelDescription & labelDescription,  
      //ProbaDescription & probaDescription);

      return new ClusteringModelOutput(modelType, nbCluster, vCriterionOutput, 
                                       likelihood, *readParameter(parameterFilename), 
                                       *readLabel(labelFilename), *readProba(probaFilename));
	}
  }

  ParameterDescription * NodeClusteringOutput::readParameter(string sFilename) {
	ValidateSchema(sFilename, IOStreamXMLFile::Parameter);
    //throw IOStreamErrorType::badLoadXML;

	DomParameter parameter;

	return parameter.readParameter(sFilename);
  }

  unique_ptr<LabelDescription>  NodeClusteringOutput::readLabel(string sFilename) {
	ValidateSchema(sFilename, IOStreamXMLFile::Data);
    //throw IOStreamErrorType::badLoadXML;

	DomLabel label;

	return label.readLabel(sFilename);
  }

  unique_ptr<ProbaDescription> NodeClusteringOutput::readProba(string sFilename) {
	ValidateSchema(sFilename, IOStreamXMLFile::Data);
    //throw IOStreamErrorType::badLoadXML;

	DomProba proba;

	return proba.readProba(sFilename);
  }

  void NodeClusteringOutput::writeOutput(ClusteringModelOutput* output, 
                                         const std::vector<CriterionName> & criterionName, string str, int64_t numOutput)
  {
	//strategy
    xmlpp::Element *outputElement = _rootOutput->add_child("Output");

	//model
    xmlpp::Element *modelElement = outputElement->add_child("Model");
    modelElement->add_child_text(ModelNameToString(output->getModelType().getModelName()));
	//nbCluster
    xmlpp::Element *nbClusterElement =  outputElement->add_child("NbCluster");
    nbClusterElement->add_child_text(std::to_string(output->getNbCluster()));
	if (output->getStrategyRunError() == NOERROR) {
      //criterion
      xmlpp::Element *listCriterionElement = outputElement->add_child("ListOutputCriterion");
      for (int64_t i = 0; i < criterionName.size() ; ++i) {
        xmlpp::Element *criterionElement = listCriterionElement->add_child("OutputCriterion");
        xmlpp::Element *nameElement = criterionElement->add_child("CriterionName");
        nameElement->add_child_text(CriterionNameToString(criterionName[i]));

        if (output->getCriterionOutput(criterionName[i]).getError() == NOERROR) {
          xmlpp::Element *valueElement = criterionElement->add_child("CriterionValue");                  
          if (IOMODE == IoMode::BINARY) {
            stringstream stream;
            int64_t tmp;
            double tmp_value = output->getCriterionOutput(criterionName[i]).getValue();
            memcpy(&tmp, &tmp_value, sizeof(tmp_value));
            stream << hex << tmp;
            string s = stream.str();
            valueElement->add_child_text(s);
          }
          else {
            valueElement->add_child_text(std::to_string(output->getCriterionOutput(criterionName[i]).getValue()));
          }
        }
        else {
          xmlpp::Element *errorElement = criterionElement->add_child("Error");
          errorElement->add_child_text(output->getStrategyRunError().what());
        }
      }

      //Likelihood
      xmlpp::Element *likelihoodElement = outputElement->add_child("Likelihood");
      if (IOMODE == IoMode::BINARY) {
        stringstream stream;
        int64_t tmp;
        double tmp_value = output->getLikelihood();
        memcpy(&tmp, &tmp_value, sizeof(tmp_value));
        stream << hex << tmp;
        string s = stream.str();
        likelihoodElement->add_child_text(s);
      }
      else {
        likelihoodElement->add_child_text(std::to_string(output->getLikelihood()));
      }
      //Parameter
      xmlpp::Element *parameterElement = outputElement->add_child("Parameter");
      string parameterFilename = str+"Param"+std::to_string(numOutput);
      parameterElement->add_child_text(parameterFilename+".mxp");

      //create parameter file (.mxp)
      DomParameter docParameter(
                                output->getParameterDescription(), parameterFilename); //.toStdString());

      //label
      xmlpp::Element *labelElement = outputElement->add_child("Label");
      string labelFilename = str + "Label" + std::to_string(numOutput);
      labelElement->add_child_text(labelFilename + ".mxd");
      //create label file (.mxd)
      DomLabel docLabel(output->getLabelDescription(), labelFilename);

      //proba
      xmlpp::Element *probaElement = outputElement->add_child("Proba");
      string probaFilename = str + "Proba" + std::to_string(numOutput);
      probaElement->add_child_text(probaFilename + ".mxd");
      //create proba file (.mxd)
      DomProba docProba(output->getProbaDescription(), probaFilename);
	}
	else {
      //error
      xmlpp::Element *errorElement = outputElement->add_child("Error");
      errorElement->add_child_text(output->getStrategyRunError().what());
	}
  }

} //end namespace
