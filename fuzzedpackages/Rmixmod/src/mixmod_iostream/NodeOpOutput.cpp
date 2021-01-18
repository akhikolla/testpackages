/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMNodeOpOutput.cpp  description
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

#include "mixmod_iostream/NodeOpOutput.h"
#include "mixmod_iostream/DomParameter.h"
#include "mixmod_iostream/DomLabel.h"
#include "mixmod_iostream/DomProba.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/Criterion/CriterionOutput.h"
#include "mixmod/Clustering/ClusteringModelOutput.h"
#include "mixmod/DiscriminantAnalysis/Learn/LearnModelOutput.h"
#include "mixmod/DiscriminantAnalysis/Predict/PredictModelOutput.h"
#include "mixmod/Kernel/IO/Label.h"
#include <stdio.h>
namespace XEM {

  NodeOpOutput::NodeOpOutput() : NodeOutput() {
  }

  NodeOpOutput::~NodeOpOutput() {
  }

  NodeOpOutput::NodeOpOutput(ClusteringOutput * output, string & s) : NodeOutput() {
	for (int64_t i = 0; i < output->getNbClusteringModelOutput(); ++i) {
      writeOutput(output->getClusteringModelOutput(i), output->getCriterionName(), s, i + 1);
	}
  }
  NodeOpOutput::NodeOpOutput(LearnOutput * output, const std::vector<CriterionName> & criterionName, string & s) : NodeOutput() {
	for (int64_t i = 0; i < output->getNbLearnModelOutput(); ++i) {
      writeOutput(output->getLearnModelOutput(i), criterionName, s, i + 1);
	}
  }
  NodeOpOutput::NodeOpOutput(PredictOutput * output, string & s) : NodeOutput() {
    auto vect = output->getPredictModelOutput();
	for (int64_t i = 0; i < vect.size(); ++i) {
      writePredictOutput(vect[i], s, i + 1);
	}
  }

  NodeOpOutput::NodeOpOutput( xmlpp::Element *rootOutput ) : NodeOutput(rootOutput) {
  }

  //read the clustering to fill XEMOutput
  template<class T>
  T * NodeOpOutput::read4Output(Input *inp) {
    PredictInput* pInput = nullptr;
    ParameterDescription *parDesc = nullptr;
    ModelType modelType;
    int64_t nbCluster;
    vector<CriterionOutput*> vCriterionOutput;
    double likelihood = 0.0;
    ///////
	std::vector<CriterionName>* vCriterion = new std::vector<CriterionName>();
	const int nbCriterion = inp->getCriterionName().size();
	for (int iCriterion = 0; iCriterion<nbCriterion; iCriterion++)
		vCriterion->push_back(inp->getCriterionName()[iCriterion]);
    ///////
	//model
    xmlpp::Element *elementModel = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("Model"));
    if(elementModel){
      ModelName modelName = StringToModelName(elementModel->get_child_text()->get_content());
      modelType = ModelType(modelName);
    } else{
        pInput = dynamic_cast<PredictInput*>(inp);
        if(!pInput) throw  IOStreamErrorType::badLoadXML; //TO DO: find a better error type...
        parDesc = new ParameterDescription(pInput->getClassificationRule());
        modelType = ModelType(*parDesc->getModelType());
        nbCluster = parDesc->getNbCluster();
    }
	//nbCluster
    xmlpp::Element *elementNbCluster = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("NbCluster"));
    if(elementNbCluster){
      nbCluster = std::stoll(elementNbCluster->get_child_text()->get_content());
    }
	//error
    xmlpp::Element *elementError = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("Error"));
	if (elementError) {
      
      Exception error = Exception(elementError->get_child_text()->get_content());

      T* moutput =  new T(modelType, nbCluster, error);
      for (auto it= vCriterion->begin(); it != vCriterion->end(); it++) {
        CriterionName criterion = *it;
        moutput->setCriterionOutput(CriterionOutput(criterion, 0.0, error));
      }
      return moutput;
	}
	else {
      //Criterion
      xmlpp::Element *elementListCriterion =  dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("ListOutputCriterion"));
      if(elementListCriterion){
        auto children = elementListCriterion->get_children();
      
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
          } else {
            //criterion Value
            xmlpp::Element* cvElt = dynamic_cast<xmlpp::Element*>(elementOutputCriterion->get_first_child("CriterionValue"));
            /*if (IOMODE == IoMode::BINARY) {
              stringstream stream;
              uint64_t tmp;
              stream << hex << cvElt->get_child_text()->get_content(); 
              stream >> tmp;
              memcpy(&criterionValue, &tmp, sizeof(tmp));
              }
              else {
              criterionValue = std_stod(cvElt->get_child_text()->get_content()); 
              }*/
            criterionValue = custom_stod(cvElt->get_child_text()->get_content()); 
            //.toElement().text().toDouble();
          }

          vCriterionOutput.push_back(new CriterionOutput(criterionName, criterionValue, *error));
          delete error;
        }
      }
      //likelihood
      xmlpp::Element *elementLikelihood = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("Likelihood"));
      if(elementLikelihood){
        likelihood = custom_stod(elementLikelihood->get_child_text()->get_content());
      }
      /*if (IOMODE == IoMode::BINARY) {
        stringstream stream;
        uint64_t tmp;
        stream << hex << elementLikelihood->get_child_text()->get_content(); 
        stream >> tmp;
        memcpy(&likelihood, &tmp, sizeof(tmp));
      }
      else {
        likelihood = std_stod(elementLikelihood->get_child_text()->get_content());
        }*/
      //Parameter
      
      string parameterFilename = "";
      xmlpp::Element *elementParameter = dynamic_cast<xmlpp::Element*>(_rootOutput->get_first_child("Parameter"));
      if(elementParameter){
        parameterFilename = elementParameter->get_child_text()->get_content();
      } /*else {
        PredictInput* pInput = dynamic_cast<PredictInput*>(inp);
        if(!pInput) throw  IOStreamErrorType::badLoadXML; //TO DO: find a better error type...
        parDesc = new ParameterDescription(pInput->getClassificationRule());
        }*/
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
      // TODO: manage memory leaks on ParameterDescription
      return new T(modelType, nbCluster, vCriterionOutput, 
                   likelihood, parameterFilename==""? *parDesc : *readParameter(parameterFilename), 
                                       *readLabel(labelFilename), *readProba(probaFilename));
	}
  }
  template ClusteringModelOutput *  NodeOpOutput::read4Output<ClusteringModelOutput>(Input *inp);
  template LearnModelOutput *  NodeOpOutput::read4Output<LearnModelOutput>(Input *inp);
  template PredictModelOutput *  NodeOpOutput::read4Output<PredictModelOutput>(Input *inp);    

  ParameterDescription * NodeOpOutput::readParameter(string sFilename) {
	ValidateSchema(sFilename, IOStreamXMLFile::Parameter);
    //throw IOStreamErrorType::badLoadXML;

	DomParameter parameter;

	return parameter.readParameter(sFilename);
  }

  unique_ptr<LabelDescription>  NodeOpOutput::readLabel(string sFilename) {
    DomLabel label;
    try{
      ValidateSchema(sFilename, IOStreamXMLFile::Label, false);
      return label.readLabel(sFilename);
    }
    catch(IOStreamErrorType &e){
      //that's because NR tests still contain obsolete Labels represented as (qualitative) Data
      //NB:  Labels as Data representation is not suitable because zeroes are not allowed as
      //qualitative data values but zeroes are valid label values in semi-supervized classification
      //TO BE REMOVED ASAP
      ValidateSchema(sFilename, IOStreamXMLFile::Data);
      return label.readLabelAsData(sFilename);
    }
  }

  unique_ptr<ProbaDescription> NodeOpOutput::readProba(string sFilename) {
	ValidateSchema(sFilename, IOStreamXMLFile::Data);
    //throw IOStreamErrorType::badLoadXML;

	DomProba proba;

	return proba.readProba(sFilename);
  }
  template <class T>
  static std::string join(T vect, int64_t sz,  std::string sep){
    std::stringstream acc;
    for(int64_t i = 0; i < sz; i++){
      if(i > 0)
        acc << sep;
      acc << std::to_string(vect[i]);
    }
    return acc.str();
  }
  template std::string join(std::vector<int64_t> vect, int64_t sz,  std::string sep);
  template std::string join(int64_t* vect, int64_t sz,  std::string sep);
  
  void NodeOpOutput::writeOutputExt(ClusteringModelOutput* output,  xmlpp::Element *outputElement, string str){}
  void NodeOpOutput::writeOutputExt(LearnModelOutput* output,  xmlpp::Element *outputElement, string str){
    /*      <xsd:group name="LearnOutputGroup">
            <xsd:sequence>
            <xsd:group ref="OutputGroup"/>
            <xsd:element name="CVLabels" minOccurs="0" type="VectOfIntType"/>
            <xsd:element name="CVClassification" minOccurs="0" type="MAPType"/>
            <xsd:element name="MapClassification" type="MAPType"/>
            <xsd:element name="MapErrorRate" type="DoubleOrHex"/>
            </xsd:sequence>
            </xsd:group>
    */


    std::vector<int64_t> labels = output->getLabelDescription()->getLabel()->getLabel();
    int64_t nbCluster = output->getNbCluster();
    if(output->getCVLabel()){
      //<xsd:element name="CVLabels" minOccurs="0" type="VectOfIntType"/>
      std::vector<int64_t> cvLabels = output->getCVLabel()->getLabel()->getLabel();      
      xmlpp::Element *cvElement = outputElement->add_child("CVLabels");
      /*for(int64_t i=0;i<cvLabels.size();i++){
        xmlpp::Element *cellElement = cvElement->add_child("cell");
        cellElement->add_child_text(std::to_string(cvLabels[i]));          
        }*/
      //cvElement->add_child_text(join(cvLabels, cvLabels.size(), " "));
      string labelFilename = str + "CVLabels";
      cvElement->add_child_text(labelFilename + ".mxl");
      //create label file (.mxl)
      LabelDescription *labelDesc = new LabelDescription(cvLabels.size(), cvLabels);
      DomLabel domLabel(labelDesc, labelFilename);
      
      //end CVLabels
      int64_t** tab = output->getCVLabel()->getLabel()->getClassificationTab(labels, nbCluster);
      // <xsd:element name="CVClassification" minOccurs="0" type="MAPType"/>
      xmlpp::Element *cvClsElement = outputElement->add_child("CVClassification");
      xmlpp::Element *mapErrElement = cvClsElement->add_child("ErrorRate");
      mapErrElement->add_child_text(custom_dtos(output->getCVLabel()->getLabel()->getErrorRate(labels)));
      xmlpp::Element *mapClassification = cvClsElement->add_child("Classification");
      for(int64_t i=0;i<nbCluster;i++){
        xmlpp::Element *rowElement = mapClassification->add_child("row");
        /*for(int64_t j=0;j<nbCluster;j++){
          xmlpp::Element *cellElement = rowElement->add_child("cell");
          cellElement->add_child_text(std::to_string(tab[i][j]));
          }*/
        rowElement->add_child_text(join(tab[i], nbCluster, " "));
      }
    }
    //<xsd:element name="MapClassification" type="MAPType"/>
    int64_t** tab = output->getLabelDescription()->getLabel()->getClassificationTab(labels, nbCluster);
    xmlpp::Element *mapElement = outputElement->add_child("MapClassification");
    xmlpp::Element *mapErrElement = mapElement->add_child("ErrorRate");
    mapErrElement->add_child_text(custom_dtos(output->getLabelDescription()->getLabel()->getErrorRate(labels)));
    xmlpp::Element *mapClassification = mapElement->add_child("Classification");    
    for(int64_t i=0;i<nbCluster;i++){
      xmlpp::Element *rowElement = mapClassification->add_child("row");
      /*for(int64_t j=0;j<nbCluster;j++){
        xmlpp::Element *cellElement = rowElement->add_child("cell");
        cellElement->add_child_text(std::to_string(tab[i][j]));
        }*/
      rowElement->add_child_text(join(tab[i], nbCluster, " "));
    }
    //<xsd:element name="MapErrorRate" type="DoubleOrHex"/>
  }
  
  template <class T>
  void NodeOpOutput::writeOutput(T* output, 
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
          /*if (IOMODE == IoMode::BINARY) {
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
            }*/
          valueElement->add_child_text(custom_dtos(output->getCriterionOutput(criterionName[i]).getValue()));
        }
        else {
          xmlpp::Element *errorElement = criterionElement->add_child("Error");
          errorElement->add_child_text(output->getStrategyRunError().what());
        }
      }
      //Likelihood
      xmlpp::Element *likelihoodElement = outputElement->add_child("Likelihood");
      /*if (IOMODE == IoMode::BINARY) {
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
        }*/
      likelihoodElement->add_child_text(custom_dtos(output->getLikelihood()));
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
      labelElement->add_child_text(labelFilename + ".mxl");
      //create label file (.mxl)
      DomLabel domLabel(output->getLabelDescription(), labelFilename);

      //proba
      xmlpp::Element *probaElement = outputElement->add_child("Proba");
      string probaFilename = str + "Proba" + std::to_string(numOutput);
      probaElement->add_child_text(probaFilename + ".mxd");
      //create proba file (.mxd)
      DomProba docProba(output->getProbaDescription(), probaFilename);
      writeOutputExt(output, outputElement, str);
      if(MASSICCC!=0){
        string entropyBasename = "entropy" + std::to_string(numOutput) + ".txt";
        string dest = PROJECT_DIRNAME + "/"+entropyBasename;
        rename(entropyBasename.c_str(), dest.c_str());
        xmlpp::Element *massicccElement =  outputElement->add_child("Custom");
        massicccElement->add_child_text(normalizeFilename(dest));
        
      }
      
	}
	else {
      //error
      xmlpp::Element *errorElement = outputElement->add_child("Error");
      errorElement->add_child_text(output->getStrategyRunError().what());
	}
  }
  template void NodeOpOutput::writeOutput<ClusteringModelOutput>(ClusteringModelOutput* output, 
                                          const std::vector<CriterionName> & criterionName, string str, int64_t numOutput);
  template void NodeOpOutput::writeOutput<LearnModelOutput>(LearnModelOutput* output, 
                                          const std::vector<CriterionName> & criterionName, string str, int64_t numOutput);

  //writePredictOutput
  void NodeOpOutput::writePredictOutput(PredictModelOutput* output, 
                                         string str, int64_t numOutput)
  {

    xmlpp::Element *outputElement = _rootOutput->add_child("Output");

	if (output->getStrategyRunError() == NOERROR) {

      //label
      xmlpp::Element *labelElement = outputElement->add_child("Label");
      string labelFilename = str + "Label" + std::to_string(numOutput);
      labelElement->add_child_text(labelFilename + ".mxl");
      //create label file (.mxl)
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
