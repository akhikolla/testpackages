/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMDomClusteringProject.cpp  description
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

#include "mixmod_iostream/DomOpProject.h"
#include "mixmod_iostream/NodeOpInput.h"
#include "mixmod_iostream/NodeOpOutput.h"
#include "mixmod/Clustering/ClusteringMain.h"

namespace XEM {

  DomOpProject::DomOpProject() : DomProject() {
  }
  
  DomOpProject::~DomOpProject() {
  }
  
  DomOpProject::DomOpProject(xmlpp::Element *root) : DomProject(root) {
  }
  
  void DomOpProject::writeMixmodXml(string& s, ClusteringMain * cMain) {
    _root->set_namespace_declaration("http://www.w3.org/2001/XMLSchema-instance", "xsi");
    _root->set_attribute("type", "Clustering", "xsi");
	//insertion of input
	NodeInput * inputNode = new NodeOpInput(cMain->getInput(), s);
	_root->import_node(inputNode->getRoot());
    
	//insertion of output
	if (cMain->getOutput() != NULL) {
      NodeOutput * outputNode = new NodeOpOutput(cMain->getOutput(), s);
      _root->import_node(outputNode->getRoot());
	}
  }
  void DomOpProject::writeMixmodXml(string& s, LearnMain * lMain) {
    _root->set_namespace_declaration("http://www.w3.org/2001/XMLSchema-instance", "xsi");
    _root->set_attribute("type", "Learn", "xsi");
	//insertion of input
	NodeInput * inputNode = new NodeOpInput(dynamic_cast<LearnInput*>(lMain->getInput()), s);
	_root->import_node(inputNode->getRoot());
    
	//insertion of output
	if (lMain->getLearnOutput() != NULL) {
      NodeOutput * outputNode = new NodeOpOutput(lMain->getLearnOutput(),lMain->getInput()->getCriterionName(), s);
      _root->import_node(outputNode->getRoot());
	}
  }

  void DomOpProject::writeMixmodXml(string& s, PredictMain * pMain) {
    _root->set_namespace_declaration("http://www.w3.org/2001/XMLSchema-instance", "xsi");
    _root->set_attribute("type", "Predict", "xsi");
	//insertion of input
	NodeInput * inputNode = new NodeOpInput(dynamic_cast<PredictInput*>(pMain->getInput()), s);
	_root->import_node(inputNode->getRoot());
    
	//insertion of output
	if (pMain->getPredictOutput() != NULL) {
      NodeOutput * outputNode = new NodeOpOutput(pMain->getPredictOutput(), s);
      _root->import_node(outputNode->getRoot());
	}


  }
  

  template<class T>
  void DomOpProject::readXmlFillIn(T *cInput) {
	//read the XML file which the case is clustering then fill the ClusteringInput or LearnInput
    
	//root in input
    xmlpp::Element* rootInput = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Input")); 
    NodeOpInput nodeInput(rootInput);

	//fill cInput from nodeInput (input node)
	nodeInput.readXmlCommand(*cInput);


  }
  template void   DomOpProject::readXmlFillIn<ClusteringInput>(ClusteringInput*);
  template void   DomOpProject::readXmlFillIn(LearnInput*);  

  PredictInput * DomOpProject::readXmlPredictInput(){
    xmlpp::Element* rootInput = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Input")); 
    NodeOpInput nodeInput(rootInput);

	//fill cInput from nodeInput (input node)
	return nodeInput.readXmlPredictInput();


  }

  static void setModelOutput(ClusteringOutput* cOutput, vector<ClusteringModelOutput*> vcmo){
    cOutput->setClusteringModelOutput(vcmo);
  }
  static void setModelOutput(LearnOutput* lOutput, vector<LearnModelOutput*> vlmo){
    lOutput->setLearnModelOutput(vlmo);
  }
  static void setModelOutput(PredictOutput* lOutput, vector<PredictModelOutput*> vlmo){
    lOutput->setPredictModelOutput(vlmo);
  }
  
  template<typename T, typename U>
  void DomOpProject::readXmlFillOut(T  * cOutput, Input *inp) {
	//read the XML file which the case is clustering then fill the ClusteringOutput

	//root in input
    xmlpp::Element* rootOutput = dynamic_cast<xmlpp::Element*>(_root->get_first_child("ListOutput"));
    xmlpp::Element* nodeElementOutput;
    vector<U *> vClusteringModelOutput;
    if (rootOutput) {
      //xmlpp::Node::NodeList
      auto children = rootOutput->get_children();
      for (auto it=children.begin(); it != children.end(); ++it){
        xmlpp::Element* nodeElementOutput = dynamic_cast<xmlpp::Element*>(*it);
        if(!nodeElementOutput) continue;
        NodeOpOutput nodeOutput(nodeElementOutput);
        U * modelOutput = nodeOutput.read4Output<U>(inp); //nodeOutput.readClustering();
        vClusteringModelOutput.push_back(modelOutput);
      }
      //cOutput->setClusteringModelOutput(vClusteringModelOutput);
      setModelOutput(cOutput, vClusteringModelOutput);
    } else {
      delete cOutput;
      cOutput = NULL;
    }
  }
  template void DomOpProject::readXmlFillOut<ClusteringOutput, ClusteringModelOutput>(ClusteringOutput*, Input *inp);
  template void DomOpProject::readXmlFillOut<LearnOutput, LearnModelOutput>(LearnOutput*, Input *inp);
  template void DomOpProject::readXmlFillOut<PredictOutput, PredictModelOutput>(PredictOutput*, Input *inp);    
}  
  /*
  void DomOpProject::readXmlFillOut(LearnOutput  * cOutput ) {}
  
} //end namespace
  */
