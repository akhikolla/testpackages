/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMNodeClusteringInput.cpp  description
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

#include "mixmod_iostream/NodeClusteringInput.h"
#include "mixmod_iostream/IOStreamUtil.h"
#include "mixmod/Kernel/Algo/EMAlgo.h"
#include <algorithm>
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod_iostream/DomLabel.h"
#include "mixmod_iostream/DomParameter.h"
#include "mixmod/Clustering/ClusteringInput.h"
#include "mixmod/Clustering/ClusteringStrategy.h"
#include "mixmod/Clustering/ClusteringStrategyInit.h"
#include "mixmod/Kernel/Algo/Algo.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/Criterion/Criterion.h"
#include "mixmod/Kernel/IO/ParameterDescription.h"
namespace XEM {

  NodeClusteringInput::NodeClusteringInput() : NodeInput() {
  }

  NodeClusteringInput::~NodeClusteringInput() {
  }

  NodeClusteringInput::NodeClusteringInput(ClusteringInput * input, string & s) 
    : NodeInput(input, s) {
	writeListModel(input);
	writeNbClusterNode(input);
	writeStrategyNode(input, s);
	writeCriterionNode(input);
	writePartitionNode(input, s);
	writeWeightsNode(input, s);
  }
  
  NodeClusteringInput::NodeClusteringInput( xmlpp::Element * rootInput ) : NodeInput(rootInput) {
    
  }

  //read the clustering to fill XEMNvInput
  void NodeClusteringInput::readClustering(ClusteringInput & input) {

	//read data from file (method from XEM::NodeInput)
	DataDescription dataDescription = readDataNode();

	//read NbCluster
	vector<int64_t> nbCluster = readNbClusterNode();

	//filling of input
	// TODO: if 'readDataNode()' above returned a pointer, we wouldn't need to (deep) copy it here.
	//       same remark for nb[Nb]Cluster (pointer...) vector.
	input.cloneInitialisation(nbCluster, dataDescription);

	readModelNode(input);
	readStrategyNode(input);
	readCriterionNode(input);
	readPartitionNode(input);
	readWeightsNode(input);
  }

  void NodeClusteringInput::writeListModel(ClusteringInput * input) {
    //model
    //XEMBinaryData * binaryData = 
    //dynamic_cast<XEMBinaryData *>((input->getDataDescription()).
    //main node
    xmlpp::Element *listModel = _rootInput->add_child("ListModel");
  
    //model type node
    xmlpp::Element *EDDAModel =  listModel->add_child("EDDAModel");
    xmlpp::Element *HDModel =  listModel->add_child("HDModel");
    xmlpp::Element *binaryModel =  listModel->add_child("BinaryModel");
    xmlpp::Element *compositeModel =  listModel->add_child("CompositeModel");
	//fill with each model name
	for (int64_t i = 0; i < input->getModelType().size(); i++) {
      ModelName modelName = input->getModelType()[i]->getModelName();
      //text for XML
      xmlpp::Element * whichModel = NULL;
      string text = ModelNameToString(modelName);
      
      if (isEDDA(modelName)) {//case EDDA
        whichModel = EDDAModel;
      }
      else if (isHD(modelName)) {//case HD
        whichModel = HDModel;
      }
      else if (isBinary(modelName)) {//case Binary
        whichModel = binaryModel;
      }
      else if (isHeterogeneous(modelName)) {
        whichModel = compositeModel;
      }
      whichModel = whichModel->add_child("Model");
      whichModel->add_child_text(text);
	}

	//selection of node non-empty node  
	//EDDAModel
    if(!EDDAModel->get_first_child()) listModel->remove_child(EDDAModel);

	//HDModel
	if (HDModel->get_first_child()) {
      //subDimensionEqual
      if (input->getModelType()[0]->_subDimensionEqual != 0) {
        xmlpp::Element *d = HDModel->add_child("D");
        d->add_child_text(std::to_string(input->getModelType()[0]->_subDimensionEqual));
      }

      //subDimensionFree
      if (input->getModelType()[0]->_nbSubDimensionFree != 0) {
        xmlpp::Element * listDk = HDModel->add_child("ListDk");
        for (int64_t i = 0; i < input->getModelType()[0]->_nbSubDimensionFree; i++) {
          xmlpp::Element *dk =  listDk->add_child("Dk");
          dk->add_child_text(std::to_string(input->getModelType()[0]->getTabSubDimensionFreeI(i)));
        }
      }
	} else {
      listModel->remove_child(HDModel);
    }

	//binaryModel
    if(!binaryModel->get_first_child()) listModel->remove_child(binaryModel);
	//binaryModel
    if(!compositeModel->get_first_child()) listModel->remove_child(compositeModel);
    if(!listModel->get_first_child()) _rootInput->remove_child(listModel);
  }

  void NodeClusteringInput::writeNbClusterNode(ClusteringInput * input) {
    
	//nbCluster
    xmlpp::Element *listCluster = _rootInput->add_child("ListNbCluster");
    vector<int64_t> vectNbCluster(input->getNbCluster());
    for (int64_t i = 0; i < vectNbCluster.size(); ++i) {
      xmlpp::Element *nbCluster = listCluster->add_child("NbCluster");
      nbCluster->add_child_text(std::to_string(vectNbCluster[i]));
    }
    //_rootInput.appendChild(listCluster);
  }

  void NodeClusteringInput::writeStrategyNode(ClusteringInput * input, string & s) {
    
	//strategy
    xmlpp::Element *strategy = _rootInput->add_child("Strategy");
    //nbTry
    xmlpp::Element *nbTry = strategy->add_child("NbTry");
    ClusteringStrategy * cStrategy = dynamic_cast<ClusteringStrategy*> (input->getStrategy());
    nbTry->add_child_text(std::to_string(cStrategy->getNbTry()));
    
	//strategy
    writeInitNode(strategy, input, s);
    xmlpp::Element *listAlgo = strategy->add_child("ListAlgo");
    for (int64_t i = 0; i < cStrategy->getNbAlgo(); ++i) {
      writeAlgoNode(listAlgo, cStrategy->getAlgo(i));
    }
  }

  // criterion
  void NodeClusteringInput::writeCriterionNode(ClusteringInput * input) {

    xmlpp::Element *listCriterion = _rootInput->add_child("ListCriterion");
    for (int64_t i = 0; i < input->getNbCriterion(); ++i) {
      xmlpp::Element *criterion = listCriterion->add_child("Criterion");
      criterion->add_child_text(CriterionNameToString(input->getCriterionName(i)));
    }
  }

  //partition
  void NodeClusteringInput::writePartitionNode(ClusteringInput * input, string & s) {

    if (input->getKnownPartition()) {
      //partition
      xmlpp::Element *partition =  _rootInput->add_child("Partition");
      //partitionFilename
      string str = s + "Partition";
      partition->add_child_text(str + ".mxl");
      DomLabel doc(input->getKnownPartition(), str);
    }
  }

  //weights
  void NodeClusteringInput::writeWeightsNode(ClusteringInput * input, string & s) {

	// HACK [bauder: I don't like that; we could use a bool e.g. "_customWeights"...]
	// Check if some weight differs from 1.0. If yes, a weights files is provided.
	int n = input->getData()->getNbSample();
	const double* weights = input->getData()->getWeight();
	for (int64_t i=0; i<n; i++) {
      if (weights[i] != 1.0) {
        //partitionFilename
        xmlpp::Element *weights = _rootInput->add_child("Weights");
        weights->add_child_text(s + "Weights.mxw");			
        break;
      }
	}
  }

  void NodeClusteringInput::writeInitNode(xmlpp::Element *strategyElement, 
		ClusteringInput * cInput, string & sFilename)
  {
	ClusteringStrategy * strategy = dynamic_cast<ClusteringStrategy*> (cInput->getStrategy());

	//Strategy Init
    xmlpp::Element *init = strategyElement->add_child("Init");
	switch (strategy->getStrategyInit()->getStrategyInitName()) {
	case RANDOM:
	case CEM_INIT:
      {
		//init
        xmlpp::Element *nbTryInInit = init->add_child("NbTry");
        nbTryInInit->add_child_text(std::to_string(strategy->getStrategyInit()->getNbTry()));
        init->set_attribute("xsi:type", StrategyInitNameToString(strategy->getStrategyInit()->getStrategyInitName()));
		break;
	}
	case SMALL_EM:
	{
		//init
      xmlpp::Element *nbTryInInit = init->add_child("NbTry");
      nbTryInInit->add_child_text(std::to_string(strategy->getStrategyInit()->getNbTry()));
      //stopRule
      xmlpp::Element *stopRule = init->add_child("StopRule");
      switch (strategy->getStrategyInit()->getStopName()) {
      case NBITERATION:
		{
          xmlpp::Element *iteration = stopRule->add_child("NbIteration");
          iteration->add_child_text(std::to_string(strategy->getStrategyInit()->getNbIteration()));
          break;
		}
      case EPSILON:
		{
          xmlpp::Element *epsilon = stopRule->add_child("Epsilon");
          epsilon->add_child_text(std::to_string(strategy->getStrategyInit()->getEpsilon()));
          break;
		}
      case NBITERATION_EPSILON:
		{
          xmlpp::Element *iteration = stopRule->add_child("NbIteration");
          iteration->add_child_text(std::to_string(strategy->getStrategyInit()->getNbIteration()));
          xmlpp::Element *epsilon = stopRule->add_child("Epsilon");
          epsilon->add_child_text(std::to_string(strategy->getStrategyInit()->getEpsilon()));
          break;
		}
      }
      //init.appendChild(stopRule);
      init->set_attribute("xsi:type", StrategyInitNameToString(strategy->getStrategyInit()->getStrategyInitName()));
      break;
	}
	case SEM_MAX:
      {
        //stopRule
        xmlpp::Element *stopRule = init->add_child("StopRule");
        xmlpp::Element *iteration = stopRule->add_child("NbIteration");
        iteration->add_child_text(std::to_string(strategy->getStrategyInit()->getNbIteration()));
        init->set_attribute("xsi:type", StrategyInitNameToString(strategy->getStrategyInit()->getStrategyInitName()));
        break;
      }
	case USER:
      {
		for (int64_t i = 0; i < strategy->getStrategyInit()->getNbInitParameter(); ++i) {
          //if (cInput->getModelType().size() == 1) {
              xmlpp::Element *filename = init->add_child("Parameter");
              string parameterFilename = sFilename + "InitParameter" + std::to_string(i+1);
              filename->add_child_text(parameterFilename + ".mxp");
              //}
        DomParameter dpar(cInput, parameterFilename);
		}
		init->set_attribute("xsi:type", "PARAMETER");

		break;
	}
	case USER_PARTITION:
      {

		for (int64_t i = 0; i < strategy->getStrategyInit()->getNbPartition(); ++i) {
          xmlpp::Element *filename = init->add_child("Partition");
          string partitionFilename = sFilename + "InitPartition"+ std::to_string(i+1) + ".mxl";
          filename->add_child_text(partitionFilename);
          string str =  sFilename + "InitPartition" + std::to_string(i + 1);
          //DomLabel doc(cInput->getKnownPartition(), str);
          DomLabel doc(cInput->getStrategy()->getStrategyInit()->getPartition(0), str);

		}
		init->set_attribute("xsi:type", "PARTITION");
	}
	default:
		break;
	}
  }

  void NodeClusteringInput::writeAlgoNode(xmlpp::Element *listAlgo, const Algo * algorithm) {

	//algo
    xmlpp::Element *algo = listAlgo->add_child("Algo");
	algo->set_attribute("xsi:type", AlgoNameToString(algorithm->getAlgoName()));
    xmlpp::Element *stopRule = algo->add_child("StopRule");
    
	//algostop
	switch (algorithm->getAlgoStopName()) {
	case NBITERATION:
      {
        xmlpp::Element *iteration = stopRule->add_child("NbIteration");
        iteration->add_child_text(std::to_string(algorithm->getNbIteration()));
        break;
      }
	case EPSILON:
      {
        xmlpp::Element *epsilon = stopRule->add_child("Epsilon");
        epsilon->add_child_text(std::to_string(algorithm->getEpsilon()));
        break;
      }
	case NBITERATION_EPSILON:
      {
        xmlpp::Element *iteration = stopRule->add_child("NbIteration");
        iteration->add_child_text(std::to_string(algorithm->getNbIteration()));
        xmlpp::Element *epsilon = stopRule->add_child("Epsilon");
        epsilon->add_child_text(std::to_string(algorithm->getEpsilon()));
        break;
      }
	}
  }

  //read the model node
  void NodeClusteringInput::readModelNode(ClusteringInput & input) {
    if (_rootInput) {
      //_root is clustering node
      xmlpp::Element *elementListModel = dynamic_cast<xmlpp::Element*>(_rootInput->get_first_child("ListModel"));
      if(!elementListModel) return;
      //fisrt model Node
      auto children = elementListModel->get_children();
      int64_t compt = 0;
      for (auto it=children.begin(); it != children.end(); ++it){
        xmlpp::Element* elementTypeModel = dynamic_cast<xmlpp::Element*>(*it);
        if(!elementTypeModel) continue; //the current node isn't an element
        //choice between 2 model types Binary or EDDA (HDModel not available in Clustering) 
        //with same following process 
        auto m_children = elementTypeModel->get_children();
        vector<string> vModelName(0);
        //cross all the model names
        for (auto it2=m_children.begin(); it2 != m_children.end(); ++it2){
          xmlpp::Element* elementModel = dynamic_cast<xmlpp::Element*>(*it2);
          if(!elementModel) continue;
          //model name of the current node
          string strModelName = elementModel->get_child_text()->get_content(); //elementModel.text().toStdString();
        
          if (find(vModelName.begin(), vModelName.end(), strModelName) == vModelName.end()) {
          
            //fill the XEMNvInput with the model current Model, set or insert according to case
            if (compt == 0) {
              input.setModelType(new ModelType(StringToModelName(strModelName)), compt);
            }
            else {
              input.insertModelType(new ModelType(StringToModelName(strModelName)), compt);
            }
          
            compt++;
            vModelName.push_back(strModelName);
          }
        }
      }
    }
  }
  void NodeClusteringInput::setInitPartition(string sFilename,ClusteringStrategy * strat){
  	//-------
	//load file in this
	//-------
    xmlpp::DomParser parser;
    parser.parse_file(sFilename);
    xmlpp::Document *doc = parser.get_document();    
    xmlpp::Element *root = doc->get_root_node();//documentElement();
     //------------------------
    //Declaration of variables
    //------------------------
    xmlpp::Element *elementNbSample, *elementNbCluster, *elementFormat, *elementType, *elementFilename;
    
    //nbSample
    elementNbSample = dynamic_cast<xmlpp::Element*>(root->get_first_child("NbSample"));
    int64_t nbSample = std::stoll(elementNbSample->get_child_text()->get_content());
    
    //nbCluster
    elementNbCluster = dynamic_cast<xmlpp::Element*>(root->get_first_child("NbCluster"));
    int64_t nbCluster = std::stoll(elementNbCluster->get_child_text()->get_content());
    
    //Format
    elementFormat = dynamic_cast<xmlpp::Element*>(root->get_first_child("Format"));
    FormatNumeric::FormatNumericFile format = 
      StringToFormatNumericFile(elementFormat->get_child_text()->get_content());

    //Type
    elementType = dynamic_cast<xmlpp::Element*>(root->get_first_child("Type"));
    TypePartition::TypePartition type = 
      StringToTypePartition(elementType->get_child_text()->get_content());

    //Parameter Filename
    elementFilename = dynamic_cast<xmlpp::Element*>(root->get_first_child("Filename"));
    string partFilename = elementFilename->get_child_text()->get_content();
    //Partition * part = new Partition(nbSample, nbCluster, pfilename);
    //Partition ** tabPartition = new Partition*[1];
    //tabPartition[0] = new Partition();
    //strat->setTabPartition(tabPartition, 1);//1 is the number of partition
	std::ifstream partitionFile(partFilename.c_str(), ios::in);
	if (! partitionFile.is_open()) {
		THROW(InputException, wrongPartitionFileName);
	}    
    Partition * part = new Partition();
    part->setDimension(nbSample, nbCluster);
    part->setPartitionFile(partFilename, type);
    partitionFile >> *part;
    strat->setInitPartition(part, 0);//0 is the place where the partition will be stocked 
  }
  void NodeClusteringInput::readStrategyNode(ClusteringInput & input) {
	if (!_rootInput) return;
    xmlpp::Element* elementStrategy = dynamic_cast<xmlpp::Element*>(_rootInput->get_first_child("Strategy"));
    if (elementStrategy) {
      //strategy
      ClusteringStrategy * strat = new ClusteringStrategy();
      auto children = elementStrategy->get_children();
      for (auto it=children.begin(); it != children.end(); ++it){
        xmlpp::Element* n = dynamic_cast<xmlpp::Element*>(*it);
        if(!n) continue;
        //nbTry
        if (n->get_name() == "NbTry") {
          strat->setNbTry(std::stoll(n->get_child_text()->get_content()));
        }
        
        //init
        readInitNode(strat, n);
        
        //algo
        readAlgoNode(strat, n);
      }
      input.setStrategy(strat);
    }
	
  }
  
  //read the InitNode
  void NodeClusteringInput::readInitNode(ClusteringStrategy * strat, xmlpp::Element *n) {

	if (n->get_name() != "Init") return;
		//initName
    strat->setStrategyInitName(
                               StringToStrategyInitName(n->get_attribute_value("type", "xsi")));
      
    switch (strat->getStrategyInit()->getStrategyInitName()) {
    case RANDOM:
    case CEM_INIT:
      {
        //nbTryInInit
        xmlpp::Element* nbTryNode = dynamic_cast<xmlpp::Element*>(n->get_first_child("NbTry"));
        if (nbTryNode) {
          strat->setNbTryInInit(std::stoll(nbTryNode->get_child_text()->get_content()));
        }
        break;
      }
    case SMALL_EM:
      {
        //nbTryInInit
        xmlpp::Element* nbTryNode = dynamic_cast<xmlpp::Element*>(n->get_first_child("NbTry"));
        if (nbTryNode) {
          strat->setNbTryInInit(std::stoll(nbTryNode->get_child_text()->get_content()));
        }
        
          //stopRule
        xmlpp::Element* stopRuleNode = dynamic_cast<xmlpp::Element*>(n->get_first_child("StopRule"));
        if (stopRuleNode) {
          xmlpp::Element* iteration = dynamic_cast<xmlpp::Element*>(stopRuleNode->get_first_child("NbIteration"));
          xmlpp::Element* epsilon = dynamic_cast<xmlpp::Element*>(stopRuleNode->get_first_child("Epsilon"));
          if (iteration && epsilon) {//NBITERATION_EPSILON
            strat->setStopNameInInit(NBITERATION_EPSILON);
            strat->setNbIterationInInit(std::stoll(iteration->get_child_text()->get_content()));//iteration.text().toLongLong());
            strat->setEpsilonInInit(std::stod(epsilon->get_child_text()->get_content())); //epsilon.text().toDouble());
            
          }
          else if (iteration) {//NBITERATION
            strat->setStopNameInInit(NBITERATION);
            strat->setNbIterationInInit(std::stoll(iteration->get_child_text()->get_content()));
            
          }
          else if (epsilon) {//EPSILON
            strat->setStopNameInInit(EPSILON);
            strat->setEpsilonInInit(std::stod(epsilon->get_child_text()->get_content()));
          }
        }
        break;
      }
    case SEM_MAX:
      {
        //stopRule
        xmlpp::Element* stopRuleNode = dynamic_cast<xmlpp::Element*>(n->get_first_child("StopRule"));        
        if (stopRuleNode) {
          xmlpp::Element* iteration = dynamic_cast<xmlpp::Element*>(stopRuleNode->get_first_child("NbIteration"));
          if (iteration) {//NBITERATION
            strat->setNbIterationInInit(std::stoll(iteration->get_child_text()->get_content()));
          }
        }
        break;
      }
    case USER_PARTITION:
      {
        xmlpp::Element* partitionNode = get_first_child_element(n);
        if (partitionNode) {
          
          //chech the .mxl
          string filename = partitionNode->get_child_text()->get_content();
          ValidateSchema(filename, IOStreamXMLFile::Partition);
          //throw IOStreamErrorType::badLoadXML;
          //Move the next part "partition filling" to XEMDomPartition
          //partition filling
          /*XEMPartition ** tabPartition = new XEMPartition*[1];
            tabPartition[0] = new XEMPartition();
            strat->setTabPartition(tabPartition, 1);//1 is the number of partition
            if (!filename.empty()){
            strat->setInitPartition(filename, 0);//0 is the place where the partition will be stocked 
            }*/

          setInitPartition(filename, strat);
        }
        break;
      }
    case USER:
      {
        //QDomElement parameterNode = n.firstChild().toElement(); //parameterNode (filename)
        xmlpp::Element* parameterNode = get_first_child_element(n);
        if (parameterNode) {
          //check the .mxp
          //string filename = parameterNode.text().toStdString();
          string filename = parameterNode->get_child_text()->get_content();
          ValidateSchema(filename, IOStreamXMLFile::Parameter);
          //throw IOStreamErrorType::badLoadXML;
        
        
        //Move to XEMDomParamter
        /*XEMParameter ** tabParameter = new XEMParameter*[1];
          tabParameter[0] = new XEMParameter();	
          strat->setTabInitParameter(tabParameter, nbParameters);
          if (!filename.empty()){
          strat->setInitParam(filename, compt);
          }*/
          
        DomParameter domParam;
        XEM::ParameterDescription * paramDesc = domParam.readParameter(filename);
        XEM::Parameter ** tabParameter = new XEM::Parameter*[1];
        tabParameter[0] = paramDesc->getParameter()->clone();
        //delete paramDesc;
        strat->setTabInitParameter(tabParameter, 1);
        }  
        break;
      }
    }
	
  }
  
  //read the AlgoNode
  void NodeClusteringInput::readAlgoNode(ClusteringStrategy * strat, xmlpp::Element *n) {

	if (n->get_name() != "ListAlgo") return;
    auto children = n->get_children();
    int64_t compt = 0;
    //while (!algoNode.isNull()) {
    for (auto it=children.begin(); it != children.end(); ++it){
      xmlpp::Element* algoNode = dynamic_cast<xmlpp::Element*>(*it);
      if(!algoNode) continue;
      switch (StringToAlgoName(algoNode->get_attribute_value("type", "xsi"))) {
      case EM:
      case CEM:
        {
          //algo
          if (compt != 0) {//after the first algo
            strat->insertAlgo(EM, compt); //Algo by default
          }
          //change AlgoName
          strat->setAlgo(
                         StringToAlgoName(algoNode->get_attribute_value("type", "xsi")), compt);
          
          //algoStopRule
          xmlpp::Element* stopRuleNode = dynamic_cast<xmlpp::Element*>(algoNode->get_first_child("StopRule"));        
          
          if (stopRuleNode) {
            xmlpp::Element* iteration = dynamic_cast<xmlpp::Element*>(stopRuleNode->get_first_child("NbIteration"));
            xmlpp::Element* epsilon = dynamic_cast<xmlpp::Element*>(stopRuleNode->get_first_child("Epsilon"));
            
            if (iteration && epsilon) {//NBITERATION_EPSILON
              strat->setAlgoStopRule(NBITERATION_EPSILON, compt);
              strat->setAlgoIteration(compt, std::stoll(iteration->get_child_text()->get_content()));
              strat->setAlgoEpsilon(compt, std::stod(epsilon->get_child_text()->get_content()));
              
            }
            else if (iteration) {//NBITERATION
              strat->setAlgoStopRule(NBITERATION, compt);
              strat->setAlgoIteration(compt, std::stoll(iteration->get_child_text()->get_content()));
              
            }
            else if (epsilon) {//EPSILON
              strat->setAlgoStopRule(EPSILON, compt);
              strat->setAlgoEpsilon(compt, std::stod(epsilon->get_child_text()->get_content()));
            }
          }
          break;
        }
      case SEM:
        {
          //algo
          if (compt != 0) {//after the first algo	   
            strat->insertAlgo(EM, compt); //Algo by default
          }
          //change AlgoName
          strat->setAlgo(
                         StringToAlgoName(algoNode->get_attribute_value("type", "xsi")), compt);
          
          //algoStopRule
          xmlpp::Element* stopRuleNode = dynamic_cast<xmlpp::Element*>(algoNode->get_first_child("StopRule"));        
          if (stopRuleNode) {
            xmlpp::Element* iteration = dynamic_cast<xmlpp::Element*>(stopRuleNode->get_first_child("NbIteration"));
            if (iteration) {
              strat->setAlgoIteration(compt, std::stoll(iteration->get_child_text()->get_content()));
            }
          }
          break;
        }
      }
      compt++;
    }
      //}
  }

  
  void NodeClusteringInput::readCriterionNode(ClusteringInput & input) {
    if (!_rootInput) return;
    xmlpp::Element *elementListCriterion = dynamic_cast<xmlpp::Element*>(_rootInput->get_first_child("ListCriterion"));
    if (!elementListCriterion) return;
    auto children = elementListCriterion->get_children();
    int64_t compt = 0;
    vector<string> vCriterionNameTmp(0);
   
    for (auto it=children.begin(); it != children.end(); ++it){
      xmlpp::Element* n = dynamic_cast<xmlpp::Element*>(*it);
      if(!n) continue;
      string sCriterionName = n->get_child_text()->get_content();
      //Criterion => insertion ou set according to the case
      if (n->get_name() == "Criterion" && 
          find(vCriterionNameTmp.begin(), vCriterionNameTmp.end(), sCriterionName) 
          == vCriterionNameTmp.end()) 
        {
          if (compt == 0) {
            //XEMNvInput has only one criterion by default
            input.setCriterion(StringtoCriterionName(sCriterionName), compt);
          }
          else {
            input.insertCriterion(StringtoCriterionName(sCriterionName), compt);
          }
          
          //filling of the tempory criterionName vector to prevent the duplicates
          vCriterionNameTmp.push_back(sCriterionName);
          
          compt++;
        }
    }
  }

  void NodeClusteringInput::readPartitionNode(ClusteringInput & input) {
	//read the partition node  
	if (!_rootInput) return;
    xmlpp::Element *elementPartition = dynamic_cast<xmlpp::Element*>(_rootInput->get_first_child("Partition"));
    string filename;
    if (elementPartition && input.getNbCluster().size() == 1) {
      filename = elementPartition->get_child_text()->get_content(); //.text().toStdString();
      ValidateSchema(filename, IOStreamXMLFile::Partition);
      //throw IOStreamErrorType::badLoadXML;
      
      // temporary [??] fix, assuming TXT format (TODO: ?!)
      // TODO: duplicated code from XEMClusteringMain * XEMIStream(...)
      xmlpp::DomParser parser;
      parser.parse_file(filename);
      xmlpp::Document *doc = parser.get_document();
      xmlpp::Element *_root = doc->get_root_node();
      xmlpp::Element *elementNbSample, *elementNbCluster, *elementFormat, *elementType, *elementFilename;

      //nbSample
      elementNbSample = dynamic_cast<xmlpp::Element*>(_root->get_first_child("NbSample"));
      int64_t nbSample = std::stoll(elementNbSample->get_child_text()->get_content());

      //nbCluster
      elementNbCluster = dynamic_cast<xmlpp::Element*>(_root->get_first_child("NbCluster"));
      int64_t nbCluster = std::stoll(elementNbCluster->get_child_text()->get_content());

      //Format
      elementFormat = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Format"));
      FormatNumeric::FormatNumericFile format = 
        StringToFormatNumericFile(elementFormat->get_child_text()->get_content());

      //Type
      elementType = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Type"));
      TypePartition::TypePartition type = 
        StringToTypePartition(elementType->get_child_text()->get_content());

      //Parameter Filename
      elementFilename = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Filename"));
      string partFilename = elementFilename->get_child_text()->get_content();


      NumericPartitionFile partitionFile;
      partitionFile._fileName = partFilename;
      partitionFile._format = format;
      partitionFile._type = type;

      input.insertKnownPartition(partitionFile);
    }
}

  void NodeClusteringInput::readWeightsNode(ClusteringInput & input) {
	//read the weights node  
	if (!_rootInput) return;
    xmlpp::Element *elementWeights = dynamic_cast<xmlpp::Element*>(_rootInput->get_first_child("Weights"));
    string filename;
    if (!elementWeights) return;
    filename = elementWeights->get_child_text()->get_content();
    ValidateSchema(filename, IOStreamXMLFile::Weights);
    //throw IOStreamErrorType::badLoadXML;
    
    // TODO: duplicated code from XEMClusteringMain * XEMIStream(...)
    xmlpp::DomParser parser;    
    parser.parse_file(filename);
    xmlpp::Document *doc = parser.get_document();
    xmlpp::Element *_root = doc->get_root_node();
    xmlpp::Element *elt = dynamic_cast<xmlpp::Element*>(_root->get_first_child("Filename"));
    string weightsFilename = elt->get_child_text()->get_content();    
    input.insertWeight(weightsFilename);


  }
      
} //end namespace

