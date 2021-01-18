/***************************************************************************
							 SRC/MIXMOD_IOSTREAM/XEMIOStreamUtil.cpp  description
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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
//#include <sys/sendfile.h>
#include <stdio.h>
#include <unistd.h>
#include "mixmod_iostream/IOStreamUtil.h"
//#include "mixmod_iostream/DomClusteringProject.h"
#include "mixmod_iostream/DomOpProject.h"
//#include "mixmod_iostream/DomDAProject.h"
#include "mixmod_iostream/NodeOpInput.h"
#include "mixmod_iostream/NodeOpOutput.h"

#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/GaussianData.h"
//#include "mixmod/Clustering/ClusteringInput.h"
#include "mixmod/Clustering/ClusteringMain.h"
//#include "mixmod/Clustering/ClusteringOutput.h"
//#include "mixmod/Clustering/ClusteringModelOutput.h"

#include "mixmod/Kernel/IO/QualitativeColumnDescription.h"
#include "mixmod/Kernel/IO/QuantitativeColumnDescription.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Clustering/ClusteringModelOutput.h"
#include "mixmod/Clustering/ClusteringOutput.h"
#include "mixmod/Kernel/IO/Input.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Parameter/Parameter.h"
#include "mixmod/Kernel/IO/Proba.h"
#include "mixmod/Kernel/IO/ProbaDescription.h"
#include "mixmod/Kernel/IO/Label.h"
#include "mixmod/Kernel/IO/LabelDescription.h"

namespace XEM {
  string PROJECT_DIRNAME = "xxx";
//IOStreamErrorTypeToString
string IOStreamErrorTypeToString(const IOStreamErrorType & errorType) {
	string res;
	switch (errorType) {
	case (IOStreamErrorType::noIOStreamError):
		res = "No error";
		break ;
	case (IOStreamErrorType::badIOStreamFormat):
		res = "Bad format";
		break ;
	case (IOStreamErrorType::badNumericFormat):
		res = "Bad numeric format";
		break ;
	case (IOStreamErrorType::badIOWriteFormat):
		res = "Bad format while writing";
		break ;
	case (IOStreamErrorType::badOpenedFile):
		res = "Bad file to open";
		break ;
	case (IOStreamErrorType::badLoadXML):
		res = "Bad XML while loading";
		break ;
	case (IOStreamErrorType::badSchema):
		res = "Bad XML Schema";
		break ;
	case (IOStreamErrorType::badXML):
		res = "Bad XML";
		break ;
	case (IOStreamErrorType::badElementInXML):
		res = "Bad element in XML";
		break ;
	case (IOStreamErrorType::badElementInDataXML):
		res = "Bad element in XML Data";
		break ;
	case (IOStreamErrorType::badColumnUsedInCreateMixmodDataFileFromUserDataFile):
		res = "Bad Column Used in create Data File";
		break ;
	case (IOStreamErrorType::errorInCreateMixmodDataFileFromUserDataFile):
		res = "Error in create Mixmod Data File";
		break ;
	case (IOStreamErrorType::notEnoughValuesInDataFile):
		res = "Error in create Mixmod Data File : not enough values in data file";
		break ;
	case (IOStreamErrorType::fileAlreadyExist):
		res = "File already exists";
		break ;
	case (IOStreamErrorType::fileDontExist):
		res = "File doesn't exist";
		break ;
	case (IOStreamErrorType::fileCantBeOpened):
		res = "File can't be opened";
		break;
	case (IOStreamErrorType::notAbsoluteFileDataName):
		res = "Data Filename must be an absolute filename";
		break ;
	}
	return res;
}

//------------//
// ISTREAM //
//------------//

// TODO [bauder]: 'bOnlyInput' should probably be removed; it may be better to create two methods:
//                 "ISTREAM_input(...)" and "ISTREAM_output(...)" (with proper names)

ClusteringMain * IStream(const string & s, IOStreamFormat format, bool bOnlyInput, IoMode iomode) {
	IOMODE = iomode; //TODO...
	ClusteringMain * res = NULL;
	// test informat
	if (format == IOStreamFormat::XML) {
		// XML format
		//-----------
		res = IStream_XML(s, bOnlyInput);
	}
	else if (format == IOStreamFormat::FLAT) {
		// FLAT format
		//-----------
		res = IStream_FLAT(s);
	}
	else {
		throw IOStreamErrorType::badIOStreamFormat;
	}

	return res;
}

//-----------------
//read the XML file
//-----------------
  ClusteringMain * IStream_XML_Clustering(const std::string& s, bool bOnlyInput, IoMode iomode) {
    IOMODE = iomode; //TODO...
    return IStream_XML(s, bOnlyInput);
  }
  ClusteringMain * IStream_XML(const std::string& s, bool bOnlyInput) {
    //take the absolute path
    const string str = s;
    ValidateSchema(str, IOStreamXMLFile::Project);
    xmlpp::DomParser parser;
    parser.parse_file(s);
    xmlpp::Document *doc = parser.get_document();
    xmlpp::Element *root = doc->get_root_node();  
    if ( root->get_name() != "Project" ) throw IOStreamErrorType::badIOStreamFormat;
    string xsitype = root->get_attribute_value("type", "xsi");
    if (xsitype.compare("Clustering") != 0) throw IOStreamErrorType::badXML;    
    //DomClusteringProject docClustering(root);
    IoModeManager ioMode = IoModeManager(root);  //set the IoMode defined by a FloatEncoding tag, if it exists    
    DomOpProject docClustering(root);
    //clusteringInput
    ClusteringInput * cInput = new ClusteringInput();

    //docClustering.readClustering(cInput);
    docClustering.readXmlFillIn<ClusteringInput>(cInput);
    //HACK: set nbVariables_binary and nbVariables_gaussian for CompositeParameter, in case of...
    //TODO: refactor...
    if (cInput->getDataType() == HeterogeneousData) {
      Global::nbVariables_binary = cInput->getData()->getBinaryData()->_pbDimension;
      Global::nbVariables_gaussian = cInput->getData()->getGaussianData()->_pbDimension;
    }
    if (cInput->getDataType() == QualitativeData || cInput->getDataType() == HeterogeneousData) {
      Global::vNbFactor.clear();
      for (int i=0; i<cInput->getData()->getBinaryData()->getPbDimension(); i++)
        Global::vNbFactor.push_back(cInput->getData()->getBinaryData()->getTabNbModality()[i]);
    }

    ClusteringOutput * cOutput = NULL;
    xmlpp::Element *listOutput = dynamic_cast<xmlpp::Element*>(root->get_first_child("ListOutput"));
    if (listOutput && !bOnlyInput) {
      //ClusteringOutput
      cOutput = new ClusteringOutput(cInput->getCriterionName());
      //docClustering.readClustering(cOutput);
      docClustering.readXmlFillOut<ClusteringOutput, ClusteringModelOutput>(cOutput, cInput);
    }

    cInput->finalize();
    return new ClusteringMain(cInput, cOutput);
}

  //////////////////////////////////////////////////////////////////::
//---------------------------------------------
//read the XML file 4 clustering with templates
//---------------------------------------------
  LearnMain * IStream_XML_Learn(const std::string& s, bool bOnlyInput, IoMode iomode) {
	IOMODE = iomode; //TODO...    
    //take the absolute path
    const string str = s;
    ValidateSchema(str, IOStreamXMLFile::Project);
    xmlpp::DomParser parser;
    parser.parse_file(s);
    xmlpp::Document *doc = parser.get_document();
    xmlpp::Element *root = doc->get_root_node();  
    if ( root->get_name() != "Project" ) throw IOStreamErrorType::badIOStreamFormat;
    string xsitype = root->get_attribute_value("type", "xsi");
    if (xsitype.compare("Learn") != 0)  throw IOStreamErrorType::badXML;    
    IoModeManager ioMode = IoModeManager(root);  //set the IoMode defined by a FloatEncoding tag, if it exists   
    DomOpProject docLearn(root);

    //learnInput
    LearnInput * lInput = new LearnInput();

    docLearn.readXmlFillIn<LearnInput>(lInput);

    //HACK: set nbVariables_binary and nbVariables_gaussian for CompositeParameter, in case of...
    //TODO: refactor...
    if (lInput->getDataType() == HeterogeneousData) {
      Global::nbVariables_binary = lInput->getData()->getBinaryData()->_pbDimension;
      Global::nbVariables_gaussian = lInput->getData()->getGaussianData()->_pbDimension;
    }
    if (lInput->getDataType() == QualitativeData || lInput->getDataType() == HeterogeneousData) {
      Global::vNbFactor.clear();
      for (int i=0; i<lInput->getData()->getBinaryData()->getPbDimension(); i++)
        Global::vNbFactor.push_back(lInput->getData()->getBinaryData()->getTabNbModality()[i]);
    }

    LearnOutput * lOutput = NULL;
    //xmlpp::Element *listOutput = dynamic_cast<xmlpp::Element*>(root->get_first_child("ListOutput"));
    //if (listOutput && !bOnlyInput) {
    if (!bOnlyInput) {
      //LearnOutput
      lOutput = new LearnOutput();//lInput->getCriterionName());
      docLearn.readXmlFillOut<LearnOutput, LearnModelOutput>(lOutput, lInput);
    }

    lInput->finalize();
    return new LearnMain(lInput, lOutput);
    
    //return NULL;
  }
  
  PredictMain * IStream_XML_Predict(const std::string& s, bool bOnlyInput, IoMode iomode) {
	IOMODE = iomode; //TODO...    
    //take the absolute path
    const string str = s;
    ValidateSchema(str, IOStreamXMLFile::Project);
    xmlpp::DomParser parser;
    parser.parse_file(s);
    xmlpp::Document *doc = parser.get_document();
    xmlpp::Element *root = doc->get_root_node();  
    if ( root->get_name() != "Project" ) throw IOStreamErrorType::badIOStreamFormat;
    string xsitype = root->get_attribute_value("type", "xsi");
    if (xsitype.compare("Predict") != 0)  throw IOStreamErrorType::badXML;    
    IoModeManager ioMode = IoModeManager(root);  //set the IoMode defined by a FloatEncoding tag, if it exists   
    DomOpProject docPredict(root);

    //predictInput
    PredictInput * pInput = docPredict.readXmlPredictInput();

    //HACK: set nbVariables_binary and nbVariables_gaussian for CompositeParameter, in case of...
    //TODO: refactor...
    if (pInput->getDataType() == HeterogeneousData) {
      Global::nbVariables_binary = pInput->getData()->getBinaryData()->_pbDimension;
      Global::nbVariables_gaussian = pInput->getData()->getGaussianData()->_pbDimension;
    }
    if (pInput->getDataType() == QualitativeData || pInput->getDataType() == HeterogeneousData) {
      Global::vNbFactor.clear();
      for (int i=0; i<pInput->getData()->getBinaryData()->getPbDimension(); i++)
        Global::vNbFactor.push_back(pInput->getData()->getBinaryData()->getTabNbModality()[i]);
    }

    PredictOutput * pOutput = NULL;
    //xmlpp::Element *listOutput = dynamic_cast<xmlpp::Element*>(root->get_first_child("ListOutput"));
    //if (listOutput && !bOnlyInput) {
    if (!bOnlyInput) {
      //PredictOutput
      pOutput = new PredictOutput();//pInput->getCriterionName());
      docPredict.readXmlFillOut<PredictOutput, PredictModelOutput>(pOutput, pInput);
    }

    pInput->finalize();
    return new PredictMain(pInput, pOutput);
    
    //return NULL;
  }


  
///to validate schema XML. The used .xsd depends on xml file 
  void ValidateSchema(const string & s, const IOStreamXMLFile & xmlFile, bool verbose) {

  bool res = false;
  string schemafile;
  string res_path = XEM_RESOURCES_PATH + string("/");
  switch (xmlFile) {
  case IOStreamXMLFile::Project:
    schemafile = res_path + "project.xsd";
    break;
  case IOStreamXMLFile::Data:
    schemafile = res_path + "data.xsd";
    break;
  case IOStreamXMLFile::Partition:
  case IOStreamXMLFile::Label:    
    schemafile = res_path + "label_or_partition.xsd";
    break;
  case IOStreamXMLFile::Parameter:
    schemafile = res_path + "parameter.xsd";
    break;
  case IOStreamXMLFile::Weights:
    schemafile = res_path + "weights.xsd";
    break;
  default:
    throw IOStreamErrorType::badSchema;
  }
  xmlpp::SchemaValidator validator(schemafile); //to change in XsdValidator later

  try {
    validator.validate(s);
  }
  catch (xmlpp::validity_error & e) {
    if(verbose)
      std::cout<< "file:"<<s<<",schema:"<<schemafile<<","<<e.what()<<std::endl;
    throw IOStreamErrorType::badXML;    
  }

}

//------------------------------------------------
// Read ClusteringInput from a 'flat format' file (txt-ascii)
//------------------------------------------------

// TODO: deprecated function, to be removed.

ClusteringMain * IStream_FLAT(const string & s) {
	ClusteringInput * input = NULL;

	// Open Stream
	ifstream fi(s.c_str(), ios::in);
	if (! fi.is_open())
		THROW(InputException, wrongInputFileName);

	int64_t i;
	bool nbFind = false;

	//--------------------//
	// Compulsory Options //
	//--------------------//
	bool binaryDataType = false;

	// Number of sample(s) //
	//---------------------//
	int64_t nbSample = 0;
	moveUntilReach(fi, "nblines");
	if (!fi.eof()) {
		fi >> nbSample;
		if (nbSample > maxNbSample) {
			THROW(InputException, nbLinesTooLarge);
		}
		else if (nbSample <= 0) {
			THROW(InputException, nbLinesTooSmall);
		}
	}
	else {
		THROW(InputException, errorNbLines);
	}

	// nbVariable //
	//-----------//
	int64_t nbVariable = 0;

	moveUntilReach(fi, "pbdimension");
	if (!fi.eof()) {
		fi >> nbVariable;
		if (nbVariable > maxPbDimension) {
			THROW(InputException, pbDimensionTooLarge);
		}
		else if (nbVariable <= 0) {
			THROW(InputException, pbDimensionTooSmall);
		}
	}
	else {
		THROW(InputException, errorPbDimension);
	}

	// Cluster(s) //
	//------------//

	/* Number of number of cluster(s) */
	int64_t nbNbCluster = 0;
	vector<int64_t> nbCluster;
	moveUntilReach(fi, "nbnbcluster");
	if (!fi.eof()) {
		fi >> nbNbCluster;
		if (nbNbCluster > maxNbNbCluster) {
			THROW(InputException, nbNbClusterTooLarge);
		}
		else if (nbNbCluster <= 0) {
			THROW(InputException, nbNbClusterTooSmall);
		}
		nbCluster.resize(nbNbCluster);
	}
	else {
		THROW(InputException, errorNbNbCluster);
	}

	/* List of number of cluster(s) */
	moveUntilReach(fi, "listnbcluster");
	if (!fi.eof()) {
		int64_t nb;
		for (i = 0; i < nbNbCluster; i++) {
			fi >> nb;
			nbCluster[i] = nb;
		}
	}
	else {
		THROW(InputException, errorListNbCluster);
	}

	// tab modality //
	//--------------//
	moveUntilReach(fi, "nbmodality");
	int64_t * tabNbModality = NULL;
	if (!fi.eof()) {
		binaryDataType = true;
		tabNbModality = new int64_t[nbVariable];
		for (i = 0; i < nbVariable; i++) {
			fi >> tabNbModality[i];
			if (tabNbModality[i] < 2)
				THROW(InputException, errorNbModality);
		}
	}

	// Data Description//
	//-----------------//
	string dataFileName = "";
	moveUntilReach(fi);
	if (fi.eof()) {
		THROW(InputException, errorDataFile);
	}
	fi >> dataFileName;

	vector<ColumnDescription *> columnDescription;
	columnDescription.resize(nbVariable);
	for (int64_t j = 0; j < nbVariable; j++) {
		if (binaryDataType) {
			columnDescription[j] = new QualitativeColumnDescription(j, tabNbModality[j]);
		}
		else {
			columnDescription[j] = new QuantitativeColumnDescription(j);
		}
	}

	DataDescription dataDescription(
			nbSample, nbVariable, columnDescription, FormatNumeric::txt, dataFileName);

	// create input with Compulsory 'inputs'
	//-------------------------------------
	input = new ClusteringInput(nbCluster, dataDescription);

	//----------------
	// Optional Inputs
	//----------------

	// init reading at the beginning of file //
	fi.clear();
	fi.seekg(0, ios::beg);
	nbFind = false;

	// Weight  //
	//---------//
	string weightFileName = "";
	moveUntilReach(fi, "weightfile");
	if (!fi.eof()) {
		fi >> weightFileName;
		Data * data = (input->getDataDescription()).getData();
		data->setWeight(weightFileName);
	}

	// Model type(s) //
	//---------------//
	nbFind = false;
	/* Number of model(s) */
	int64_t nbModelName = 1;
	moveUntilReach(fi, "nbmodel");
	if (!fi.eof()) {
		nbFind = true;
		fi >> nbModelName;

		if (nbModelName > maxNbModel) {
			THROW(InputException, nbModelTypeTooLarge);
		}
		else if (nbModelName <= 0) {
			THROW(InputException, nbModelTypeTooSmall);
		}

		// listmodel
		moveUntilReach(fi, "listmodel");
		if ( (!fi.eof()) && (nbFind) ) {
			ModelName modelName;
			string strModelName;
			for (int64_t i = 0; i < nbModelName; ++i) {
				fi >> strModelName;
				modelName = StringToModelName(strModelName);
				if (isHD(modelName)) {
					THROW(InputException, wrongModelName);
				}
				else {
					if (i == 0) {
  					ModelType temp = ModelType(modelName);
  					input->setModelType(&temp, i);
  				}
  				else {
  					ModelType temp = ModelType(modelName);
  					input->insertModelType(&temp, i);
					}
				}
			}
		}
	}

	// Criterion(ia) //
	//---------------//
	nbFind = false;

	/* Number of criterion(ia) */
	int64_t nbCriterion = 1;
	moveUntilReach(fi, "nbcriterion");
	if (!fi.eof()) {
		nbFind = true;
		fi >> nbCriterion;

		if (nbCriterion > maxNbCriterion) {
			THROW(InputException, nbCriterionTooLarge);
		}
		else if (nbCriterion <= 0) {
			THROW(InputException, nbCriterionTooSmall);
		}

		moveUntilReach(fi, "listcriterion");
		if ( (!fi.eof()) && (nbFind) ) {
			string strCriterionName;
			CriterionName criterionName;
			for (int64_t i = 0; i < nbCriterion; i++) {
				fi >> strCriterionName;
				criterionName = StringtoCriterionName(strCriterionName);
				if (i == 0) {
					input->setCriterion(criterionName, i);
				}
				else {
					input->insertCriterion(criterionName, i);
				}
			}
		}
		else if ( (!fi.eof()) && (!nbFind) ) {
			THROW(InputException, errorNbCriterion);
		}
		else if ( (fi.eof()) && (nbFind) ) {
			THROW(InputException, errorListCriterion);
		}
	}

	// KnownLabel //
	//------------//

	nbFind = false;
	moveUntilReach(fi, "partitionfile");
	if (!fi.eof()) {
		if (nbNbCluster != 1) {
			THROW(InputException, knownPartitionNotAvailable);
		}
		string partitionFileName;
		fi >> partitionFileName;

    xmlpp::DomParser parser;
    parser.parse_file(partitionFileName);
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

    string partitionFilename = elementFilename->get_child_text()->get_content();

    NumericPartitionFile partitionFile;
    partitionFile._fileName = partitionFilename;
    partitionFile._format = format;
    partitionFile._type = type;

		input->insertKnownPartition(partitionFile);
	}

	/* Strategy */
	// nb of strategies
	moveUntilReach(fi, "nbstrategy");
	int64_t nbStrategy = 1;
	if (!fi.eof()) {
		fi >> nbStrategy;
		if (nbStrategy > 1) {
			THROW(InputException, nbStrategyTypeTooLarge);
		}
		else if (nbStrategy <= 0) {
			THROW(InputException, nbStrategyTypeTooSmall);
		}
		ClusteringStrategy * strat = new ClusteringStrategy();
		int64_t * tabNbCluster = new int64_t[nbNbCluster];
		for (int64_t k = 0; k < nbNbCluster; k++) {
			tabNbCluster[k] = nbCluster[k];
		}
		ModelType * modelType = input->getModelType()[0];
		Data * data = (input->getDataDescription()).getData();
		strat->input_FLAT_FORMAT(fi, data, nbNbCluster, tabNbCluster, modelType);
		((ClusteringInput*) input)->setStrategy(strat);
		delete [] tabNbCluster;
	}

	input->finalize();

	fi.close();

	ClusteringMain * temp = new ClusteringMain( input );
	return temp;
}

//------------//
// OSTREAM //
//------------//

void OStream(const string & s, IOStreamFormat format, ClusteringMain * cMain, IoMode iomode) {
	IOMODE = iomode; //TODO...
	//test format
	if (format == IOStreamFormat::XML) {
		// XML format
		//----------

		OStream_XML(s, cMain); 
	}
	else if (format == IOStreamFormat::FLAT) {
		// FLAT format
		//-----------
		OStream_Clustering_FLAT(cMain);

	}
	else {
		throw IOStreamErrorType::badIOWriteFormat;
	}
}
  string getBaseName(const string& s){
    std::size_t i = s.find_last_of("/");
    return s.substr(i+1);
  }
  string getDirName(const string& s){
    std::size_t i = s.find_last_of("/");
    return s.substr(0, i);
  }
  vector<string> getPathElements(const string& s){
    vector<string> res(2);
    std::size_t i = s.find_last_of("/");
    res[0] = s.substr(0, i);
    res[1] = s.substr(i+1);
    return res;
  }
  string normalizeFilename(const string& s){
    if(PATHMODE == PathMode::ABSOLUTE) return s;
    vector<string> vs = getPathElements(s);
    if(vs[0].compare(PROJECT_DIRNAME)!=0) return s;
    return vs[1];
    
  }
  string getAbsolutePath(const string& s){
    if(PATHMODE == PathMode::ABSOLUTE) return s;
    return PROJECT_DIRNAME+"/"+s;
  }
  template<class T>
  void OStream_XML(const string & s, T * cMain, IoMode iomode, PathMode pathMode) {
    IOMODE = iomode; //TODO...
    string dataFilename = string(s);
    
    //fill project with our input and output
    DomOpProject project ;
    PATHMODE = pathMode;
    vector<string> elts = getPathElements(s);
    PROJECT_DIRNAME = elts[0];
    string prefix = (pathMode == PathMode::ABSOLUTE) ? dataFilename : elts[1];
    project.writeMixmodXml(prefix, cMain);
    string filename = dataFilename + ".mixmod";
    removeIfExists(filename);
    //creates the .mixmod with project
    project.write_to_file(filename);
  }
  template void   OStream_XML(const string & s, ClusteringMain* cMain, IoMode iomode, PathMode pathMode);
  template void   OStream_XML(const string & s, LearnMain* cMain, IoMode iomode, PathMode pathMode);
  template void   OStream_XML(const string & s, PredictMain* cMain, IoMode iomode, PathMode pathMode);    
  
void OStream_DiscriminantAnalysis_XML(const string & s, ClusteringMain * cMain) {
}

//------------------------------------------------------------
// write output Data in FLAT format (txt) for clustering study
//------------------------------------------------------------
/*
Note : this function works only with a ClusteringMain which has been computed (not with one created from a XML file

The following files are created :
- complete.txt
- XposteriorProbabilities.txt (where X=BIC, ICL, NEC)
- Xlabel.txt
- Xlikelihood.txt (where X=BIC, ICL, NEC)
- XError.txt (where X=BIC, ICL, NEC)
- errorMixmod.txt
- errorModel.txt
 */

// TODO: deprecated function, to be removed.

void OStream_Clustering_FLAT(ClusteringMain * cMain) {

	int64_t iEst, iCrit;

	ClusteringOutput * cOutput =  cMain->getOutput();
	ClusteringModelOutput *  currentCMOutput = NULL;
	//Model * currentEstimation = NULL;

	// Test if CMain has Estimations (has been computed)
	currentCMOutput = cOutput->getClusteringModelOutput(0);
	try {
		//    currentEstimation = currentCMOutput->getModel();
	}
	catch (Exception & e) {
		THROW(OtherException, nonImplementedMethod);
	}

	//-------------
	// complete.txt
	//-------------
	ofstream completeStream("complete.txt", ios::out);
	completeStream.setf(ios::fixed, ios::floatfield);
	if (! completeStream.is_open()) {
		THROW(InputException, errorOpenFile);
	}

	completeStream << "---------------------------------\n";
	completeStream << "     MIXMOD  Complete Output     \n";
	completeStream << "---------------------------------\n";
	completeStream << "\n";
	completeStream << "Number of samples : " 
			<< currentCMOutput->getLabelDescription()->getNbSample() << endl << endl;
	for (iEst = 0; iEst < cOutput->getNbClusteringModelOutput(); iEst++) {
		currentCMOutput = cOutput->getClusteringModelOutput(iEst) ;
		//    currentEstimation = currentCMOutput->getModel();
		//dynamic_cast<ClusteringStrategy*>(cMain->getInput()->getStrategy())->edit(completeStream);
		completeStream << "\n";
		completeStream << "\t\tNumber of Clusters : " << currentCMOutput->getNbCluster() << endl;
		completeStream << "\t\t------------------" << endl << endl;
		currentCMOutput->getModelType().edit(completeStream);
		if (currentCMOutput->getStrategyRunError() == NOERROR) {
			for (iCrit = 0; iCrit < cMain->getInput()->getCriterionName().size() ; iCrit++) {
				currentCMOutput->getCriterionOutput(cMain->getInput()
						->getCriterionName(iCrit)).editTypeAndValue(completeStream);
			}
			// TODO change to currentCMOutput->gettParameterDescription()->edit();
			currentCMOutput->getParameterDescription()->getParameter()->edit(completeStream, true);
		}
		else {
			completeStream << "\t\t\tError" << endl << endl;
		}
	}
	completeStream.close();

	//-------------
	// numericComplete.txt
	//-------------
	ofstream numericCompleteStream("numericComplete.txt", ios::out);
	numericCompleteStream.setf(ios::fixed, ios::floatfield);
	if (! numericCompleteStream.is_open()) {
		THROW(InputException, errorOpenFile);
	}

	for (iEst = 0; iEst < cOutput->getNbClusteringModelOutput(); iEst++) {
		currentCMOutput = cOutput->getClusteringModelOutput(iEst) ;
		//    currentEstimation = currentCMOutput->getModel();

		if (currentCMOutput->getStrategyRunError() == NOERROR) {
			for (iCrit = 0; iCrit < cMain->getInput()->getCriterionName().size() ; iCrit++) {
				currentCMOutput->getCriterionOutput(cMain->getInput()
						->getCriterionName(iCrit)).editValue(numericCompleteStream);
			}
			// TODO change to currentCMOutput->gettParameterDescription()->edit();
			currentCMOutput->getParameterDescription()->getParameter()
					->edit(numericCompleteStream, false);
		}
	}
	numericCompleteStream.close();

	//ClusteringOutput originalOutput(*cOutput);

	/* ------------------------------------------------
  - XposteriorProbabilities.txt (where X=BIC, ICL, NEC)
  - XnumericLikelihood.txt (where X=BIC, ICL, NEC)
  - XError.txt (where X=BIC, ICL, NEC)
  ------------------------------------------------------*/

	// get all criterionName
	ClusteringInput * cInput = dynamic_cast<ClusteringInput*> (cMain->getInput());
	vector<CriterionName> criterionNameInInput = cInput->getCriterionName();

	vector<CriterionName> criterionName;
	criterionName.push_back(BIC);
	criterionName.push_back(ICL);
	criterionName.push_back(NEC);

	CriterionName iCritInCriterionNameInInput;
	// loop over criterion name
	for (iCrit = 0; iCrit < criterionName.size(); iCrit++) {
		// Is this criterion used in Input ?
		bool used = false;
		for ( int i = 0; i < criterionNameInInput.size(); i++) {
			if (criterionNameInInput[i] == criterionName[iCrit]) {
				used = true;
				iCritInCriterionNameInInput = criterionNameInInput[i];
			}
		}

		//-----------------------------------
		// create XposteriorProbabilities.txt
		//-----------------------------------
		string criterionString = CriterionNameToString(criterionName[iCrit]);
		string fileName = criterionString + "posteriorProbabilities.txt";
		ofstream postProbaStream(fileName.c_str(), ios::out);
		postProbaStream.setf(ios::fixed, ios::floatfield);
		if (! postProbaStream.is_open()) {
			THROW(InputException, errorOpenFile);
		}

		if (used) {
			cOutput->sort(criterionName[iCrit]);
			currentCMOutput = cOutput->getClusteringModelOutput(0); // the best
			if (currentCMOutput->getStrategyRunError() == NOERROR) {
				ProbaDescription * probaDescription = currentCMOutput->getProbaDescription();
				if (probaDescription) {
					probaDescription->getProba()->edit(postProbaStream);
				}
			}
			else {
				// error in strategyRun
				//TODO
			}

		}
		postProbaStream.close();

		//-----------------------------------
		// create Xlabel.txt
		//-----------------------------------
		fileName = criterionString + "label.txt";
		ofstream labelStream(fileName.c_str(), ios::out);
		labelStream.setf(ios::fixed, ios::floatfield);
		if (! labelStream.is_open()) {
			THROW(InputException, errorOpenFile);
		}

		if (used) {
			if (currentCMOutput->getStrategyRunError() == NOERROR) {
				LabelDescription * labelDescription = currentCMOutput->getLabelDescription();
				if (labelDescription) {
					const_cast<Label *> (labelDescription->getLabel())->edit(labelStream);
				}
			}
			else {
				// error in strategyRun
				//TODO
			}

		}
		labelStream.close();

		//-----------------------
		// create XnumericLikelihood.txt
		//-----------------------
		fileName = criterionString + "numericLikelihood.txt";
		ofstream LLStream(fileName.c_str(), ios::out);
		LLStream.setf(ios::fixed, ios::floatfield);
		if (! LLStream.is_open()) {
			THROW(InputException, errorOpenFile);
		}

		if (used) {
			if (currentCMOutput->getStrategyRunError() == NOERROR) {

				//Model * currentModel = currentCMOutput->getModel();

				LLStream << currentCMOutput->getParameterDescription()->getParameter()
						->getFreeParameter() << endl;
				LLStream << currentCMOutput->getParameterDescription()->getParameter()
						// false : to not compute fik because already done
						->getModel()->getLogLikelihood(false) << endl;
				LLStream << currentCMOutput->getParameterDescription()->getParameter()
						->getModel()->getCompletedLogLikelihood() << endl;
				LLStream << currentCMOutput->getParameterDescription()->getParameter()
						->getModel()->getEntropy() << endl;
			}
			else {
				// error in strategyRun
				//TODO
			}

		}
		postProbaStream.close();

		//------------------
		// create XError.txt
		//------------------
		fileName = criterionString + "Error.txt";
		ofstream errorStream(fileName.c_str(), ios::out);
		errorStream.setf(ios::fixed, ios::floatfield);
		if (! errorStream.is_open()) {
			THROW(InputException, errorOpenFile);
		}

		if (used) {
			CriterionOutput currentCriterionOutput;
			for (iEst = 0; iEst < cOutput->getNbClusteringModelOutput(); iEst++) {
				currentCMOutput = cOutput->getClusteringModelOutput(iEst) ;
				//        currentEstimation = currentCMOutput->getModel();
				if (currentCMOutput->getStrategyRunError() == NOERROR) {
					currentCriterionOutput = 
							currentCMOutput->getCriterionOutput(iCritInCriterionNameInInput);
					errorStream << currentCriterionOutput.getError().what() << endl;
				}
				else {
					errorStream << currentCMOutput->getStrategyRunError().what() << endl;
				}
			}
		}
		errorStream.close();
	} // end if iCrit

	//----------------------
	// create errorModel.txt
	//----------------------
	ofstream errorModelStream("errorModel.txt", ios::out);
	errorModelStream.setf(ios::fixed, ios::floatfield);
	if (! errorModelStream.is_open()) {
		THROW(InputException, errorOpenFile);
	}
	for (iEst = 0; iEst < cOutput->getNbClusteringModelOutput(); iEst++) {
		currentCMOutput = cOutput->getClusteringModelOutput(iEst);
		//    currentEstimation = currentCMOutput->getModel();
		errorModelStream << currentCMOutput->getStrategyRunError().what() << endl;
	}
	errorModelStream.close();

	/*
	 //-----------------------
	 // create errorMixmod.txt
	 //-----------------------
	 ofstream errorMixmodStream("errorMixmod.txt", ios::out);
	 errorMixmodStream.setf(ios::fixed, ios::floatfield);
	 if (! errorMixmodStream.is_open()){
	   throw errorOpenFile;
	 }
	 errorMixmodStream<<cOutput->getError()<<endl;
	 errorMixmodStream.close();
	 */
}
  // Tools for floats
  double custom_stod(string s){
    //cout<<"custom_stod:"<<s<<":"<<endl;
    //if (IOMODE != IoMode::BINARY || s.find('.')!=string::npos) {
    if (IOMODE != IoMode::BINARY) {      
      //cout<<"custom_stod_HR:"<<s<<":"<<endl;
      return std::stod(s);
    }
    double result;      
    stringstream stream;
    uint64_t tmp;
    stream << hex << s; 
    stream >> tmp;
    memcpy(&result, &tmp, sizeof(tmp));
    return result;
  }
  double std_stod(string s){
    //cout<<"custom_stod:"<<s<<":"<<endl;
    return std::stod(s);
  }
  string custom_dtos(double d){
    if (IOMODE != IoMode::BINARY) {
      return std::to_string(d);
    }
    stringstream stream;
    int64_t tmp;
    double tmp_value = d;
    memcpy(&tmp, &tmp_value, sizeof(tmp_value));
    stream << hex << tmp;
    return stream.str();
  }
//VariableTypeToString
string ColumnTypeToString(const IOStreamColumnType & columnType) {
	string res;
	switch (columnType) {
	case IOStreamColumnType::Quantitative:
		res = "Quantitative";
		break;
	case IOStreamColumnType::Qualitative:
		res = "Qualitative";
		break;
	case IOStreamColumnType::Individual:
		res = "Individual";
		break;
	case IOStreamColumnType::Weight:
		res = "Weight";
		break;
	}
	return res;
}

//StringToVariableType
IOStreamColumnType StringToColumnType(const string & strVariableType) {
	IOStreamColumnType res;
	if (strVariableType.compare("Quantitative") == 0)
		res = IOStreamColumnType::Quantitative;
	if (strVariableType.compare("Qualitative") == 0)
		res = IOStreamColumnType::Qualitative;
	if (strVariableType.compare("Individual") == 0)
		res = IOStreamColumnType::Individual;
	if (strVariableType.compare("Weight") == 0)
		res = IOStreamColumnType::Weight;
	if (strVariableType.compare("Unused") == 0)
			res = IOStreamColumnType::Unused;
	return res;
}

//-------------------------------------
// CreateMixmodDataFileFromUserDataFile
// Example
// - userDataFileName="/home/user/foo.dat"
// - mixmodDataFileName="/home/user/mixmod/mixmod_foo.dat"
//-------------------------------------

void createMixmodDataFileFromUserDataFile(string userDataFileName, string mixmodDataFileName, 
		vector<bool> & columnUsed, FormatNumeric::FormatNumericFile format, 
		bool individualNameMustBeGenerated, int64_t & nbSample) 
{
	int64_t nbSampleOut = 0;

	//------
	// verif
	//------
	// file format
	// 	if (format==FormatNumeric::txt){
	// 		string st="file "+userDataFileName+" >res.txt";
	// 		int res = system(st.c_str());
	// 		ifstream formatIstream("res.txt", ios::in);
	// 		string contenu;
	// 		getline(formatIstream, contenu);
	//     if (!strstr(contenu.c_str(), "ASCII")){
	//       cout<<"not ascii"<<endl;
	//     }
	//     else{
	//       cout<<"ascii"<<endl;
	//     }
	// 		formatIstream.close();
	// 	}

	// column
	if (columnUsed.size() < 1) {
		throw IOStreamErrorType::badColumnUsedInCreateMixmodDataFileFromUserDataFile;
	}

	int64_t i;
	int64_t cptTrue = 0;
	for (i = 0; i < columnUsed.size(); i++) {
		if (columnUsed[i]) {
			cptTrue++;
		}
	}
	if (cptTrue == 0) {
		throw IOStreamErrorType::badColumnUsedInCreateMixmodDataFileFromUserDataFile;
	}

	ifstream dataIstream((userDataFileName).c_str(), ios::in);
	if (!dataIstream) {
		throw IOStreamErrorType::fileDontExist;
	}
	if (! dataIstream.is_open()) {
		throw IOStreamErrorType::fileCantBeOpened;
	}

	//test if file exists
	fstream fs(mixmodDataFileName.c_str(), ios_base::in); // attempt open for read
	if (fs) {//ok, file exists. close and reopen in write mode
		fs.close();
		throw IOStreamErrorType::fileAlreadyExist;
	}

	//  int64_t columnUsedSize = columnUsed.size();
	//cout<<"columnUsedSize="<<columnUsedSize<<endl;
	//cout<<"cptTrue="<<cptTrue<<endl;

	ofstream dataOstream((mixmodDataFileName).c_str(), ios::out);
	//------------
	//Read/write
	//------------
	if (format == FormatNumeric::txt) {
		vector<string> strings;
		string tmp;
		//strings.reserve(columnUsedSize);
		bool end = false;

		do {
			int64_t cptIstream = 0;
			int64_t cptOstream = 0;
			// can we read columnUsedSize string ?
			while (!dataIstream.eof() && cptOstream < cptTrue) {
				//cout<<endl<<"Entree dans le while, cptIstream="<<cptIstream<<endl;
				//cout<<"Entree dans le while, cptOstream="<<cptOstream<<endl;
				dataIstream >> tmp;
				//cout<<"lu :"<<tmp<<endl;
				if (columnUsed[cptIstream]) {
					//cout<<"mis dans strings - size :"<<strings.size()<<endl;
					strings.push_back(tmp);
					cptOstream++;
				}
				//else{cout<<"pas mis dans strings"<<endl;}
				cptIstream++;
				char c = dataIstream.peek();
				if (c == '\r') {
					//cout<<"rr, on avance"<<endl;
					c = dataIstream.get();
				}
			}
			if (cptOstream < cptTrue) {
				end = true;
			}
			else {
				// on verifie la taille de strings
				if (strings.size() == cptTrue) {
					//cout<<"write in oStream"<<endl;
					// write in dataOstream
					// 1rst : individualNameMustBeGenerated ?
					if (individualNameMustBeGenerated) {
						dataOstream << "Individual_" << nbSampleOut + 1 << "\t";
					}
					for (i = 0; i < cptTrue; i++) {
						dataOstream << strings[0] << "\t";
						strings.erase(strings.begin());
					}
					//cout<<"apres avoir depil�, size:"<<strings.size()<<endl;
					nbSampleOut++;
					dataOstream << "\n";
				}
				else {
					//cout<<"erreur size"<<endl;
					throw IOStreamErrorType::notEnoughValuesInDataFile;
				}
			}

		}
		while (!dataIstream.eof() && !end);

		dataOstream << endl;

		if (nbSample != 0 && nbSampleOut < nbSample) {
			throw IOStreamErrorType::notEnoughValuesInDataFile;
		}
		if (nbSample == 0) {
			nbSample = nbSampleOut;
		}
	}
	else {
		throw IOStreamErrorType::badNumericFormat;
	}

	//close Stream
	dataOstream.close();
	dataIstream.close();
}


int Global::nbVariables_binary = 0;
int Global::nbVariables_gaussian = 0;
std::vector<int64_t> Global::vNbFactor;

  //CPOLI
  void removeIfExists(const std::string& filename){
    struct stat st_;
    if (stat(filename.c_str(), &st_) != -1){
      remove(filename.c_str());
    }


  }
  //CPOLI
  xmlpp::Element *get_first_child_element(xmlpp::Node *parent){
      auto children = parent->get_children();
      for (auto it=children.begin(); it != children.end(); ++it){
        xmlpp::Element* elt = dynamic_cast<xmlpp::Element*>(*it);
        if(elt) return elt;
      }
      return NULL;
  }
  /*  
  void copyFile(const std::string& src, const std::string& dest){
    struct stat stat_src;
    off_t offset = 0;
    int src_fd = open (src.c_str(), O_RDONLY);
    fstat (src_fd, &stat_src);
    int dest_fd = open (dest.c_str(), O_WRONLY | O_CREAT, stat_src.st_mode);
    sendfile (dest_fd, src_fd, &offset, stat_src.st_size);
    close (dest_fd);
    close (src_fd);
    }*/
  
  ProjectType getProjectType(string filename){
    ValidateSchema(filename, IOStreamXMLFile::Project);
    xmlpp::DomParser parser;
    parser.parse_file(filename);
    xmlpp::Document *doc = parser.get_document();
    xmlpp::Element *root = doc->get_root_node();  
    if ( root->get_name() != "Project" ) throw IOStreamErrorType::badIOStreamFormat;
    string xsitype = root->get_attribute_value("type", "xsi");
    if (xsitype.compare("Clustering") == 0) return ProjectType::Clustering;
    if (xsitype.compare("Learn") == 0) return ProjectType::Learn;
    if (xsitype.compare("Predict") == 0) return ProjectType::Predict;
    return ProjectType::Unknown;            
  }

  IoModeManager::IoModeManager(xmlpp::Element* root){
    //xmlpp::Element* elt = dynamic_cast<xmlpp::Element*>(root->get_first_child("FloatEncoding"));
    /*if(elt){
      previous_ = IOMODE;      
      string encoding = elt->get_child_text()->get_content();
      if(encoding.compare("HexBinary")==0){
        IOMODE = IoMode::BINARY;
      } else {
        IOMODE = IoMode::NUMERIC;
      }
      } */
  }

  
  IoModeManager::IoModeManager(IoMode mode){
    //previous_ = IOMODE;
    //IOMODE = mode;
  }
  IoModeManager::~IoModeManager(){
    //IOMODE = previous_;
  }

} //end namespace


