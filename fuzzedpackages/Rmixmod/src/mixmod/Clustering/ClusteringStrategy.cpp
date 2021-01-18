/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringStrategy.cpp  description
    copyright            : (C) MIXMOD Team - 2001-2016
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

#include "mixmod/Clustering/ClusteringStrategy.h"
#include "mixmod/Clustering/ClusteringStrategyInit.h"
#include "mixmod/Kernel/Algo/EMAlgo.h"
#include "mixmod/Kernel/Algo/CEMAlgo.h"
#include "mixmod/Kernel/Algo/SEMAlgo.h"
#include "mixmod/Kernel/Parameter/GaussianParameter.h"
#include "mixmod/Kernel/Parameter/GaussianSphericalParameter.h"
#include "mixmod/Kernel/Parameter/GaussianDiagParameter.h"
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Kernel/Parameter/GaussianHDDAParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkjhParameter.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/IO/Partition.h"
#include "mixmod/Kernel/Model/Model.h"

namespace XEM {

//-----------
//Constructor
//-----------
ClusteringStrategy::ClusteringStrategy() {
	_nbTry = defaultNbTryInStrategy;
	_strategyInit = new ClusteringStrategyInit();
	_nbAlgo  = defaultNbAlgo;
	_tabAlgo.reserve(_nbAlgo);
	//  _tabAlgo = new XEMAlgo * [_nbAlgo];
	for (int64_t i = 0; i < _nbAlgo; i++) {
		_tabAlgo.push_back(createDefaultClusteringAlgo());
	}
}

ClusteringStrategy::ClusteringStrategy(const ClusteringStrategy & strategy) {
	_nbTry = strategy.getNbTry();
	_strategyInit = new ClusteringStrategyInit(*(strategy.getStrategyInit()));
	_nbAlgo = strategy.getNbAlgo();
	//  _tabAlgo = new XEMAlgo * [_nbAlgo];
	//  XEMAlgo ** tabA = strategy.getTabAlgo();
	std::vector<Algo*> tabA = strategy.getTabAlgo();
	for (int64_t i = 0; i < _nbAlgo; i++) {
		//_tabAlgo[i] = tabA[i]->clone();
		_tabAlgo.push_back(tabA[i]->clone());
	}
}

ClusteringStrategy* ClusteringStrategy::clone() {
	return new ClusteringStrategy(*this);
}

//----------
//Destructor
//----------
ClusteringStrategy::~ClusteringStrategy() {

	for (unsigned int i = 0; i < _tabAlgo.size(); i++) {
		delete _tabAlgo[i];
	}
	if (_strategyInit) delete _strategyInit;
}

// setAlgoEpsilon
void ClusteringStrategy::setAlgoEpsilon(int64_t position, double epsilonValue) {
	_tabAlgo[position]->setEpsilon(epsilonValue);
}

// setAlgoStopRuleTypeValue
void ClusteringStrategy::setAlgoStopRule(AlgoStopName stopName, int64_t position) {
	_tabAlgo[position]->setAlgoStopName(stopName);
}

void ClusteringStrategy::setAlgoIteration( int64_t position, int64_t nbIterationValue) {
	_tabAlgo[position]->setNbIteration(nbIterationValue);
}

// setAlgo
void ClusteringStrategy::setAlgo(AlgoName algoName, int64_t position) {
	if (_tabAlgo[position] != NULL) {
		delete _tabAlgo[position];
	}
	switch (algoName) {
	case EM:
		_tabAlgo[position] = new EMAlgo();
		break;
	case CEM:
		_tabAlgo[position] = new CEMAlgo();
		break;
	case SEM:
		_tabAlgo[position] = new SEMAlgo();
		break;
	default:
		THROW(OtherException, internalMixmodError);
	}
}

// setAlgo
void ClusteringStrategy::addAlgo(AlgoName algoName) {

	switch (algoName) {
	case EM:
		_tabAlgo.push_back(new EMAlgo());
		break;
	case CEM:
		_tabAlgo.push_back(new CEMAlgo());
		break;
	case SEM:
		_tabAlgo.push_back(new SEMAlgo());
		break;
	default:
		THROW(OtherException, internalMixmodError);
	}
	_nbAlgo++;
}

// set init parameter
void ClusteringStrategy::setInitParam(std::string & paramFileName, int64_t position) {
	_strategyInit->setInitParam(paramFileName, position);
}

void ClusteringStrategy::setTabInitParameter(
		Parameter ** tabInitParameter, int64_t nbInitParameter)
{
	_strategyInit->setTabInitParameter(tabInitParameter, nbInitParameter);
}

// set init partition
void ClusteringStrategy::setInitPartition(std::string & partitionFileName, int64_t position) {
	_strategyInit->setPartition(partitionFileName, position);
}

// set init partition
void ClusteringStrategy::setInitPartition(Partition * part, int64_t position) {
	_strategyInit->setPartition(part, position);
}

void ClusteringStrategy::setTabPartition(Partition ** tabPartition, int64_t nbPartition) {
	_strategyInit->setTabPartition(tabPartition, nbPartition);
}

// insert algo
void ClusteringStrategy::insertAlgo(AlgoName algoName, int64_t position) {
	switch (algoName) {
	case EM:
		_tabAlgo.insert(_tabAlgo.begin() + position, new EMAlgo());
		break;
	case CEM:
		_tabAlgo.insert(_tabAlgo.begin() + position, new CEMAlgo());
		break;
	case SEM:
		_tabAlgo.insert(_tabAlgo.begin() + position, new SEMAlgo());
		break;
	default:
		THROW(OtherException, internalMixmodError);
	}
	_nbAlgo++;
}

// removeAlgo
void ClusteringStrategy::removeAlgo( unsigned int position) {
	if ( position < _tabAlgo.size() ) {
      if(_tabAlgo[position]) delete _tabAlgo[position];
      _tabAlgo.erase(_tabAlgo.begin() + position);
      _nbAlgo--;
	}
}

//--------------------
// setStrategyInit
//--------------------
void ClusteringStrategy::setStrategyInit(ClusteringStrategyInit * iStrategyInit) {
	_strategyInit = iStrategyInit; //copy constructor
}

void ClusteringStrategy::setStrategyInit(StrategyInitName initName, Data *& data,
		int64_t nbNbCluster, int64_t * tabNbCluster, ModelType * modelType)
{
	int64_t nbSample    = data->_nbSample;
	int64_t pbDimension  = data->_pbDimension;

	// TODO [bauder]: Why do we initialize with an empty filename ?
	//                Maybe add a function parameter ?
	std::string fileName = "";

	Parameter ** tabInitParameter = NULL;
	Partition ** tabInitPartition = NULL;

	switch (initName) {
	case RANDOM:
	case CEM_INIT:
	case SEM_MAX:
	case SMALL_EM:
		_strategyInit->setStrategyInitName(initName);
		break;

	case  USER:
		_strategyInit->setStrategyInitName(initName);
		tabInitParameter = new Parameter * [nbNbCluster];

		for (int64_t k = 0; k < nbNbCluster; k++) {
			if (isEDDA(modelType->_nameModel)) {
				tabInitParameter[k] = new GaussianGeneralParameter(
						tabNbCluster[k], pbDimension, modelType, fileName);
			}
			else if (getModelGenre(modelType->_nameModel) == QualitativeModel) {
				int64_t * tabNbModality = (data->getBinaryData())->getTabNbModality();
				tabInitParameter[k] = new BinaryEkjhParameter(
						tabNbCluster[k], pbDimension, modelType, tabNbModality, fileName);
			}
			else if (isHD(modelType->_nameModel)) {
				tabInitParameter[k] = new GaussianHDDAParameter(
						tabNbCluster[k], pbDimension, modelType, fileName);
			}
			else THROW(OtherException, internalMixmodError);
		}
		_strategyInit->setTabInitParameter(tabInitParameter, nbNbCluster);
		break;

	case USER_PARTITION:
		_strategyInit->setStrategyInitName(initName);
		tabInitPartition = new Partition * [nbNbCluster];
		for (int64_t k = 0; k < nbNbCluster; k++) {
			NumericPartitionFile partitionFile;
			partitionFile._fileName = fileName;
			partitionFile._format = FormatNumeric::defaultFormatNumericFile;
			partitionFile._type = TypePartition::defaultTypePartition;
			tabInitPartition[k] = new Partition(nbSample, tabNbCluster[k], partitionFile);
		}
		_strategyInit->setTabPartition(tabInitPartition, nbNbCluster);
		break;
	}
}

//-------------
// setNbTry
//-------------
void ClusteringStrategy::setNbTry(int64_t nbTry) {
	if ((_strategyInit->getStrategyInitName() == USER) ||
			(_strategyInit->getStrategyInitName() == USER_PARTITION)) {
		THROW(InputException, badSetNbTry);
	}
	if (nbTry < minNbTryInStrategy) {
		THROW(InputException, nbTryInStrategyTooSmall);
	}
	else
		if (nbTry > maxNbTryInStrategy) {
		THROW(InputException, nbTryInStrategyTooLarge);
	}
	else {
		_nbTry = nbTry;
	}
}

//---
//run
//---
void ClusteringStrategy::run(Model *& model) const {
	//cout<<"XEMClusteringStrategy Init, nbTry="<<_nbTry<<endl;
	if (_nbTry == 1) {
      oneTry(model, true);
	}
	else {
		Model * bestModel = model->clone();
	  oneTry(bestModel);

		// first tries : try until one model has "noerror"
		int iTry = 1;
		while ( !((bestModel->getErrorType()) == NOERROR) && ( iTry < _nbTry ) ) {
			delete bestModel;
			bestModel = model->clone();
			oneTry(bestModel);
			++iTry;
		}

		if ( (bestModel->getErrorType()) == NOERROR ) {
			double bestLLorCLL = bestModel->getCompletedLogLikelihoodOrLogLikelihood();
			// other tries
			for (int64_t i = iTry; i < _nbTry; i++) {
				Model * currentModel = model->clone();
				oneTry(currentModel);
				if ( (currentModel->getErrorType()) == NOERROR ) {
					double lastLLorCLL = currentModel->getCompletedLogLikelihoodOrLogLikelihood();
					if (lastLLorCLL > bestLLorCLL) {
						delete bestModel;
						bestModel = currentModel->clone();
						bestLLorCLL = currentModel->getCompletedLogLikelihoodOrLogLikelihood();
					}
				}
				delete currentModel;
			}
		}
		else{//there is no try with NOERROR
      THROW(OtherException, AllTriesGotErros);
		}
		delete model;
		model = bestModel;
		/*    oneTry(model);
			XEMModel * bestModel = new XEMModel(model);
			double bestLLorCLL = model->getCompletedLogLikelihoodOrLogLikelihood();
			// others tries
			for (int64_t i=1; i<_nbTry; i++){
			  oneTry(model);
			  double lastLLorCLL = model->getCompletedLogLikelihoodOrLogLikelihood();
			  cout<<"Try n°="<< i << " lastLL=" << lastLLorCLL << " bestLL=" << bestLLorCLL
		      << " nbCluster=" << model->getNbCluster() << " modelName="
		      << XEMModelNameToString(model->getModelType()->getModelName()) << endl;
			  if (lastLLorCLL > bestLLorCLL){
				delete bestModel;
				bestModel = new XEMModel(model);
				bestLLorCLL = model->getCompletedLogLikelihoodOrLogLikelihood();
			  }
			}
			model = bestModel;
			delete bestModel;
		 */
	}
}

//------
//oneTry
//------
  void ClusteringStrategy::oneTry(Model *& model, bool doThrow) const {
  try{
	//1 init model
	//------------
	switch (_strategyInit->getStrategyInitName()) {

	case RANDOM:
		model->initRANDOM(_strategyInit->getNbTry());
		break;

	case USER:
	{
		// get initPartition
		int64_t nbCluster = model->getNbCluster();
		int64_t index = 0;
		bool ok = false;
		int64_t nbInitParameter = _strategyInit->getNbInitParameter();
		while (ok == false && index < nbInitParameter) {
			int64_t nbClusterOfInitParameter =
					_strategyInit->getInitParameter(index)->getNbCluster();
			if (nbCluster == nbClusterOfInitParameter) {
				ok = true;
			}
			else {
				index++;
			}
		}
		if (!ok)
			THROW(OtherException, internalMixmodError);
		Parameter * initParameter = _strategyInit->getInitParameter(index);
		model->initUSER(initParameter);
	}
		break;

	case USER_PARTITION:
	{
		// get initPartition
		int64_t nbCluster = model->getNbCluster();
		int64_t index = 0;
		bool ok = false;
		int64_t nbPartition = _strategyInit->getNbPartition();
		while (ok == false && index < nbPartition) {
			int64_t nbClusterOfInitPartition = _strategyInit->getPartition(index)->_nbCluster;
			if (nbCluster == nbClusterOfInitPartition) {
				ok = true;
			}
			else {
				index++;
			}
		}
		if (!ok)
			THROW(OtherException, internalMixmodError);
		Partition * initPartition = _strategyInit->getPartition(index);
		int64_t nbTyInInit = _strategyInit->getNbTry();
		model->initUSER_PARTITION(initPartition, nbTyInInit);
	}
		break;

	case SMALL_EM:
		_strategyInit->initSMALL_EM(model);
		break;

	case CEM_INIT:
		_strategyInit->initCEM_INIT(model);
		break;

	case SEM_MAX:
		_strategyInit->initSEM_MAX(model);
		break;

	default:
		THROW(InputException, wrongStrategyInitName);
	}


	// 2. Algo(s)
	//-----------
	model->setAlgoName(UNKNOWN_ALGO_NAME);

	// model->getParameter()->edit();

	if (DEBUG > 0) {
		cout << "After initialization :" << endl;
		model->editDebugInformation();
	}

	// runs algos
	// define a number of errors iterator
	int nbErrorInAlgo = 0;
	//  _tabAlgo[0]->run(model);
	for (int64_t i = 0; i < _nbAlgo ; i++) {
		try {
			_tabAlgo[i]->run(model);
		}
		catch (Exception&errorType) {

			if (VERBOSE == 1)
				errorType.run();

			++nbErrorInAlgo;
			// set error for that model
			if ( nbErrorInAlgo == _nbAlgo ) {
				model->setError(errorType);
				throw; //TODO FL
			}
		}
	}
  }//fin try
  //s'il y a eu un throw, normalement le model n'est pas noerror
  // il ne faut pas que le throw aille plus haut car il y a éventuellement plusieurs nbTry
  catch(Exception & error){
    //cout<<"erreur dans un oneRun, on continue"<<endl;
      //Nothing to do
    // except when nbTry==1
    if(doThrow) THROW(OtherException, AllTriesGotErros); 
  }
}

//-------
// verify
//-------
bool ClusteringStrategy::verify() {

	bool res = true;

	// Test
	//-----
	// nbAlogType > 0
	if (_nbAlgo < 1 || _tabAlgo.empty()) {
		res = false;
		THROW(InputException, nbAlgoTypeTooSmall);
	}
	if (_nbTry < minNbTryInStrategy) {
		res = false;
		THROW(InputException, nbTryInStrategyTooSmall);
	}
	if (_nbTry > maxNbTryInStrategy) {
		res = false;
		THROW(InputException, nbTryInStrategyTooLarge);
	}

	if (res) {
		res = _strategyInit->verify();
	}

	return res;
}

//---------------
// Input strategy
// TODO XEMInput : a enlever
//---------------
void ClusteringStrategy::input_FLAT_FORMAT(std::ifstream & fi, Data *& data,
		int64_t nbNbCluster, int64_t * tabNbCluster, ModelType * modelType)
{
	int64_t j;
	std::string keyWord = "";
	bool alreadyRead = false;
	std::string a = "";

	// nbTry
	//------
	fi >> keyWord;
	ConvertBigtoLowString(keyWord);
	if (keyWord.compare("nbtry") == 0) {
		int64_t nbTry;
		fi >> nbTry;
		setNbTry(nbTry);
	}

	// StrategyInit
	//-----------------
	_strategyInit->input(fi, data, nbNbCluster, tabNbCluster, modelType, alreadyRead);

	// Algos
	//-------
	/* Number of algorithms */
	moveUntilReach(fi, "nbAlgorithm");
	if (!fi.eof()) {
		for (j = 0; j < _nbAlgo; j++) {
			delete _tabAlgo[j];
		}
		//    delete[] _tabAlgo;

		fi >> _nbAlgo;
		if (_nbAlgo > maxNbAlgo) {
			THROW(InputException, nbAlgoTooLarge);
		}
		else if (_nbAlgo <= 0) {
			THROW(InputException, nbAlgoTooSmall);
		}

		//    _tabAlgo = new XEMAlgo * [_nbAlgo];
		_tabAlgo.resize(_nbAlgo);

		for (j = 0; j < _nbAlgo; j++) {
			fi >> keyWord;
			ConvertBigtoLowString(keyWord);

			if (keyWord.compare("algorithm") == 0) {
				// tabAlgoType._type
				fi >> a;

				if (a.compare("CEM") == 0) {
					_tabAlgo[j] = new CEMAlgo();
				}
				else if (a.compare("EM") == 0) {
					_tabAlgo[j] = new EMAlgo();
				}
				else if (a.compare("SEM") == 0) {
					_tabAlgo[j] = new SEMAlgo();
				}
				else {
					THROW(InputException, wrongAlgoType);
				}

				fi >> keyWord;
				ConvertBigtoLowString(keyWord);

				if (keyWord.compare("stoprule") == 0) {
					fi >> a;

					if (a.compare("NBITERATION") == 0) {
						_tabAlgo[j]->setAlgoStopName(NBITERATION);
					}

					else if (a.compare("EPSILON") == 0) {
						_tabAlgo[j]->setAlgoStopName(EPSILON);
					}

					else if (a.compare("NBITERATION_EPSILON") == 0) {
						_tabAlgo[j]->setAlgoStopName(NBITERATION_EPSILON);
					}
					else {
						THROW(InputException, wrongAlgoStopName);
					}
					// nbIteration && epsilon
					fi >> keyWord;
					ConvertBigtoLowString(keyWord);

					if (keyWord.compare("stoprulevalue") == 0) {
						if (_tabAlgo[j]->getAlgoStopName() == NBITERATION) {
							int64_t nbIteration;
							fi >> nbIteration;
							_tabAlgo[j]->setNbIteration(nbIteration);
							//_tabAlgo[j]->setEpsilon(minEpsilon);
						}
						else if (_tabAlgo[j]->getAlgoStopName() == EPSILON) {
							double epsilon;
							fi >> epsilon;
							_tabAlgo[j]->setEpsilon(epsilon);
							// _tabAlgo[j]->setNbIteration(maxNbIteration);
						}
						else if (_tabAlgo[j]->getAlgoStopName() == NBITERATION_EPSILON) {
							int64_t nbIteration;
							double epsilon;
							fi >> nbIteration;
							_tabAlgo[j]->setNbIteration(nbIteration);
							fi >> epsilon;
							_tabAlgo[j]->setEpsilon(epsilon);
						}

					}// end if stopRuleValue
					else {
						THROW(InputException, errorStopRuleValue);
					}

				}// end if StopRule
				else {
					THROW(InputException, errorStopRule);
				}

			}// end if algorithm
			else {
				THROW(InputException, errorAlgo);
			}

		}// end for j<nbAlgo*

	}// end if NbAlgo
	else {
		THROW(InputException, errorNbAlgo);
	}
}

void ClusteringStrategy::edit(std::ostream & out) {

	out << "Strategy : " << endl;

	out << "\tInitial start parameters method : " << endl;
	out << *_strategyInit << endl;

	out << "\tNumber of tries : " << _nbTry << endl;

	out << "\tNumber of algorithms in the strategy : " << _nbAlgo << endl;

	for (int64_t i = 0; i < _nbAlgo; i++) {
		out << "\tAlgorithm " << i + 1 << endl;
		_tabAlgo[i]->edit(out);
	}
	//out<<"\tNumber of strategy repetitions : "<<_nbStrategyTry<<endl;
}

//-----------
// ostream <<
//-----------
std::ostream & operator << (std::ostream & fo, ClusteringStrategy & strategy) {
	// nbTry
	fo << "nbTry : " << strategy._nbTry << endl;

	fo << "init : " << endl;
	fo << *(strategy._strategyInit) << endl;

	// nbAlgo
	fo << "nbAlgo : " << strategy._nbAlgo << endl;
	for (int64_t j = 0; j < strategy._nbAlgo; j++) {
		Algo * curAlgo = strategy._tabAlgo[j];
		fo << "Algo n " << j + 1 << " : " << endl;
		fo << (*curAlgo);
		fo << endl;
	}
	return fo;
}

const int64_t ClusteringStrategy::getNbTryInInit() const {
	return _strategyInit->getNbTry();
}

const int64_t ClusteringStrategy::getNbIterationInInit() const {
	return _strategyInit->getNbIteration();
}

const double ClusteringStrategy::getEpsilonInInit() const {
	return _strategyInit->getEpsilon();
}

void ClusteringStrategy::setNbTryInInit(int64_t nbTry) {
	_strategyInit->setNbTry(nbTry);
}

void ClusteringStrategy::setStrategyInitName(StrategyInitName initName) {
	_strategyInit->setStrategyInitName(initName);
}

void ClusteringStrategy::setNbIterationInInit(int64_t nbIteration) {
	_strategyInit->setNbIteration(nbIteration);
}

void ClusteringStrategy::setEpsilonInInit(double epsilon) {
	_strategyInit->setEpsilon(epsilon);
}

const AlgoStopName ClusteringStrategy::getStopNameInInit() const {
	return (_strategyInit->getStopName());
}

void ClusteringStrategy::setStopNameInInit(AlgoStopName stopName) {
	_strategyInit->setStopName(stopName);
}

}
