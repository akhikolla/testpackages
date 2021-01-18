/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringStrategyInit.cpp  description
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

#include "mixmod/Clustering/ClusteringStrategyInit.h"
#include "mixmod/Kernel/Parameter/GaussianSphericalParameter.h"
#include "mixmod/Kernel/Parameter/GaussianDiagParameter.h"
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Kernel/Parameter/GaussianHDDAParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkjhParameter.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/IO/Partition.h"

namespace XEM {

//------------
// Constructor
//------------
ClusteringStrategyInit::ClusteringStrategyInit() {
	_strategyInitName   = defaultStrategyInitName;
	_nbInitParameter    = 0;
	_tabInitParameter   = NULL;
	_nbPartition        = 0;
	_tabPartition       = NULL;
	_deleteTabParameter = false;

	_nbTry = defaultNbTryInInit;
	_nbIteration = defaultNbIterationInInit;
	_epsilon = defaultEpsilonInInit;
	setStopName(NBITERATION_EPSILON);
}

//------------------
// Copy constructor
//------------------
ClusteringStrategyInit::ClusteringStrategyInit(const ClusteringStrategyInit & strategyInit) {
	_strategyInitName = strategyInit.getStrategyInitName();
	_nbInitParameter  = strategyInit.getNbInitParameter();
	_nbPartition      = strategyInit.getNbPartition();
	_tabPartition     = NULL;
	
	if (_nbPartition != 0) {
		_tabPartition = new Partition*[_nbPartition];
		const Partition ** iPartition = strategyInit.getTabPartition();
		for (int64_t i = 0; i < _nbPartition; i++) {
			_tabPartition[i] = new Partition(*(iPartition[i])); //copy constructor
		}
	}

	_nbInitParameter = strategyInit.getNbInitParameter();
	_tabInitParameter = NULL;
	if (_nbInitParameter != 0) {
		_tabInitParameter = new Parameter*[_nbInitParameter];
		const Parameter ** iParameter = strategyInit.getTabInitParameter();
		for (int64_t i = 0; i < _nbInitParameter; i++) {
			_tabInitParameter[i] = (iParameter[i])->clone();
		}
	}

	_deleteTabParameter = false; // TODO ?

	_nbTry = strategyInit.getNbTry();
	_nbIteration = strategyInit.getNbIteration();
	_epsilon = strategyInit.getEpsilon();
	_stopName = strategyInit.getStopName();
}

//-----------
// Destructor
//-----------
ClusteringStrategyInit::~ClusteringStrategyInit() {
	if (_tabInitParameter && _deleteTabParameter) {
		for (unsigned int i = 0; i < _nbInitParameter; i++) {
			delete _tabInitParameter[i];
		}
		delete [] _tabInitParameter;
		_tabInitParameter = NULL;
	}

	if (_tabPartition) {
		for (unsigned  int i = 0; i < _nbPartition; i++) {
			delete _tabPartition[i];
			_tabPartition[i] = NULL;
		}
		delete [] _tabPartition;
		_tabPartition = NULL;
	}
}

//---------------------
// setStrategyInitName 
//---------------------
void ClusteringStrategyInit::setStrategyInitName(StrategyInitName initName) {

	// delete ?
	int64_t i;
	if (_tabInitParameter && _deleteTabParameter) {
		for (i = 0; i < _nbInitParameter; i++) {
			delete _tabInitParameter[i];
		}
		delete [] _tabInitParameter;
		_tabInitParameter = NULL;
	}

	if (_tabPartition) {
		for (i = 0; i < _nbPartition; i++) {
			delete _tabPartition[i];
			_tabPartition[i] = NULL;
		}
		delete [] _tabPartition;
		_tabPartition = NULL;
	}

	_strategyInitName = initName;
	_nbInitParameter   = 0;
	_tabInitParameter  = NULL;
	_nbPartition           = 0;
	_tabPartition          = NULL;
	_deleteTabParameter = false;

	_nbTry = defaultNbTryInInit;
	_nbIteration = defaultNbIterationInInit;
	if (_strategyInitName == SEM_MAX) {
		_nbIteration = defaultNbIterationInInitForSemMax;
		setStopName(NBITERATION);
	}
	if (_strategyInitName == USER || _strategyInitName == USER_PARTITION) {
		_nbTry = 1;
	}
	_epsilon = defaultEpsilonInInit;

}

//-------------------
// set init parameter
//-------------------
void ClusteringStrategyInit::setInitParam(std::string & paramFileName, int64_t position) {
	std::ifstream paramFile(paramFileName.c_str(), ios::in);
	if (! paramFile.is_open()) {
		THROW(InputException, wrongParamFileName);
	}
	if (_tabInitParameter != NULL) {
		_tabInitParameter[position]->input(paramFile);
		paramFile.close();
	}
	else {
		THROW(OtherException, internalMixmodError);
	}
}

//----------------
// setTabInitParam
//----------------
void ClusteringStrategyInit::setTabInitParameter(
		Parameter ** tabInitParameter, int64_t nbInitParameter) 
{
	int64_t i;
	if (_tabInitParameter && _deleteTabParameter) {
		for (i = 0; i < _nbInitParameter; i++) {
			delete _tabInitParameter[i];
		}
		delete [] _tabInitParameter;
		_tabInitParameter = NULL;
	}

	_tabInitParameter = tabInitParameter;
	_nbInitParameter = nbInitParameter;
}

//------------------  
//set Init Partition
//------------------
void ClusteringStrategyInit::setPartition(Partition * part, int64_t position) {

	if (position < 0) {
		THROW(OtherException, internalMixmodError);
	}

	if (part != NULL) {
		if (position < _nbPartition) {
			delete _tabPartition[position];
			_tabPartition[position] = part;
		}
		else {
			if (position == 0) { // this is the only situation that will be maintained
				_nbPartition = 1;
				_tabPartition = new Partition*[1];
				_tabPartition[0] = part;
			}
			else {
				THROW(InputException, badInitPart);
			}
		}
	}
	else {
		THROW(OtherException, internalMixmodError);
	}
}

//------------------
//set Init Partition
//----------------
void ClusteringStrategyInit::setPartition(std::string & partitionFileName, int64_t position) {
	std::ifstream partitionFile(partitionFileName.c_str(), ios::in);
	if (! partitionFile.is_open()) {
		THROW(InputException, wrongPartitionFileName);
	}
	try {
		if (position < _nbPartition) {
			delete _tabPartition[position];
			partitionFile >> (*_tabPartition[position]);
		}
		else {
			if (position == 0) { // this is the only situation that will be maintained
				_nbPartition = 1;
				_tabPartition = new Partition*[1];
				partitionFile >> (*_tabPartition[0]);
			}
			else {
				THROW(InputException, badInitPart);
			}
		}
		partitionFile.close();
	}
	catch (...) {
		THROW(InputException, badInitPart);
	}
}

//----------------
// setTabPartition
//----------------
void ClusteringStrategyInit::setTabPartition(Partition ** tabPartition, int64_t nbPartition) {
	int64_t i;
	if (_tabPartition) {
		for (i = 0; i < _nbPartition; i++) {
			delete _tabPartition[i];
			_tabPartition[i] = NULL;
		}
		delete [] _tabPartition;
		_tabPartition = NULL;
	}

	_tabPartition = tabPartition;
	_nbPartition = nbPartition;
}

//------------  
// setStopName
//-------------
void ClusteringStrategyInit::setStopName(AlgoStopName stopName) {
	if (_strategyInitName == SMALL_EM) {
		_stopName = stopName;
	}
	else if (_strategyInitName == SEM_MAX && stopName == NBITERATION) {
		_stopName = NBITERATION;
	}
	else {
		THROW(InputException, badSetStopNameInInit);
	}
}

//---------
// setNbTry
//---------
void ClusteringStrategyInit::setNbTry(int64_t nbTry) {
	if (_strategyInitName == SMALL_EM || _strategyInitName == CEM_INIT 
			|| _strategyInitName == RANDOM)
	{
		if (nbTry > maxNbTryInInit) {
			THROW(InputException, nbTryInInitTooLarge);
		}
		else if (nbTry < minNbTryInInit) {
			THROW(InputException, nbTryInInitTooSmall);
		}
		else {
			_nbTry = nbTry;
		}
	}
	else {
		THROW(InputException, badSetNbTryInInit);
	}
}

//----------------
// set NbIteration
//----------------
void ClusteringStrategyInit::setNbIteration(int64_t nbIteration) {
	if (_strategyInitName == SMALL_EM || _strategyInitName == SEM_MAX) {
		if (nbIteration > maxNbIterationInInit) {
			THROW(InputException, nbIterationInInitTooLarge);
		}
		else if (nbIteration < minNbIterationInInit) {
			THROW(InputException, nbIterationInInitTooSmall);
		}
		else {
			_nbIteration = nbIteration;
		}
	}
	else {
		THROW(InputException, badSetNbIterationInInit);
	}
}

//-----------
// setEpsilon
//-----------
void ClusteringStrategyInit::setEpsilon(double epsilon) {
	if (_strategyInitName == SMALL_EM) {
		if (epsilon > maxEpsilonInInit) {
			THROW(InputException, epsilonInInitTooLarge);
		}
		else if (epsilon < minEpsilonInInit) {
			THROW(InputException, epsilonInInitTooSmall);
		}
		else {
			_epsilon = epsilon;
		}
	}
	else {
		THROW(InputException, badSetEpsilonInInit);
	}
}

//-------
// verify
//-------
bool ClusteringStrategyInit::verify() const {
	bool res = true;
	//if init is USER or PARTITION, nbTry must be 1
	if ((_strategyInitName == USER_PARTITION || _strategyInitName == USER_PARTITION) 
			&& _nbTry != 1) 
	{
		res = false;
		THROW(InputException, wrongNbStrategyTryValue);
	}

	if (_strategyInitName == USER && _nbInitParameter != 1) {
		res = false;
		THROW(InputException, badNbParameterInInit);
	}

	if (_strategyInitName == USER_PARTITION && _nbPartition != 1) {
		res = false;
		THROW(InputException, badNbPartitionInInit);
	}

	return res;
}

// input
void ClusteringStrategyInit::input(std::ifstream & fi, Data *& data, int64_t nbNbCluster, 
		int64_t * tabNbCluster, ModelType * modelType, bool & alreadyRead) 
{
	std::string keyWord = "";
	std::string a = "";
	int64_t k;
	int64_t pbDimension  = data->_pbDimension;
	int64_t nbSample    = data->_nbSample;

	moveUntilReach(fi, "inittype");
	if (!fi.eof()) {

		fi >> a;
		if (a.compare("RANDOM") == 0) {
			setStrategyInitName(RANDOM);
			// nbTryInInit
			fi >> keyWord;
			ConvertBigtoLowString(keyWord);
			if (keyWord.compare("nbtryininit") == 0) {
				int64_t nbTry;
				fi >> nbTry;
				setNbTry(nbTry);
			}
			else {
				alreadyRead = true;
			}
		}

			// USER init type //
		else if (a.compare("USER") == 0) {
			setStrategyInitName(USER);

			// tabInitFileName
			fi >> keyWord;
			ConvertBigtoLowString(keyWord);

			if (keyWord.compare("initfile") == 0) {
				// _strategyInit->_nbInitParameter  = nbNbCluster;
				Parameter ** tabInitParameter = new Parameter * [nbNbCluster];
				std::string * tabFileName = new std::string[nbNbCluster];
				for (k = 0; k < nbNbCluster; k++) {
					tabFileName[k] = "";
				}
				readTabFileName(fi, nbNbCluster, tabFileName, keyWord);
				alreadyRead = true;

				for (k = 0; k < nbNbCluster; k++) {
					if (isEDDA(modelType->_nameModel)) {
						tabInitParameter[k] = new GaussianGeneralParameter(
								tabNbCluster[k], pbDimension, modelType, tabFileName[k]);
                        tabInitParameter[k]->setFilename(tabFileName[k]); //to be factorized when isHeterogeneous is implemented
					}
					else if (isBinary(modelType->_nameModel)) {
						int64_t * tabNbModality = (data->getBinaryData())->getTabNbModality();
						tabInitParameter[k] = new BinaryEkjhParameter(
								tabNbCluster[k], pbDimension, modelType, 
								tabNbModality, tabFileName[k]);
                        tabInitParameter[k]->setFilename(tabFileName[k]); //to be factorized when isHeterogeneous is implemented                       
					}
					else if (isHD(modelType->_nameModel)) {
						tabInitParameter[k] = new GaussianHDDAParameter(
								tabNbCluster[k], pbDimension, modelType, tabFileName[k]);
                        tabInitParameter[k]->setFilename(tabFileName[k]); //to be factorized when isHeterogeneous is implemented                        
					}
					else if (isHeterogeneous(modelType->_nameModel)) {
						//TODO
//						tabInitParameter[k] = new CompositeParameter(
//								tabNbCluster[k], pbDimension, modelType, (data->getBinaryData())->getTabNbModality());
					}
					else THROW(OtherException, internalMixmodError);
				}
				setTabInitParameter(tabInitParameter, nbNbCluster);
				delete[] tabFileName;
			}
			else {
				THROW(InputException, errorInitFile);
			}
		}

			// USER_PARTITION init type //
		else if (a.compare("USER_PARTITION") == 0) {
			setStrategyInitName(USER_PARTITION);

			//label
			fi >> keyWord;
			ConvertBigtoLowString(keyWord);
			if (keyWord.compare("initfile") == 0) {
				//_strategyInit->_nbPartition  = nbNbCluster;
				Partition ** tabPartition = new Partition * [nbNbCluster];
				std::string * tabFileName = new std::string[nbNbCluster];
				for (k = 0; k < nbNbCluster; k++) {
					tabFileName[k] = "";
				}
				readTabFileName(fi, nbNbCluster, tabFileName, keyWord);
				alreadyRead = true;

				for (k = 0; k < nbNbCluster; k++) {
					NumericPartitionFile partitionFile;
					partitionFile._fileName = tabFileName[k];
					partitionFile._format = FormatNumeric::defaultFormatNumericFile;
					partitionFile._type = TypePartition::defaultTypePartition;
					tabPartition[k] = new Partition(nbSample, tabNbCluster[k], partitionFile);
				}
				setTabPartition(tabPartition, nbNbCluster);
				delete[] tabFileName;
			}
			else {
				THROW(InputException, errorInitFile);
			}
		}

		else if (a.compare("SMALL_EM") == 0) {
			setStrategyInitName(SMALL_EM);
			// nbTry
			fi >> keyWord;
			ConvertBigtoLowString(keyWord);
			if (keyWord.compare("nbtryininit") == 0) {
				int64_t nbTry;
				fi >> nbTry;
				setNbTry(nbTry);
			}
			else {
				alreadyRead = true;
			}

			bool nbIteration = false;
			bool epsilon = false;
			//nbIterationInInit
			if (!alreadyRead) {
				fi >> keyWord;
			}
			ConvertBigtoLowString(keyWord);
			if (keyWord.compare("nbiterationininit") == 0) {
				nbIteration = true;
				int64_t nbIteration;
				fi >> nbIteration;
				setNbIteration(nbIteration);
				alreadyRead = false;
			}
			else {
				alreadyRead = true;
			}

			//epsilonInInit
			if (!alreadyRead) {
				fi >> keyWord;
				alreadyRead = false;
			}
			ConvertBigtoLowString(keyWord);
			if (keyWord.compare("epsilonininit") == 0) {
				epsilon = true;
				double epsilon;
				fi >> epsilon;
				setEpsilon(epsilon);
				alreadyRead = false;
			}
			else {
				alreadyRead = true;
			}
			if (nbIteration && epsilon) {
				setStopName(NBITERATION_EPSILON);
			}
			if (nbIteration && !epsilon) {
				setStopName(NBITERATION);
			}
			if (!nbIteration && epsilon) {
				setStopName(EPSILON);
			}
			if (!nbIteration && !epsilon) {
				setStopName(NBITERATION_EPSILON);
			}
		}

		else if (a.compare("CEM_INIT") == 0) {
			setStrategyInitName(CEM_INIT);
			fi >> keyWord;
			ConvertBigtoLowString(keyWord);
			if (keyWord.compare("nbtryininit") == 0) {
				int64_t nbTry;
				fi >> nbTry;
				setNbTry(nbTry);
			}
			else {
				alreadyRead = true;
			}
		}

		else if (a.compare("SEM_MAX") == 0) {
			setStrategyInitName(SEM_MAX);
			fi >> keyWord;
			ConvertBigtoLowString(keyWord);
			if (keyWord.compare("nbiterationininit") == 0) {
				int64_t nbIteration;
				fi >> nbIteration;
				setNbIteration(nbIteration);
			}
			else {
				alreadyRead = true;
			}
		}

		else {
			THROW(InputException, wrongStrategyInitName);
		}
	}
}

// ostream <<
//-----------
std::ostream & operator << (std::ostream & fo, ClusteringStrategyInit & strategyInit) {
	//strategyInitType
	std::string init = StrategyInitNameToString(strategyInit._strategyInitName);
	fo << "\t strategyInitName : " << init << endl;

	switch (strategyInit._strategyInitName) {
	case SMALL_EM:
		fo << "\t nbTryInInit : " << strategyInit._nbTry << endl;
		fo << "\t nbIterationInInit : " << strategyInit._nbIteration << endl;
		fo << "\t epsilonInInit : " << strategyInit._epsilon << endl;
		break;

	case RANDOM:
		fo << "\t nbTryInInit : " << strategyInit._nbTry << endl;
		break;

	case CEM_INIT:
		fo << "\t nbTryInInit : " << strategyInit._nbTry << endl;
		break;

	case SEM_MAX:
		fo << "\t nbIterationInInit : " << strategyInit._nbIteration << endl;
		break;

	case USER:
		// none or
		/*int64_t nbInitParameter = strategyInit._nbInitParameter;
		fo<<"\t nbInitParameter : "<<nbInitParameter<<endl;
		for (int64_t p=0; p<nbInitParameter; p++){
			(strategyInit.getTabInitParameter())[p]->edit();
		}
		fo<<endl;	*/
		break;

	case USER_PARTITION:
		// none
		break;
	}

	/*
int64_t nbPartition = strategyInit._nbPartition;
fo<<"\t nbLabel : "<<nbPartition<<endl;
for (int64_t p=0; p<nbPartition; p++){
  const XEMPartition * part = strategyInit._tabPartition[p];
  fo<<*(part)<<endl;
}*/
	return fo;
}

/*------------------------------------------------------
 initSMALL_EM
 ------------
 
 updated in this method : 
 - _parameter
 - _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk 
   (because an Estep is called to choose the bestParameter)
 Note : _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk wil be 're'computed in the following EStep
 So only _parameter have to be updated in this method
 
 -------------------------------------------------------*/
void ClusteringStrategyInit::initSMALL_EM(Model*& model) {
	// cout<<"init SMALL_EM, nbTryInInit="<<strategyInit->getNbTry()<<",
	// nbIteration = "<<strategyInit->getNbIteration()<<", epsilon = "<<strategyInit->getEpsilon()<<endl;
	model->setAlgoName(EM);
	double logLikelihood, bestLogLikelihood;
	// get pointer to Parameter
	Parameter * bestParameter = model->getParameter()->clone();
	int64_t  i, nbRunOfSmallEMOk = 0;
	bestLogLikelihood = 0.0;
	try {
		for (i = 0; i < _nbTry; i++) {
			nbRunOfSmallEMOk++;
			try {
				// one run of small EM
				model->getParameter()->reset();
				oneRunOfSmallEM(model, logLikelihood);
				//cout<<"sortie de oneRunOfSmallEM, LL = "<<logLikelihood<<endl;
				if ((nbRunOfSmallEMOk == 1) || (logLikelihood > bestLogLikelihood)) {
					bestLogLikelihood = logLikelihood;
					bestParameter->recopy(model->getParameter());
					//cout<<"best LL dans SMALL_EM : " <<bestLogLikelihood<<endl;
				}
			}
			catch (Exception& errorType) {
				nbRunOfSmallEMOk--;
			}
		}
	}
	catch (Exception&errorType) {
		throw;
	}

	if (nbRunOfSmallEMOk == 0) {
		THROW(NumericException, SMALL_EM_error);
	}

	// set best parameter
	model->setParameter(bestParameter);
	model->getParameter()->setModel(model);

	//delete bestParameter;
	//  cout<<"fin de init SMALL_EM, nb d'essais effectues="<<i<<endl;
}

//---------------------
// one run if small EM
//--------------------
void ClusteringStrategyInit::oneRunOfSmallEM(Model*& model, double & logLikelihood) {
	double lastLogLikelihood, eps;
	eps = 1000;
	model->initRANDOM(1);
	model->Estep();
	model->Mstep();
	logLikelihood = model->getLogLikelihood(true);  // true : to compute fik
	int64_t  nbIteration = 1;
	bool continueAgain = true;
	while (continueAgain) {
		//    cout<<"while de oneRunOfSmallEM, nbIteration = "<<nbIteration<<endl;
		//(nbIteration < strategyInit->getNbIteration()) && (eps > strategyInit->getEpsilon())){
		lastLogLikelihood = logLikelihood;
		model->Estep();
		model->Mstep();
		nbIteration++;
		// update continueAgain
		switch (_stopName) {
		case NBITERATION:
			continueAgain = (nbIteration < _nbIteration);
			break;
		case EPSILON:
			logLikelihood = model->getLogLikelihood(true);  // true : to compute fik
			eps = fabs(logLikelihood - lastLogLikelihood);
			// on ajoute un test pour ne pas faire trop d'iterations quand meme ....
			continueAgain = (eps > _epsilon && (nbIteration < maxNbIterationInInit));
			//continueAgain = (eps > strategyInit->getEpsilon());
			break;
		case NBITERATION_EPSILON:
			logLikelihood = model->getLogLikelihood(true);  // true : to compute fi
			eps = fabs(logLikelihood - lastLogLikelihood);
			continueAgain = ((eps > _epsilon) && (nbIteration < _nbIteration));
			break;
		default: THROW(OtherException, internalMixmodError);
		}
	}
	if (_stopName == NBITERATION) { // logLikelihood is an output
		logLikelihood = model->getLogLikelihood(true);  // true : to compute fi
	}
	//cout<<"Fin de oneRunOfSmallEM, nb d'iterations effectuees = "<<nbIteration<<", logLikelihood = "<<logLikelihood<<endl;
}

/*---------------------------------------------------
 initCEM_INIT
 ---------------------------------------------------
 updated in this method : 
 - _parameter
 - _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk 
   (because an Estep and a CStep are called to choose the bestParameter)
 Note : _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk wil be 're'computed in the following EStep
 So only _parameter have to be updated in this method
 
 ---------------------------------------------------*/
void ClusteringStrategyInit::initCEM_INIT(Model*& model) {
	//cout<<"init CEM, nbTryInInit="<<strategyInit->getNbTry()<<endl;
	model->setAlgoName(CEM);
	int64_t i;
	double cLogLikelihood = 0.0, oldLogLikelihood = 0.0;

	// get pointer to Parameter
	Parameter * bestParameter = model->getParameter()->clone();
	int64_t nbRunOfCEMOk = 0;
	double bestCLogLikelihood = 0.0;

	for (i = 0; i < _nbTry; i++) {
		nbRunOfCEMOk++;
		try {
			model->getParameter()->reset(); // reset to default values
			model->initRANDOM(1);
			model->setAlgoName(CEM);
			int64_t nbIter = 0;
			bool fin = false;
			while (!fin && nbIter <= maxNbIterationInCEM_INIT) {
				model->Estep();
				model->Cstep();
				model->Mstep();
				nbIter++;
				if (nbIter == 1) {
					oldLogLikelihood = model->getCompletedLogLikelihood();
				}
				else {
					cLogLikelihood = model->getCompletedLogLikelihood();
					if (cLogLikelihood == oldLogLikelihood) {
						fin = true;
					}
					else {
						oldLogLikelihood = cLogLikelihood;
					}
				}
			}
			//cout<<"dans init CEM, nb d'iterations effectuées : "<<nbIter<<endl;
			// Compute log-likelihood 
			cLogLikelihood = model->getCompletedLogLikelihood();
			// Comparaison of log-likelihood between step p and p-1 
			if ((nbRunOfCEMOk == 1) || (cLogLikelihood > bestCLogLikelihood)) {
				bestCLogLikelihood = cLogLikelihood;
				bestParameter->recopy(model->getParameter());
			}
			//cout<<"nbIter : "<<nbIter<<endl;
		}
		catch (...) {
			nbRunOfCEMOk--;
		}
	}

	if (nbRunOfCEMOk == 0) {
		// set best parameter
		model->setParameter(bestParameter);
		model->getParameter()->setModel(model);
		THROW(NumericException, CEM_INIT_error);
	}

	// set Best parameter
	model->setParameter(bestParameter);
	model->getParameter()->setModel(model);
	//cout<<"fin de init CEM, nb d'essais effectues="<<i<<endl;
}

/*---------------------------------------------------------
 Initialization by SEM
 ---------------------
 
 updated in this method : 
 - _parameter
 - _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk 
   (because an Estep and a SStep are called to choose the bestParameter)
 Note : _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk wil be 're'computed in the following EStep
 So, only _parameter have to be updated in this method
 
 -------------------------------------------------------*/
void ClusteringStrategyInit::initSEM_MAX(Model*& model) {
	//cout<<"init SEM_MAX, nbTryInInit="<<strategyInit->getNbIteration()<<endl;
	model->setAlgoName(SEM);
	int64_t  j;
	double logLikelihood, bestLogLikelihood;
	// get pointer to Parameter
	Parameter * bestParameter = model->getParameter()->clone();
	int64_t  nbRunOfSEMMAXOk = 0;
	bestLogLikelihood = 0.0;
	//  int64_t  bestIndex=0;

	for (j = 0; j < _nbIteration; j++) {
		nbRunOfSEMMAXOk++;
		try {
			model->getParameter()->reset();
			model->initRANDOM(1);
			model->Estep();
			model->Sstep();
			model->Mstep();
			// Compute log-likelihood 
			logLikelihood = model->getLogLikelihood(true);  // true : to compute fik
			if ((nbRunOfSEMMAXOk == 1) || (logLikelihood > bestLogLikelihood)) {
				bestLogLikelihood = logLikelihood;
				bestParameter->recopy(model->getParameter());
				//        bestIndex = j;
			}
		}
		catch (...) {
			nbRunOfSEMMAXOk--;
		}
	}

	if (nbRunOfSEMMAXOk == 0) {
		THROW(NumericException, SEM_MAX_error);
	}

	// set Best parameter
	model->setParameter(bestParameter);
	model->getParameter()->setModel(model);
	//cout<<"fin de init SEM_MAX, nb d'iterations effectuees="<<j<<" meilleure solution : "<<bestIndex<<endl;
}

}
