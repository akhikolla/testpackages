/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/GaussianHDDAParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/GaussianHDDAParameter.h"
#include "mixmod/Kernel/Parameter/GaussianParameter.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/GaussianSample.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Utilities/Util.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Utilities/Random.h"
#include "mixmod/Matrix/Matrix.h"
#include "mixmod/Matrix/DiagMatrix.h"
#include "mixmod/Matrix/SphericalMatrix.h"
#include "mixmod/Matrix/SymmetricMatrix.h"
#include "mixmod/Matrix/GeneralMatrix.h"

namespace XEM {

/****************/
/* Constructors */
/****************/

GaussianHDDAParameter::GaussianHDDAParameter() : GaussianParameter() {
	THROW(OtherException, wrongConstructorType);
}

//-------------------------------------------------------------------------------------
// constructor called by XEMModel
//-------------------------------------------------------------------------------------
GaussianHDDAParameter::GaussianHDDAParameter(Model * iModel, ModelType * iModelType)
: GaussianParameter(iModel, iModelType)
{
	int64_t k;

	_tabAkj = new double*[_nbCluster];
	_tabBk = new double[_nbCluster];
	_tabShape = new DiagMatrix*[_nbCluster];
	_tabQk = new GeneralMatrix*[_nbCluster];
	_W = new SymmetricMatrix(_pbDimension); //Id
	_tabDk = new int64_t [_nbCluster];
	_tabGammak = NULL;
	_Gammak = NULL;

	for (k = 0; k < _nbCluster; k++) {
		_tabShape[k] = new DiagMatrix(_pbDimension); //Id
		_tabQk[k] = new GeneralMatrix(_pbDimension); //Id
		_tabWk[k] = new SymmetricMatrix(_pbDimension); // Id
		_tabDk[k] = 0;
	}
	__storeDim = _pbDimension * (_pbDimension + 1) / 2;

	if ((iModelType->_tabSubDimensionFree != NULL) &&
			isFreeSubDimension(iModelType->_nameModel)) {
		for (k = 0; k < _nbCluster; k++) {
			_tabDk[k] = iModelType->_tabSubDimensionFree[k];
		}
	}
	else if ((iModelType->_subDimensionEqual != 0) &&
			!isFreeSubDimension(iModelType->_nameModel)) {
		for (k = 0; k < _nbCluster; k++) {
			_tabDk[k] = iModelType->_subDimensionEqual;
		}
	}

	for (k = 0; k < _nbCluster; k++) {
		_tabAkj[k] = new double[_tabDk[k]];
		for (int64_t j = 0; j < _tabDk[k]; j++) {
			_tabAkj[k][j] = 1.0;
		}
		_tabBk[k] = 1.0; // /(_pbDimension - _tabDk[k]);
	}
}

//constructeur avec une initialisation USER
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
GaussianHDDAParameter::GaussianHDDAParameter(int64_t iNbCluster, int64_t iPbDimension,
		ModelType * iModelType, std::string & iFileName)
: GaussianParameter(iNbCluster, iPbDimension, iModelType)
{
	int64_t k;
	_tabAkj = new double*[_nbCluster];
	_tabBk = new double[_nbCluster];
	_tabDk = new int64_t [_nbCluster];
	_tabGammak = NULL;
	_Gammak = NULL;
	__storeDim = _pbDimension * (_pbDimension + 1) / 2;

	_tabShape = new DiagMatrix*[_nbCluster];
	_tabQk = new GeneralMatrix*[_nbCluster];

	for (k = 0; k < _nbCluster; k++) {
		_tabShape[k] = new DiagMatrix(_pbDimension); //Id
		_tabQk[k] = new GeneralMatrix(_pbDimension); //Id
		_tabWk[k] = new SymmetricMatrix(_pbDimension); //Id
		_tabAkj[k] = NULL;
	}
	_W = new SymmetricMatrix(_pbDimension); //Id

	// read parameters in file iFileName//
	if (iFileName.compare("") != 0) {
		std::ifstream paramFile(iFileName.c_str(), ios::in);
		if (!paramFile.is_open()) {
			THROW(InputException, wrongParamFileName);
		}
		input(paramFile);
		paramFile.close();
	}
}

//constructeur par copie
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
GaussianHDDAParameter::GaussianHDDAParameter(const GaussianHDDAParameter * iParameter)
: GaussianParameter(iParameter)
{
	int64_t k;
	__storeDim = _pbDimension * (_pbDimension + 1) / 2;
	int64_t * iTabD = iParameter->getTabD();
	double ** iTabA = iParameter->getTabA();
	double *iTabB = iParameter->getTabB();

	_tabShape = new DiagMatrix*[_nbCluster];
	_tabQk = new GeneralMatrix*[_nbCluster];
	_tabDk = new int64_t [_nbCluster];
	_tabAkj = new double*[_nbCluster];
	_tabBk = new double[_nbCluster];

	DiagMatrix ** iTabShape = iParameter->getTabShape();
	GeneralMatrix ** iTabQ = iParameter->getTabQ();
	Matrix ** iTabWk = iParameter->getTabWk();

	_tabGammak = NULL;
	_Gammak = NULL;
	_W = new SymmetricMatrix(_pbDimension); //Id
	(* _W) = iParameter->getW();

	recopyTab(iTabD, _tabDk, _nbCluster);
	recopyTab(iTabB, _tabBk, _nbCluster);

	for (k = 0; k < _nbCluster; k++) {
		_tabAkj[k] = new double[_tabDk[k]];
		recopyTab(iTabA[k], _tabAkj[k], _tabDk[k]);
		_tabShape[k] = new DiagMatrix(iTabShape[k]);
		_tabQk[k] = new GeneralMatrix(iTabQ[k]); // copy constructor
		_tabWk[k] = new SymmetricMatrix(_pbDimension); //Id
		(* _tabWk[k]) = iTabWk[k];
	}
}

/**************/
/* Destructor */
/**************/
GaussianHDDAParameter::~GaussianHDDAParameter() {
	int64_t k;
	if (_tabShape) {
		for (k = 0; k < _nbCluster; k++) {
			delete _tabShape[k];
			_tabShape[k] = NULL;
		}
		delete[] _tabShape;
		_tabShape = NULL;
	}

	if (_tabQk) {
		for (k = 0; k < _nbCluster; k++) {
			delete _tabQk[k];
			_tabQk[k] = NULL;
		}
		delete[] _tabQk;
		_tabQk = NULL;
	}

	if (_tabAkj) {
		for (k = 0; k < _nbCluster; k++) {
			delete[] _tabAkj[k];
			_tabAkj[k] = NULL;
		}
		delete[] _tabAkj;
		_tabAkj = NULL;
	}

	if (_tabBk) {
		delete[] _tabBk;
		_tabBk = NULL;
	}

	if (_tabDk) {
		delete[] _tabDk;
		_tabDk = NULL;
	}

	if (_Gammak) {
		for (k = 0; k < _nbCluster; k++) {
			delete[] _Gammak[k];
			_Gammak[k] = NULL;
		}
		delete[] _Gammak;
		_Gammak = NULL;
	}

	if (_tabGammak) {
		for (k = 0; k < _nbCluster; k++) {
			delete _tabGammak[k];
		}
		delete[] _tabGammak;
		_tabGammak = NULL;
	}
}

//------------------------
// reset to default values
//------------------------
void GaussianHDDAParameter::reset() {
	THROW(OtherException, internalMixmodError);
	// faire ensuite : XEMGaussianParameter::reset();
}

/*********/
/* clone */
/*********/
Parameter * GaussianHDDAParameter::clone() const {
	GaussianHDDAParameter * newParam = new GaussianHDDAParameter(this);
	return (newParam);
}

/*********/
/* MStep */
/********/
void GaussianHDDAParameter::MStep() {
	GaussianParameter::MStep();

	switch (_modelType->_nameModel) {
	case (Gaussian_HD_pk_AkjBkQkDk):
	case (Gaussian_HD_p_AkjBkQkDk):
	case (Gaussian_HD_pk_AkjBkQkD):
	case (Gaussian_HD_p_AkjBkQkD):
		computeAkjBkQk();
		break;

	case (Gaussian_HD_pk_AkBkQkDk):
	case (Gaussian_HD_p_AkBkQkDk):
	case (Gaussian_HD_pk_AkBkQkD):
	case (Gaussian_HD_p_AkBkQkD):
		computeAkBkQk();
		break;

	case (Gaussian_HD_p_AkjBQkD):
	case (Gaussian_HD_pk_AkjBQkD):
		computeAkjBQk();
		break;

	case (Gaussian_HD_p_AkBQkD):
	case (Gaussian_HD_pk_AkBQkD):
		computeAkBQk();
		break;

	case (Gaussian_HD_p_AjBkQkD):
	case (Gaussian_HD_pk_AjBkQkD):
		computeAjBkQk();
		break;

	case (Gaussian_HD_p_AjBQkD):
	case (Gaussian_HD_pk_AjBQkD):
		computeAjBQk();
		break;
	default: THROW(OtherException, internalMixmodError);
	}
}

//**********/
/* getPdf */
/**********/
// returns the density of a sample using the cost function property Cost = -2LL
double GaussianHDDAParameter::getPdf(int64_t iSample, int64_t kCluster)const {

	GaussianData * data = _model->getGaussianData();
	double * ligne = (data->getYStore())[iSample];
	double normPdf = 0, K = 0;

	Parameter * parameter = _model->getGaussianParameter();
	GaussianParameter* gparameter = (GaussianParameter*) parameter;
	double ** tabMean = gparameter->getTabMean();
	double* tabProportion = gparameter->getTabProportion();
	double * xiMoinsMuk = new double[_pbDimension];
	double* tabShapek_store = new double[_pbDimension];

	//----------------calcul le produit Q*t(Q)---------------
	SymmetricMatrix* Pk = new SymmetricMatrix(_pbDimension); //Id
	Pk->compute_as_M_tM(_tabQk[kCluster], _tabDk[kCluster]);

	//---------------calcul du produit A = Qi_L^-1_t(Qi)------------------
	SymmetricMatrix* A = new SymmetricMatrix(_pbDimension); //id

	double sum_lambda = 0.0;
	for (int64_t j = 0; j < _tabDk[kCluster]; j++) {
		tabShapek_store[j] = 1 / _tabAkj[kCluster][j];
		sum_lambda += log(_tabAkj[kCluster][j]);
	}
	for (int64_t j = _tabDk[kCluster]; j < _pbDimension; j++) {
		tabShapek_store[j] = 0.0;
	}
	A->compute_as_O_S_O(_tabQk[kCluster], tabShapek_store);

	double constante = sum_lambda + (_pbDimension - _tabDk[kCluster]) * log(_tabBk[kCluster])
			- 2 * log(tabProportion[kCluster]) + _pbDimension * log(2 * XEMPI);

	//-----------soustraction des moyennes
	for (int64_t j = 0; j < _pbDimension; j++) {
		xiMoinsMuk[j] = ligne[j] - tabMean[kCluster][j];
	}

	//-----------produit Qt(Q)(x-muk)----------------------
	SymmetricMatrix * Pi = new SymmetricMatrix(_pbDimension); //Id
	Pi->compute_as_M_V(Pk, xiMoinsMuk);
	double * storePi = Pi->getStore();

	//----------------norme(muk-Pi(x))_A = t(muk-Pi(x))_A_(muk-Pi(x))-----------
	double normeA = A->norme(xiMoinsMuk);

	//----------------norme(x-Pi(x))------------------------------
	double norme = 0.0;
	for (int64_t i = 0; i < _pbDimension; i++) {
		storePi[i] += tabMean[kCluster][i];
		norme += (ligne[i] - storePi[i])*(ligne[i] - storePi[i]);
	}

	//----------------calcul de normPdf---------
	K = normeA + 1.0 / _tabBk[kCluster] * norme + constante;
	//normPdf = exp(-1 / 2 * K);
	normPdf = exp(-0.5 * K);


	delete Pk;
	delete A;
	delete Pi;
	delete[] xiMoinsMuk;
	xiMoinsMuk = NULL;
	delete[] tabShapek_store;
	tabShapek_store = NULL;

	return normPdf;
}

/*************/
/* getAllPdf */
/*************/
// returns the density of all the data (tabFik)
void GaussianHDDAParameter::getAllPdf(double ** tabFik, double * tabProportion)const {
	//cout<<"XEMGaussianHDDAParameter::getAllPdf"<<endl;
	/*cout<<"_tabBk : "<<_tabBk[0]<<endl;
	cout<<"_tabBk : "<<_tabBk[1]<<endl;*/
	double ** Cost = computeCost(_tabQk);
	int64_t nbSample = _model->getNbSample();

	for (int64_t i = 0; i < nbSample; i++) {
		for (int64_t k = 0; k < _nbCluster; k++) {
			tabFik[i][k] = exp(-0.5 * Cost[k][i]);
			//cout<<" Cost[k][i] :  "<<Cost[k][i]<<endl;
			//cout<<"tabFik[i][k] :  "<<tabFik[i][k]<<endl;
		}
	}

	for (int64_t k = 0; k < _nbCluster; k++) {
		delete[] Cost[k];
		Cost[k] = NULL;
	}
	delete[] Cost;
	Cost = NULL;
}

/**********/
/* getPdf */
/*********/
// returns the density of a sample using the cost function property Cost = -2LL
double GaussianHDDAParameter::getPdf(Sample * x, int64_t kCluster)const {

	double * ligne = ((GaussianSample*) x)->getTabValue();
	double normPdf = 0, K = 0;
	Parameter * parameter = _model->getGaussianParameter();
	GaussianParameter* gparameter = (GaussianParameter*) parameter;
	double ** tabMean = gparameter->getTabMean();
	double* tabProportion = gparameter->getTabProportion();
	double * xiMoinsMuk = new double[_pbDimension];
	double* tabShapek_store = new double[_pbDimension];

	//----------------calcul le produit Q*t(Q)---------------
	SymmetricMatrix* Pk = new SymmetricMatrix(_pbDimension); //Id
	Pk->compute_as_M_tM(_tabQk[kCluster], _tabDk[kCluster]);

	//---------------calcul du produit A = Qi_L^-1_t(Qi)------------------
	SymmetricMatrix* A = new SymmetricMatrix(_pbDimension); //Id

	double sum_lambda = 0.0;
	for (int64_t j = 0; j < _tabDk[kCluster]; j++) {
		tabShapek_store[j] = 1 / _tabAkj[kCluster][j];
		sum_lambda += log(_tabAkj[kCluster][j]);
	}
	for (int64_t j = _tabDk[kCluster]; j < _pbDimension; j++) {
		tabShapek_store[j] = 0.0;
	}

	A->compute_as_O_S_O(_tabQk[kCluster], tabShapek_store);

	double constante = sum_lambda + (_pbDimension - _tabDk[kCluster]) * log(_tabBk[kCluster])
			- 2 * log(tabProportion[kCluster]) + _pbDimension * log(2 * XEMPI);

	//-----------soustraction des moyennes
	for (int64_t j = 0; j < _pbDimension; j++) {
		xiMoinsMuk[j] = ligne[j] - tabMean[kCluster][j];
	}

	//-----------produit Qt(Q)(x-muk)----------------------
	SymmetricMatrix * Pi = new SymmetricMatrix(_pbDimension); //Id
	Pi->compute_as_M_V(Pk, xiMoinsMuk);
	double * storePi = Pi->getStore();

	//----------------norme(muk-Pi(x))_A = t(muk-Pi(x))_A_(muk-Pi(x))-----------
	double normeA = A->norme(xiMoinsMuk);

	//----------------norme(x-Pi(x))------------------------------
	double norme = 0.0;
	for (int64_t i = 0; i < _pbDimension; i++) {
		storePi[i] += tabMean[kCluster][i];
		norme += (ligne[i] - storePi[i])*(ligne[i] - storePi[i]);
	}

	//----------------calcul de normPdf---------
	K = normeA + 1.0 / _tabBk[kCluster] * norme + constante;
	//normPdf = exp(-1 / 2 * K);
	normPdf = exp(-0.5*K);

	delete Pk;
	delete A;
	delete Pi;
	delete[] xiMoinsMuk;
	xiMoinsMuk = NULL;
	delete[] tabShapek_store;
	tabShapek_store = NULL;

	return normPdf;
}

/*********/
/* input */
/*********/
// reads the input data file for HDDA models
void GaussianHDDAParameter::input(std::ifstream & fi) {
	int64_t j, k;
  for (k = 0; k < _nbCluster; k++) {
    // Proportions //
		_tabProportion[k] = getDoubleFromStream(fi);

    // Center (mean) //
    for (j = 0; j < _pbDimension; j++)
      _tabMean[k][j] = getDoubleFromStream(fi);

    // Sub Dimension  //
    fi >> _tabDk[k];
    if (_tabAkj[k]) {
      //cout<<_tabAkj[k][0]<<endl;
			//cout<<"ok"<<endl;
			delete [] _tabAkj[k];
			_tabAkj[k] = NULL;
		}
		_tabAkj[k] = new double[_tabDk[k]];

		// Parameters Akj //
		for (j = 0; j < _tabDk[k]; j++) {
			fi >> _tabAkj[k][j];
		}

		// Parameters Bk //
		fi >> _tabBk[k];
		// Orientation matrix //
		_tabQk[k]->input(fi, _tabDk[k]); // virtual method
	} // end for k
}

/*****************/
/* computeTabWkW */
/*****************/
//computes W_k and W but also _tabGammak and Gammak
void GaussianHDDAParameter::computeTabWkW() {
	double* tabNk = _model->getTabNk();
	double ** tabCik = _model->getTabCik();
	int64_t nbSample = _model->getNbSample();
	GaussianData * data = _model->getGaussianData();
	double * weight = data->_weight;
	int64_t i;
	double ** matrix = data->getYStore(); // to store x_i i=1,...,n
	int64_t j, k, l, dimStoreG;

	k = 0;
	while (k < _nbCluster) {
		if (tabNk[k] < _pbDimension) {
			_tabGammak = new SymmetricMatrix*[_nbCluster];
			k = _nbCluster;
		}
		else {
			k++;
		}
	}

	GaussianParameter::computeTabWkW();

	for (k = 0; k < _nbCluster; k++) {
		if ((tabNk[k] < _pbDimension) && (_tabDk[k] < (tabNk[k] + 1))) {
			l = 0;
			double test = floor(tabNk[k]);
			if (test == tabNk[k]) {
				_Gammak = new double*[_nbCluster];
				int64_t nk = (int64_t) tabNk[k];
				_tabGammak[k] = new SymmetricMatrix(nk); //Id
				dimStoreG = nk * _pbDimension;
				_Gammak[k] = new double[dimStoreG];
				for (i = 0; i < nbSample; i++) {
					if (tabCik[i][k] == 1) {
						for (j = 0; j < _pbDimension; j++) {
							//multiplier par le poids ?
							_Gammak[k][l] = matrix[i][j] * weight[i] - _tabMean[k][j];
							// _Gammak[k][l] = matrix[i][j] - _tabMean[k][j];
							l += 1;
						}
					}
				}
				//----------------calcul le produit Gamma*t(Gamma)---------------
				*(_tabGammak[k]) = 0.0;
				_tabGammak[k]->compute_M_tM(_Gammak[k], dimStoreG);
			}
			else {
				THROW(NumericException, tabNkNotInteger);
			} // end test
		}
	} // fin for k
}

//-------------------------------------------
// initialize attributes before an InitRandom
//-------------------------------------------
void GaussianHDDAParameter::initForInitRANDOM() {
	THROW(OtherException, internalMixmodError);
}

//---------------------------
// initForInitUSER_PARTITION
//--------------------------
void GaussianHDDAParameter::initForInitUSER_PARTITION(int64_t & nbInitializedCluster,
		bool * tabNotInitializedCluster, Partition * initPartition)
{
	computeTabMeanInitUSER_PARTITION(
			nbInitializedCluster, tabNotInitializedCluster, initPartition);
	// initialization of _tabAkj, ;... :

	DiagMatrix * matrixDataVar = new DiagMatrix(_pbDimension, 0.0);
	computeGlobalDiagDataVariance(matrixDataVar);
	matrixDataVar->sortDiagMatrix();

	double * store = matrixDataVar->getStore();
	double sum_lambda = 0.0;

	for (int64_t k = 0; k < _nbCluster; k++) {
		*_tabQk[k] = 1.0;
	}

	for (int64_t j = 0; j < _tabDk[0]; j++) {
		_tabAkj[0][j] = store[j];
		sum_lambda += _tabAkj[0][j];
	}

	double trace = matrixDataVar->computeTrace();

	_tabBk[0] = 1.0 / (_pbDimension - _tabDk[0])*(trace - sum_lambda);

	for (int64_t k = 1; k < _nbCluster; k++) {
		for (int64_t j = 0; j < _tabDk[k]; j++) {
			_tabAkj[k][j] = store[j];
		}
		_tabBk[k] = _tabBk[0];
	}

	if (nbInitializedCluster != _nbCluster) {
		/* ca ne devrait pas arriver puisque on a verifie dans XEMOldInput que :
		- la partition doit etre 'complete' (y compris que toutes les classes
		  sont representees) quand on a l'algo M
		- que l'on autorise que M ou MAP pour les modeles HD
		- que MAP => USER (et donc pas 'USER_PARTITION' avec MAP)
		 */
		THROW(OtherException, internalMixmodError);
	}
	delete matrixDataVar;
}

/**********************/
/* initUSER_PARTITION */
/**********************/
/*void XEMGaussianHDDAParameter:: initUSER_PARTITION(){
	computeTabMeanInitUSER_PARTITION();

   _model->computeNk();
   computeTabWkW();
   if (_tabDk[0]==0){
	 computeTabDk();
   }

   switch(_modelType){
	 case (Gaussian_HD_pk_AkjBkQkDk) :
	 case (Gaussian_HD_p_AkjBkQkDk) :
	 case (Gaussian_HD_pk_AkjBkQkD) :
	 case (Gaussian_HD_p_AkjBkQkD) :
	   computeAkjBkQk();
	   break;

	 case (Gaussian_HD_pk_AkBkQkDk) :
	 case (Gaussian_HD_p_AkBkQkDk) :
	 case (Gaussian_HD_pk_AkBkQkD) :
	 case (Gaussian_HD_p_AkBkQkD) :
	   computeAkBkQk();
	   break;

	 case (Gaussian_HD_p_AkjBQkD):
	 case (Gaussian_HD_pk_AkjBQkD):
	   computeAkjBQk();
	   break;

	 case (Gaussian_HD_p_AkBQkD):
	 case (Gaussian_HD_pk_AkBQkD):
	   computeAkBQk();
	   break;

	 case (Gaussian_HD_p_AjBkQkD):
	 case (Gaussian_HD_pk_AjBkQkD):
	   computeAjBkQk();
	   break;

	 case (Gaussian_HD_p_AjBQkD):
	 case (Gaussian_HD_pk_AjBQkD):
	   computeAjBQk();
	   break;
   }
}*/

/************/
/* initUSER */
/***********/
void GaussianHDDAParameter::initUSER(Parameter * iParam) {
	GaussianHDDAParameter * param = (GaussianHDDAParameter *) iParam;

	double ** iTabMean = param->getTabMean();
	double * iTabProportion = param->getTabProportion();
	Matrix ** iTabWk = param->getTabWk();
	int64_t * iTabDk = param->getTabD();
	double * iTabBk = param->getTabB();
	double ** iTabAkj = param->getTabA();
	GeneralMatrix** iTabQk = param->getTabQ();

	int64_t k;
	recopyTab(iTabBk, _tabBk, _nbCluster);
	for (k = 0; k < _nbCluster; k++) {
		if (_tabDk[k] != iTabDk[k]) {
			THROW(InputException, differentSubDimensionsWithMAP);
		}
		else {
			recopyTab(iTabMean[k], _tabMean[k], _pbDimension);
			recopyTab(iTabAkj[k], _tabAkj[k], _tabDk[k]);
			// recopy _tabWk
			(*_tabWk[k]) = iTabWk[k];

			if (!hasFreeProportion(_modelType->_nameModel)) {
				_tabProportion[k] = 1.0 / _nbCluster;
			}
			else {
				_tabProportion[k] = iTabProportion[k];
			}
			(*_tabQk[k]) = iTabQk[k];
		}
	}
}

/***********/
/* recopy  */
/***********/
void GaussianHDDAParameter::recopy(Parameter * otherParameter) {

	GaussianParameter * iParameter = otherParameter->getGaussianParameter();

	recopyTab(iParameter->getTabMean(), _tabMean, _nbCluster, _pbDimension);
	//_W->recopy(iParameter->getW());
	(* _W) = iParameter->getW();

	int64_t k;
	Matrix ** iTabWk = iParameter->getTabWk();
	for (k = 0; k < _nbCluster; k++) {
		//_tabWk[k]->recopy(iTabWk[k]);
		(* _tabWk[k]) = iTabWk[k];
	}
}

//-------------------
//getLogLikelihoodOne
//-------------------
double GaussianHDDAParameter::getLogLikelihoodOne()const {
	//cout<<"XEMGaussianHDDAParameter::getLogLikelihoodOne"<<endl;
	/* Compute log-likelihood for one cluster
	   useful for NEC criterion */
	/* Initialization */
	/*int64_t nbSample = _model->getNbSample();
	int64_t i;
	XEMGaussianData * data = _model->getGaussianData();
	int64_t j, n, l;
	double logLikelihoodOne;         // Log-likelihood for k=1
	double * Mean = new double[_pbDimension] ;
	double ** y  = data->_yStore;
	double * yi;
   // XEMGeneralMatrix * Sigma = new XEMGeneralMatrix(_pbDimension);
	XEMSymmetricMatrix * W     = new XEMSymmetricMatrix(_pbDimension);
	 *W = 0.0;
	double norme, logDet, detDiagW ;
	double * weight = data->_weight;

	//  Mean Estimator (empirical estimator)
	double totalWeight = data->_weightTotal;
	computeMeanOne(Mean,weight,y,nbSample,totalWeight);
	weight = data->_weight;
	  // Compute the Cluster Scattering Matrix W
	int64_tp; // parcours
	double * xiMoinsMuk = data->getTmpTabOfSizePbDimension();
	double tmp;
	  for(i=0; i<nbSample ; i++){
		 yi = y[i];
		 for(p=0 ; p<_pbDimension ; p++){
			xiMoinsMuk[p] = yi[p] - Mean[p];
		 }
		 W->add(xiMoinsMuk, weight[i]);
	  }

	//Compute determinant of diag(W)
   // logDet   = W->detDiag(minDeterminantDiagWValueError);  // virtual
	detDiagW = powAndCheckIfNotNull(logDet ,1.0/_pbDimension);

	Sigma->divise(W, totalWeight); // virtual

	// inverse of Sigma
	XEMGeneralMatrix * SigmaMoins1 = new XEMGeneralMatrix(_pbDimension);
	SigmaMoins1->inverse(Sigma);// virtual

	double detSigma  = Sigma->determinant(minDeterminantSigmaValueError); // virtual

	// Compute the log-likelihood for one cluster (k=1)
	logLikelihoodOne = 0.0;
	for (i=0; i<nbSample; i++){
	  yi = y[i];
	  for(p=0; p<_pbDimension; p++){
		xiMoinsMuk[p] = yi[p] - Mean[p];
	  }

	  norme             = SigmaMoins1->norme(xiMoinsMuk);   // virtual
	  logLikelihoodOne += norme * weight[i];
	}

   logLikelihoodOne += totalWeight * ( data->getPbDimensionLog2Pi() + log(detSigma)) ;
   logLikelihoodOne *= -0.5 ;

   delete W;
   delete Sigma;
   delete SigmaMoins1;

   delete[] Mean;

   return logLikelihoodOne;*/
	return 0.0;
}

//Calcule la log-vraismeblance avec la fonction de coût
/*************************/
/* computeLoglikelihoodK */
/*************************/
double* GaussianHDDAParameter::computeLoglikelihoodK(double**K) {
	int64_t i, k;
	int64_t nbSample = _model->getNbSample();
	int64_t ** tabZikKnown = _model->getTabZikKnown();

	double* L = new double[_nbCluster];
	for (k = 0; k < _nbCluster; k++) {
		L[k] = 0.0;
	}

	for (i = 0; i < nbSample; i++) {
		for (k = 0; k < _nbCluster; k++) {
			if (tabZikKnown[i][k] == 1) {
				L[k] += K[k][i];
			}
		}
	}
	for (k = 0; k < _nbCluster; k++) {
		L[k] = -L[k] / 2.0;
	}
	return L;
}

/********************/
/* getFreeParameter */
/********************/
int64_t GaussianHDDAParameter::getFreeParameter()const {
	int64_t nbParameter; // Number of parameters
	int64_t k = _nbCluster; // Sample size
	int64_t p = _pbDimension; // Sample dimension

	int64_t roF = k * p + k - 1;
	int64_t roR = k*p;
	int64_t toEqual;
	int64_t toFree = 0;
	int64_t dMean = 0;

	switch (_modelType->_nameModel) {
	case (Gaussian_HD_pk_AkjBkQkDk):
		for (int64_t cls = 0; cls < k; cls++) {
			toFree += _tabDk[cls]*(p - (_tabDk[cls] + 1) / 2);
			dMean += _tabDk[cls];
		}
		toFree /= k;
		dMean /= k;
		nbParameter = roF + k * (toFree + dMean + 2);
		break;

	case (Gaussian_HD_pk_AkBkQkDk):
		for (int64_t cls = 0; cls < k; cls++) {
			toFree += _tabDk[cls]*(p - (_tabDk[cls] + 1) / 2);
			dMean += _tabDk[cls];
		}
		toFree /= k;
		dMean /= k;
		nbParameter = roF + k * (toFree + 3);
		break;

	case (Gaussian_HD_pk_AkjBkQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roF + k * (toEqual + _tabDk[0] + 1) + 1;
		break;
	case (Gaussian_HD_p_AkjBkQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roR + k * (toEqual + _tabDk[0] + 1) + 1;
		break;

	case (Gaussian_HD_pk_AkjBQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roF + k * (toEqual + _tabDk[0]) + 2;
		break;
	case (Gaussian_HD_p_AkjBQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roR + k * (toEqual + _tabDk[0]) + 3;
		break;

	case (Gaussian_HD_pk_AjBkQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roF + k * (toEqual + 1) + 1;
		break;
	case (Gaussian_HD_p_AjBkQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roR + k * (toEqual + 1) + 1;
		break;

	case (Gaussian_HD_pk_AjBQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roF + k * toEqual + 2;
		break;
	case (Gaussian_HD_p_AjBQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roR + k * toEqual + 2;
		break;

	case (Gaussian_HD_pk_AkBkQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roF + k * (toEqual + 2) + 1;
		break;
	case (Gaussian_HD_p_AkBkQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roR + k * (toEqual + 2) + 1;
		break;

	case (Gaussian_HD_pk_AkBQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roF + k * (toEqual + 1) + 2;
		break;
	case (Gaussian_HD_p_AkBQkD):
		toEqual = _tabDk[0]*(p - (_tabDk[0] + 1) / 2);
		nbParameter = roR + k * (toEqual + 1) + 2;
		break;

	case (Gaussian_HD_p_AkjBkQkDk):
		for (int64_t cls = 0; cls < k; cls++) {
			toFree += _tabDk[cls]*(p - (_tabDk[cls] + 1) / 2);
			dMean += _tabDk[cls];
		}
		toFree /= k;
		dMean /= k;
		nbParameter = roR + k * (toFree + dMean + 2);
		break;

	case (Gaussian_HD_p_AkBkQkDk):
		for (int64_t cls = 0; cls < k; cls++) {
			toFree += _tabDk[cls]*(p - (_tabDk[cls] + 1) / 2);
			dMean += _tabDk[cls];
		}
		toFree /= k;
		dMean /= k;
		nbParameter = roR + k * (toFree + 3);
		break;

	default:
		THROW(OtherException, internalMixmodError);
		break;
	}
	return nbParameter;
}

/********/
/* edit */ //debug
/********/
void GaussianHDDAParameter::edit() {
	int64_t k;
	for (k = 0; k < _nbCluster; k++) {
		cout << "\tcomponent : " << k << endl;
		cout << "\t\tproportion : " << _tabProportion[k] << endl;
		editTab(_tabMean + k, 1, _pbDimension, cout, " ", "\t\tmean : ");
		cout << "\tSub dimension : " << _tabDk[k] << endl;
		editTab(_tabAkj + k, 1, _tabDk[k], cout, " ", "\t\t\tParameters Akj : ");
		cout << "\t\t\tParameter Bk : " << _tabBk[k] << endl;
		cout << "\t\tOrientation : " << endl;
		_tabQk[k]->edit(cout, "\t\t\t", " ", _tabDk[k]);

		cout << "\t\tWk : " << endl;
		_tabWk[k]->edit(cout, "\t\t\t");
	}

	cout << "\tW : " << endl;
	_W->edit(cout, "\t\t");
}

/********/
/* edit */
/********/
void GaussianHDDAParameter::edit(std::ofstream & oFile, bool text) {
	int64_t k;

	if (text) {
		for (k = 0; k < _nbCluster; k++) {
			oFile << "\t\t\tComponent " << k + 1 << endl;
			oFile << "\t\t\t---------" << endl;
			oFile << "\t\t\tMixing proportion : " << _tabProportion[k] << endl;
			editTab(_tabMean + k, 1, _pbDimension, oFile, " ", "\t\t\tMean : ");
			oFile << "\t\t\tSub Dimension  : " << _tabDk[k] << endl;
			editTab(_tabAkj + k, 1, _tabDk[k], oFile, " ", "\t\t\tParameters Akj : ");
			oFile << "\t\t\tParameter Bk : " << _tabBk[k] << endl;
			oFile << "\t\t\tOrientation matrix : " << endl;
			_tabQk[k]->edit(oFile, "\t\t\t\t\t", " ", _tabDk[k]);
			oFile << endl;
		}
		oFile << endl;
	}
  else {
    for (k = 0; k < _nbCluster; k++) {
			putDoubleInStream(oFile, _tabProportion[k]);
			editTab(_tabMean + k, 1, _pbDimension, oFile, " ", "");
			oFile << _tabDk[k] << endl;
			editTab(_tabAkj + k, 1, _tabDk[k], oFile, " ", "");
			oFile << _tabBk[k] << endl;
			_tabQk[k]->edit(oFile, "", " ", _tabDk[k]);
			oFile << endl;

		}
		oFile << endl;
	}
}

/****************/
/* computeTabDk */
/****************/
//-------------------------------------------------------------------
//compute parameters for each model
//-------------------------------------------------------------------
void GaussianHDDAParameter::computeTabDk() {
	THROW(InputException, ungivenSubDimension);

	/*XEMSymmetricMatrix* W_k;
	 double* tabNk = _model->getTabNk();
	 int64_tk;
  //   double* tabShapek_store = new double[_pbDimension];
	 for (k=0;k<_nbCluster;k++){
	   if (tabNk[k]<_pbDimension){
		 int64_tnk = (int) tabNk[k];
		 XEMGeneralMatrix * tabQk = new XEMGeneralMatrix(nk);
		  tabQk->resetToZero();
		  _tabQk[k]->resetToZero();
		 //W_k = (XEMSymmetricMatrix* )(_tabGammak[k]);
		 //W_k = _tabGammak[k];
		 //W_k -> computeSVD(_tabShape[k], tabQk);
		 _tabGammak[k] -> computeSVD(_tabShape[k], tabQk);
		 _tabQk[k]->multiply(_Gammak[k], nk, tabQk);
		 delete tabQk;
		 delete _tabGammak[k];
		 _tabGammak[k] = NULL;
	   }
	   else{
		 //W_k = (XEMSymmetricMatrix* )(_tabWk[k]);
		 //W_k = _tabWk[k];
		 //W_k -> computeSVD(_tabShape[k], _tabQk[k]);
		  _tabWk[k] -> computeSVD(_tabShape[k], _tabQk[k]);
	   }
	   double * tabShapek_store = _tabShape[k]->getStore();

		double max_shape = (tabShapek_store[0] - tabShapek_store[1])/tabNk[k];
		int64_ti = 2;
		double seuil = 0.1; // to be decided
		while (i<_pbDimension){
		  double diff = (tabShapek_store[i-1] - tabShapek_store[i])/tabNk[k];

		  int64_tres = 0;
		   if( diff < seuil*max_shape){
			  int64_tr = i+1;
			   while (r < _pbDimension){
				double diff2 = (tabShapek_store[r-1] - tabShapek_store[r])/tabNk[k];
				  if (diff2<seuil*max_shape){
					res+=1;
					r++;
				  }
				  else{
					r = _pbDimension;
				  }
			   } // end while
			   if (res==(_pbDimension-1-i)){
				  _tabDk[k] = i-1;
				  i = _pbDimension;
			   }
			   else{
				  i++;
			   }
		   } // end if diff
		   else{
			i++;
		   }
		} //end while
	 }

	 if ((_modelType!=Gaussian_HD_pk_AkjBkQkDk)&(_modelType!=Gaussian_HD_pk_AkBkQkDk)
	    &(_modelType!=Gaussian_HD_p_AkjBkQkDk) & (_modelType!=Gaussian_HD_p_AkBkQkDk)){
	   int64_td_mean = 0;
	   for (k=0;k<_nbCluster;k++){
		 d_mean += _tabDk[k];
	   }
	   d_mean /= _nbCluster;
	   for (k=0;k<_nbCluster;k++){
		_tabDk[k] = d_mean;
	   }
	 }

  //delete tabShapek_store;*/
}

/*void XEMGaussianHDDAParameter::computeTabDk(){
  double* tabNk = _model->getTabNk();
  int64_tnbSample = _model->getNbSample();
  int64_tk;

  for (k=0;k<_nbCluster;k++){
	_tabDk[k] = 1;
  }

  double nbFreeParameter;// = getFreeParameter();

   switch (_modelType){
	 case Gaussian_HD_pk_AkjBkQkDk:
	 case Gaussian_HD_p_AkjBkQkDk:
	 case Gaussian_HD_pk_AkjBkQkD:
	 case Gaussian_HD_p_AkjBkQkD:
	   computeAkjBkQk();
	   break;
	 case Gaussian_HD_pk_AkBkQkDk:
	 case Gaussian_HD_p_AkBkQkDk:
	 case Gaussian_HD_pk_AkBkQkD:
	 case Gaussian_HD_p_AkBkQkD:
	   computeAkBkQk();
	   break;
	 case Gaussian_HD_pk_AkjBQkD:
	 case Gaussian_HD_p_AkjBQkD:
	   computeAkjBQk();
	   break;
	 case Gaussian_HD_pk_AkBQkD:
	 case Gaussian_HD_p_AkBQkD:
	   computeAkBQk();
	   break;
	 case Gaussian_HD_pk_AjBQkD:
	 case Gaussian_HD_p_AjBQkD:
	   computeAjBQk();
	   break;
	 case Gaussian_HD_pk_AjBkQkD:
	 case Gaussian_HD_p_AjBkQkD:
	   computeAjBkQk();
	   break;
   } // end switch

  double** Cost1 = computeCost(_tabQk);
  double* Lk_min = computeLoglikelihoodK(Cost1);
  double * BIC_min = new double[_nbCluster];
  double* BIC = new double[_nbCluster];
  int64_t* tabD_min = new int[_nbCluster];
  double* Lk;
  double** Cost;

  for (int64_tk=0;k<_nbCluster;k++){
	nbFreeParameter = (_pbDimension+1.0-1.0/_nbCluster)
             + _tabDk[k]*(_pbDimension-(_tabDk[k]+1.0)/2.0)+_tabDk[k]+2.0;
	BIC_min[k] =  (-2*Lk_min[k] + nbFreeParameter * log(tabNk[k])) / (tabNk[k]);
   // cout<<"d  :  "<<1<<"  BIC_min["<<k<<"] :  "<<BIC_min[k]<<endl;
	tabD_min[k] = _tabDk[k];
  }

for (int64_ti=2;i<_pbDimension;i++){

   for (k=0;k<_nbCluster;k++){
	_tabDk[k] = i;
   }
	switch (_modelType){
	 case Gaussian_HD_pk_AkjBkQkDk:
	 case Gaussian_HD_p_AkjBkQkDk:
	 case Gaussian_HD_pk_AkjBkQkD:
	 case Gaussian_HD_p_AkjBkQkD:
	   computeAkjBkQk();
	   break;
	 case Gaussian_HD_pk_AkBkQkDk:
	 case Gaussian_HD_p_AkBkQkDk:
	 case Gaussian_HD_pk_AkBkQkD:
	 case Gaussian_HD_p_AkBkQkD:
	   computeAkBkQk();
	   break;
	 case Gaussian_HD_pk_AkjBQkD:
	 case Gaussian_HD_p_AkjBQkD:
	   computeAkjBQk();
	   break;
	 case Gaussian_HD_pk_AkBQkD:
	 case Gaussian_HD_p_AkBQkD:
	   computeAkBQk();
	   break;
	 case Gaussian_HD_pk_AjBQkD:
	 case Gaussian_HD_p_AjBQkD:
	   computeAjBQk();
	   break;
	case Gaussian_HD_pk_AjBkQkD:
	case Gaussian_HD_p_AjBkQkD:
	  computeAjBkQk();
	  break;

	} // end switch

	Cost = computeCost(_tabQk);
	Lk = computeLoglikelihoodK(Cost);
	for (int64_tk=0;k<_nbCluster;k++){
	 nbFreeParameter = (_pbDimension+1.0-1.0/_nbCluster)
           + _tabDk[k]*(_pbDimension-(_tabDk[k]+1.0)/2.0)+_tabDk[k]+2.0;
	 BIC[k] = (-2*Lk[k] + nbFreeParameter * log(tabNk[k])) / (tabNk[k]);
	 //cout<<"d  : "<<i<<"  BIC["<<k<<"] :  "<<BIC[k]<<endl;
	 if (BIC[k] < BIC_min[k]){
	   tabD_min[k] = i;
	   BIC_min[k] = BIC[k];
	 }
	}
	delete[] Lk;
	Lk = NULL;

   for (int64_tcls=0;cls<_nbCluster;cls++){
	 delete[] Cost[cls];
	 Cost[cls] = NULL;
   }
   delete[] Cost;
   Cost = NULL;
  } // end for i

  if ((_modelType==Gaussian_HD_pk_AkjBkQkDk)||(_modelType==Gaussian_HD_pk_AkBkQkDk)
         ||(_modelType==Gaussian_HD_p_AkjBkQkDk) || (_modelType==Gaussian_HD_p_AkBkQkDk)){
	 for (int64_tk=0;k<_nbCluster;k++){
	   _tabDk[k] = tabD_min[k];
	 }
  }
  else{
   int64_t d_mean = 0;
   for (k=0;k<_nbCluster;k++){
	d_mean += tabD_min[k];
   }
   d_mean /= _nbCluster;
   for (k=0;k<_nbCluster;k++){
	_tabDk[k] = d_mean;
   }
  }

  delete[] Lk_min;
  Lk_min = NULL;
 // delete[] Lk;
 // Lk = NULL;
 delete[] BIC_min;
  BIC_min = NULL;
  delete[] BIC;
  BIC = NULL;
  delete[] tabD_min;
  tabD_min = NULL;

  for (k=0;k<_nbCluster;k++){
   delete[] Cost1[k];
   Cost1[k] = NULL;
  }
  delete[] Cost1;
  Cost1 = NULL;
} // end computetabD
 */

/***************/
/* computeCost */
/***************/
double** GaussianHDDAParameter::computeCost(GeneralMatrix ** tabQ)const {
	double ** K = new double*[_nbCluster];
	int64_t j;
	Parameter * parameter = _model->getGaussianParameter();
	GaussianParameter* gparameter = (GaussianParameter*) parameter;
	double ** tabMean = gparameter->getTabMean();
	double* tabProportion = gparameter->getTabProportion();
	int64_t nbSample = _model->getNbSample();
	GaussianData * data = (GaussianData *) (_model->getData());
	double ** matrix = data->getYStore();
	double * xiMoinsMuk = new double[_pbDimension];
	SymmetricMatrix* Pk = new SymmetricMatrix(_pbDimension); //Id
	SymmetricMatrix* A = new SymmetricMatrix(_pbDimension); //Id
	SymmetricMatrix * Pi = new SymmetricMatrix(_pbDimension); //id

	for (int64_t classe = 0; classe < _nbCluster; classe++) {
		double* tabShapek_store = new double[_pbDimension];
		K[classe] = new double[nbSample];
		//----------------calcul le produit Q*t(Q)---------------

		Pk->compute_as_M_tM(tabQ[classe], _tabDk[classe]);
		//---------------calcul du produit A = Qi_L^-1_t(Qi)------------------
		double sum_lambda = 0.0;
		for (j = 0; j < _tabDk[classe]; j++) {
			tabShapek_store[j] = 1.0 / _tabAkj[classe][j];
			sum_lambda += log(_tabAkj[classe][j]);
		}
		for (j = _tabDk[classe]; j < _pbDimension; j++) {
			tabShapek_store[j] = 0.0;
		}
		A->compute_as_O_S_O(tabQ[classe], tabShapek_store);
		double constante = sum_lambda + (_pbDimension - _tabDk[classe]) * log(_tabBk[classe])
				- 2 * log(tabProportion[classe]) + _pbDimension * log(2 * XEMPI);

		for (int64_t sample = 0; sample < nbSample; sample++) {

			//-----------soustraction des moyennes
			for (j = 0; j < _pbDimension; j++) {
				xiMoinsMuk[j] = matrix[sample][j] - tabMean[classe][j];
			}
			//-----------produit Pi = Qt(Q)(x-muk) (resultat vecteur)----------------------

			Pi->compute_as_M_V(Pk, xiMoinsMuk);
			double * storePi = Pi->getStore();
			//----------------norme(muk-Pi(x))_A = t(muk-Pi(x))_A_(muk-Pi(x))-----------

			double normeA = A->norme(xiMoinsMuk);

			//----------------norme(x-Pi(x))------------------------------
			double norme = 0.0;
			for (int64_t i = 0; i < _pbDimension; i++) {
				storePi[i] += tabMean[classe][i];
				norme += (matrix[sample][i] - storePi[i])*(matrix[sample][i] - storePi[i]);
			}

			//----------------calcul de K(x)---------
			K[classe][sample] = normeA + 1.0 / _tabBk[classe] * norme + constante;
		}//fin sample

		delete[] tabShapek_store;
		tabShapek_store = NULL;
	}//fin classe

	delete Pk;
	delete A;

	delete Pi;
	//delete[] storePi;

	delete[] xiMoinsMuk;
	xiMoinsMuk = NULL;

	return K;
}

/******************/
/* computeAkjBkQk */
/*****************/
// Dans mixmod, les matrices Wk ne sont pas divisees par nk alors que dans
//le these de Charles elles le sont. L'expression des parametres dans le code
//est differente de celle de la these car elle rectifie cette difference
void GaussianHDDAParameter::computeAkjBkQk() {
	Matrix* W_k;
	double* tabNk = _model->getTabNk();

	for (int64_t k = 0; k < _nbCluster; k++) {
		double sum_lambda = 0.0;
		if (tabNk[k] < _pbDimension) {
			int64_t nk = (int64_t) tabNk[k];
			GeneralMatrix * tabQk = new GeneralMatrix(nk); //Id
			W_k = _tabGammak[k];
			W_k -> computeSVD(_tabShape[k], tabQk);
			_tabQk[k]->multiply(_Gammak[k], nk, tabQk);
			delete tabQk;
			delete _tabGammak[k];
			_tabGammak[k] = NULL;
		}
		else {
			W_k = _tabWk[k];
			W_k -> computeSVD(_tabShape[k], _tabQk[k]);
		}
		double * storeShapek = _tabShape[k]->getStore();
		for (int64_t j = 0; j < _tabDk[k]; j++) {
			_tabAkj[k][j] = storeShapek[j] / tabNk[k];
			sum_lambda += _tabAkj[k][j];
			//cout<<"tabA :  "<<_tabAkj[k][j]<<endl;
		}
		double trace = W_k->computeTrace();
		_tabBk[k] = 1.0 / (_pbDimension - _tabDk[k])*(trace / tabNk[k] - sum_lambda);
		//cout<<"tabB :  "<<_tabBk[k]<<endl;
	}
}

/*****************/
/* computeAkjBQk */
/*****************/
void GaussianHDDAParameter::computeAkjBQk() {
	DiagMatrix* tabShapeW = new DiagMatrix(_pbDimension); //Id
	GeneralMatrix* tabQW = new GeneralMatrix(_pbDimension); //Id
	double* tabNk = _model->getTabNk();

	_W -> computeSVD(tabShapeW, tabQW);
	double trace = _W->computeTrace() / _model->getNbSample();
	double somme = 0.0;
	for (int64_t k = 0; k < _nbCluster; k++) {
		double sum_lambda = 0.0;
		if (tabNk[k] < _pbDimension) {
			int64_t nk = (int64_t) tabNk[k];
			GeneralMatrix * tabQk = new GeneralMatrix(nk); //Id
			_tabGammak[k]-> computeSVD(_tabShape[k], tabQk);
			_tabQk[k]->multiply(_Gammak[k], nk, tabQk);
			delete tabQk;
		}
		else {
			_tabWk[k] -> computeSVD(_tabShape[k], _tabQk[k]);
		}
		double * storeShapek = _tabShape[k]->getStore();
		for (int64_t j = 0; j < _tabDk[k]; j++) {
			_tabAkj[k][j] = storeShapek[j] / tabNk[k];
			sum_lambda += _tabAkj[k][j];
			//cout<<"tabA :  "<<_tabA[k][j]<<endl;
		}
		somme += tabNk[k] * sum_lambda;
	}
	somme = somme / _model->getNbSample();
	for (int64_t k = 0; k < _nbCluster; k++) {
		_tabBk[k] = 1.0 / (_pbDimension - _tabDk[k])*(trace - somme);
		//cout<<"tabB  :  "<<_tabB[k]<<endl;
	}
	delete tabShapeW;
	delete tabQW;
}

/******************/
/* computeAjBkQk */
/*****************/
void GaussianHDDAParameter::computeAjBkQk() {
	Matrix * W_k;
	DiagMatrix* tabShapeW = new DiagMatrix(_pbDimension); //Id
	GeneralMatrix* tabQW = new GeneralMatrix(_pbDimension); //Id
	double sum_lambda;
	double* tabNk = _model->getTabNk();
	_W -> computeSVD(tabShapeW, tabQW);
	double * storeShape = tabShapeW->getStore();
	for (int64_t k = 0; k < _nbCluster; k++) {
		sum_lambda = 0.0;
		if (tabNk[k] < _pbDimension) {
			int64_t nk = (int64_t) tabNk[k];
			GeneralMatrix * tabQk = new GeneralMatrix(nk); //Id
			W_k = _tabGammak[k];
			W_k -> computeSVD(_tabShape[k], tabQk);
			_tabQk[k]->multiply(_Gammak[k], nk, tabQk);
			delete tabQk;
		}
		else {
			W_k = _tabWk[k];
			W_k->computeSVD(_tabShape[k], _tabQk[k]);
		}
		double * storeShapek = _tabShape[k]->getStore();
		for (int64_t j = 0; j < _tabDk[k]; j++) {
			_tabAkj[k][j] = storeShape[j] / _model->getNbSample();
			sum_lambda += storeShapek[j] / tabNk[k];
			//cout<<"tabAij :  "<<_tabA[k][j]<<endl;
		}
		double trace = W_k->computeTrace() / tabNk[k];
		_tabBk[k] = 1.0 / (_pbDimension - _tabDk[k])*(trace - sum_lambda);
		//cout<<"tabB :  "<<_tabB[k]<<endl;
	}
	delete tabShapeW;
	delete tabQW;
}

/****************/
/* computeAjBQk */
/***************/
void GaussianHDDAParameter::computeAjBQk() {
	double somme = 0.0;
	DiagMatrix* tabShapeW = new DiagMatrix(_pbDimension); //Id
	GeneralMatrix* tabQW = new GeneralMatrix(_pbDimension); //Id
	double* tabNk = _model->getTabNk();

	double trace = _W->computeTrace() / _model->getNbSample();
	_W -> computeSVD(tabShapeW, tabQW);
	double * storeShape = tabShapeW->getStore();
	for (int64_t k = 0; k < _nbCluster; k++) {
		double sum_lambda = 0.0;
		if (tabNk[k] < _pbDimension) {
			int64_t nk = (int64_t) tabNk[k];
			GeneralMatrix * tabQk = new GeneralMatrix(nk); //Id
			_tabGammak[k] -> computeSVD(_tabShape[k], tabQk);
			_tabQk[k]->multiply(_Gammak[k], nk, tabQk);
			delete tabQk;
		}
		else {
			_tabWk[k] -> computeSVD(_tabShape[k], _tabQk[k]);
		}
		double * storeShapek = _tabShape[k]->getStore();
		for (int64_t j = 0; j < _tabDk[k]; j++) {
			_tabAkj[k][j] = storeShape[j] / _model->getNbSample();
			sum_lambda += storeShapek[j];
			//cout<<"tabA :  "<<_tabA[k][j]<<endl;
		}
		somme += sum_lambda;
	}
	somme = somme / _model->getNbSample();
	for (int64_t k = 0; k < _nbCluster; k++) {
		_tabBk[k] = 1.0 / (_pbDimension - _tabDk[k])*(trace - somme);
		// cout<<"tabB :  "<<_tabB[k]<<endl;
	}
	delete tabShapeW;
	delete tabQW;
}

/*****************/
/* computeAkBkQk */
/****************/
void GaussianHDDAParameter::computeAkBkQk() {
	double sum_lambda;
	double* tabNk = _model->getTabNk();
	Matrix * W_k;
	for (int64_t k = 0; k < _nbCluster; k++) {
		sum_lambda = 0.0;
		if (tabNk[k] < _pbDimension) {
			int64_t nk = (int64_t) tabNk[k];
			GeneralMatrix * tabQk = new GeneralMatrix(nk); //Id
			W_k = _tabGammak[k];
			W_k -> computeSVD(_tabShape[k], tabQk);
			_tabQk[k]->multiply(_Gammak[k], nk, tabQk);
			delete tabQk;
		}
		else {
			W_k = _tabWk[k];
			W_k -> computeSVD(_tabShape[k], _tabQk[k]);
		}
		double * storeShapek = _tabShape[k]->getStore();
		for (int64_t j = 0; j < _tabDk[k]; j++) {
			sum_lambda += storeShapek[j] / tabNk[k];
		}
		for (int64_t j = 0; j < _tabDk[k]; j++) {
			_tabAkj[k][j] = 1.0 / _tabDk[k] * sum_lambda;
			//cout<<"tabA : "<<_tabA[k][j]<<endl;
		}
		double trace = W_k->computeTrace() / tabNk[k];
		_tabBk[k] = 1.0 / (_pbDimension - _tabDk[k]) * (trace - sum_lambda);
		//cout<<"tabB : "<<_tabB[k]<<endl;
	}
}

/****************/
/* computeAkBQk */
/****************/
void GaussianHDDAParameter::computeAkBQk() {
	double sum_lambda;
	double somme = 0.0;
	DiagMatrix* tabShapeW = new DiagMatrix(_pbDimension); //Id
	GeneralMatrix* tabQW = new GeneralMatrix(_pbDimension); //Id
	double* tabNk = _model->getTabNk();

	double trace = _W->computeTrace() / _model->getNbSample();
	_W -> computeSVD(tabShapeW, tabQW);
	for (int64_t k = 0; k < _nbCluster; k++) {
		sum_lambda = 0.0;
		if (tabNk[k] < _pbDimension) {
			int64_t nk = (int64_t) tabNk[k];
			GeneralMatrix * tabQk = new GeneralMatrix(nk); //Id
			_tabGammak[k] -> computeSVD(_tabShape[k], tabQk);
			_tabQk[k]->multiply(_Gammak[k], nk, tabQk);
			delete tabQk;
		}
		else {
			_tabWk[k] -> computeSVD(_tabShape[k], _tabQk[k]);
		}
		double * storeShapek = _tabShape[k]->getStore();
		for (int64_t j = 0; j < _tabDk[k]; j++) {
			sum_lambda += storeShapek[j] / tabNk[k];
		}
		for (int64_t j = 0; j < _tabDk[k]; j++) {
			_tabAkj[k][j] = 1.0 / _tabDk[k] * sum_lambda;
			//cout<<"tabA :  "<<_tabA[k][j]<<endl;
		}
		somme += tabNk[k] * sum_lambda;
	}
	somme = somme / _model->getNbSample();
	for (int64_t k = 0; k < _nbCluster; k++) {
		_tabBk[k] = 1.0 / (_pbDimension - _tabDk[k])*(trace - somme);
		//cout<<"tabB :  "<<_tabB[k]<<endl;
	}
	delete tabShapeW;
	delete tabQW;
}

}
