/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/BinaryEkjParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/BinaryEkjParameter.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/BinarySample.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Utilities/Random.h"

namespace XEM {

//--------------------
// Default constructor
//--------------------
BinaryEkjParameter::BinaryEkjParameter() {
	THROW(OtherException, wrongConstructorType);
}

//-------------------------------
// Constructor called by XEMModel
//-------------------------------
BinaryEkjParameter::BinaryEkjParameter(
		Model * iModel, ModelType * iModelType, int64_t * tabNbModality) 
: BinaryParameter(iModel, iModelType, tabNbModality) 
{
	_scatter = new double*[_nbCluster];
	for (int64_t k = 0; k < _nbCluster; k++) {
		_scatter[k] = new double[_pbDimension];
		for (int64_t j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = 0.0;
		}
	}
}

//-----------------
// copy Constructor
//-----------------
BinaryEkjParameter::BinaryEkjParameter(const BinaryEkjParameter * iParameter) 
: BinaryParameter(iParameter) 
{
	_scatter = new double*[_nbCluster];
	for (int64_t k = 0; k < _nbCluster; k++) {
		_scatter[k] = new double[_pbDimension];
	}
	double ** iScatter = iParameter->getScatter();
	for (int64_t k = 0; k < _nbCluster; k++) {
		for (int64_t j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = iScatter[k][j];
		}
	}
}

//---------
// clone 
//---------
Parameter * BinaryEkjParameter::clone() const {
	BinaryEkjParameter * newParam = new BinaryEkjParameter(this);
	return (newParam);
}

//-----------
// Destructor
//-----------
BinaryEkjParameter::~BinaryEkjParameter() {
	if (_scatter) {
		for (int64_t k = 0; k < _nbCluster; k++) {
			delete [] _scatter[k];
		}
	}
	delete [] _scatter;
	_scatter = NULL;
}

//---------------------
/// Comparison operator
//---------------------
bool BinaryEkjParameter::operator ==(const BinaryEkjParameter & param) const {
	if (!BinaryParameter::operator==(param)) return false;
	for (int64_t k = 0; k < _nbCluster; k++) {
		for (int64_t j = 0; j < _pbDimension; j++) {
			if (_scatter[k][j] != param.getScatter()[k][j]) return false;
		}
	}
	return true;
}

//------------------------
// reset to default values
//------------------------
void BinaryEkjParameter::reset() {
	int64_t k, j;
	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = 0.0;
		}
	}
	BinaryParameter::reset();
}

//-----------
// getFreeParameter
//-----------
int64_t BinaryEkjParameter::getFreeParameter() const {
	int64_t nbFreeParameter = _pbDimension * _nbCluster;
	if (_freeProportion) {
		nbFreeParameter += _nbCluster - 1;
	}
	return nbFreeParameter;
}

//-------
// getPdf
//-------
double BinaryEkjParameter::getPdf(int64_t iSample, int64_t kCluster) const {
	//cout<<" XEMBinaryEkjParameter::getPdf"<<endl;
	int64_t j;
	double bernPdf = 1.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		//cout<<"curSample :  "<<curSample->getDataValue(j)<<endl;
		//cout<<" _tabCenter[kCluster][j] :  "<< _tabCenter[kCluster][j]<<endl;
		// iSample have major modality ?//
		if (curSample->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter[kCluster][j];
		}
		else {
			bernPdf *= _scatter[kCluster][j] / (_tabNbModality[j] - 1.0);
		}
	}
	return bernPdf;
}

//----------
// getLogPdf
//----------
long double BinaryEkjParameter::getLogPdf(int64_t iSample, int64_t kCluster) const {
	//cout<<" XEMBinaryEkjParameter::getPdf"<<endl;
	int64_t j;
	double bernPdf = 0.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		//cout<<"curSample :  "<<curSample->getDataValue(j)<<endl;
		//cout<<" _tabCenter[kCluster][j] :  "<< _tabCenter[kCluster][j]<<endl;
		// iSample have major modality ?//
		if (curSample->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf += log(1.0 - _scatter[kCluster][j]);
		}
		else {
			bernPdf += log(_scatter[kCluster][j] / (_tabNbModality[j] - 1.0));
		}
	}
	return bernPdf;
}

//-------
// getPdf
//-------
/* Compute normal probability density function
	   for x vector and kCluster th cluster
 */
double BinaryEkjParameter::getPdf(Sample * x, int64_t kCluster) const {
	int64_t j;
	double bernPdf = 1.0;
	BinarySample * binaryX = x->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ? //
		if (binaryX->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter[kCluster][j];
		}
		else {
			bernPdf *= _scatter[kCluster][j] / (_tabNbModality[j] - 1.0);
		}
	}
	return bernPdf;
}

//--------------------
// getlogLikelihoodOne (one cluster)
//--------------------
double BinaryEkjParameter::getLogLikelihoodOne() const {
	int64_t i;
	int64_t j;
	double logLikelihoodOne = 0.0, pdf;//, * Scatter;
	//Scatter = new double[_pbDimension];
	std::unique_ptr<double[]> Scatter(new double[_pbDimension]);    
	//int64_t * Center = new int64_t[_pbDimension];
	std::unique_ptr<int64_t[]> Center(new int64_t[_pbDimension]);    
	//double * tabNbSampleInMajorModality = new double[_pbDimension];
	std::unique_ptr<double[]> tabNbSampleInMajorModality(new double[_pbDimension]);    
	int64_t nbSample = _model->getNbSample();
	BinaryData * data = _model->getBinaryData();

	// Compute Center fo One cluster //
	getTabCenterIfOneCluster(Center.get(), tabNbSampleInMajorModality.get());

	// Compute Scatter for One cluster //
	for (j = 0; j < _pbDimension; j++)
		//Scatter[j] = 1- (tabNbSampleInMajorModality[j] / data->_weightTotal);
		Scatter[j] = 1 - ((tabNbSampleInMajorModality[j] 
				+ 1. / _tabNbModality[j]) / (data->_weightTotal + 1));

	// Compute the log-likelihood for one cluster (k=1) //
	//--------------------------------------------------//
	for (i = 0; i < nbSample; i++) {
      pdf = computePdfOneCluster(data->_matrix[i], Center.get(), Scatter.get(), _tabNbModality);
		logLikelihoodOne += log(pdf) * data->_weight[i];
	}

	//delete[] Center;
	//delete[] Scatter;
	//delete[] tabNbSampleInMajorModality;

	return logLikelihoodOne;
}

//----------------
// Compute scatter 
//----------------
void BinaryEkjParameter::computeScatter() {
	int64_t j, k;
	int64_t i;
	double ekj; // nb d'individus de la classe k prenant la modalite maj sur la variable j
	double * tabNk = _model->getTabNk();
	double ** tabCik = _model->getTabCik();

	BinaryData * data = _model->getBinaryData();
	Sample ** dataMatrix = data->getDataMatrix();
	BinarySample * curSample;
	int64_t nbSample = _model->getNbSample();

	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			ekj = 0.0;
			for (i = 0; i < nbSample; i++) {
				curSample = dataMatrix[i]->getBinarySample();
				if (curSample->getDataValue(j) == _tabCenter[k][j]) {
					ekj += (tabCik[i][k] * data->_weight[i]);
				}
			}
			//_scatter[k][j] = 1 - (ekj / tabNk[k]);
			_scatter[k][j] = 1 - ((ekj + 1. / _tabNbModality[j]) / (tabNk[k] + 1));
		} // end for j
	} // end for k
}

//--------------------------------------------------
// Compute scatter(s)  as if there was only one cluster
//---------------------------------------------------
/*void XEMBinaryEkjParameter::computeScatterIfOneCluster(double totalWeight, 
      double * tabNbSampleInMajorModality, double ** tabNbSamplePerModality){
  for (int64_t k=0; k<_nbCluster; k++){
	for (int64_t j=0; j<_pbDimension; j++){
	  _scatter[k][j] = 1 - (tabNbSampleInMajorModality[j] / totalWeight);
	}
  }
}*/

//--------------------------
// Compute random scatter(s)
//--------------------------
void BinaryEkjParameter::computeRandomScatter() {
	for (int64_t k = 0; k < _nbCluster; k++) {
		// tirage d'une valeur comprise entre 0 et 1./_tabNbModality[j] 
		for (int64_t j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = rnd() / _tabNbModality[j];
		}
	}
}

//---------------
//recopy scatter from iParam 
//---------------
// Note : iParam must be a XEMBinaryEkjParameter*
void BinaryEkjParameter::recopyScatter(Parameter * iParam) {
	if (typeid (*iParam) != typeid (*this)) {
		THROW(OtherException, badBinaryParameterClass);
	}
	double ** iScatter = ((BinaryEkjParameter*) iParam)->getScatter();
	for (int64_t k = 0; k < _nbCluster; k++) {
		for (int64_t j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = iScatter[k][j];
		}
	}
}

//---------------
//create scatter from Scatter Ekjh 
//---------------
void BinaryEkjParameter::createScatter(double *** scatter) {
	int64_t k, j, h;
	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			h = _tabCenter[k][j];
			_scatter[k][j] = scatter[k][j][h - 1];
		}
	}
}

//------------
// editScatter (for debug)
//------------
void BinaryEkjParameter::editScatter(int64_t k) {
	int64_t j, h;
	for (j = 0; j < _pbDimension; j++) {
		for (h = 1; h <= _tabNbModality[j]; h++) {
			if (h == _tabCenter[k][j]) {
				cout << "\t" << _scatter[k][j];
			}
			else {
				cout << "\t" << _scatter[k][j] / (_tabNbModality[j] - 1);
			}
		}
		cout << endl;
	}
}

// editScatter 
//------------
void BinaryEkjParameter::editScatter(std::ofstream & oFile, int64_t k, bool text) {
	int64_t j, h;
	if (text) {
		oFile << "\t\t\tScattering : \n";
	}
	for (j = 0; j < _pbDimension; j++) {
		if (text) {
			oFile << "\t\t\t\t\t";
			;
		}
		for (h = 1; h <= _tabNbModality[j]; h++) {
			if (h == _tabCenter[k][j])
				putDoubleInStream(oFile, _scatter[k][j], "  ");
      else
				putDoubleInStream(oFile, _scatter[k][j] / (_tabNbModality[j] - 1), "  ");
		}
		oFile << endl;
	}
}

// Read Scatter in input file
//---------------------------
void BinaryEkjParameter::inputScatter(std::ifstream & fi, int64_t k) {
	THROW(OtherException, internalMixmodError);
}

// Read Scatter in input containers
//---------------------------
void BinaryEkjParameter::inputScatter(double *** scatters) {
	THROW(OtherException, internalMixmodError);
}

double *** BinaryEkjParameter::scatterToArray() const {
	int64_t k, j, h;
	double *** tabScatter = new double**[_nbCluster];
	for (k = 0; k < _nbCluster; k++) {
		tabScatter[k] = new double*[_pbDimension];
		for (j = 0; j < _pbDimension; j++) {
			tabScatter[k][j] = new double[_tabNbModality[j]];
			for (h = 1; h <= _tabNbModality[j]; h++) {
				if (h == _tabCenter[k][j]) {
					tabScatter[k][j][h - 1] = _scatter[k][j];
				}
				else {
					tabScatter[k][j][h - 1] = _scatter[k][j] / (_tabNbModality[j] - 1);
				}
			}
		}
	}
	return tabScatter;
}

}
