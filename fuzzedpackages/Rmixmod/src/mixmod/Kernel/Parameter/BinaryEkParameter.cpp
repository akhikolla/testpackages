/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/BinaryEkParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/BinaryEkParameter.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/BinarySample.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Utilities/Random.h"

namespace XEM {

//--------------------
// Default constructor
//--------------------

BinaryEkParameter::BinaryEkParameter() {
	THROW(OtherException, wrongConstructorType);
}

//-------------------------------
// Constructor called by XEMModel
//-------------------------------
BinaryEkParameter::BinaryEkParameter(
		Model * iModel, ModelType * iModelType, int64_t * tabNbModality) 
: BinaryParameter(iModel, iModelType, tabNbModality) 
{
	_scatter = new double[_nbCluster];
	int64_t k;
	for (k = 0; k < _nbCluster; k++) {
		_scatter[k] = 0.0;
	}
}

//-----------------
// copy Constructor
//-----------------
BinaryEkParameter::BinaryEkParameter(const BinaryEkParameter * iParameter) 
: BinaryParameter(iParameter) 
{
	_scatter = new double[_nbCluster];
	double * iScatter = iParameter->getScatter();
	for (int64_t k = 0; k < _nbCluster; k++) {
		_scatter[k] = iScatter[k];
	}
}

//---------
// clone 
//---------
Parameter * BinaryEkParameter::clone() const {
	BinaryEkParameter * newParam = new BinaryEkParameter(this);
	return (newParam);
}

//-----------
// Destructor
//-----------
BinaryEkParameter::~BinaryEkParameter() {
	if (_scatter) {
		delete [] _scatter;
	}
}

//---------------------
/// Comparison operator
//---------------------
bool BinaryEkParameter::operator ==(const BinaryEkParameter & param) const {
	if (!BinaryParameter::operator==(param)) return false;
	for (int64_t k = 0; k < _nbCluster; k++) {
		if (_scatter[k] != param.getScatter()[k]) return false;
	}
	return true;
}

//------------------------
// reset to default values
//------------------------
void BinaryEkParameter::reset() {
	int64_t k;
	for (k = 0; k < _nbCluster; k++) {
		_scatter[k] = 0.0;
	}
	BinaryParameter::reset();
}

//-----------
// getFreeParameter
//-----------
int64_t BinaryEkParameter::getFreeParameter() const {
	int64_t nbFreeParameter = _nbCluster;
	if (_freeProportion) {
		nbFreeParameter += _nbCluster - 1;
	}
	return nbFreeParameter;
}

//-------
// getPdf
//-------
double BinaryEkParameter::getPdf(int64_t iSample, int64_t kCluster) const {
	int64_t j;
	double bernPdf = 1.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ?//
		if (curSample->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter[kCluster];
		}
		else {
			bernPdf *= _scatter[kCluster] / (_tabNbModality[j] - 1.0);
		}
	}
	return bernPdf;
}

//----------
// getLogPdf
//----------
long double BinaryEkParameter::getLogPdf(int64_t iSample, int64_t kCluster) const {
	int64_t j;
	double bernPdf = 0.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ?//
		if (curSample->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf += log(1.0 - _scatter[kCluster]);
		}
		else {
			bernPdf += log(_scatter[kCluster] / (_tabNbModality[j] - 1.0));
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
double BinaryEkParameter::getPdf(Sample * x, int64_t kCluster) const {
	int64_t j;
	double bernPdf = 1.0;
	BinarySample * binaryX = x->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ? //
		if (binaryX->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter[kCluster];
		}
		else {
			bernPdf *= _scatter[kCluster] / (_tabNbModality[j] - 1.0);
		}
	}
	return bernPdf;
}

//--------------------
// getlogLikelihoodOne (one cluster)
//--------------------
double BinaryEkParameter::getLogLikelihoodOne() const {
	int64_t i;
	int64_t j;
	double logLikelihoodOne = 0.0, value, pdf, Scatter;
	//int64_t * Center = new int64_t[_pbDimension];
    std::unique_ptr<int64_t[]> Center(new int64_t[_pbDimension]);
	//double * tabNbSampleInMajorModality = new double[_pbDimension];
    std::unique_ptr<double[]> tabNbSampleInMajorModality(new double[_pbDimension]);
	int64_t nbSample = _model->getNbSample();
	BinaryData * data = _model->getBinaryData();

	// Compute Center fo One cluster //
	getTabCenterIfOneCluster(Center.get(), tabNbSampleInMajorModality.get());

	// Compute Scatter for One cluster //
	value = 0.0;
	for (j = 0; j < _pbDimension; j++) {
		//value += tabNbSampleInMajorModality[j];
		value += tabNbSampleInMajorModality[j] + 1. / _tabNbModality[j];
	}
	//Scatter = 1- ( value /(data->_weightTotal * _pbDimension));
	Scatter = 1 - (value / ((data->_weightTotal + 1) * _pbDimension));

	// Compute the log-likelihood for one cluster (k=1) //
	//--------------------------------------------------//
	for (i = 0; i < nbSample; i++) {
      pdf = computePdfOneCluster(data->_matrix[i], Center.get(), Scatter, _tabNbModality);
		logLikelihoodOne += log(pdf) * data->_weight[i];
	}

	//delete[] Center;
	//delete[] tabNbSampleInMajorModality;

	return logLikelihoodOne;
}

//----------------
// Compute scatter 
//----------------
void BinaryEkParameter::computeScatter() {
	int64_t j, k;
	int64_t i;
	double ek; // nb d'individus de la classe k prenant la modalite 
	           // majoritaire (quelle que soit la variable)
	double * tabNk = _model->getTabNk();
	double ** tabCik = _model->getTabCik();

	BinaryData * data = _model->getBinaryData();
	Sample ** dataMatrix = data->getDataMatrix();
	BinarySample * curSample = NULL;
	int64_t nbSample = _model->getNbSample();

	for (k = 0; k < _nbCluster; k++) {
		ek = 0.0;
		for (j = 0; j < _pbDimension; j++) {
			for (i = 0; i < nbSample; i++) {
				curSample = dataMatrix[i]->getBinarySample();
				if (curSample->getDataValue(j) == _tabCenter[k][j]) {
					ek += (tabCik[i][k] * data->_weight[i]);
				}
			}// end for i
			ek += 1. / _tabNbModality[j];
		} // end for j
		_scatter[k] = 1 - (ek / ((tabNk[k] + 1) * _pbDimension));
		//_scatter[k] = 1 - (ek / ((tabNk[k] )*_pbDimension));
	} // end for k
}

//--------------------------
// Compute random scatter(s)   
//--------------------------
void BinaryEkParameter::computeRandomScatter() {
	int64_t minNbModality = _tabNbModality[0];
	for (int64_t j = 1; j < _pbDimension; j++) {
		if (_tabNbModality[j] < minNbModality) {
			minNbModality = _tabNbModality[j];
		}
	}

	// tirage d'une valeur comprise entre 0 et 1/minNbModality
	for (int64_t k = 0; k < _nbCluster; k++) {
		_scatter[k] = rnd() / minNbModality;
	}
}

//---------------
//recopy scatter from iParam 
//---------------
// Note : iParam must be a XEMBinaryEkParameter*
void BinaryEkParameter::recopyScatter(Parameter * iParam) {
	if (typeid (*iParam) != typeid (*this)) {
		THROW(OtherException, badBinaryParameterClass);
	}
	double * iScatter = ((BinaryEkParameter*) iParam)->getScatter();
	for (int64_t k = 0; k < _nbCluster; k++) {
		_scatter[k] = iScatter[k];
	}
}

//---------------
//create scatter from Scatter Ekjh 
//---------------
//on fait une moyenne
void BinaryEkParameter::createScatter(double *** scatter) {
	int64_t k, j, h;
	for (k = 0; k < _nbCluster; k++) {
		_scatter[k] = 0.0;
		for (j = 0; j < _pbDimension; j++) {
			h = _tabCenter[k][j];
			_scatter[k] += scatter[k][j][h - 1];
		}
		_scatter[k] /= _pbDimension;
	}
}

//------------
// editScatter (for debug)
//------------
void BinaryEkParameter::editScatter(int64_t k) {
	int64_t j, h;
	for (j = 0; j < _pbDimension; j++) {
		for (h = 1; h <= _tabNbModality[j]; h++) {
			if (h == _tabCenter[k][j]) {
				cout << "\t" << _scatter[k];
			}
			else {
				cout << "\t" << _scatter[k] / (_tabNbModality[j] - 1);
			}
		}
		cout << endl;
	}
}

// editScatter 
//------------
void BinaryEkParameter::editScatter(std::ofstream & oFile, int64_t k, bool text) {
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
				putDoubleInStream(oFile, _scatter[k], "  ");
      else
				putDoubleInStream(oFile, _scatter[k] / (_tabNbModality[j] - 1), "  ");
    }
    oFile << endl;
  }
}

// Read Scatter in input file
//---------------------------
void BinaryEkParameter::inputScatter(std::ifstream & fi, int64_t k) {
	THROW(OtherException, internalMixmodError);
}

// Read Scatter in input containers
//---------------------------
void BinaryEkParameter::inputScatter(double *** scatters) {
	THROW(OtherException, internalMixmodError);
}

double *** BinaryEkParameter::scatterToArray() const {
	int64_t k, j, h;
	double *** tabScatter = new double**[_nbCluster];
	for (k = 0; k < _nbCluster; k++) {
		tabScatter[k] = new double*[_pbDimension];
		for (j = 0; j < _pbDimension; j++) {
			tabScatter[k][j] = new double[_tabNbModality[j]];
			for (h = 1; h <= _tabNbModality[j]; h++) {
				if (h == _tabCenter[k][j]) {
					tabScatter[k][j][h - 1] = _scatter[k];
				}
				else {
					tabScatter[k][j][h - 1] = _scatter[k] / (_tabNbModality[j] - 1);
				}
			}
		}
	}
	return tabScatter;
}

}
