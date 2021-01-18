/***************************************************************************
                             SRC/mixmod/Kernel/IO/GaussianData.cpp  description
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

#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/GaussianSample.h"
#include "mixmod/Kernel/IO/DataDescription.h"
#include "mixmod/Kernel/IO/WeightColumnDescription.h"
#include "mixmod/Kernel/IO/QuantitativeColumnDescription.h"
#include "mixmod/Kernel/IO/IndividualColumnDescription.h"
#include <clocale> //for setlocale()

namespace XEM {

//-----------
//Constructor
//-----------
GaussianData::GaussianData() {
}

//------------
// Constructor
//------------
GaussianData::GaussianData(const GaussianData & iData) : Data(iData) {
	int64_t i;
	Sample ** matrix = iData._matrix;

	_matrix = new Sample*[_nbSample];
	_yStore = new double* [_nbSample];

	for (i = 0; i < _nbSample; i++) {
		_matrix[i] = new GaussianSample(matrix[i]->getGaussianSample());
		_yStore[i] = (_matrix[i]->getGaussianSample())->getTabValue();
	}
	_Inv2PiPow = iData.getInv2PiPow();
	_pbDimensionLog2Pi = iData.getPbDimensionLog2Pi();
	_halfPbDimensionLog2Pi = _pbDimensionLog2Pi / 2.0;
	__tmpTabOfSizePbDimension = new double[_pbDimension];
	_deleteSamples = true;
}

//------------
// Constructor
//------------
GaussianData::GaussianData(int64_t nbSample, int64_t pbDimension) : Data(nbSample, pbDimension) {
	int64_t i;
	_Inv2PiPow = 1.0 / pow(2.0 * XEMPI, pbDimension / 2.0);
	_pbDimensionLog2Pi = pbDimension * log(2.0 * XEMPI);
	_halfPbDimensionLog2Pi = _pbDimensionLog2Pi / 2.0;
	__tmpTabOfSizePbDimension = new double[_pbDimension];

	_matrix = new Sample*[_nbSample];
	_yStore = new double * [_nbSample];

	for (i = 0; i < _nbSample; i++) {
		_weight[i] = 1.0;
		_matrix[i] = new GaussianSample(_pbDimension);
		_yStore[i] = ((_matrix[i])->getGaussianSample())->getTabValue();
	}
	_weightTotal = _nbSample;
}

//------------
// Constructor
//------------
GaussianData::GaussianData(int64_t nbSample, int64_t pbDimension, double ** matrix) 
: Data(nbSample, pbDimension) 
{
	int64_t i;

	if (matrix == NULL) THROW(OtherException, internalMixmodError);

	_Inv2PiPow = 1.0 / pow(2.0 * XEMPI, pbDimension / 2.0);
	_pbDimensionLog2Pi = pbDimension * log(2.0 * XEMPI);
	_halfPbDimensionLog2Pi = _pbDimensionLog2Pi / 2.0;
	__tmpTabOfSizePbDimension = new double[_pbDimension];

	_matrix = new Sample*[_nbSample];
	_yStore = new double * [_nbSample];

	for (i = 0; i < _nbSample; i++) {
		_weight[i] = 1.0;
		_matrix[i] = new GaussianSample(_pbDimension, matrix[i]);
		_yStore[i] = ((_matrix[i])->getGaussianSample())->getTabValue();
	}
	_deleteSamples = true;
	_weightTotal = _nbSample;
}

//------------
// Constructor
//------------
GaussianData::GaussianData(int64_t nbSample, int64_t pbDimension, const std::string & dataFileName) 
: Data(nbSample, pbDimension) 
{
	int64_t i;

	_Inv2PiPow = 1.0 / pow(2.0 * XEMPI, pbDimension / 2.0);
	_pbDimensionLog2Pi = pbDimension * log(2.0 * XEMPI);
	_halfPbDimensionLog2Pi = _pbDimensionLog2Pi / 2.0;
	__tmpTabOfSizePbDimension = new double[_pbDimension];

	_matrix = new Sample*[_nbSample];
	_yStore = new double * [_nbSample];

	for (i = 0; i < _nbSample; i++) {
		_matrix[i] = new GaussianSample(_pbDimension);
		_yStore[i] = ((_matrix[i])->getGaussianSample())->getTabValue();
	}

	std::ifstream dataStream((dataFileName).c_str(), ios::in);
	if (!dataStream.is_open()) {
		THROW(InputException, wrongDataFileName);
	}
	input(dataStream);
	dataStream.close();
	_deleteSamples = true;
	_fileNameData = dataFileName;
}

//------------
// Constructor for dataReduce
//------------
GaussianData::GaussianData(int64_t nbSample, int64_t pbDimension, 
		double weightTotal, Sample **& matrix, double * weight) 
: Data(nbSample, pbDimension, weightTotal, weight) 
{
	// 1/ (2 * pi)^(d/2)
	_Inv2PiPow = 1.0 / pow(2.0 * XEMPI, pbDimension / 2.0);
	_pbDimensionLog2Pi = pbDimension * log(2.0 * XEMPI);
	_halfPbDimensionLog2Pi = _pbDimensionLog2Pi / 2.0;
	_pbDimensionLog2Pi = pbDimension * log(2.0 * XEMPI);
	__tmpTabOfSizePbDimension = new double[_pbDimension];

	_matrix = matrix;
	_yStore = new double *[nbSample];
	int64_t i;

	for (i = 0; i < _nbSample; i++) {
		_yStore[i] = ((_matrix[i])->getGaussianSample())->getTabValue();
	}
	_deleteSamples = true;
}

//------------
// Constructor (used in DCV context)
//------------
GaussianData::GaussianData(int64_t nbSample, int64_t pbDimension, 
		Data * originalData, CVBlock & block) 
: Data(nbSample, pbDimension) 
{
	GaussianData * origData = (GaussianData *) (originalData);
	Sample ** origMatrix = origData->_matrix;

	// 1/ (2 * pi)^(d/2)
	_Inv2PiPow = 1.0 / pow(2.0 * XEMPI, pbDimension / 2.0);
	_pbDimensionLog2Pi = pbDimension * log(2.0 * XEMPI);
	_halfPbDimensionLog2Pi = _pbDimensionLog2Pi / 2.0;
	__tmpTabOfSizePbDimension = new double[_pbDimension];
	_deleteSamples = false;


	_weightTotal = block._weightTotal;
	_matrix = new Sample*[_nbSample];
	for (int64_t i = 0; i < _nbSample; i++) {
		_matrix[i] = origMatrix[block._tabWeightedIndividual[i].val];
		//cout<<"ind : "<<block._tabWeightedIndividual[i].val;
		_weight[i] = block._tabWeightedIndividual[i].weight;
		//cout<<" - weight : "<<block._tabWeightedIndividual[i].weight<<endl;
	}

	_yStore = new double *[nbSample];
	for (int64_t j = 0; j < _nbSample; j++) {
		_yStore[j] = ((_matrix[j])->getGaussianSample())->getTabValue();
	}
}

//----------
//Destructor
//----------
GaussianData::~GaussianData() {
	int64_t i;
	if (_matrix) {

		if (_deleteSamples) {
			for (i = 0; i < _nbSample; i++) {
				delete _matrix[i];
				_matrix[i] = NULL;
			}
		}

		delete[] _matrix;
		_matrix = NULL;
	}
	if (_yStore) {
		delete[] _yStore;
		_yStore = NULL;
	}
	if (__tmpTabOfSizePbDimension) {
		delete[] __tmpTabOfSizePbDimension;
		__tmpTabOfSizePbDimension = NULL;
	}
}

//-----------
// Clone data
//-----------
Data * GaussianData::clone() const {
	GaussianData * newData = new GaussianData(*this);
	return (newData);
}

//------------------
// Clone data matrix
//------------------
Sample ** GaussianData::cloneMatrix() {
	int64_t i;

	Sample ** newMatrix = new Sample*[_nbSample];
	for (i = 0; i < _nbSample; i++) {
		newMatrix[i] = new GaussianSample(_matrix[i]->getGaussianSample());
	}

	return newMatrix;
}

//------
// input
//------
void GaussianData::input(std::ifstream & fi) {

	int64_t j;
	int64_t i;
	double ** p_y;
	double * p_y_i;

	p_y = _yStore;
	for (i = 0; i < _nbSample; i++) {
		p_y_i = *p_y;
		for (j = 0; j < _pbDimension; j++) {
			if (fi.eof()) {
				THROW(InputException, endDataFileReach);
			}
			fi >> p_y_i[j];
		}
		_weight[i] = 1.0;
		p_y++;
	}
	_weightTotal = _nbSample;
}

//------
// input
//------
void GaussianData::input(const DataDescription & dataDescription) {
  //double* curSampleValue = new double[_pbDimension];
  std::unique_ptr<double[]> curSampleValue(new double[_pbDimension]);    
	_weightTotal = 0.0;

	_fileNameData = dataDescription.getFileName();
	std::ifstream fi(_fileNameData.c_str(), ios::in);
	if (!fi.is_open())
		THROW(InputException, wrongDataFileName);

	// Guess separator: ',', ' ' or '\t'
	char sep = '*';
	do
		sep = fi.get();
	while (sep != ' ' && sep != '\t' && sep != ',');
	fi.seekg(0);

	//Dubious needed HACK ??! Otherwise atof() truncate numbers to their integer part.
	//                        (At least on some systems).
	setlocale(LC_NUMERIC, "C");

	string line;
	for (int64_t i = 0; i < _nbSample; i++) {
		getline(fi, line);
		istringstream s(line);
		string atom;
		int64_t gaussianVariableIndex = 0;
		for (int64_t j = 0; j < dataDescription.getNbColumn(); j++) {
			if (s.eof())
				THROW(InputException, endDataFileReach);
			do
				// allow any number of separators between items
				getline(s, atom, sep);
			while (atom.empty());

			if (typeid (*(dataDescription.getColumnDescription(j)))
				== typeid (QuantitativeColumnDescription))
			{
				curSampleValue[gaussianVariableIndex] = atof(atom.c_str());
				_yStore[i][gaussianVariableIndex] = curSampleValue[gaussianVariableIndex];
				gaussianVariableIndex++;
			}
			else {
				if (typeid (*(dataDescription.getColumnDescription(j)))
						== typeid (WeightColumnDescription))
				{
					_weight[i] = atof(atom.c_str());
				}
				else {
					//ignore everything else
				}
			}
		}
		_matrix[i]->getGaussianSample()->setDataTabValue(curSampleValue.get());
		_weightTotal += _weight[i];
	}

	//delete [] curSampleValue;
}

// output
//-------
void GaussianData::output(std::ostream & fo) {
	
	if (VERBOSE == 1) {
		cout << "Sample size: " << _nbSample << endl;
		cout << "  Dimension: " << _pbDimension << endl;
	}
	
	editTab(_yStore, _nbSample, _pbDimension, fo, " ", "");
}

bool GaussianData::verify()const {
	bool res = Data::verify();

	// others verif  ?
	return res;
}

}
