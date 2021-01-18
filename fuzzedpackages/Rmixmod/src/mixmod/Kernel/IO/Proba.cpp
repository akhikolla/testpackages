/***************************************************************************
                             SRC/mixmod/Kernel/IO/Proba.cpp  description
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

#include "mixmod/Kernel/IO/Proba.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/BinaryModel.h"
#include "mixmod/Kernel/Model/ModelType.h"

namespace XEM {

//------------
// Constructor
//------------
Proba::Proba() {
	_nbCluster = 0;
	_nbSample = 0;
}

//------------
// Constructor
//------------
Proba::Proba(int64_t nbSample, int64_t nbCluster) {
	_nbCluster = nbCluster;
	_nbSample = nbSample;
	_proba.resize(_nbSample);
	for (int64_t i = 0; i < _nbSample; i++) {
		_proba[i].resize(_nbCluster);
	}
}

//------------
// Constructor
//------------
Proba::Proba(Model * model) {
	_nbCluster = model->getNbCluster();

	if (model == NULL) {
		THROW(OtherException, internalMixmodError);
	}

	// compute tabProba
	//-----------------
	double ** tabProba = NULL;
	ModelType * modelType = model->getModelType();
	bool binary = isBinary(modelType->_nameModel);
	
	if (!binary || (binary && !DATA_REDUCE)) {
		_nbSample = model->getNbSample();
		tabProba = copyTab(model->getPostProba(), _nbSample, _nbCluster); // copy
	}
	
	else {
		//binary + DATA_REDUCE case [currently deactivated ; see Util.h]
		const std::vector<int64_t> & correspondenceOriginDataToReduceData = 
				dynamic_cast<BinaryModel*> (model)->getCorrespondenceOriginDataToReduceData();
		_nbSample = correspondenceOriginDataToReduceData.size();
		tabProba = new double*[_nbSample];
		for (int64_t i = 0; i < _nbSample; i++) {
			tabProba[i] = new double[_nbCluster];
		}
		int64_t nbSampleOfDataReduce = model->getNbSample();
		double ** tabPostProbaReduce = NULL;
		// copy
		tabPostProbaReduce = copyTab(model->getPostProba(), nbSampleOfDataReduce, _nbCluster);
		//editTab<double>(tabPostProbaReduce,nbSampleOfDataReduce, _nbCluster);
		// convert labelReduce, partitionReduce, postProbaReduce to label, partition, postProba
		for (int64_t i = 0; i < _nbSample; i++) {
			for (int64_t k = 0; k < _nbCluster; k++) {
				int64_t index = correspondenceOriginDataToReduceData[i];
				tabProba[i][k] = tabPostProbaReduce[index][k];
			}
		}

		//delete
		for (int64_t i = 0; i < nbSampleOfDataReduce; i++) {
			delete [] tabPostProbaReduce[i];
		}
		delete [] tabPostProbaReduce;
	}

	// compute _proba
	recopyTabToVector(tabProba, _proba, _nbSample, _nbCluster);
	for (int64_t i = 0; i < _nbSample; i++) {
		delete[] tabProba[i];
	}
	delete [] tabProba;
}

//------------
// Constructor
//------------
Proba::Proba(const Proba & iProba) {
	_nbSample = iProba.getNbSample();
	_nbCluster = iProba.getNbCluster();
	_proba = iProba.getProba();
}

//-----------
// Destructor
//-----------
Proba::~Proba() {
}

//--------------------
/// Comparison operator
//--------------------
bool Proba::operator ==(const Proba & proba) const {
	if (_nbSample != proba.getNbSample()) return false;
	if (_nbCluster != proba.getNbCluster()) return false;
	for (int64_t i = 0; i < _nbSample; i++) {
		for (int64_t k = 0; k < _nbCluster; k++) {
			if (_proba[i][k] != proba.getProba()[i][k]) return false;
		}
	}
	return true;
}

//----------
// editProba
//----------
void Proba::edit(std::ostream & stream) {
	stream.setf(ios::fixed, ios::floatfield);
	for (int64_t i = 0; i < _nbSample; i++) {
		for (int64_t k = 0; k < _nbCluster; k++)
			putDoubleInStream(stream, _proba[i][k], "\t");
		stream << endl;
	}
}

//---------
// getProba
//---------
double ** Proba::getTabProba() const {
	double ** res;
	recopyVectorToTab(_proba, res);
	return res;
}

// -----------
//input stream
// -----------
void Proba::input(std::ifstream & flux) {
	int64_t i = 0;
	int64_t k = 0;

	while (i < _nbSample && !flux.eof()) {
		k = 0;
		while (k < _nbCluster && !flux.eof()) {
			_proba[i][k] = getDoubleFromStream(flux);
			k++;
		}
		if (!flux.eof() && k != _nbCluster) {
			THROW(InputException, notEnoughValuesInProbaInput);
		}
		i++;
	}
	if (!flux.eof() && i != _nbSample) {
		THROW(InputException, notEnoughValuesInProbaInput);
	}
}

}
