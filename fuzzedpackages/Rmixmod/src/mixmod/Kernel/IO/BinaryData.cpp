/***************************************************************************
                             SRC/mixmod/Kernel/IO/BinaryData.cpp  description
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

#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/BinarySample.h"
#include "mixmod/Kernel/IO/Partition.h"
#include "mixmod/Kernel/IO/DataDescription.h"
#include "mixmod/Kernel/IO/WeightColumnDescription.h"
#include "mixmod/Kernel/IO/QualitativeColumnDescription.h"
#include "mixmod/Kernel/IO/IndividualColumnDescription.h"
#include "mixmod/Utilities/Util.h"
#include <list>

namespace XEM {

//------------
// Constructor
//------------
BinaryData::BinaryData() {
	_reducedData = NULL;
}

//------------
// Constructor
//------------
BinaryData::BinaryData(const BinaryData & iData) : Data(iData) {
	int64_t i;
	_reducedData = NULL;
	Sample ** matrix = iData._matrix;

	_matrix = new Sample*[_nbSample];
	for (i = 0; i < _nbSample; i++)
		_matrix[i] = new BinarySample(matrix[i]->getBinarySample());

	_tabNbModality = new int64_t[_pbDimension];
	recopyTab<int64_t>(iData._tabNbModality, _tabNbModality, _pbDimension);
	/*for (j=0; j<_pbDimension; j++)
	  _tabNbModality[j] = iData->_tabNbModality[j];*/
}

//------------
// Constructor
//------------
BinaryData::BinaryData(int64_t nbSample, int64_t pbDimension, std::vector<int64_t> nbModality) 
: Data(nbSample, pbDimension) 
{
	_reducedData = NULL;
	int64_t j;
	int64_t i;

	_matrix = new Sample*[_nbSample];
	for (i = 0; i < _nbSample; i++) {
		_matrix[i] = new BinarySample(_pbDimension);
	}

	_tabNbModality = new int64_t[_pbDimension];
	for (j = 0; j < _pbDimension; j++) {
		_tabNbModality[j] = nbModality[j];
	}
}

//------------
// Constructor
//------------
BinaryData::BinaryData(int64_t nbSample, int64_t pbDimension, 
		std::vector<int64_t> nbModality, int64_t ** matrix) 
: Data(nbSample, pbDimension) 
{
	_reducedData = NULL;
	int64_t j;
	int64_t i;

	_matrix = new Sample*[_nbSample];
	for (i = 0; i < _nbSample; i++) {
		_matrix[i] = new BinarySample(_pbDimension, matrix[i]);
	}

	_tabNbModality = new int64_t[_pbDimension];
	for (j = 0; j < _pbDimension; j++) {
		_tabNbModality[j] = nbModality[j];
	}
}

//------------
// Constructor
//------------
BinaryData::BinaryData(int64_t nbSample, int64_t pbDimension,
		const std::string & dataFileName, int64_t * tabNbModality) 
: Data(nbSample, pbDimension) 
{
	_reducedData = NULL;
	int64_t j;
	int64_t i;

	_matrix = new Sample*[_nbSample];
	for (i = 0; i < _nbSample; i++) {
		_matrix[i] = new BinarySample(_pbDimension);
	}

	_tabNbModality = new int64_t[_pbDimension];
	for (j = 0; j < _pbDimension; j++)
		_tabNbModality[j] = tabNbModality[j];

	std::ifstream dataFileStream(dataFileName.c_str(), ios::in);
	if (!dataFileStream.is_open()) {
		dataFileStream.close();
		THROW(InputException, wrongDataFileName);
	}
	input(dataFileStream);
	dataFileStream.close();
	_fileNameData = dataFileName;
}

//------------
// Constructor 
//------------
BinaryData::BinaryData(int64_t nbSample, int64_t pbDimension, int64_t * tabNbModality, 
		double weightTotal, Sample **& matrix, double * weight) 
: Data(nbSample, pbDimension, weightTotal, weight) 
{
	_reducedData = NULL;
	_matrix = matrix;
	_tabNbModality = new int64_t[_pbDimension];
	for (int64_t i = 0; i < _pbDimension; i++) {
		_tabNbModality[i] = tabNbModality[i];
	}
}

//------------
// Constructor (used in DCV context)
//------------
// originalData contains all the data set
BinaryData::BinaryData(int64_t nbSample, int64_t pbDimension, 
		Data * originalData, CVBlock & block) 
: Data(nbSample, pbDimension) 
{
	_reducedData = NULL;
	BinaryData * origData = (BinaryData *) (originalData);
	Sample ** origMatrix = origData->_matrix;

	_tabNbModality = new int64_t[_pbDimension];
	for (int64_t j = 0; j < _pbDimension; j++)
		_tabNbModality[j] = origData->_tabNbModality[j];

	_weightTotal = block._weightTotal;
	_matrix = new Sample*[_nbSample];

	for (int64_t i = 0; i < _nbSample; i++) {
		_matrix[i] = new BinarySample(pbDimension, 
				(origMatrix[block._tabWeightedIndividual[i].val]->getBinarySample())->getTabValue());
		//cout<<"ind : "<<block._tabWeightedIndividual[i].val;

		_weight[i] = block._tabWeightedIndividual[i].weight;
		//cout<<" - weight : "<<block._tabWeightedIndividual[i].weight<<endl;
	}
}

//----------
//Destructor
//----------
BinaryData::~BinaryData() {
	int64_t i;

	if (_matrix) {
		for (i = 0; i < _nbSample; i++) {
			delete _matrix[i];
			//_matrix[i] = NULL;
		}
		delete[] _matrix;
		_matrix = NULL;
	}

	if (_tabNbModality) {
		delete[] _tabNbModality;
		_tabNbModality = NULL;
	}

	if (_reducedData) {
		delete _reducedData;
		_reducedData = NULL;
	}
}

//-----------
// Clone data
//-----------
Data * BinaryData::clone() const {
	BinaryData * newData = new BinaryData(*this);
	return (newData);
}

//------------------
// Clone data matrix
//------------------
Sample ** BinaryData::cloneMatrix() {
	int64_t i;

	Sample ** newMatrix = new Sample*[_nbSample];
	for (i = 0; i < _nbSample; i++) {
		newMatrix[i] = new BinarySample(_matrix[i]->getBinarySample());
	}

	return newMatrix;
}

//------
// input
//------
void BinaryData::input(std::ifstream & fi) {
	int64_t i;
	int64_t j;
	int64_t * curSampleValue = new int64_t[_pbDimension];

	for (i = 0; i < _nbSample; i++) {
		for (j = 0; j < _pbDimension; j++) {
			if (fi.eof()) {
				THROW(InputException, endDataFileReach);
			}
			fi >> curSampleValue[j];

			// Test data value //
			if (curSampleValue[j] > _tabNbModality[j] || curSampleValue[j] <= 0) {
				THROW(InputException, wrongValueInMultinomialCase);
			}
		}// end j
		(_matrix[i]->getBinarySample())->setDataTabValue(curSampleValue);
		_weight[i] = 1;
	}
	_weightTotal = _nbSample;

	delete [] curSampleValue;
}

//------
// input
//------
void BinaryData::input(const DataDescription & dataDescription) {
	int64_t* curSampleValue = new int64_t[_pbDimension];
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

	string line;
	for (int64_t i = 0; i < _nbSample; i++) {
		getline(fi, line);
		istringstream s(line);
		string atom;
		int64_t binaryVariableIndex = 0;
    for (int64_t j = 0; j < dataDescription.getNbColumn(); j++) {
			if (s.eof())
				THROW(InputException, endDataFileReach);
			do
				// allow any number of separators between items
				getline(s, atom, sep);
			while (atom.empty());

			if (typeid (*(dataDescription.getColumnDescription(j))) 
				== typeid (QualitativeColumnDescription)) 
			{
				int64_t binValue = atoi(atom.c_str());
				if (binValue > _tabNbModality[binaryVariableIndex] || binValue <= 0)
					THROW(InputException, wrongValueInMultinomialCase);
				curSampleValue[binaryVariableIndex] = binValue;
				binaryVariableIndex++;
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
		_matrix[i]->getBinarySample()->setDataTabValue(curSampleValue);
		_weightTotal += _weight[i];
	}

	delete [] curSampleValue;
}

//-------
// output
//-------
void BinaryData::output(std::ostream & fo) {
	fo << "Sample size: " << _nbSample;
	fo << "  Dimension: " << _pbDimension;
	fo << " values : " << endl;
	for (int64_t i = 0; i < _nbSample; i++) {
		int64_t * values = getDataTabValue(i);
		for (int64_t j = 0; j < _pbDimension; j++) {
			fo << values[j] << " ";
		}
		fo << " - weight : " << _weight[i];
		fo << endl;
	}
}

//-------------------
// Create reduce data
//
// On passe d'une écriture sur n chiffres (tous d'une base a priori différente : les modalités)
// à une base 10 pour reconnaître les individus identiques et on repasse à la base d'origine
// avec - potentiellement - beaucoup moins d'individus (et des poids).
// En pratique, cela permet d'accélérer réellement les calculs.
//--------------------
Data * BinaryData::reduceData(std::vector<int64_t> & correspondenceOriginDataToReduceData, 
		Partition * knownPartition, Partition * initPartition, 
		Partition *& oKnownPartition, Partition *& oInitPartition) 
{
	int64_t rOld, value;
	int64_t * tabBaseMj = new int64_t[_nbSample];
	int64_t j, k, sizeTabFactor;
	int64_t idxDim;
	int64_t idxSample, sizeList, i;
	double * weight;
	//Sample ** data;

	correspondenceOriginDataToReduceData.resize(_nbSample);

	// test for int64_t limit
	double maxTabFactor = 1;
	for (j = 1; j < _pbDimension; j++) {
		maxTabFactor = maxTabFactor * _tabNbModality[j - 1];
	}
	if (knownPartition && initPartition) {
		maxTabFactor = maxTabFactor * _tabNbModality[_pbDimension - 1];
		maxTabFactor = maxTabFactor * (knownPartition->_nbCluster + 1);
	}
	if ((knownPartition && !initPartition) || (!knownPartition && initPartition)) {
		maxTabFactor = maxTabFactor * _tabNbModality[_pbDimension - 1];
	}
	if (maxTabFactor > int64_t_max) {
		THROW(NumericException, int64_t_max_error);
	}

	list<TWeightedIndividual*> listDiffIndiv;
	list<TWeightedIndividual*>::iterator listIterator;
	list<TWeightedIndividual*>::iterator listBegin;
	list<TWeightedIndividual*>::iterator listEnd;

	// tab of multiplicator value to set sample in base mj //
	//-----------------------------------------------------//
	sizeTabFactor = _pbDimension;
	if (knownPartition)
		sizeTabFactor++;
	if (initPartition)
		sizeTabFactor++;
	int64_t * tabFactor = new int64_t[sizeTabFactor];
	tabFactor[0] = 1;
	for (j = 1; j < _pbDimension; j++) {
		tabFactor[j] = tabFactor[j - 1] * _tabNbModality[j - 1];
		//cout<<"_tabNbModality["<<j-1<<"] : "<<_tabNbModality[j-1]<<endl;
		//cout<<"tabFactor["<<j<<"] : "<<tabFactor[j]<<endl;
	}
	if (knownPartition && initPartition) {
		tabFactor[_pbDimension] = tabFactor[_pbDimension - 1] * _tabNbModality[_pbDimension - 1];
		//cout<<"tabFactor["<<_pbDimension<<"] : "<<tabFactor[_pbDimension]<<endl;
		tabFactor[_pbDimension + 1] = tabFactor[_pbDimension]* (knownPartition->_nbCluster + 1);
		//cout<<"tabFactor["<<_pbDimension+1<<"] : "<<tabFactor[_pbDimension+1]<<endl;
	}
	if ((knownPartition && !initPartition) || (!knownPartition && initPartition)) {
		tabFactor[_pbDimension] = tabFactor[_pbDimension - 1] * _tabNbModality[_pbDimension - 1];
		//cout<<"tabFactor["<<_pbDimension<<"] : "<<tabFactor[_pbDimension]<<endl;
	}

	sizeList = 0;

	// compute all samples in base mj //
	// create list of differents samples in base mj //
	//----------------------------------------------//
	for (i = 0; i < _nbSample; i++) {

		value = 0;
		for (j = 0; j < _pbDimension; j++)
			value += ((_matrix[i]->getBinarySample())->getDataValue(j) - 1) * tabFactor[j];

		if (knownPartition && initPartition) {
			value += (knownPartition->getGroupNumber(i) + 1) * tabFactor[_pbDimension];
			value += (initPartition->getGroupNumber(i) + 1) * tabFactor[_pbDimension + 1];
		}
		if (knownPartition && !initPartition)
			value += (knownPartition->getGroupNumber(i) + 1) * tabFactor[_pbDimension];
		if (!knownPartition && initPartition)
			value += (initPartition->getGroupNumber(i) + 1) * tabFactor[_pbDimension];

		//cout<<"value de l'indvidu "<<i+1<<" : "<<value<<endl;
		tabBaseMj[i] = value;

		// Search if sample already exist in list
		listBegin = listDiffIndiv.begin();
		listEnd = listDiffIndiv.end();

		listIterator = listBegin;
		while ((listIterator != listEnd) && (value > (*listIterator)->val))
			listIterator++;

		// list empty
		if (listBegin == listEnd) {
			TWeightedIndividual * elem = new TWeightedIndividual;
			elem->val = value;
			elem->weight = _weight[i];

			listDiffIndiv.push_front(elem);
			sizeList++;
		}

		else {
			// if element in end of list
			if (listIterator == listEnd) {
				listIterator--;
				if ((*listIterator)->val == value)
					(*listIterator)->weight += _weight[i];
				else {
					TWeightedIndividual * elem = new TWeightedIndividual;
					elem->val = value;
					elem->weight = _weight[i];

					listDiffIndiv.push_back(elem);
					sizeList++;
				}
			}

			// element in begin or in middle of list
			else {
				if ((*listIterator)->val == value)
					(*listIterator)->weight += _weight[i];
				else {
					TWeightedIndividual * elem = new TWeightedIndividual;
					elem->val = value;
					elem->weight = _weight[i];

					if (listIterator == listDiffIndiv.begin())
						listDiffIndiv.push_front(elem);
					else
						listDiffIndiv.insert(listIterator, 1, elem);

					sizeList++;
				}
			}
		}
	} // end for i

	// Create reduce matrix with reduce weight //
	//-----------------------------------------//
	weight = new double[sizeList];
	Sample ** data = new Sample *[sizeList];
	idxSample = 0;
	listBegin = listDiffIndiv.begin();
	listEnd = listDiffIndiv.end();
	for (listIterator = listBegin; listIterator != listEnd; listIterator++) {

		data[idxSample] = new BinarySample(_pbDimension);
		weight[idxSample] = (*listIterator)->weight;

		idxDim = _pbDimension - 1;
		rOld = (*listIterator)->val;
		// knownLabel and initLabel included in computation of base Mj //
		if (knownPartition && initPartition) {
			rOld = rOld % tabFactor[_pbDimension + 1];
			rOld = rOld % tabFactor[_pbDimension];
		}
		// knownLabel or initLabel included in computation of base Mj //
		if ((knownPartition && !initPartition) || (!knownPartition && initPartition))
			rOld = rOld % tabFactor[_pbDimension];

		BinarySample * curDataSample = data[idxSample]->getBinarySample();
		while (idxDim >= 0) {
			curDataSample->setDataValue(idxDim, 1 + rOld / tabFactor[idxDim]);
			rOld = rOld % tabFactor[idxDim];
			idxDim--;
		}
		idxSample++;

	}

	// Set array of correspondence //
	for (i = 0; i < _nbSample; i++) {
		listIterator = listBegin;
		idxSample = 0;
		while ((listIterator != listEnd) && (tabBaseMj[i] != (*listIterator)->val)) {
			listIterator++;
			idxSample++;
		}
		correspondenceOriginDataToReduceData[i] = idxSample;
		//cout<<correspondcenceOriginDataToReduceData[i]<<endl;
	}

	// Set reduce labels //
	int64_t nbCluster;
	if (initPartition) {
		nbCluster = initPartition->_nbCluster;
		oInitPartition = new Partition();
		oInitPartition->setDimension(sizeList, nbCluster);
		oInitPartition->_tabValue = new int64_t*[sizeList];
		for (i = 0; i < sizeList; i++) {
			oInitPartition->_tabValue[i] = new int64_t[nbCluster];
		}
		// l'algo suivant peut etre optimise (on fait eventuellement plusieurs fois la meme chose)
		for (i = 0; i < _nbSample; i++) {
			for (k = 0; k < nbCluster; k++)
				oInitPartition->_tabValue[correspondenceOriginDataToReduceData[i]][k] = 
						initPartition->_tabValue[i][k];
		}
		//oInitLabel->_complete = true;
	}

	if (knownPartition) {
		nbCluster = knownPartition->_nbCluster;
		oKnownPartition = new Partition();
		oKnownPartition->setDimension(sizeList, nbCluster);
		oKnownPartition->_tabValue = new int64_t*[sizeList];
		for (i = 0; i < sizeList; i++)
			oKnownPartition->_tabValue[i] = new int64_t[nbCluster];
		// l'algo suivant peut etre optimise(on fait eventuellement plusieurs fois la meme chose)
		for (i = 0; i < _nbSample; i++) {
			for (k = 0; k < nbCluster; k++)
				oKnownPartition->_tabValue[correspondenceOriginDataToReduceData[i]][k] = 
						knownPartition->_tabValue[i][k];
		}
		//oKnownLabel->_complete = true;
	}

	///////////////////////////////////////////////
	/*cout<<endl<<"nb d'individus differents : "<<sizeList<<endl;
	 cout<<"individu | knownLabel | initLabel | Poids  "<<endl<<endl;
	 for (listIterator = listDiffIndiv.begin(); listIterator != listDiffIndiv.end(); listIterator++)
	 cout<<(*listIterator)->val<<"   "<<(*listIterator)->weight<<endl;

	 for (i=0; i<sizeList; i++){
	 for (j=0; j<_pbDimension; j++){
	 cout<<" "<<((BinarySample*)data[i])->getDataValue(j)<<" ";
   }
	 cout<<"  | ";
	 if (knownPartition){
	 cout<<oKnownPartition->getGroupNumber(i)+1;
   }// 0 if unknown
	 cout<<"  | ";
	 if (initPartition){
	 cout<<oInitPartition->getGroupNumber(i)+1;
   }// 0 if unknown
	 cout<<"  | ";
	 cout<<weight[i]<<" ";
	 cout<<endl;
   }
   
	 for (i=0; i<_nbSample; i++){
	 cout<<"correspondcenceOriginDataToReduceData["<<i<<"] : "<<correspondcenceOriginDataToReduceData[i]<<endl;
   }*/
	//////////////////////////////////////////////

	if (_reducedData) {
		delete _reducedData;
	}

	_reducedData = new BinaryData(
			sizeList, _pbDimension, _tabNbModality, _weightTotal, data, weight);

	delete[] tabFactor;
	delete[] tabBaseMj;
	listIterator = listDiffIndiv.begin();

	while (!listDiffIndiv.empty()) {
		TWeightedIndividual * pelem = *listDiffIndiv.begin();
		listDiffIndiv.pop_front();
		delete pelem;
	}

	delete[] weight;

	return _reducedData;
}

//-------
// Verify
//-------
bool BinaryData::verify()const {
	bool res = Data::verify();

	// others verif  ?
	// les valeurs sont bien comprises entre 1 et _tabModality[i] , ...
	return res;
}

//------------------
// Return tab value
//-----------------
int64_t * BinaryData::getDataTabValue(int64_t idxSample) const {
	return (_matrix[idxSample]->getBinarySample())->getTabValue();
}

}
