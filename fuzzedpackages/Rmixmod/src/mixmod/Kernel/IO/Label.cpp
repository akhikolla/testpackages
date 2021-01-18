/***************************************************************************
                             SRC/mixmod/Kernel/IO/Label.cpp  description
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
#include "mixmod/Kernel/IO/Label.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/BinaryModel.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/IO/IndividualColumnDescription.h"
#include <algorithm>

namespace XEM {

//------------
// Constructor
//------------
Label::Label() {
	_nbSample = 0;
}

//------------
// Constructor
//------------
Label::Label(int64_t nbSample) {
	_nbSample = nbSample;
	_label.resize(_nbSample);
}

//------------
// Constructor
//------------
Label::Label(Model * model) {

	if (model == NULL) {
		THROW(OtherException, internalMixmodError);
	}

	// compute tabLabel
	//-----------------
	//int64_t * tabLabel_p = NULL;
    std::unique_ptr<int64_t[]>  tabLabel;
	int64_t nbCluster = model->getNbCluster();
	ModelType * modelType = model->getModelType();
	bool binary = isBinary(modelType->_nameModel);

	if (!binary || (binary && !DATA_REDUCE)) {
		_nbSample = model->getNbSample();
		//int64_t ** tabPartition = new int64_t*[_nbSample];
        std::unique_ptr<int64_t*[], TabDeleter<int64_t>>  tabPartition(new int64_t*[_nbSample], TabDeleter<int64_t>(_nbSample));
		for (int64_t i = 0; i < _nbSample; i++) {
			tabPartition[i] = new int64_t[nbCluster];
		}
		//tabLabel = new int64_t[_nbSample];
        tabLabel.reset(new int64_t[_nbSample]); //provides exception-safe deletion
		model->getLabelAndPartitionByMAPOrKnownPartition(tabLabel.get(), tabPartition.get());
		//for (int64_t i = 0; i < _nbSample; i++) {
		//	delete [] tabPartition[i];
		//}
		//delete [] tabPartition;
	}

	else {
		//binary case
		const vector<int64_t> & correspondenceOriginDataToReduceData =
				dynamic_cast<BinaryModel*> (model)->getCorrespondenceOriginDataToReduceData();
		_nbSample = correspondenceOriginDataToReduceData.size();
		//tabLabel = new int64_t[_nbSample];
        tabLabel.reset(new int64_t[_nbSample]); //provides exception-safe deletion
		//label et partition on reduceData
		int64_t nbSampleOfDataReduce = model->getNbSample();
		//int64_t * tabLabelReduce = new int64_t[nbSampleOfDataReduce];
		std::unique_ptr<int64_t[]> tabLabelReduce(new int64_t[nbSampleOfDataReduce]);        
		//int64_t ** tabPartitionReduce = new int64_t*[nbSampleOfDataReduce];
        std::unique_ptr<int64_t*[], TabDeleter<int64_t>>  tabPartitionReduce(new int64_t*[nbSampleOfDataReduce], TabDeleter<int64_t>(nbSampleOfDataReduce));
		for (int64_t i = 0; i < nbSampleOfDataReduce; i++) {
			tabPartitionReduce[i] = new int64_t[nbCluster];
		}
		model->getLabelAndPartitionByMAPOrKnownPartition(tabLabelReduce.get(), tabPartitionReduce.get());

		//  double ** tabPostProbaReduce = NULL;
		//  tabPostProbaReduce = copyTab(estimation->getModel()->getPostProba(),
		//  nbSampleOfDataReduce, nbCluster); // copy

		// convert labelReduce, partitionReduce, postProbaReduce to label, partition, postProba
		for (int64_t i = 0; i < _nbSample; i++) {
			tabLabel[i] = tabLabelReduce[correspondenceOriginDataToReduceData[i]];
		}

		//delete //deletion is made by unique_ptr
		//for (int64_t i = 0; i < nbSampleOfDataReduce; i++) {
		//	delete [] tabPartitionReduce[i];
		//}
		//delete [] tabPartitionReduce;
		//delete[] tabLabelReduce;
	}

	// compute _label
	recopyTabToVector(tabLabel.get(), _label, _nbSample);
	//delete [] tabLabel; //deletion done by unique_ptr
}

//------------
// Constructor
//------------
Label::Label(const Label & iLabel) {
	_nbSample = iLabel.getNbSample();
	_label = iLabel.getLabel();
}

//-----------
// Destructor
//-----------
Label::~Label() {
}

//--------------------
/// Comparison operator
//--------------------
bool Label::operator ==(const Label & label) const {
	if (_nbSample != label.getNbSample()) return false;
	for (int64_t i = 0; i < _nbSample; i++) {
		if (_label[i] != label.getLabel()[i]) return false;
	}
	return true;
}

//----------
// editProba
//----------
void Label::edit(std::ostream & stream) const {
	stream.setf(ios::fixed, ios::floatfield);
	for (int64_t i = 0; i < _nbSample; i++) {
		stream << _label[i] << endl;
	}
}

//---------
// getProba
//---------
int64_t * Label::getTabLabel() const {
	int64_t * res;
	recopyVectorToTab(_label, res);
	return res;
}

//---------
// get Error Rate
//---------
const double Label::getErrorRate(std::vector<int64_t> const & label) const {
	if (_nbSample != (int64_t) label.size()) {
		THROW(InputException, badNumberOfValuesInLabelInput);
	}

	double missClass = 0.0;
	for (int64_t i = 0; i < _nbSample; i++) {
		if (_label[i] != label[i]) ++missClass;
	}
	return missClass / _nbSample;
}

//---------
// get getClassificationTab
//---------
int64_t** Label::getClassificationTab(std::vector<int64_t> const & label, int64_t nbCluster) const {
	if (_nbSample != (int64_t) label.size()) {
		THROW(InputException, badNumberOfValuesInLabelInput);
	}

	// memory allocation
	int64_t** classTab = new int64_t*[nbCluster];
	for (unsigned int i = 0; i < nbCluster; i++) {
		classTab[i] = new int64_t[nbCluster];
	}
	// initialization
	for (unsigned int i = 0; i < nbCluster; i++)
		for (unsigned int j = 0; j < nbCluster; j++)
			classTab[i][j] = 0;

	// loop over labels
	for (int64_t i = 0; i < _nbSample; i++) {
	  //cout<<_label[i]<<endl;
      if(label[i] > 0){
		++classTab[_label[i] - 1][label[i] - 1];
      }
	}

	return classTab;
}

// -----------
//input stream
// read labels between 1 and nbCluster
// -----------
void Label::input(std::ifstream & flux, int64_t nbCluster) {
	int64_t i = 0;
	int64_t read;

	while (i < _nbSample && !flux.eof()) {
		flux >> read;
		if (read >= 1 && read <= nbCluster) {
			_label[i] = read;
		}
		else {
			THROW(InputException, badValueInLabelInput);
		}
		i++;
	}

	if (!flux.eof() && i != _nbSample) {
		THROW(InputException, badNumberOfValuesInLabelInput);
	}
}

// -----------
//input stream
// read labels between 1 and nbCluster
// -----------
void Label::input(const LabelDescription & labelDescription) {
	int64_t i = 0;
	int64_t readLabel;

	std::string labelFilename = labelDescription.getFileName();

	_nbSample = labelDescription.getNbSample();

	std::ifstream fi((labelFilename).c_str(), ios::in);
	if (!fi.is_open()) {
		THROW(InputException, wrongDataFileName);
	}

	while (i < _nbSample && !fi.eof()) {
		for (int64_t j = 0; j < labelDescription.getNbColumn(); ++j) {
			if (fi.eof()) {
				THROW(InputException, endDataFileReach);
			}
			if (typeid (*(labelDescription.getColumnDescription(j)))
					== typeid (IndividualColumnDescription))
			{
				std::string stringTmp;
				fi >> stringTmp;
				//cout<<stringTmp<<endl;
			}
			else {
				fi >> readLabel;
				_label.push_back(readLabel);
			}
		}
		++i;
	}
	if (!fi.eof() && i != _nbSample) {
		THROW(InputException, badNumberOfValuesInLabelInput);
	}
}

}
