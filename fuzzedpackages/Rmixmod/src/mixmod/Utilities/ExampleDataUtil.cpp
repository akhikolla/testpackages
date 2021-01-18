/***************************************************************************
                             SRC/mixmod/Utilities/ExampleDataUtil.cpp  description
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
#include "mixmod/Utilities/ExampleDataUtil.h"
#include "mixmod/Utilities/Util.h"

namespace XEM {

/**
 * @brief Structure to store raw data parsed from text files.
 * @var nbSample samples count (integer)
 * @var bData integer 2D array to store binary data
 * @var gData double 2D array to store continuous data
 * @var bNbColumn binary data columns count
 * @var gNbColumn continuous data columns count
 * @var modality (pointer to) modalities for binary data.
 * @var labels optional pointer to a vector of labels (used in discriminant analysis case).
 */
struct GenericData {
	int64_t nbSample;
	int64_t** bData;
	double** gData;
	int64_t bNbColumn;
	int64_t gNbColumn;
	vector<int64_t>* modality;
	vector<int64_t>* labels;
	GenericData (int64_t nbSample, int64_t** bData, double** gData, int64_t bNbColumn,
	             int64_t gNbColumn, vector<int64_t>* modality, vector<int64_t>* labels) {
		this->nbSample = nbSample;
		this->bData = bData;
		this->gData = gData;
		this->bNbColumn = bNbColumn;
		this->gNbColumn = gNbColumn;
		this->modality = modality;
		this->labels = labels;
	}
	~GenericData() {
		if (bData)
			DeleteData (bData, nbSample);
		if (gData)
			DeleteData (gData, nbSample);
		if (modality)
			delete modality;
		if (labels)
			delete labels;
	}
};

GenericData* readGenericData (string fileName) {
	fstream file;
	file.open (fileName.c_str(), std::ios::in);
	if (!file.is_open()) {
		THROW (InputException, wrongDataFileName);
	}

	//first line of input file describes columns types (and optionally modalities counts)
	std::vector<char> columnType;
	std::string csvLine;
	std::getline (file, csvLine);
	istringstream csvStream (csvLine);
	int64_t bNbColumn = 0;
	int64_t gNbColumn = 0;
	string val;
	vector<int64_t>* modality = new vector<int64_t>();
	vector<int64_t>* labels = 0;
	while (std::getline (csvStream, val, ',')) {
		switch (val[0]) {
		case 'c':
		case 'C':
			// continuous (gaussian)
			columnType.push_back ('C');
			gNbColumn++;
			break;
		case 'b':
		case 'B':
			// binary data (categories)
			columnType.push_back ('B');
			modality->push_back (val.length() > 1 ? atoi (val.substr (1).c_str()) : 0);
			bNbColumn++;
			break;
		case 'l':
		case 'L':
			// 'L' or 'l' indicates label (at most one occurrence)
			columnType.push_back ('L');
			labels = new vector<int64_t>();
			break;
		}
	}

	// count samples
	int64_t nbSample = 0;
	while (std::getline (file, csvLine)) {
		if (csvLine.length() <= 1)
			break; //skip new lines at end of file
		nbSample++;
	}
	file.close();

	int64_t** bData = 0;
	if (bNbColumn > 0) {
		bData = new int64_t*[nbSample];
		for (int64_t i = 0; i < nbSample; ++i) {
			bData[i] = new int64_t[bNbColumn];
		}
	}
	double** gData = 0;
	if (gNbColumn > 0) {
		gData = new double*[nbSample];
		for (int64_t i = 0; i < nbSample; ++i) {
			gData[i] = new double[gNbColumn];
		}
	}

	if (labels)
		labels->resize (nbSample);

	file.open (fileName.c_str(), std::ios::in);
	std::getline (file, csvLine); //read away columns descriptors line
	int64_t row = 0;
	while (std::getline (file, csvLine)) {
		std::istringstream csvStream (csvLine);
		int64_t column = 0, bColumn = 0, gColumn = 0;
		while (std::getline (csvStream, val, ',')) {
			switch (columnType[column]) {
			case 'C':
				// continuous (gaussian)
				gData[row][gColumn] = atof (val.c_str());
				gColumn++;
				break;
			case 'B': {
				// binary (categories represented by integers)
				int64_t binValue = atoi (val.c_str());
				bData[row][bColumn] = binValue;
				if (binValue > (*modality) [bColumn]) {
					(*modality) [bColumn] = binValue;
				}
				bColumn++;
				break;
			}
			case 'L':
				// label (at most one occurrence)
				(*labels) [row] = atoi (val.c_str());
				break;
			}
			column++;
		}
		row++;
	}

	return new GenericData (nbSample, bData, gData, bNbColumn, gNbColumn, modality, labels);
}

DataDescription* getXEMDataDescriptionFromGeneric (GenericData* genericData) {
	// infer data type from genericData fields
	DataType dataType;
	if (genericData->bData && genericData->gData)
		dataType = HeterogeneousData;
	else if (genericData->bData)
		dataType = QualitativeData;
	else if (genericData->gData)
		dataType = QuantitativeData;
	
	// call appropriate data constructor
	DataDescription* dataDescription = 0;
	switch (dataType) {
	case QualitativeData: {
		BinaryData* data = new BinaryData (genericData->nbSample, genericData->bNbColumn,
		        * (genericData->modality), genericData->bData);
		dataDescription = new DataDescription (data);
		delete data;
		break;
	}
	case QuantitativeData: {
		GaussianData* data = new GaussianData (genericData->nbSample, genericData->gNbColumn,
		        genericData->gData);
		dataDescription = new DataDescription (data);
		delete data;
		break;
	}
	case HeterogeneousData: {
		BinaryData* bData = new BinaryData (genericData->nbSample, genericData->bNbColumn,
		        * (genericData->modality), genericData->bData);
		GaussianData* gData = new GaussianData (genericData->nbSample, genericData->gNbColumn,
		        genericData->gData);
		CompositeData* data = new CompositeData (bData, gData);
		dataDescription = new DataDescription (data);
		delete data;
		break;
	}
	}
	return dataDescription;
}

ClusteringInput* getClusteringInput (string fileName, const vector<int64_t>& nbCluster) {
	// parse informations from file in a generic way
	GenericData* genericData = readGenericData (fileName);
	// deduce 'specialized' data description from these informations
	DataDescription* dataDescription = getXEMDataDescriptionFromGeneric (genericData);
	delete genericData;

	// NOTE: any found label is ignored (same file format clustering & learn & predict)
	ClusteringInput* cInput = new ClusteringInput (nbCluster, *dataDescription);
	delete dataDescription;

	return cInput;
}

LearnInput* getLearnInput (string fileName) {
	// parse informations from file in a generic way
	GenericData* genericData = readGenericData (fileName);

	// deduce 'specialized' data description from these informations
	DataDescription* dataDescription = getXEMDataDescriptionFromGeneric (genericData);

	LabelDescription* labelDescription = new LabelDescription (genericData->nbSample, * (genericData->labels));
    delete genericData;

	LearnInput* lInput =  new LearnInput (dataDescription, labelDescription);
	delete dataDescription;
	delete labelDescription;

	return lInput;
}

PredictInput* getPredictInput (string fileName, LearnModelOutput* lOutput) {
	// parse informations from file in a generic way
	GenericData* genericData = readGenericData (fileName);

	// deduce 'specialized' data description from these informations
	DataDescription* dataDescription = getXEMDataDescriptionFromGeneric (genericData);
	delete genericData;

	ParameterDescription* paramPredict = new ParameterDescription (lOutput);
	PredictInput* pInput = new PredictInput (dataDescription, paramPredict);
 	delete dataDescription;
// 	delete paramPredict; //TODO: fix memory leak
                         // HINT: data[parameter]Description removal has a side effect: remove lOutput (?!)

	return pInput;
}

}
