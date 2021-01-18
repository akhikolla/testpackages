/***************************************************************************
                             SRC/mixmod/Kernel/IO/Input.h  description
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
#ifndef XEMINPUT_H
#define XEMINPUT_H

#include "mixmod/Utilities/Util.h"
#include "mixmod/Kernel/IO/DataDescription.h"

namespace XEM {

// pre-declaration
class LabelDescription;
class Partition;

/** 
 \class XEMInput
 Main class for Mixmod input
 @author F. Langrognet
		@date 2011
		@brief XEMInput classe
 */
class Input {

public:

	/// Default Constructor
	Input();

	/// Copy Constructor
	Input(const Input & input);

	/// Initialisation constructor
	Input(const std::vector<int64_t> & iNbCluster, const DataDescription & iDataDescription);

	/// clone initialisation
	void cloneInitialisation(const std::vector<int64_t> & iNbCluster, const DataDescription & iDataDescription);

	// Constructor used in DCV context
	Input(Input * originalInput, CVBlock & learningBlock);

	/// Destructor
	virtual ~Input();

	
	//----------------------------
	// Modifiers
	//----------------------------

	
	//------ Sample  ----//
	/// getNbSample 
	int64_t getNbSample() const;

	
	//------ Dimension  ----//
	/// getPbDimension
	int64_t getPbDimension() const;

	
	//------ nbCluster  ----//
	/// get all NbCluster
	std::vector<int64_t> getNbCluster() const;

	/// get ith NbCluster
	int64_t getNbCluster(int64_t index) const;

	/// get NbCluster size
	int64_t getNbClusterSize() const;

	
	//------ Criterion  ----//
	/// get All Criterion Name
	std::vector <CriterionName> const & getCriterionName() const;

	///getNbCriterionName
	int64_t getNbCriterion() const;

	///getCriterionName[i]
	CriterionName getCriterionName(unsigned int index) const;

	/// setCriterionName
	virtual void setCriterion(std::vector<CriterionName> const & criterionName) = 0;

	/// setCriterionName
	virtual void setCriterion(const CriterionName criterion, unsigned int index) = 0;

	///insertCriterionName[i]
	virtual void insertCriterion(const CriterionName criterion, unsigned int index) = 0;

	// add new criterion
	virtual void addCriterion(const CriterionName criterion) = 0;

	/// removeCriterionName
	void removeCriterion(unsigned int index);

	
	// ----- ModelType ----//
	/// get ModelType vector
	std::vector<ModelType*> getModelType() const;

	/// setModelType
  /*virtual*/ void setModelType(const ModelType * modelType, unsigned int index);

	/// insertModelType
  /*virtual*/ void insertModelType(const ModelType * modelType, unsigned int index);

	/// add new model type (at the end)
	/*virtual*/ void addModelType(const ModelType * modelType);

	/// removeModelType
	void removeModelType(unsigned int index);
	
	/// add new model (modelName -> modelType)
	/*virtual*/ void addModel(ModelName const modelName);

	/// setModel (modelName -> modelType)
	/*virtual*/ void setModel(std::vector<ModelName> const & modelName);


	/// setSubDimensionEqual
	// void setSubDimensionEqual(int64_t modelTypePosition, int64_t subDimensionValue);

	/// setSubDimensionFreel
	// void setSubDimensionFree(int64_t modelTypePosition, int64_t subDimensionValue, int64_t dimensionPosition);

	///setWeight();
	void setWeight(std::string weightFileName);

	///setWeight();
	void setWeight(double* weight);

	///removeWeight();
	void removeWeight();

	///insertWeight();
	void insertWeight(std::string weightFileName);

	///removeWeight();
	//void removeWeight();
	
	
	// ----- KnownPartition ----//
	/// getKnownPartition
	Partition * getKnownPartition() const;

	/// setKnownPartition
	void setKnownPartition(std::string iFileName);

	/// insertKnownPartition
	void insertKnownPartition(NumericPartitionFile partitionFile);

	/// removeKnownPartition
	void removeKnownPartition();

	
	//------- KnownLabel -----------//
	///getKnownLabelDescription
	const LabelDescription * getKnownLabelDescription() const;
	LabelDescription * getKnownLabelDescriptionNC() const;    

	/// setKnownLabelDescription
	void setKnownLabelDescription(LabelDescription & labeldescription);
	void setKnownLabelDescription(LabelDescription *labeldescription);    

	/// removeLabel
	void removeKnownLabelDescription();

	/// isBinaryData
	const DataType getDataType() const;

	// finalize 
	void finalize();

	/// get Data Description
	const DataDescription & getDataDescription() const;

	/// get Data 
	Data * getData() const;

	bool isFinalized() const;

	// print input
	virtual void edit(std::ostream & out) const;

protected:
	
	//---------
	/// verification of inputs validity  
	virtual bool verif();

	// Criterion
	std::vector <CriterionName> _criterionName;

	// ModelType
	std::vector <ModelType*> _modelType;

	/** a input object must be finalized (verif, ...)
	 */
	bool _finalized;

private:

	/// Data Description
	DataDescription _dataDescription;

	// Known Partition
	Partition * _knownPartition;
	LabelDescription * _knownLabelDescription;

	// Number of samples (no reduced data)
	int64_t _nbSample;

	// Problem dimension
	int64_t _pbDimension;

	// nbCluster
	std::vector<int64_t> _nbCluster;
};

inline int64_t Input::getNbSample() const {
	return _nbSample;
}

inline int64_t Input::getPbDimension() const {
	return _pbDimension;
}

inline std::vector<int64_t> Input::getNbCluster() const {
	return _nbCluster;
}

inline int64_t Input::getNbCluster(int64_t index) const {
	return _nbCluster[index];
}

inline int64_t Input::getNbClusterSize() const {
	return _nbCluster.size();
}

inline int64_t Input::getNbCriterion() const {
	return _criterionName.size();
}

inline std::vector<ModelType*> Input::getModelType() const {
	return _modelType;
}

inline Partition * Input::getKnownPartition() const {
	return _knownPartition;
}

inline const LabelDescription * Input::getKnownLabelDescription() const {
	return _knownLabelDescription;
}
inline LabelDescription * Input::getKnownLabelDescriptionNC() const{
	return _knownLabelDescription;
}

inline const DataType Input::getDataType() const {
	return _dataDescription.getDataType();
}

inline bool Input::isFinalized() const {
	return _finalized;
}

inline const DataDescription & Input::getDataDescription() const {
	return _dataDescription;
}

inline Data * Input::getData() const {
	return _dataDescription.getData();
}

inline std::vector<CriterionName> const & Input::getCriterionName() const {
	return _criterionName;
}

}

#endif
