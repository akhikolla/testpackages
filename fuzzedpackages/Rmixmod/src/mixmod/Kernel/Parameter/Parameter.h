/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/Parameter.h  description
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
#ifndef XEMParameter_H
#define XEMParameter_H

#include "mixmod/Utilities/Util.h"
#include "mixmod/Kernel/Model/ModelType.h"

namespace XEM {

// pre-declaration
class Model;
class Sample;
class Partition;
class GaussianParameter;
class BinaryParameter;

/**
  @brief Base class for XEMParameter(s)
  @author F Langrognet
 */

class Parameter {

public:

	/// Default constructor
	Parameter();

	/// Constructor
	// called by GaussianParameter or BinaryParameter (called by StrategyType)
	// if USER initialisation
	Parameter(int64_t iNbCluster, int64_t iPbDimension, ModelType * iModelType);

	/// Constructor
	// called by XEMParameter (called by XEMModel)
	Parameter(Model * iModel, ModelType * iModelType);

	/// Constructor ccopy)
	Parameter(const Parameter * iParameter);

	/// Destructor
	virtual ~Parameter();

	virtual GaussianParameter* getGaussianParameter() {
		return (GaussianParameter*)this;
	}

	virtual BinaryParameter* getBinaryParameter() {
		return (BinaryParameter*)this;
	}
	/// Comparison operator
	virtual bool operator ==(const Parameter & param) const;

	/// reset to default values
	virtual void reset() = 0;

	//void initRandomForMeanOrCenter();

	/** @brief Selector
		@return A copy of the model
	 */
	virtual Parameter * clone()const = 0;


	//--------
	// compute
	//--------

	/** Compute normal probability density function
		 for iSample the sample and kCluster th cluster
	 */
	virtual double getPdf(int64_t iSample, int64_t kCluster) const = 0;

	// compute normal probability density function
	// for all i=1,..,n and k=1,..,K
	virtual void getAllPdf(double ** tabFik, double * tabProportion) const = 0;

	/** Compute normal probability density function
		 for x vector and kCluster th cluster
	 */
	virtual double getPdf(Sample * x, int64_t kCluster) const = 0;

	// computeTabProportion
	void computeTabProportion();

	/** @brief Selector
		@return The number of free parameters
	 */
	virtual int64_t getFreeParameter() const = 0;

	/// get loglikelihood with one cluster
	virtual double getLogLikelihoodOne() const = 0;

	/// compute Tik for xi (i=0 -> _nbSample-1) when underflow
	virtual void computeTikUnderflow(int64_t i, double ** tabTik);


	//---------------
	// initialization
	//---------------

	/// init user
	virtual void initUSER(Parameter * iParam) = 0;

	int64_t generateRandomIndex(bool * tabIndividualCanBeUsedForInitRandom, 
			double * weight, double totalWeight);

	/// initialize attributes before an InitRandom  
	virtual void initForInitRANDOM() = 0;

	virtual void updateForInitRANDOMorUSER_PARTITION(
			Sample ** tabSampleForInit, bool * tabClusterToInitialze) = 0;

	/// initialize attributes for init USER_PARTITION
	/// outputs :
	/// -  nbInitializedCluster
	/// - tabNotInitializedCluster (array of size _nbCluster)

	virtual void initForInitUSER_PARTITION(int64_t & nbInitializedCluster, 
			bool * tabNotInitializedCluster, Partition * initPartition) = 0;


	//-----------
	// Algorithms
	//-----------

	/// Maximum a posteriori step method
	//virtual void MAPStep() = 0;

	/// Maximization step method
	virtual void MStep() = 0;

	
	//-------
	// in/out
	//-------

	// edit (for debug)
	virtual void edit() = 0;

	/// edit
	virtual void edit(std::ofstream & oFile, bool text = false) = 0;

	/// input
	virtual void input(std::ifstream & fi) = 0;


	//-------
	// select
	//-------

	/// get TabProportion
	double * getTabProportion() const;

	/// get nbCluster
	int64_t getNbCluster() const;

	/// get pbDimension
	int64_t getPbDimension() const;

	/// getFreeProportion
	bool getFreeProportion() const;

	/// getModel
	Model * getModel() const;

	/// getModelType
	ModelType * getModelType() const;

	/// setModel , made it virtual so that composite paramter class can  override it.
	virtual void setModel(Model * iModel);

	///set modeltype
	void setModelType(ModelType * iModeltype);

	///getFilename
	const std::string & getFilename() const;

	///setFilename
	void setFilename(const std::string & filename);

	/// recopie sans faire construction / destruction
	// utilisé par SMALL_EM, CEM_INIT
	virtual void recopy(Parameter * otherParameter) = 0;

	virtual void updateForCV(Model * originalModel, CVBlock & CVBlock) = 0;

	///get Format
	FormatNumeric::FormatNumericFile getFormat()const;

	///set FormatNumeric
	void setFormat(const FormatNumeric::FormatNumericFile format);

protected:

	/// Number of classes
	int64_t _nbCluster;

	/// problem dimension
	int64_t _pbDimension;

	/// Table of proportion of each cluster
	double * _tabProportion;

	/// free proportion ?
	bool _freeProportion;

	/// model
	Model * _model; // not agregated

	///modelType
	ModelType * _modelType;

private:

	///filename
	std::string _filename;

	FormatNumeric::FormatNumericFile _format;
};

//---------------
// inline methods
//---------------

inline int64_t Parameter::getNbCluster() const {
	return _nbCluster;
}

inline int64_t Parameter::getPbDimension() const {
	return _pbDimension;
}

inline double * Parameter::getTabProportion() const {
	return _tabProportion;
}

inline bool Parameter::getFreeProportion() const {
	return _freeProportion;
}

inline Model * Parameter::getModel() const {
	return _model;
}

inline void Parameter::setModel(Model * iModel) {
	_model = iModel;
}

inline void Parameter::setModelType(ModelType * iModeltype) {
	_modelType = iModeltype;
}

inline ModelType * Parameter::getModelType() const {
	return _modelType;
}

inline const std::string & Parameter::getFilename() const {
	return _filename;
}

inline void Parameter::setFilename(const std::string & filename) {
	_filename = filename;
}

inline FormatNumeric::FormatNumericFile Parameter::getFormat() const {
	return _format;
}

inline void Parameter::setFormat(const FormatNumeric::FormatNumericFile format) {
	_format = format;
}

}

#endif
