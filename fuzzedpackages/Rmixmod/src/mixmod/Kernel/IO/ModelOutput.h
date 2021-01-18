/***************************************************************************
                             SRC/mixmod/Kernel/IO/ModelOutput.h  description
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
#ifndef XEMMODELOUTPUT_H
#define XEMMODELOUTPUT_H

#include "mixmod/Utilities/Util.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/Criterion/CriterionOutput.h"
#include "mixmod/Kernel/Parameter/Parameter.h"
#include "mixmod/Kernel/IO/ParameterDescription.h"

namespace XEM {

// pre-declaration
class ProbaDescription;
class LabelDescription;

/* Note :
A XEMModelOutput object could be created :
- with an XEMEstimation (after calculation)
- or without XEMEstimation. In this case, XEMModelOutput is created from XML mixmod file which contains input and output information. 
 */

/** 
 \class XEMModelOutput
 @author F. Langrognet
		@date 2011
		@brief XEMModelOutput class
 */
class ModelOutput {

public:

	/// Default Constructor
	ModelOutput();

	/// Copy Constructor
	ModelOutput(const ModelOutput & modelOutput);

	/// Initialization Constructor 1
	ModelOutput(Model * estimation);

	/// Initialization Constructor 2
	ModelOutput(ModelType & modelType, int64_t nbCluster, std::vector<CriterionOutput*> & criterionOutput, double likelihood, ParameterDescription & parameterDescription, LabelDescription & labelDescription, ProbaDescription & probaDescription);

	///Initialization Constructor 3
	ModelOutput(ModelType & modelType, int64_t nbCluster, Exception& error);

	/// Destructor
	virtual ~ModelOutput();

	/// Comparison operator
	bool operator ==(const ModelOutput & modelOutput) const;

	// --- get --- ///
	//-------------//
	ModelType getModelType() const;

	int64_t getNbCluster() const;

	ParameterDescription * getParameterDescription() const;

	LabelDescription * getLabelDescription() const;

	ProbaDescription * getProbaDescription() const;

	Exception & getStrategyRunError() const;

	Model * getModel() const;

	double getLikelihood() const;

	CriterionOutput const & getCriterionOutput(CriterionName criterionName) const;
	CriterionOutput const & getCriterionOutput(const int index) const;
	CriterionOutput & getCriterionOutput(CriterionName criterionName);

	// set criterion output
	void setCriterionOutput(CriterionOutput const & criterionOutput);

protected:

	// criterion output
	CriterionOutput _criterionOutput[maxNbCriterion];

	// type of the model
	ModelType _modelType;

	// the number of cluster
	int64_t _nbCluster;

	// parameter description for that model
	ParameterDescription * _parameterDescription;

	// labels for the model
	LabelDescription * _labelDescription;

	// the probabilities of the model
	ProbaDescription * _probaDescription;

	// the model likelihood
	double _likelihood;

	// the error
	Exception * _strategyRunError;
};

inline ModelType ModelOutput::getModelType() const {
	return _modelType;
}

inline int64_t ModelOutput::getNbCluster() const {
	return _nbCluster;
}

inline ParameterDescription * ModelOutput::getParameterDescription() const {
	return _parameterDescription;
}

inline LabelDescription * ModelOutput::getLabelDescription() const {
	return _labelDescription;
}

inline ProbaDescription * ModelOutput::getProbaDescription() const {
	return _probaDescription;
}

inline Exception & ModelOutput::getStrategyRunError() const {
	return *_strategyRunError;
}

inline double ModelOutput::getLikelihood() const {
	return _likelihood;
}

inline CriterionOutput const & ModelOutput::getCriterionOutput(CriterionName criterionName) const {
	return _criterionOutput[criterionName];
}

inline CriterionOutput const & ModelOutput::getCriterionOutput(const int index) const {
	return _criterionOutput[index];
}

inline CriterionOutput & ModelOutput::getCriterionOutput(CriterionName criterionName) {
	return _criterionOutput[criterionName];
}

// Define structure to sort by criterion name
struct SortByCriterion {

	// Constructor
	SortByCriterion(CriterionName criterionName) : _criterionName(criterionName) {
	}

	// Destructor
	~SortByCriterion() {
	}

	// operator()
	inline bool operator ()(const ModelOutput * m1, const ModelOutput * m2) const {
		Exception& error1 = (m1->getCriterionOutput(_criterionName).getError());
		Exception& error2 = (m2->getCriterionOutput(_criterionName).getError());
		if (error1 != NOERROR && error2 != NOERROR) return false;
		if (error1 != NOERROR) return false;
		if (error2 != NOERROR) return true;
		const double value1 = m1->getCriterionOutput(_criterionName).getValue();
		const double value2 = m2->getCriterionOutput(_criterionName).getValue();
		if (value1 == value2) {
			return m1->getParameterDescription()->getParameter()->getFreeParameter() < m2->getParameterDescription()->getParameter()->getFreeParameter();
		}
		else {
			return value1 < value2;
		}
	}

private:
	
	// criterion name
	CriterionName _criterionName;
};

}

#endif
