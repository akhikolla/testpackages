/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Predict/PredictInput.h  description
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
#ifndef XEMPREDICTINPUT_H
#define XEMPREDICTINPUT_H

#include "mixmod/Kernel/IO/Input.h"

namespace XEM {

// pre-declaration
class PredictStrategy;
class Parameter;
class ParameterDescription;

/** 
 \class XEMPredictInput
 Main class for Predict Input (2nd step of discriminant analysis)
 @author F. Langrognet - R Lebret
		@date 2012
		@brief XEMPredictInput class derived from XEMInput
 */
class PredictInput : public Input {

public:

	/// Default Constructor
	PredictInput();

	/// Copy Constructor
	PredictInput(const PredictInput & CInput);

	/// Initialisation constructor
	PredictInput(DataDescription * predictData, ParameterDescription * classificationRule);

	/// Destructor
	virtual ~PredictInput();

	// accessor
	Parameter * getClassificationRule() const;

	/// setCriterionName
	virtual void setCriterion(std::vector<CriterionName> const & criterionName);

	/// setCriterion
	virtual void setCriterion(const CriterionName criterionName, unsigned int index);

	///insertCriterion
	virtual void insertCriterion(const CriterionName criterionName, unsigned int index);

	///addCriterion
	virtual void addCriterion(const CriterionName criterionName);

	/// getCriterion
	virtual CriterionName getCriterionName(unsigned int index) const;

	// remove criterion
	virtual void removeCriterion(unsigned int index);

	// set model type
	virtual void setModelType(const ModelType * modelType, unsigned int index);

	// insert model type
	virtual void insertModelType(const ModelType * modelType, unsigned int index);

	// remove model type
	virtual void removeModelType(unsigned int index);

	// add model type
	virtual void addModel(ModelName const modelName );

protected:
	
	/// pointer to a classification rule
	Parameter * _classificationRule;
	/// verif
	virtual bool verif();
};

// accessor
inline Parameter * PredictInput::getClassificationRule() const {
	return _classificationRule;
}

}

#endif
