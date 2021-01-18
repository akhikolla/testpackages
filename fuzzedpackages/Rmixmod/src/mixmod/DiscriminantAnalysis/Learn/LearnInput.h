/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Learn/LearnInput.h  description
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
#ifndef XEMLEARNINPUT_H
#define XEMLEARNINPUT_H

#include "mixmod/Kernel/IO/Input.h"

namespace XEM {

/** 
 \class XEMLearnInput
 Main class for Learn Input (1rst step of discriminant analysis)
 @author F. Langrognet - R Lebret
		@date 2012
		@brief XEMLearnInput class derived from XEMInput
 */
class LearnInput : public Input {

public:

	/// Default Constructor
	LearnInput();

	/// Copy Constructor
	LearnInput(const LearnInput & CInput);

	/// Initialisation constructor
	LearnInput(DataDescription * learnData, LabelDescription * knownLabelDescription);

	/// Destructor
	virtual ~LearnInput();

	// Accessors
	const int64_t getNbCVBlock() const;

	/// setCriterionName
	virtual void setCriterion(std::vector<CriterionName> const & criterionName);

	/// setCriterion
	virtual void setCriterion(const CriterionName criterionName, unsigned int index);

	///insertCriterion
	virtual void insertCriterion(const CriterionName criterionName, unsigned int index);

	///addCriterion
	virtual void addCriterion(const CriterionName criterionName);

	/// set the number of CV blocks
	void setNbCVBlock(int64_t nbCVBlock);

protected:
	
	/// verif
	virtual bool verif();

private:
	
	// number of CV blocks
	int64_t _nbCVBlock;
};

inline const int64_t LearnInput::getNbCVBlock() const {
	return _nbCVBlock;
}

}

#endif
