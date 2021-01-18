/***************************************************************************
                             SRC/mixmod/Kernel/Criterion/CVCriterion.h  description
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
#ifndef XEMCVCRITERION_H
#define XEMCVCRITERION_H

#include "mixmod/Kernel/Criterion/Criterion.h"

namespace XEM {

// pre-declaration
class LabelDescription;

/**
  @brief Derived class of XEMCriterion for CV Criterion
  @author F Langrognet 
 */

class CVCriterion : public Criterion {

public:

	/// Default constructor
	CVCriterion(Model * model, const int64_t nbCVBlock);

	/// Destructor
	virtual ~CVCriterion();

	/// Run method
	virtual void run(CriterionOutput & output);

	/// Accessor
	std::vector<int64_t> & getCVLabel();

private:

	void createCVBlocks();

	/// Table of XEMCVBlock : testBlock
	CVBlock* _tabCVBlock;

	// pointer to a XEMLabelDescription
	std::vector<int64_t> _cvLabel;

	// number of CV Blocks  
	int64_t _nbCVBlock;

	// initialisation method of cv blocks
	CVinitBlocks _CVinitBlocks;
};

inline std::vector<int64_t> & CVCriterion::getCVLabel() {
	return _cvLabel;
}

}

#endif
