/***************************************************************************
                             SRC/mixmod/Kernel/Criterion/ICLCriterion.cpp  description
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
#include "mixmod/Kernel/Criterion/ICLCriterion.h"
#include "mixmod/Kernel/Criterion/CriterionOutput.h"
#include "mixmod/Kernel/Model/Model.h"

namespace XEM {

//------------
// Constructor
//------------
ICLCriterion::ICLCriterion(Model * model) : Criterion(model) {
}

//-----------
// Destructor
//-----------
ICLCriterion::~ICLCriterion() {
}

//---
//run
//---
void ICLCriterion::run(CriterionOutput & output) {
	
	/* Compute ICL (Integrated Completed Likelihood) Criterion */

	/* new version (mixmodLib 3.0 - 2013-10-03)
  ICL = BIC + 2.Sum(Sum(tik.ln(tik))
 */
	

	// initialize value
	double value = 0.0;
	// initialize error
	Exception * error = &NOERROR;

	try {
		double bic_value = 0.0;
		const double loglikelihood = _model->getLogLikelihood(false);
		const int64_t freeParameter = _model->getFreeParameter();
		const double logN = _model->getLogN();
		bic_value = (-2 * loglikelihood)+(freeParameter * logN);
		double entropy = _model->getEntropy();
		value = bic_value + 2.0 * entropy; 
	}
	catch (Exception&e) {
		// add name to criterion output
		output.setCriterionName(ICL);
		// add error to criterion output
		output.setError(e);
		throw;
	}
	// add name to criterion output
	output.setCriterionName(ICL);
	// add value to criterion output
	output.setValue(value);
	// add error to criterion output
	output.setError(*error);
}

}//end of namespace XEM
