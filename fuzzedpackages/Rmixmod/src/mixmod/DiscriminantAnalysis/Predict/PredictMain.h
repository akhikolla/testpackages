/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Predict/PredictMain.h  description
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
#ifndef XEMPREDICTMAIN_H
#define XEMPREDICTMAIN_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class Input;
class PredictInput;
class PredictOutput;
class Model;

/** 
 \class XEMPredictMain
 Main class for Predict treatment (2nd step of discriminant analysis)
 @author F. Langrognet - R Lebret
		@date 2012
		@brief XEMPredictMain class
 */
class PredictMain {

public:

	/// Default constructor
	PredictMain();

	/// Constructor
	PredictMain(PredictInput * input,  PredictOutput * output = NULL);

	/// Destructor
	virtual ~PredictMain();


	// Accessor
	PredictOutput * getPredictOutput() const;

	// get Input
	Input * getInput();

	/// Run method
	void run(IoMode iomode = IoMode::NUMERIC, int verbose = 0, int massiccc = 0);
	void setOutputNull();

private:
	
	// input container
	PredictInput * _input;
	// output container
	PredictOutput * _output;
};

inline PredictOutput * PredictMain::getPredictOutput() const {
	return _output;
}

inline void PredictMain::setOutputNull() {
	_output = NULL;
}

}

#endif
