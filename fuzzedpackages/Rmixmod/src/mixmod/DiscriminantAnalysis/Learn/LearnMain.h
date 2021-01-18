/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Learn/LearnMain.h  description
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
#ifndef XEMLEARNMAIN_H
#define XEMLEARNMAIN_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class Input;
class LearnInput;
class LearnOutput;
class Model;

/** 
 \class XEMLearnMain
 Main class for Learn Treatment (1rst step of discriminant analysis)
 @author F. Langrognet - R Lebret
		@date 2012
		@brief XEMLearnMain class
 */
class LearnMain {

public:

	/// Default constructor
	LearnMain();

	/// Constructor
	LearnMain(LearnInput * input,  LearnOutput * output = NULL);

	/// Destructor
	virtual ~LearnMain();


	// Accessor
	LearnOutput * getLearnOutput() const;

	// get Input
	Input * getInput();

	/// Run method
	void run(int seed = -1, IoMode iomode = IoMode::NUMERIC, int verbose = 0, int massiccc = 0);
	void setOutputNull();

private:
	
	// input container
	LearnInput * _input;
	// output container
	LearnOutput * _output;
};

inline LearnOutput * LearnMain::getLearnOutput() const {
	return _output;
}

inline void LearnMain::setOutputNull() {
	_output = NULL;
}

}

#endif
