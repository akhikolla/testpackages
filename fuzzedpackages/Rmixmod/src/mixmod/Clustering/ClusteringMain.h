/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringMain.h  description
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

#ifndef XEMCLUSTERINGMMAIN_H
#define XEMCLUSTERINGMMAIN_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

class Model;
class ClusteringInput;
class ClusteringOutput;

/** 
 * @class ClusteringMain
 * @brief Main class for Clustering using mixture model(s)
 * @author F. Langrognet
 */
class ClusteringMain {

public:

	/// Invalid default constructor: clustering input and output must be provided.
	ClusteringMain();

	/// Build a ClusteringMain object from ClusteringInput and ClusteringOutput
	ClusteringMain(ClusteringInput * cInput,  ClusteringOutput * output = NULL);

	/// Destructor
	virtual ~ClusteringMain();

	/// Run clustering(s) task(s) described by _input. Fill _output.
	void run(int seed = -1, IoMode iomode = IoMode::NUMERIC, int verbose = 0, int massiccc = 0);

	/// Return pointer to input
	ClusteringInput * getInput() const;

	/// Return pointer to output
	ClusteringOutput * getOutput() const;

	// Set output pointer to null
	void setOutputNull();

private:

	/// Input
	ClusteringInput * _input;
	
	/// Output
	ClusteringOutput * _output;
};

inline  ClusteringInput * ClusteringMain::getInput() const {
	if (_input)
		return _input;
	else
		THROW(OtherException, nullPointerError);
}

inline ClusteringOutput * ClusteringMain::getOutput() const {
	return _output;
}

inline void ClusteringMain::setOutputNull() {
	_output = NULL;
}

}

#endif
