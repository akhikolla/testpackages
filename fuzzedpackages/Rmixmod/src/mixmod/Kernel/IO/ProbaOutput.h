/***************************************************************************
                             SRC/mixmod/Kernel/IO/ProbaOutput.h  description
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
#ifndef XEMProbaOutput_H
#define XEMProbaOutput_H

/** @brief Base class for Label(s)
	@author F Langrognet
 */

#include <iostream>
#include <stdint.h>

namespace XEM {

// pre-declaration
class Model;

class ProbaOutput {

public:

	/// Default constructor
	ProbaOutput();

	/// Constructor
	ProbaOutput(Model * model);

	ProbaOutput(ProbaOutput * iProbaOutput);

	/// Destructor
	virtual ~ProbaOutput();

	/// setCVLabel
	void setCVLabel(int64_t * CVLabel);

	/// editPartition
	void editPartition(std::ofstream & oFile);

	/// editLabel
	void editLabel(std::ofstream & oFile);

	/// editLabel
	void editLabel();

	/// editPostProba
	void editPostProba(std::ofstream & oFile);

	/// editCVLabel
	void editCVLabel(std::ofstream & oFile);

	/// Selector
	int64_t * getTabLabel() const;

	/// Selector
	int64_t * getTabCVLabel() const;

	/// Selector
	int64_t ** getTabPartition() const;

	/// Selector
	double ** getTabPostProba() const;
	/// Selector
	int64_t getNbSample() const;
	int64_t getNbCluster() const;
	///clone
	ProbaOutput * clone();

private:

	/// Number of samples
	int64_t _nbSample;

	/// Number of cluster
	int64_t _nbCluster;

	/// Vector of sample label (dim :_nbSample)
	int64_t * _tabLabel;

	/// Yes if CV criterion
	bool _CVLabelAvailable;

	/// Table of sample label from cross validation (dim :_nbSample)
	int64_t * _tabCVLabel;

	/// Matrix of partition  _partition(i,j)=1 if sample i in class j
	/// dim :  _nbSample * _nbCluster

	//TODO RD : utiliser XEMPartition plutôt
	int64_t ** _tabPartition;

	/// Matrix of posterior probabilities
	/// dim : _nbSample *_nbCluster
	double ** _tabPostProba;
};

//---------------
// inline methods
//---------------

inline int64_t * ProbaOutput::getTabLabel() const {
	return _tabLabel;
}

inline int64_t * ProbaOutput::getTabCVLabel() const {
	return _tabCVLabel;
}

inline int64_t ** ProbaOutput::getTabPartition() const {
	return _tabPartition;
}

inline double ** ProbaOutput::getTabPostProba() const {
	return _tabPostProba;
}

inline int64_t ProbaOutput::getNbSample() const {
	return _nbSample;
}

inline int64_t ProbaOutput::getNbCluster() const {
	return _nbCluster;
}

}

#endif
