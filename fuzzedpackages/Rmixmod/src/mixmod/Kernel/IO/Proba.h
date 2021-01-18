/***************************************************************************
                             SRC/mixmod/Kernel/IO/Proba.h  description
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
#ifndef XEMProba_H
#define XEMProba_H

/** @brief Base class for Label(s)
	@author F Langrognet
 */

#include <vector>
#include <iostream>
#include <stdint.h>

namespace XEM {

// pre-declaration
class Model;

class Proba {

public:

	/// Default constructor
	Proba();

	/// Constructor
	Proba(int64_t nbSample, int64_t nbCluster);

	/// Constructor
	Proba(Model * model);

	Proba(const Proba & iProba);

	/// Destructor
	virtual ~Proba();

	/// Comparison operator
	bool operator ==(const Proba & proba) const;

	/// editProba
	void edit(std::ostream & stream);

	/// getProba
	double ** getTabProba() const;

	/// getProba
	std::vector<std::vector<double> > getProba() const;

	/// set Proba
	void setProba(double ** proba, int64_t nbSample, int64_t nbCluster);
	void setProba(std::vector<std::vector<double> > proba);

	/// Selector
	int64_t getNbSample() const;
	int64_t getNbCluster() const;

	///input stream
	void input(std::ifstream & flux);

private:

	/// Number of samples
	int64_t _nbSample;

	/// Number of cluster
	int64_t _nbCluster;

	/// dim : _nbSample *_nbCluster
	std::vector<std::vector<double> > _proba;

};

inline std::vector<std::vector<double> > Proba::getProba() const {
	return _proba;
}

inline int64_t Proba::getNbSample() const {
	return _nbSample;
}

inline int64_t Proba::getNbCluster() const {
	return _nbCluster;
}

}

#endif
