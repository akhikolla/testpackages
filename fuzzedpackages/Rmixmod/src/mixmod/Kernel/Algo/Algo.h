/***************************************************************************
                             SRC/mixmod/Kernel/Algo/Algo.h  description
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
#ifndef XEMALGO_H
#define XEMALGO_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class Model;

/**
		@brief Base class for Algorithm(s)
		@author F Langrognet
 */

class Algo {

public:

	/// Default constructor
	Algo();

	/// copy constructor
	Algo(const Algo & algo);

	/// Constructor
	Algo(AlgoStopName algoStopName, double epsilon, int64_t nbIteration);

	/// Destructor
	virtual ~Algo();

	/// clone
	virtual Algo * clone() = 0;

	/// Run method
	virtual void run(Model *& model) = 0;

	//void edit(std::ofstream & oFile);
	void edit(std::ostream & out);

	virtual const AlgoStopName getAlgoStopName() const;

	virtual const AlgoName getAlgoName() const = 0;

	virtual void setAlgoStopName(AlgoStopName algoStopName);

	virtual void setNbIteration(int64_t nbIteration);

	virtual const int64_t getNbIteration() const;

	virtual void setEpsilon(double epsilon);

	virtual const double getEpsilon() const;

	friend std::ostream & operator <<(std::ostream & fo, Algo & algo);

protected:

	/// Type of stopping rule of the algorithm
	AlgoStopName _algoStopName;

	/// Number of iterations
	int64_t _nbIteration;

	/// Current iteration number
	int64_t _indexIteration;

	/// Value of Epsilon (default 1.e-4)
	double _epsilon;

	/** @brief Selector
	@return     1 if algorithm not reached else 0
	 */
	bool continueAgain();

	double _xml_old;

	double _xml;

#if SAVE_ALL_MODELS
	Model ** _tabModel; // aggregate
#endif
};

inline void Algo::setAlgoStopName(AlgoStopName algoStopName) {
	_algoStopName = algoStopName;
}

inline void Algo::setNbIteration(int64_t nbIteration) {
	if (nbIteration < minNbIteration) {
		THROW(InputException, nbIterationTooSmall);
	}
	else if (nbIteration > maxNbIteration) {
		THROW(InputException, nbIterationTooLarge);
	}
	else {
		_nbIteration = nbIteration;
	}
}

inline const int64_t Algo::getNbIteration() const {
	return _nbIteration;
}

inline const double Algo::getEpsilon() const {
	return _epsilon;
}

inline const AlgoStopName Algo::getAlgoStopName() const {
	return _algoStopName;
}

// others functions
Algo * createDefaultClusteringAlgo();
}

#endif
