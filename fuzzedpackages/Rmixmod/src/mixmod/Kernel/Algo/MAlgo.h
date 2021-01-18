/***************************************************************************
                             SRC/mixmod/Kernel/Algo/MAlgo.h  description
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
#ifndef XEMMALGO_H
#define XEMMALGO_H

#include "mixmod/Kernel/Algo/Algo.h"

namespace XEM {

/**
  @brief Derived class of XEMAlgo for M Algorithm
  @author F Langrognet
 */

class MAlgo : public Algo {

public:

	/// Default constructor
	MAlgo();

	/// copy constructor
	MAlgo(const MAlgo & mAlgo);

	/// Constructor
	MAlgo(AlgoStopName algoStopName, double epsilon, int64_t nbIteration);

	/// Destructor
	virtual ~MAlgo();

	/// clone
	virtual Algo * clone();

	/// Run method
	virtual void run(Model *& model);

	virtual const AlgoName getAlgoName() const;

	virtual void setEpsilon(double epsilon);

	virtual const double getEpsilon() const;
};

inline const AlgoName MAlgo::getAlgoName() const {
	return M;
}

inline const double MAlgo::getEpsilon() const {
	THROW(OtherException, internalMixmodError);
}

inline void MAlgo::setEpsilon(double eps) {
	THROW(OtherException, internalMixmodError);
}

}

#endif
