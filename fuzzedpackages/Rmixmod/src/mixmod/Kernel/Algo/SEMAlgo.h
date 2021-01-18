/***************************************************************************
                             SRC/mixmod/Kernel/Algo/SEMAlgo.h  description
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
#ifndef XEMSEMALGO_H
#define XEMSEMALGO_H

#include "mixmod/Kernel/Algo/Algo.h"

namespace XEM {

/**
  @brief Derived class of XEMAlgo for SEM Algorithm(s)
  @author F Langrognet
 */

class SEMAlgo : public Algo {

public:

	/// Default constructor
	SEMAlgo();

	/// Copy constructor
	SEMAlgo(const SEMAlgo & semAlgo);

	/// Constructor
	SEMAlgo(AlgoStopName algoStopName, int64_t nbIteration);

	/// Destructor
	virtual ~SEMAlgo();

	/// clone
	virtual Algo * clone();

	/// Run method
	virtual void run(Model *& model);

	virtual const AlgoName getAlgoName() const;

	virtual void setNbIteration(int64_t nbIteration);

	virtual void setEpsilon(double epsilon);

	virtual const double getEpsilon() const;
};

inline const AlgoName SEMAlgo::getAlgoName() const {
	return SEM;
}

inline const double SEMAlgo::getEpsilon() const {
	THROW(OtherException, internalMixmodError);
}

inline void SEMAlgo::setEpsilon(double eps) {
	THROW(InputException, wrongAlgoStopName);
}

}

#endif
