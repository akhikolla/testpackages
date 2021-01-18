/***************************************************************************
                             SRC/mixmod/Kernel/Algo/EMAlgo.h  description
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
#ifndef XEMEMALGO_H
#define XEMEMALGO_H

#include "mixmod/Kernel/Algo/Algo.h"

namespace XEM {

/**
  @brief Derived class of XEMAlgo for EM Algorithm(s)
  @author F Langrognet & A Echenim
 */

class EMAlgo : public Algo {

public:

	/// Default constructor
	EMAlgo();

	/// Copy constructor
	EMAlgo(const EMAlgo & emAlgo);

	/// Constructor
	EMAlgo(AlgoStopName algoStopName, double epsilon, int64_t nbIteration);

	/// Destructor
	virtual ~EMAlgo();

	/// clone
	virtual Algo * clone();

	/// Run method
	virtual void run(Model *& model);

	virtual const AlgoName getAlgoName() const;
};

inline const AlgoName EMAlgo::getAlgoName() const {
	return EM;
}

}

#endif
