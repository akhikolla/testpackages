/***************************************************************************
                             SRC/mixmod/Kernel/Algo/MAPAlgo.h  description
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
#ifndef XEMMAPALGO_H
#define XEMMAPALGO_H

#include "mixmod/Kernel/Algo/Algo.h"

namespace XEM {

/**
  @brief Derived class of XEMAlgo for MAP Algorithm
  @author F Langrognet
 */

class MAPAlgo : public Algo {

public:

	/// Default constructor
	MAPAlgo();

	// copy constructor
	MAPAlgo(const MAPAlgo & mapAlgo);

	/// Constructor
	MAPAlgo(AlgoStopName algoStopName, double epsilon, int64_t nbIteration);

	/// Destructor
	virtual ~MAPAlgo();


	/// clone
	virtual Algo * clone();

	/// Run method
	virtual void run(Model *& model);

	virtual const AlgoName getAlgoName() const;

	virtual void setEpsilon(double epsilon);

	virtual const double getEpsilon() const;

protected:

	/// Number of clusters
	int64_t _nbCluster;
};

inline const AlgoName MAPAlgo::getAlgoName() const {
	return MAP;
}

inline const double MAPAlgo::getEpsilon() const {
	THROW(OtherException, internalMixmodError);
}

inline void MAPAlgo::setEpsilon(double eps) {
	THROW(OtherException, internalMixmodError);
}

}

#endif
