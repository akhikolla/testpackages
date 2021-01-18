/***************************************************************************
                             SRC/mixmod/Kernel/Algo/SEMAlgo.cpp  description
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

#include "mixmod/Kernel/Algo/SEMAlgo.h"
#include "mixmod/Kernel/Model/Model.h"

namespace XEM {

//------------
// Constructor
//------------

SEMAlgo::SEMAlgo() {
	_algoStopName = NBITERATION;
}

/// Copy constructor
SEMAlgo::SEMAlgo(const SEMAlgo & semAlgo) : Algo(semAlgo) {
}

//------------
// Constructor
//------------
SEMAlgo::SEMAlgo(AlgoStopName algoStopName, int64_t nbIteration)
: Algo(algoStopName, defaultEpsilon, nbIteration) {
}

//----------
//Destructor
//----------
SEMAlgo::~SEMAlgo() {
}

// clone
//------
Algo * SEMAlgo::clone() {
	return (new SEMAlgo(*this));
}

//----------------
// setNbIteration
//----------------
void SEMAlgo::setNbIteration(int64_t nbIteration) {
	if (nbIteration < minNbIterationForSEM) {
		THROW(InputException, nbIterationTooSmall);
	}
	else if (nbIteration > maxNbIteration) {
		THROW(InputException, nbIterationTooLarge);
	}
	else {
		_nbIteration = nbIteration;
	}
}

//---
//run
//---
void SEMAlgo::run(Model *& model) {
	model->setAlgoName(SEM);

	if (DEBUG > 0) {
		cout << "Debut de l'Algo SEM :" << endl;
		model->editDebugInformation();
	}

	_indexIteration = 1;
	// 1rst SEM : to initialize bestModel and bestLL
	model->Estep(); // E Step
	model->Sstep(); // S Step
	model->Mstep(); // M Step
	
	if (DEBUG > 0) {
		cout << "Apres la 1ere iteration de SEM :" << endl;
		model->editDebugInformation();
	}
	
	//Model * bestModel = new Model(model);
    std::unique_ptr<Model>  bestModel(new Model(model));
	double bestLL = bestModel->getLogLikelihood(true); // true : to update fik

	// others iterations
	_indexIteration++;
	double lastLL;
	while (_indexIteration <= _nbIteration) {
		model->Estep(); // E Step
		model->Sstep(); // S Step
		model->Mstep(); // M Step
		
		if (DEBUG > 0) {
			cout << "Apres la " << _indexIteration << " eme iteration de SEM :" << endl;
			model->editDebugInformation();
		}
		
		// select Best Model
		lastLL = model->getLogLikelihood(true); // true : to compute fik
		if (lastLL > bestLL) {

			if (DEBUG > 0)
				cout << "BestModel prend la valeur de LastModel" << endl;

			bestLL = lastLL;
			//delete bestModel;
			//bestModel = new Model(model);
            bestModel.reset(new Model(model));
		}
		_indexIteration++;
	}

	bestModel->Estep(); // E step to update Tik and fik
	// update model (to output)
	model = bestModel.release();
}

}
