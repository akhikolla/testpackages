/***************************************************************************
                             SRC/mixmod/Utilities/ExampleDataUtil.h  description
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
/**
 * @file ExampleDataUtil
 * @brief Alternative functions to build Clustering and Discriminant Analysis inputs from files
 * @author Benjamin Auder
 */

#ifndef XEM_EXAMPLEDATAUTIL_H
#define XEM_EXAMPLEDATAUTIL_H

#include "mixmod/Kernel/IO/Data.h"
#include "mixmod/Kernel/IO/CompositeData.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/DataDescription.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Clustering/ClusteringInput.h"
#include "mixmod/DiscriminantAnalysis/Learn/LearnInput.h"
#include "mixmod/DiscriminantAnalysis/Predict/PredictInput.h"
#include "mixmod/DiscriminantAnalysis/Learn/LearnModelOutput.h"
#include "mixmod/Kernel/Model/ModelType.h"

namespace XEM {

// @param fileName Input file name. It should be CSV file with first row containing type for corresponding column.
// Type can be 'C' for continuous (gaussian) or 'B' for binary (categories).
ClusteringInput* getClusteringInput(string fileName, const vector<int64_t>& nbCluster);
LearnInput* getLearnInput(string fileName);
PredictInput* getPredictInput(string fileName, LearnModelOutput* lOutput);

template<class T>
inline void DeleteData(T ** data, int nbSample){
	for (int i=0; i<nbSample; i++)
		delete [] data[i];
	delete [] data;
}

}

#endif
