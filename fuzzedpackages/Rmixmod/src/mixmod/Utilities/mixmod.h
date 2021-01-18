/***************************************************************************
                             SRC/mixmod/Utilities/mixmod.h  description
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

#ifndef MIXMOD_H
#define MIXMOD_H

#ifndef WANT_STREAM
#define WANT_STREAM
#endif

#ifndef WANT_MATH
#define WANT_MATH
#endif
#ifndef XEMmathLib
#define XEMmathLib 1 // default is Eigen
#endif
// mixmod includes
#include "mixmod/Clustering/ClusteringInput.h"
#include "mixmod/Clustering/ClusteringOutput.h"
#include "mixmod/Clustering/ClusteringModelOutput.h"
#include "mixmod/Clustering/ClusteringMain.h"
#include "mixmod/Clustering/ClusteringStrategy.h"

#include "mixmod/DiscriminantAnalysis/Learn/LearnInput.h"
#include "mixmod/DiscriminantAnalysis/Learn/LearnOutput.h"
#include "mixmod/DiscriminantAnalysis/Learn/LearnModelOutput.h"
#include "mixmod/DiscriminantAnalysis/Learn/LearnMain.h"

#include "mixmod/DiscriminantAnalysis/Predict/PredictInput.h"
#include "mixmod/DiscriminantAnalysis/Predict/PredictOutput.h"
#include "mixmod/DiscriminantAnalysis/Predict/PredictModelOutput.h"
#include "mixmod/DiscriminantAnalysis/Predict/PredictMain.h"

#include "mixmod/Utilities/Random.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/CompositeData.h"
#include "mixmod/Kernel/Parameter/Parameter.h"
#include "mixmod/Kernel/IO/Label.h"
#include "mixmod/Kernel/IO/Proba.h"
#include "mixmod/Kernel/IO/Partition.h"
#include "mixmod/Kernel/Criterion/CriterionOutput.h"

#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/IO/ParameterDescription.h"
#include "mixmod/Kernel/IO/ProbaDescription.h"

#include "mixmod/Utilities/Error.h"

#include "mixmod/Kernel/Parameter/BinaryEjParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkjhParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkjParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEParameter.h"
#include "mixmod/Kernel/Parameter/BinaryParameter.h"

#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Parameter/GaussianDiagParameter.h"
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Kernel/Parameter/GaussianSphericalParameter.h"
#include "mixmod/Kernel/Parameter/GaussianHDDAParameter.h"

#endif
