#ifndef __SteepDescent_H_
#define __SteepDescent_H_

#include <valarray>
#include <math.h>
#include "omxState.h"
#include "omxFitFunction.h"
#include "omxExportBackendState.h"
#include "Compute.h"
#include "matrix.h"

void omxSD(GradientOptimizerContext &);

#endif
