#ifndef MANIFOLD_FACTORY_H
#define MANIFOLD_FACTORY_H

#include "ManifoldOptimException.h"
#include "Rcpp.h"
#include "Manifold.h"
#include "Stiefel.h"
#include "Sphere.h"
#include "ProductManifold.h"
#include "L2Sphere.h"
#include "OrthGroup.h"
#include "CpxNStQOrth.h"
#include "Euclidean.h"
#include "PreShapeCurves.h"
#include "SPDVector.h"
#include "LowRank.h"
#include "Grassmann.h"
#include "SPDManifold.h"
#include "LinearOPE.h"
#include "Oblique.h"
#include "ProductElement.h"
#include "def.h"

using namespace ROPTLIB;

class ManifoldFactory
{
public:
	static Manifold* GetManifold(const std::string& name, integer n, integer m, integer p = 1);
};

#endif

