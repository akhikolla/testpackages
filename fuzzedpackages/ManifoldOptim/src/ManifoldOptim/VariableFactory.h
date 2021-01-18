#ifndef VARIABLE_FACTORY_H
#define VARIABLE_FACTORY_H

#include "ManifoldOptimException.h"
#include "StieVariable.h"
#include "SphereVariable.h"
#include "SPDTVariable.h"
#include "L2SphereVariable.h"
#include "OrthGroupVariable.h"
#include "CSOVariable.h"
#include "EucVariable.h"
#include "PSCVariable.h"
#include "LowRankVariable.h"
#include "GrassVariable.h"
#include "SPDVariable.h"
#include "SPDVariable.h"
#include "ObliqueVariable.h"
#include "def.h"

using namespace ROPTLIB;

class VariableFactory
{
public:
	static Element* GetVariable(const std::string& name, integer n, integer m, integer p = 1);
};

#endif


