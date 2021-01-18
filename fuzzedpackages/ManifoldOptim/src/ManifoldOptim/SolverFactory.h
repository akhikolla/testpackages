#ifndef SOLVER_FACTORY_H
#define SOLVER_FACTORY_H

#include "ManifoldOptimException.h"
#include "Solvers.h"
#include "LRBFGS.h"
#include "LRTRSR1.h"
#include "RBFGS.h"
#include "RBroydenFamily.h"
#include "RCG.h"
#include "RNewton.h"
#include "RSD.h"
#include "RTRNewton.h"
#include "RTRSD.h"
#include "RTRSR1.h"
#include "RWRBFGS.h"
#include "Problem.h"
#include "def.h"

using namespace ROPTLIB;

class SolverFactory
{
public:
	static Solvers* GetSolver(const std::string& solverName, const Problem* prob,
		const Variable* X_init, LinearOPE* H_init);
};

#endif

