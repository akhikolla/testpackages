#include "SolverFactory.h"

Solvers* SolverFactory::GetSolver(const std::string& solverName, const Problem* prob,
	const Variable* X_init, LinearOPE* H_init)
{
	Solvers* solver;

	if (solverName == "LRBFGS") {
		solver = static_cast<Solvers*>(new LRBFGS(prob, X_init));
	} else if (solverName == "LRTRSR1") {
		solver = static_cast<Solvers*>(new LRTRSR1(prob, X_init));
	} else if (solverName == "MRankAdaptive") {
		throw ManifoldOptimException("MRankAdaptive solver currently not supported");
	} else if (solverName == "RBFGS") {
		solver = static_cast<Solvers*>(new RBFGS(prob, X_init, H_init));
	} else if (solverName == "RBroydenFamily") {
		solver = static_cast<Solvers*>(new RBroydenFamily(prob, X_init, H_init));
	} else if (solverName == "RCG") {
		solver = static_cast<Solvers*>(new RCG(prob, X_init));
	} else if (solverName == "RNewton") {
		solver = static_cast<Solvers*>(new RNewton(prob, X_init));
	} else if (solverName == "RSD") {
		solver = static_cast<Solvers*>(new RSD(prob, X_init));
	} else if (solverName == "RTRNewton") {
		solver = static_cast<Solvers*>(new RTRNewton(prob, X_init));
	} else if (solverName == "RTRSD") {
		solver = static_cast<Solvers*>(new RTRSD(prob, X_init));
	} else if (solverName == "RTRSR1") {
		solver = static_cast<Solvers*>(new RTRSR1(prob, X_init, H_init));
	} else if (solverName == "RWRBFGS") {
		solver = static_cast<Solvers*>(new RWRBFGS(prob, X_init, H_init));
	} else {
		throw ManifoldOptimException("Invalid solver specified");
	}

	return solver;
}

