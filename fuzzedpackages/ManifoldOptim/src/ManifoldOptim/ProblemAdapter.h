#ifndef PROBLEM_ADAPTER_H
#define PROBLEM_ADAPTER_H

#include <RcppArmadillo.h>
#include <iostream>
#include "Problem.h"
#include "Util.h"
#include "ManifoldOptimProblem.h"

class ProblemAdapter : public Problem
{
public:
	ProblemAdapter(ManifoldOptimProblem* up);
	virtual ~ProblemAdapter();

	double f(Variable* x) const;
	void EucGrad(Variable* x, Vector* egf) const;
	void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;
	bool GetUsedHessian() const { return m_usedHessian; };

private:
	mutable ManifoldOptimProblem* m_up;
	mutable bool m_usedHessian;
};

#endif
