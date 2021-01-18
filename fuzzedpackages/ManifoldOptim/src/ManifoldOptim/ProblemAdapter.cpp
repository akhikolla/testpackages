#include "ProblemAdapter.h"

ProblemAdapter::ProblemAdapter(ManifoldOptimProblem* up)
: m_up(up), m_usedHessian(false)
{
}

ProblemAdapter::~ProblemAdapter(){}

double ProblemAdapter::f(Variable* x) const
{
	const arma::vec& X = ToArmaVec(x);
	return m_up->objFun(X);
}

void ProblemAdapter::EucGrad(Variable* x, Vector* egf) const
{
	const arma::vec& X = ToArmaVec(x);
	const arma::mat& G = m_up->gradFun(X);
	CopyFrom(egf, G);
}

void ProblemAdapter::EucHessianEta(Variable* x, Vector* etax, Vector* exix) const
{
	m_usedHessian = true;
	const arma::vec& X = ToArmaVec(x);
	const arma::vec& eta = ToArmaVec(etax);
	const arma::vec& hessEta = m_up->hessEtaFun(X, eta);
	CopyFrom(exix, hessEta);
}
