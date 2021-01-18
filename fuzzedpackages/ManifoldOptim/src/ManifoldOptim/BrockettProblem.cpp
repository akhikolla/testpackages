#include <RcppArmadillo.h>
#include "ManifoldOptimProblem.h"
#include "ManifoldOptimException.h"

using namespace Rcpp;  // This is needed for RCPP_MODULE

class BrockettProblem : public ManifoldOptimProblem
{
public:
	BrockettProblem(const arma::mat& B, const arma::mat& D)
	: ManifoldOptimProblem(), m_B(B), m_D(D)
	{
	}

	virtual ~BrockettProblem() { }

	double objFun(const arma::vec& x) const
	{
		arma::mat X;
		tx(X, x);
		return arma::trace(X.t() * m_B * X * m_D);
	}

	arma::vec gradFun(const arma::vec& x) const
	{
		arma::mat X;
		tx(X, x);
		return reshape(2 * m_B * X * m_D, x.n_elem, 1);
	}

	// Only need to declare this if we want to expose to R via Module
	arma::vec hessEtaFun(const arma::vec& x, const arma::vec& eta) const
	{
		return ManifoldOptimProblem::hessEtaFun(x, eta);
	}

	void tx(arma::mat& X, const arma::vec& x) const
	{
		X = x;
		X.reshape(m_B.n_rows, m_D.n_rows);
	}

	const arma::mat& GetB() const
	{
		return m_B;
	}

	const arma::mat& GetD() const
	{
		return m_D;
	}

private:
	arma::mat m_B;
	arma::mat m_D;
};

RCPP_MODULE(Brockett_module) {
	class_<BrockettProblem>("BrockettProblem")
	.constructor<arma::mat,arma::mat>()
	.method("objFun", &BrockettProblem::objFun)
	.method("gradFun", &BrockettProblem::gradFun)
	.method("GetB", &BrockettProblem::GetB)
	.method("GetD", &BrockettProblem::GetD)
	.method("hessEtaFun", &BrockettProblem::hessEtaFun)
	;
}

