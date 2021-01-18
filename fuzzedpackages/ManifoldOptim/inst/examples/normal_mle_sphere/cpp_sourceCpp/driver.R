library(ManifoldOptim)
library(mvtnorm)

set.seed(1234)

Rcpp::sourceCpp(code = '
//[[Rcpp::depends(RcppArmadillo,ManifoldOptim)]]
#include <RcppArmadillo.h>
#include <ManifoldOptim.h>

using namespace Rcpp;  // This is needed for RCPP_MODULE

class NormalMLEProblem : public ManifoldOptimProblem
{
public:
	NormalMLEProblem(const arma::mat& Y)
	: ManifoldOptimProblem(), m_Y(Y)
	{
	}

	virtual ~NormalMLEProblem() { }

	double objFun(const arma::vec& x) const
	{
		size_t n = m_Y.n_rows;
		size_t p = m_Y.n_cols;

		arma::vec mu(p);
		arma::mat Sigma(p,p);
		tx(mu, Sigma, x);

		double ll = -0.5 * n * p * log(2*arma::datum::pi);
		ll -= 0.5 * n * arma::log_det(Sigma).real();
		for (size_t i = 0; i < n; i++) {
			arma::vec q_i = trans(m_Y.row(i)) - mu;
			ll -= 0.5 * arma::dot(q_i, arma::solve(Sigma, q_i));
		}
		return -ll;
	}

	// Only need to declare this if we want to expose to R via Module
	arma::vec gradFun(const arma::vec& x) const
	{
		return ManifoldOptimProblem::gradFun(x);
	}

	// Only need to declare this if we want to expose to R via Module
	arma::vec hessEtaFun(const arma::vec& x, const arma::vec& eta) const
	{
		return ManifoldOptimProblem::hessEtaFun(x, eta);
	}

	void tx(arma::vec& o_mu, arma::mat& o_Sigma, const arma::vec& x) const
	{
		size_t n = m_Y.n_rows;
		size_t p = m_Y.n_cols;

		o_mu = x(arma::span(0, p-1));
		o_Sigma = x(arma::span(p, p*(p+1)-1));
		o_Sigma.reshape(p, p);
	}

	const arma::mat& GetY() const
	{
		return m_Y;
	}

private:
	arma::mat m_Y;
};

RCPP_MODULE(NormalMLE_module) {
	class_<NormalMLEProblem>("NormalMLEProblem")
	.constructor<arma::mat>()
	.method("objFun", &NormalMLEProblem::objFun)
	.method("gradFun", &NormalMLEProblem::gradFun)
	.method("hessEtaFun", &NormalMLEProblem::hessEtaFun)
	.method("GetY", &NormalMLEProblem::GetY)
	;
}
')

# Try to estimate jointly: Sigma in SPD manifold and mu in sphere manifold
n <- 400
p <- 3
mu.true <- rep(1/sqrt(p), p)
Sigma.true <- diag(2,p) + 0.1
y <- rmvnorm(n, mean = mu.true, sigma = Sigma.true)

tx <- function(x) {
	idx.mu <- 1:p
	idx.S <- 1:p^2 + p
	mu <- x[idx.mu]
	S <- matrix(x[idx.S], p, p)
	# S <- (S + t(S)) / 2
	list(mu = mu, Sigma = S)
}

# Take the problem above (defined in R) and make a Problem object for it.
# The Problem can be passed to the optimizer in C++
prob <- new(NormalMLEProblem, y)

mu0 <- diag(1, p)[,1]
Sigma0 <- diag(1, p)
x0 <- c(mu0, as.numeric(Sigma0))

# Test the obj, grad, and Hess fn
eta <- diag(1, p^2 + p, 1)
prob$objFun(x0)
head(prob$gradFun(x0))
head(prob$hessEtaFun(x0, eta))

mani.defn <- get.product.defn(get.sphere.defn(p), get.spd.defn(p))
mani.params <- get.manifold.params(IsCheckParams = TRUE)
solver.params <- get.solver.params(DEBUG = 0, Tolerance = 1e-4,
	Max_Iteration = 1000, IsCheckParams = TRUE, IsCheckGradHess = FALSE)
deriv.params <- get.deriv.params(EpsNumericalGrad = 1e-8, EpsNumericalHessEta = 1e-4)

res <- manifold.optim(prob, mani.defn, method = "LRBFGS", 
	mani.params = mani.params, solver.params = solver.params,
	deriv.params = deriv.params, x0 = x0)
print(res)
print(tx(res$xopt))
