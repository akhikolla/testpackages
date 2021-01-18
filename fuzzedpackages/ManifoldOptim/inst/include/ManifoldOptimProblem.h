#ifndef MANIFOLD_OPTIM_PROBLEM_H
#define MANIFOLD_OPTIM_PROBLEM_H

#include <RcppArmadillo.h>
#include <cstring>
#include <assert.h>
#include <iostream>
#include "ManifoldOptimException.h"

class ManifoldOptimProblem
{

public:
	ManifoldOptimProblem()
	: m_epsNumericalGrad(1e-6), m_epsNumericalHessEta(1e-4), m_usedNumericalHessian(false)
	{
	}
	virtual ~ManifoldOptimProblem() {};

	void SetEpsNumericalGrad(double eps)
	{
		if (eps <= 0) {
			throw ManifoldOptimException("In SetEpsNumericalGrad, eps must be positive");
		}
		m_epsNumericalGrad = eps;
	}

	void SetEpsNumericalHessEta(double eps)
	{
		if (eps <= 0) {
			throw ManifoldOptimException("In SetEpsNumericalHessEta, eps must be positive");
		}
		m_epsNumericalHessEta = eps;
	}

	bool GetUsedNumericalHessian() const { return m_usedNumericalHessian; };

	virtual double objFun(const arma::vec& x) const = 0;

	/*
	 * Calculate a numerical approximation to the Jacobian (Euclidean Gradient).
	 * Reference: http://www.maths.lth.se/na/courses/FMN081/FMN081-06/lecture7.pdf
	 * Note that there are better approximations if we are willing to evaluate the
	 * objective more times. See for example:
	 * http://www.boost.org/doc/libs/1_65_1/libs/multiprecision/doc/html/boost_multiprecision/tut/floats/fp_eg/nd.html
	**/
	virtual arma::vec gradFun(const arma::vec& x) const
	{
		double eps = m_epsNumericalGrad;
		double fx = objFun(x);
		arma::vec x_eps = x;
		arma::vec grad(x.n_elem);
	
		for (size_t i = 0; i < x.n_elem; ++i) {
			x_eps(i) += eps;
			double fp = objFun(x_eps);
			grad(i) = (fp - fx) / eps;
			x_eps(i) = x(i);
		}
	
		return grad;
	}

	/*
	 * Calculate a numerical approximation to the Hessian.
	 * Reference: http://objectmix.com/fortran/730003-how-calculate-hessian-matrix.html
	 * As noted in the gradient function, there are better
	 * approximation formulas available when we can get to it.
	**/
	virtual arma::vec hessEtaFun(const arma::vec& x, const arma::vec& eta) const
	{
		m_usedNumericalHessian = true;

		if (x.n_elem != eta.n_elem) {
			Rcpp::stop("eta must be same length as x");
		}

		double eps = m_epsNumericalHessEta;
		arma::vec Heta = arma::zeros(x.n_elem);

		arma::vec x_00 = x;
		arma::vec x_01 = x;
		arma::vec x_10 = x;
		arma::vec x_11 = x;

		for (size_t i = 0; i < x.n_elem; ++i) {
			for (size_t j = 0; j < x.n_elem; ++j) {
				x_00(i) -= eps;
				x_00(j) -= eps;
				x_01(i) -= eps;
				x_01(j) += eps;
				x_10(i) += eps;
				x_10(j) -= eps;
				x_11(i) += eps;
				x_11(j) += eps;

				double H_ij = (objFun(x_11) - objFun(x_10) - objFun(x_01) +
					objFun(x_00)) / (4*eps*eps);
				Heta(i) += H_ij * eta(j);

				x_00(i) = x(i);
				x_00(j) = x(j);
				x_01(i) = x(i);
				x_01(j) = x(j);
				x_10(i) = x(i);
				x_10(j) = x(j);
				x_11(i) = x(i);
				x_11(j) = x(j);
			}
		}

		return Heta;
	}

protected:
	double m_epsNumericalGrad;
	double m_epsNumericalHessEta;
	mutable bool m_usedNumericalHessian;
};

#endif

