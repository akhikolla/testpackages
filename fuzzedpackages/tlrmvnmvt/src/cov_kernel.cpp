#include <RcppEigen.h>
#include <Rmath.h>
#include <functional>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;
using namespace Eigen;

std::function<double(double)> exp_kernel(double beta)
{
        return [beta](double dist){return exp(-dist/beta);};
}

std::function<double(double)> matern_kernel(double rho, double nu)
{
	double c1 = pow(2.0, 1.0 - nu) / tgamma(nu);
	return [c1, rho, nu](double dist){return dist > 0.0 ? c1 * 
		pow(dist/rho, nu) * boost::math::cyl_bessel_k(nu, dist/rho) : 1.0;};
}
