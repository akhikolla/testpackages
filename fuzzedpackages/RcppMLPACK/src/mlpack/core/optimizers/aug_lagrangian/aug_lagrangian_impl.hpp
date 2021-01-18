/**
 * @file aug_lagrangian_impl.hpp
 * @author Ryan Curtin
 *
 * Implementation of AugLagrangian class (Augmented Lagrangian optimization
 * method).
 *
 * This file is part of MLPACK 1.0.10.
 *
 * MLPACK is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * MLPACK is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details (LICENSE.txt).
 *
 * You should have received a copy of the GNU General Public License along with
 * MLPACK.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __MLPACK_CORE_OPTIMIZERS_AUG_LAGRANGIAN_AUG_LAGRANGIAN_IMPL_HPP
#define __MLPACK_CORE_OPTIMIZERS_AUG_LAGRANGIAN_AUG_LAGRANGIAN_IMPL_HPP

#include <mlpack/core/optimizers/lbfgs/lbfgs.hpp>
#include "aug_lagrangian_function.hpp"

namespace mlpack {
namespace optimization {

template<typename LagrangianFunction>
AugLagrangian<LagrangianFunction>::AugLagrangian(LagrangianFunction& function) :
    function(function),
    augfunc(function),
    lbfgsInternal(augfunc),
    lbfgs(lbfgsInternal)
{
  lbfgs.MaxIterations() = 1000;
}

template<typename LagrangianFunction>
AugLagrangian<LagrangianFunction>::AugLagrangian(
    AugLagrangianFunction<LagrangianFunction>& augfunc,
    L_BFGSType& lbfgs) :
    function(augfunc.Function()),
    augfunc(augfunc),
    lbfgs(lbfgs)
{
  // Nothing to do.  lbfgsInternal isn't used in this case.
}

// This overload just sets the lambda and sigma and calls the other overload.
template<typename LagrangianFunction>
bool AugLagrangian<LagrangianFunction>::Optimize(arma::mat& coordinates,
                                                 const arma::vec& initLambda,
                                                 const double initSigma,
                                                 const size_t maxIterations)
{
  augfunc.Lambda() = initLambda;
  augfunc.Sigma() = initSigma;

  return Optimize(coordinates, maxIterations);
}

// Convert the object to a string.
template<typename LagrangianFunction>
std::string AugLagrangian<LagrangianFunction>::ToString() const
{
  std::ostringstream convert;
  convert << "AugLagrangian [" << this << "]" << std::endl;
  convert << "  Function:" << std::endl;
  convert << mlpack::util::Indent(function.ToString(), 2);
  convert << "  L-BFGS optimizer:" << std::endl;
  convert << mlpack::util::Indent(lbfgs.ToString(), 2);
  return convert.str();
}

template<typename LagrangianFunction>
bool AugLagrangian<LagrangianFunction>::Optimize(arma::mat& coordinates,
                                                 const size_t maxIterations)
{
  // Ensure that we update lambda immediately.
  double penaltyThreshold = DBL_MAX;

  // Track the last objective to compare for convergence.
  double lastObjective = function.Evaluate(coordinates);

  // Then, calculate the current penalty.
  double penalty = 0;
  for (size_t i = 0; i < function.NumConstraints(); i++)
    penalty += std::pow(function.EvaluateConstraint(i, coordinates), 2);

  Rcpp::Rcout << "Penalty is " << penalty << " (threshold " << penaltyThreshold
      << ")." << std::endl;

  // The odd comparison allows user to pass maxIterations = 0 (i.e. no limit on
  // number of iterations).
  size_t it;
  for (it = 0; it != (maxIterations - 1); it++)
  {
    Rcpp::Rcout << "AugLagrangian on iteration " << it
        << ", starting with objective "  << lastObjective << "." << std::endl;

 //   Rcpp::Rcout << coordinates << std::endl;

//    Rcpp::Rcout << trans(coordinates) * coordinates << std::endl;

    if (!lbfgs.Optimize(coordinates))
      Rcpp::Rcout << "L-BFGS reported an error during optimization."
          << std::endl;

    // Check if we are done with the entire optimization (the threshold we are
    // comparing with is arbitrary).
    if (std::abs(lastObjective - function.Evaluate(coordinates)) < 1e-10 &&
        augfunc.Sigma() > 500000)
      return true;

    lastObjective = function.Evaluate(coordinates);

    // Assuming that the optimization has converged to a new set of coordinates,
    // we now update either lambda or sigma.  We update sigma if the penalty
    // term is too high, and we update lambda otherwise.

    // First, calculate the current penalty.
    double penalty = 0;
    for (size_t i = 0; i < function.NumConstraints(); i++)
    {
      penalty += std::pow(function.EvaluateConstraint(i, coordinates), 2);
//      Rcpp::Rcout << "Constraint " << i << " is " <<
//          function.EvaluateConstraint(i, coordinates) << std::endl;
    }

    Rcpp::Rcout << "Penalty is " << penalty << " (threshold "
        << penaltyThreshold << ")." << std::endl;

    for (size_t i = 0; i < function.NumConstraints(); ++i)
    {
//      arma::mat tmpgrad;
//      function.GradientConstraint(i, coordinates, tmpgrad);
//      Rcpp::Rcout << "Gradient of constraint " << i << " is " << std::endl;
//      Rcpp::Rcout << tmpgrad << std::endl;
    }

    if (penalty < penaltyThreshold) // We update lambda.
    {
      // We use the update: lambda_{k + 1} = lambda_k - sigma * c(coordinates),
      // but we have to write a loop to do this for each constraint.
      for (size_t i = 0; i < function.NumConstraints(); i++)
        augfunc.Lambda()[i] -= augfunc.Sigma() *
            function.EvaluateConstraint(i, coordinates);

      // We also update the penalty threshold to be a factor of the current
      // penalty.  TODO: this factor should be a parameter (from CLI).  The
      // value of 0.25 is taken from Burer and Monteiro (2002).
      penaltyThreshold = 0.25 * penalty;
      Rcpp::Rcout << "Lagrange multiplier estimates updated." << std::endl;
    }
    else
    {
      // We multiply sigma by a constant value.  TODO: this factor should be a
      // parameter (from CLI).  The value of 10 is taken from Burer and Monteiro
      // (2002).
      augfunc.Sigma() *= 10;
      Rcpp::Rcout << "Updated sigma to " << augfunc.Sigma() << "." << std::endl;
    }
  }

  return false;
}

}; // namespace optimization
}; // namespace mlpack

#endif // __MLPACK_CORE_OPTIMIZERS_AUG_LAGRANGIAN_AUG_LAGRANGIAN_IMPL_HPP

