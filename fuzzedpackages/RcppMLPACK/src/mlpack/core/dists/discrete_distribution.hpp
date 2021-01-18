/**
 * @file discrete_distribution.hpp
 * @author Ryan Curtin
 *
 * Implementation of the discrete distribution, where each discrete observation
 * has a given probability.
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
#ifndef __MLPACK_METHODS_HMM_DISTRIBUTIONS_DISCRETE_DISTRIBUTION_HPP
#define __MLPACK_METHODS_HMM_DISTRIBUTIONS_DISCRETE_DISTRIBUTION_HPP

#include <mlpack/core.hpp>

namespace mlpack {
namespace distribution /** Probability distributions. */ {

/**
 * A discrete distribution where the only observations are discrete
 * observations.  This is useful (for example) with discrete Hidden Markov
 * Models, where observations are non-negative integers representing specific
 * emissions.
 *
 * No bounds checking is performed for observations, so if an invalid
 * observation is passed (i.e. observation > numObservations), a crash will
 * probably occur.
 *
 * This distribution only supports one-dimensional observations, so when passing
 * an arma::vec as an observation, it should only have one dimension
 * (vec.n_rows == 1).  Any additional dimensions will simply be ignored.
 *
 * @note
 * This class, like every other class in MLPACK, uses arma::vec to represent
 * observations.  While a discrete distribution only has positive integers
 * (size_t) as observations, these can be converted to doubles (which is what
 * arma::vec holds).  This distribution internally converts those doubles back
 * into size_t before comparisons.
 * @endnote
 */
class DiscreteDistribution
{
 public:
  /**
   * Default constructor, which creates a distribution that has no observations.
   */
  DiscreteDistribution() { /* nothing to do */ }

  /**
   * Define the discrete distribution as having numObservations possible
   * observations.  The probability in each state will be set to (1 /
   * numObservations).
   *
   * @param numObservations Number of possible observations this distribution
   *    can have.
   */
  DiscreteDistribution(const size_t numObservations) :
      probabilities(arma::ones<arma::vec>(numObservations) / numObservations)
  { /* nothing to do */ }

  /**
   * Define the discrete distribution as having the given probabilities for each
   * observation.
   *
   * @param probabilities Probabilities of each possible observation.
   */
  DiscreteDistribution(const arma::vec& probabilities)
  {
    // We must be sure that our distribution is normalized.
    double sum = accu(probabilities);
    if (sum > 0)
      this->probabilities = probabilities / sum;
    else
    {
      this->probabilities.set_size(probabilities.n_elem);
      this->probabilities.fill(1 / probabilities.n_elem);
    }
  }

  /**
   * Get the dimensionality of the distribution.
   */
  size_t Dimensionality() const { return 1; }

  /**
   * Return the probability of the given observation.  If the observation is
   * greater than the number of possible observations, then a crash will
   * probably occur -- bounds checking is not performed.
   *
   * @param observation Observation to return the probability of.
   * @return Probability of the given observation.
   */
  double Probability(const arma::vec& observation) const
  {
    // Adding 0.5 helps ensure that we cast the floating point to a size_t
    // correctly.
    const size_t obs = size_t(observation[0] + 0.5);

    // Ensure that the observation is within the bounds.
    if (obs >= probabilities.n_elem)
    {
      Rcpp::Rcout << "DiscreteDistribution::Probability(): received observation "
          << obs << "; observation must be in [0, " << probabilities.n_elem
          << "] for this distribution." << std::endl;
    }

    return probabilities(obs);
  }

  /**
   * Return a randomly generated observation (one-dimensional vector; one
   * observation) according to the probability distribution defined by this
   * object.
   *
   * @return Random observation.
   */
  arma::vec Random() const;

  /**
   * Estimate the probability distribution directly from the given observations.
   * If any of the observations is greater than numObservations, a crash is
   * likely to occur.
   *
   * @param observations List of observations.
   */
  void Estimate(const arma::mat& observations);

  /**
   * Estimate the probability distribution from the given observations, taking
   * into account the probability of each observation actually being from this
   * distribution.
   *
   * @param observations List of observations.
   * @param probabilities List of probabilities that each observation is
   *    actually from this distribution.
   */
  void Estimate(const arma::mat& observations,
                const arma::vec& probabilities);

  //! Return the vector of probabilities.
  const arma::vec& Probabilities() const { return probabilities; }
  //! Modify the vector of probabilities.
  arma::vec& Probabilities() { return probabilities; }

  /*
   * Returns a string representation of this object.
   */
  std::string ToString() const;

 private:
  arma::vec probabilities;
};

}; // namespace distribution
}; // namespace mlpack

#endif
