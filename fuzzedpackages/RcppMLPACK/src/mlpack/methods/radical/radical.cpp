/**
 * @file radical.cpp
 * @author Nishant Mehta
 *
 * Implementation of Radical class
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

#include "radical.hpp"

using namespace std;
using namespace arma;
using namespace mlpack;
using namespace mlpack::radical;

// Set the parameters to RADICAL.
Radical::Radical(const double noiseStdDev,
                 const size_t replicates,
                 const size_t angles,
                 const size_t sweeps,
                 const size_t m) :
    noiseStdDev(noiseStdDev),
    replicates(replicates),
    angles(angles),
    sweeps(sweeps),
    m(m)
{
  // Nothing to do here.
}

void Radical::CopyAndPerturb(mat& xNew, const mat& x) const
{
  //Timer::Start("radical_copy_and_perturb");
  xNew = repmat(x, replicates, 1) + noiseStdDev * randn(replicates * x.n_rows,
      x.n_cols);
  //Timer::Stop("radical_copy_and_perturb");
}


double Radical::Vasicek(vec& z) const
{
  z = sort(z);

  // Apparently slower.
  /*
  vec logs = log(z.subvec(m, z.n_elem - 1) - z.subvec(0, z.n_elem - 1 - m));
  //vec val = sum(log(z.subvec(m, nPoints - 1) - z.subvec(0, nPoints - 1 - m)));
  return (double) sum(logs);
  */

  // Apparently faster.
  double sum = 0;
  uword range = z.n_elem - m;
  for (uword i = 0; i < range; i++)
  {
    sum += log(z(i + m) - z(i));
  }

  return sum;
}


double Radical::DoRadical2D(const mat& matX)
{
  CopyAndPerturb(perturbed, matX);

  mat::fixed<2, 2> matJacobi;

  vec values(angles);

  for (size_t i = 0; i < angles; i++)
  {
    const double theta = (i / (double) angles) * M_PI / 2.0;
    const double cosTheta = cos(theta);
    const double sinTheta = sin(theta);

    matJacobi(0, 0) = cosTheta;
    matJacobi(1, 0) = -sinTheta;
    matJacobi(0, 1) = sinTheta;
    matJacobi(1, 1) = cosTheta;

    candidate = perturbed * matJacobi;
    vec candidateY1 = candidate.unsafe_col(0);
    vec candidateY2 = candidate.unsafe_col(1);

    values(i) = Vasicek(candidateY1) + Vasicek(candidateY2);
  }

  uword indOpt;
  values.min(indOpt); // we ignore the return value; we don't care about it
  return (indOpt / (double) angles) * M_PI / 2.0;
}


void Radical::DoRadical(const mat& matXT, mat& matY, mat& matW)
{
  // matX is nPoints by nDims (although less intuitive than columns being
  // points, and although this is the transpose of the ICA literature, this
  // choice is for computational efficiency when repeatedly generating
  // two-dimensional coordinate projections for Radical2D).
  //Timer::Start("radical_transpose_data");
  mat matX = trans(matXT);
  //Timer::Stop("radical_transpose_data");

  // If m was not specified, initialize m as recommended in
  // (Learned-Miller and Fisher, 2003).
  if (m < 1)
    m = floor(sqrt((double) matX.n_rows));

  const size_t nDims = matX.n_cols;
  const size_t nPoints = matX.n_rows;

  //Timer::Start("radical_whiten_data");
  mat matXWhitened;
  mat matWhitening;
  WhitenFeatureMajorMatrix(matX, matY, matWhitening);
  //Timer::Stop("radical_whiten_data");
  // matY is now the whitened form of matX.

  // In the RADICAL code, they do not copy and perturb initially, although the
  // paper does.  We follow the code as it should match their reported results
  // and likely does a better job bouncing out of local optima.
  //GeneratePerturbedX(X, X);

  // Initialize the unmixing matrix to the whitening matrix.
  //Timer::Start("radical_do_radical");
  matW = matWhitening;

  mat matYSubspace(nPoints, 2);

  mat matJ = eye(nDims, nDims);

  for (size_t sweepNum = 0; sweepNum < sweeps; sweepNum++)
  {
    Rcpp::Rcout << "RADICAL: sweep " << sweepNum << "." << std::endl;

    for (size_t i = 0; i < nDims - 1; i++)
    {
      for (size_t j = i + 1; j < nDims; j++)
      {
        Rcpp::Rcout << "RADICAL 2D on dimensions " << i << " and " << j << "."
            << std::endl;

        matYSubspace.col(0) = matY.col(i);
        matYSubspace.col(1) = matY.col(j);

        const double thetaOpt = DoRadical2D(matYSubspace);

        const double cosThetaOpt = cos(thetaOpt);
        const double sinThetaOpt = sin(thetaOpt);

        // Set elements of transformation matrix.
        matJ(i, i) = cosThetaOpt;
        matJ(j, i) = -sinThetaOpt;
        matJ(i, j) = sinThetaOpt;
        matJ(j, j) = cosThetaOpt;

        matY *= matJ;

        // Unset elements of transformation matrix, so matJ = eye() again.
        matJ(i, i) = 1;
        matJ(j, i) = 0;
        matJ(i, j) = 0;
        matJ(j, j) = 1;
      }
    }
  }
  //Timer::Stop("radical_do_radical");

  // The final transposes provide W and Y in the typical form from the ICA
  // literature.
  //Timer::Start("radical_transpose_data");
  matW = trans(matW);
  matY = trans(matY);
  //Timer::Stop("radical_transpose_data");
}

void mlpack::radical::WhitenFeatureMajorMatrix(const mat& matX,
                                               mat& matXWhitened,
                                               mat& matWhitening)
{
  mat matU, matV;
  vec s;
  svd(matU, s, matV, cov(matX));
  matWhitening = matU * diagmat(1 / sqrt(s)) * trans(matV);
  matXWhitened = matX * matWhitening;
}

// return a string of this object.
std::string Radical::ToString() const
{
  std::ostringstream convert;
  convert << "Radical  [" << this << "]" << std::endl;
  convert << "  StdDev of Noise: " << noiseStdDev << std::endl;
  convert << "  Number of Replicates: " << replicates << std::endl;
  convert << "  Number of Angles: " << angles << std::endl;
  convert << "  M value: " << m << std::endl;
  return convert.str();
}
