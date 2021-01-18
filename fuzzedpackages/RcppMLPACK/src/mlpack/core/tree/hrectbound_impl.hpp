/**
 * @file hrectbound_impl.hpp
 *
 * Implementation of hyper-rectangle bound policy class.
 * Template parameter Power is the metric to use; use 2 for Euclidean (L2).
 *
 * @experimental
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
#ifndef __MLPACK_CORE_TREE_HRECTBOUND_IMPL_HPP
#define __MLPACK_CORE_TREE_HRECTBOUND_IMPL_HPP

#include <math.h>

// In case it has not been included yet.
#include "hrectbound.hpp"

namespace mlpack {
namespace bound {

/**
 * Empty constructor.
 */
template<int Power, bool TakeRoot>
HRectBound<Power, TakeRoot>::HRectBound() :
    dim(0),
    bounds(NULL),
    minWidth(0)
{ /* Nothing to do. */ }

/**
 * Initializes to specified dimensionality with each dimension the empty
 * set.
 */
template<int Power, bool TakeRoot>
HRectBound<Power, TakeRoot>::HRectBound(const size_t dimension) :
    dim(dimension),
    bounds(new math::Range[dim]),
    minWidth(0)
{ /* Nothing to do. */ }

/***
 * Copy constructor necessary to prevent memory leaks.
 */
template<int Power, bool TakeRoot>
HRectBound<Power, TakeRoot>::HRectBound(const HRectBound& other) :
    dim(other.Dim()),
    bounds(new math::Range[dim]),
    minWidth(other.MinWidth())
{
  // Copy other bounds over.
  for (size_t i = 0; i < dim; i++)
    bounds[i] = other[i];
}

/***
 * Same as the copy constructor.
 */
template<int Power, bool TakeRoot>
HRectBound<Power, TakeRoot>& HRectBound<Power, TakeRoot>::operator=(
    const HRectBound& other)
{
  if (dim != other.Dim())
  {
    // Reallocation is necessary.
    if (bounds)
      delete[] bounds;

    dim = other.Dim();
    bounds = new math::Range[dim];
  }

  // Now copy each of the bound values.
  for (size_t i = 0; i < dim; i++)
    bounds[i] = other[i];

  minWidth = other.MinWidth();

  return *this;
}

/**
 * Destructor: clean up memory.
 */
template<int Power, bool TakeRoot>
HRectBound<Power, TakeRoot>::~HRectBound()
{
  if (bounds)
    delete[] bounds;
}

/**
 * Resets all dimensions to the empty set.
 */
template<int Power, bool TakeRoot>
void HRectBound<Power, TakeRoot>::Clear()
{
  for (size_t i = 0; i < dim; i++)
    bounds[i] = math::Range();
  minWidth = 0;
}

/***
 * Calculates the centroid of the range, placing it into the given vector.
 *
 * @param centroid Vector which the centroid will be written to.
 */
template<int Power, bool TakeRoot>
void HRectBound<Power, TakeRoot>::Centroid(arma::vec& centroid) const
{
  // Set size correctly if necessary.
  if (!(centroid.n_elem == dim))
    centroid.set_size(dim);

  for (size_t i = 0; i < dim; i++)
    centroid(i) = bounds[i].Mid();
}

/**
 * Calculates minimum bound-to-point squared distance.
 */
template<int Power, bool TakeRoot>
template<typename VecType>
double HRectBound<Power, TakeRoot>::MinDistance(
    const VecType& point,
    typename boost::enable_if<IsVector<VecType> >* /* junk */) const
{
  //Log::Assert(point.n_elem == dim);

  double sum = 0;

  double lower, higher;
  for (size_t d = 0; d < dim; d++)
  {
    lower = bounds[d].Lo() - point[d];
    higher = point[d] - bounds[d].Hi();

    // Since only one of 'lower' or 'higher' is negative, if we add each's
    // absolute value to itself and then sum those two, our result is the
    // nonnegative half of the equation times two; then we raise to power Power.
    sum += pow((lower + fabs(lower)) + (higher + fabs(higher)), (double) Power);
  }

  // Now take the Power'th root (but make sure our result is squared if it needs
  // to be); then cancel out the constant of 2 (which may have been squared now)
  // that was introduced earlier.  The compiler should optimize out the if
  // statement entirely.
  if (TakeRoot)
    return pow(sum, 1.0 / (double) Power) / 2.0;
  else
    return sum / pow(2.0, Power);
}

/**
 * Calculates minimum bound-to-bound squared distance.
 */
template<int Power, bool TakeRoot>
double HRectBound<Power, TakeRoot>::MinDistance(const HRectBound& other) const
{
  //Log::Assert(dim == other.dim);

  double sum = 0;
  const math::Range* mbound = bounds;
  const math::Range* obound = other.bounds;

  double lower, higher;
  for (size_t d = 0; d < dim; d++)
  {
    lower = obound->Lo() - mbound->Hi();
    higher = mbound->Lo() - obound->Hi();
    // We invoke the following:
    //   x + fabs(x) = max(x * 2, 0)
    //   (x * 2)^2 / 4 = x^2
    sum += pow((lower + fabs(lower)) + (higher + fabs(higher)), (double) Power);

    // Move bound pointers.
    mbound++;
    obound++;
  }

  // The compiler should optimize out this if statement entirely.
  if (TakeRoot)
    return pow(sum, 1.0 / (double) Power) / 2.0;
  else
    return sum / pow(2.0, Power);
}

/**
 * Calculates maximum bound-to-point squared distance.
 */
template<int Power, bool TakeRoot>
template<typename VecType>
double HRectBound<Power, TakeRoot>::MaxDistance(
    const VecType& point,
    typename boost::enable_if<IsVector<VecType> >* /* junk */) const
{
  double sum = 0;

  //Log::Assert(point.n_elem == dim);

  for (size_t d = 0; d < dim; d++)
  {
    double v = std::max(fabs(point[d] - bounds[d].Lo()),
        fabs(bounds[d].Hi() - point[d]));
    sum += pow(v, (double) Power);
  }

  // The compiler should optimize out this if statement entirely.
  if (TakeRoot)
    return pow(sum, 1.0 / (double) Power);
  else
    return sum;
}

/**
 * Computes maximum distance.
 */
template<int Power, bool TakeRoot>
double HRectBound<Power, TakeRoot>::MaxDistance(const HRectBound& other) const
{
  double sum = 0;

  //Log::Assert(dim == other.dim);

  double v;
  for (size_t d = 0; d < dim; d++)
  {
    v = std::max(fabs(other.bounds[d].Hi() - bounds[d].Lo()),
        fabs(bounds[d].Hi() - other.bounds[d].Lo()));
    sum += pow(v, (double) Power); // v is non-negative.
  }

  // The compiler should optimize out this if statement entirely.
  if (TakeRoot)
    return pow(sum, 1.0 / (double) Power);
  else
    return sum;
}

/**
 * Calculates minimum and maximum bound-to-bound squared distance.
 */
template<int Power, bool TakeRoot>
math::Range HRectBound<Power, TakeRoot>::RangeDistance(const HRectBound& other)
    const
{
  double loSum = 0;
  double hiSum = 0;

  //Log::Assert(dim == other.dim);

  double v1, v2, vLo, vHi;
  for (size_t d = 0; d < dim; d++)
  {
    v1 = other.bounds[d].Lo() - bounds[d].Hi();
    v2 = bounds[d].Lo() - other.bounds[d].Hi();
    // One of v1 or v2 is negative.
    if (v1 >= v2)
    {
      vHi = -v2; // Make it nonnegative.
      vLo = (v1 > 0) ? v1 : 0; // Force to be 0 if negative.
    }
    else
    {
      vHi = -v1; // Make it nonnegative.
      vLo = (v2 > 0) ? v2 : 0; // Force to be 0 if negative.
    }

    loSum += pow(vLo, (double) Power);
    hiSum += pow(vHi, (double) Power);
  }

  if (TakeRoot)
    return math::Range(pow(loSum, 1.0 / (double) Power),
                       pow(hiSum, 1.0 / (double) Power));
  else
    return math::Range(loSum, hiSum);
}

/**
 * Calculates minimum and maximum bound-to-point squared distance.
 */
template<int Power, bool TakeRoot>
template<typename VecType>
math::Range HRectBound<Power, TakeRoot>::RangeDistance(
    const VecType& point,
    typename boost::enable_if<IsVector<VecType> >* /* junk */) const
{
  double loSum = 0;
  double hiSum = 0;

  //Log::Assert(point.n_elem == dim);

  double v1, v2, vLo, vHi;
  for (size_t d = 0; d < dim; d++)
  {
    v1 = bounds[d].Lo() - point[d]; // Negative if point[d] > lo.
    v2 = point[d] - bounds[d].Hi(); // Negative if point[d] < hi.
    // One of v1 or v2 (or both) is negative.
    if (v1 >= 0) // point[d] <= bounds_[d].Lo().
    {
      vHi = -v2; // v2 will be larger but must be negated.
      vLo = v1;
    }
    else // point[d] is between lo and hi, or greater than hi.
    {
      if (v2 >= 0)
      {
        vHi = -v1; // v1 will be larger, but must be negated.
        vLo = v2;
      }
      else
      {
        vHi = -std::min(v1, v2); // Both are negative, but we need the larger.
        vLo = 0;
      }
    }

    loSum += pow(vLo, (double) Power);
    hiSum += pow(vHi, (double) Power);
  }

  if (TakeRoot)
    return math::Range(pow(loSum, 1.0 / (double) Power),
                       pow(hiSum, 1.0 / (double) Power));
  else
    return math::Range(loSum, hiSum);
}

/**
 * Expands this region to include a new point.
 */
template<int Power, bool TakeRoot>
template<typename MatType>
HRectBound<Power, TakeRoot>& HRectBound<Power, TakeRoot>::operator|=(
    const MatType& data)
{
  //Log::Assert(data.n_rows == dim);

  arma::vec mins(min(data, 1));
  arma::vec maxs(max(data, 1));

  minWidth = DBL_MAX;
  for (size_t i = 0; i < dim; i++)
  {
    bounds[i] |= math::Range(mins[i], maxs[i]);
    const double width = bounds[i].Width();
    if (width < minWidth)
      minWidth = width;
  }

  return *this;
}

/**
 * Expands this region to encompass another bound.
 */
template<int Power, bool TakeRoot>
HRectBound<Power, TakeRoot>& HRectBound<Power, TakeRoot>::operator|=(
    const HRectBound& other)
{
  assert(other.dim == dim);

  minWidth = DBL_MAX;
  for (size_t i = 0; i < dim; i++)
  {
    bounds[i] |= other.bounds[i];
    const double width = bounds[i].Width();
    if (width < minWidth)
      minWidth = width;
  }

  return *this;
}

/**
 * Determines if a point is within this bound.
 */
template<int Power, bool TakeRoot>
template<typename VecType>
bool HRectBound<Power, TakeRoot>::Contains(const VecType& point) const
{
  for (size_t i = 0; i < point.n_elem; i++)
  {
    if (!bounds[i].Contains(point(i)))
      return false;
  }

  return true;
}

/**
 * Returns the diameter of the hyperrectangle (that is, the longest diagonal).
 */
template<int Power, bool TakeRoot>
double HRectBound<Power, TakeRoot>::Diameter() const
{
  double d = 0;
  for (size_t i = 0; i < dim; ++i)
    d += std::pow(bounds[i].Hi() - bounds[i].Lo(), (double) Power);

  if (TakeRoot)
    return std::pow(d, 1.0 / (double) Power);
  else
    return d;
}

/**
 * Returns a string representation of this object.
 */
template<int Power, bool TakeRoot>
std::string HRectBound<Power, TakeRoot>::ToString() const
{
  std::ostringstream convert;
  convert << "HRectBound [" << this << "]" << std::endl;
  convert << "  Power: " << Power << std::endl;
  convert << "  TakeRoot: " << (TakeRoot ? "true" : "false") << std::endl;
  convert << "  Dimensionality: " << dim << std::endl;
  convert << "  Bounds: " << std::endl;
  for (size_t i = 0; i < dim; ++i)
    convert << util::Indent(bounds[i].ToString()) << std::endl;
  convert << "  Minimum width: " << minWidth << std::endl;

  return convert.str();
}

}; // namespace bound
}; // namespace mlpack

#endif // __MLPACK_CORE_TREE_HRECTBOUND_IMPL_HPP
