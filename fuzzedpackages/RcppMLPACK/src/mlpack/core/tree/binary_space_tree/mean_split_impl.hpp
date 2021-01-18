/**
 * @file mean_split_impl.hpp
 * @author Yash Vadalia
 * @author Ryan Curtin
 *
 * Implementation of class(MeanSplit) to split a binary space partition tree.
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
#ifndef __MLPACK_CORE_TREE_BINARY_SPACE_TREE_MEAN_SPLIT_IMPL_HPP
#define __MLPACK_CORE_TREE_BINARY_SPACE_TREE_MEAN_SPLIT_IMPL_HPP

#include "mean_split.hpp"

namespace mlpack {
namespace tree {

template<typename BoundType, typename MatType>
bool MeanSplit<BoundType, MatType>::SplitNode(const BoundType& bound,
                                              MatType& data,
                                              const size_t begin,
                                              const size_t count,
                                              size_t& splitDimension,
                                              size_t& splitCol)
{
  splitDimension = data.n_rows; // Indicate invalid.
  double maxWidth = -1;

  // Find the split dimension.
  for (size_t d = 0; d < data.n_rows; d++)
  {
    double width = bound[d].Width();

    if (width > maxWidth)
    {
      maxWidth = width;
      splitDimension = d;
    }
  }

  if (maxWidth == 0) // All these points are the same.  We can't split.
    return false;

  // Split in the middle of that dimension.
  double splitVal = bound[splitDimension].Mid();

  // Perform the actual splitting.  This will order the dataset such that points
  // with value in dimension splitDimension less than or equal to splitVal are
  // on the left of splitCol, and points with value in dimension splitDimension
  // greater than splitVal are on the right side of splitCol.
  splitCol = PerformSplit(data, begin, count, splitDimension, splitVal);

  return true;
}

template<typename BoundType, typename MatType>
bool MeanSplit<BoundType, MatType>::SplitNode(const BoundType& bound,
                                              MatType& data,
                                              const size_t begin,
                                              const size_t count,
                                              size_t& splitDimension,
                                              size_t& splitCol,
                                              std::vector<size_t>& oldFromNew)
{
  splitDimension = data.n_rows; // Indicate invalid.
  double maxWidth = -1;

  // Find the split dimension.
  for (size_t d = 0; d < data.n_rows; d++)
  {
    double width = bound[d].Width();

    if (width > maxWidth)
    {
      maxWidth = width;
      splitDimension = d;
    }
  }

  if (maxWidth == 0) // All these points are the same.  We can't split.
    return false;

  // Split in the middle of that dimension.
  double splitVal = bound[splitDimension].Mid();

  // Perform the actual splitting.  This will order the dataset such that points
  // with value in dimension splitDimension less than or equal to splitVal are
  // on the left of splitCol, and points with value in dimension splitDimension
  // greater than splitVal are on the right side of splitCol.
  splitCol = PerformSplit(data, begin, count, splitDimension, splitVal,
      oldFromNew);

  return true;
}

template<typename BoundType, typename MatType>
size_t MeanSplit<BoundType, MatType>::
    PerformSplit(MatType& data,
                 const size_t begin,
                 const size_t count,
                 const size_t splitDimension,
                 const double splitVal)
{
  // This method modifies the input dataset.  We loop both from the left and
  // right sides of the points contained in this node.  The points less than
  // splitVal should be on the left side of the matrix, and the points greater
  // than splitVal should be on the right side of the matrix.
  size_t left = begin;
  size_t right = begin + count - 1;

  // First half-iteration of the loop is out here because the termination
  // condition is in the middle.
  while ((data(splitDimension, left) < splitVal) && (left <= right))
    left++;
  while ((data(splitDimension, right) >= splitVal) && (left <= right))
    right--;

  while (left <= right)
  {
    // Swap columns.
    data.swap_cols(left, right);

    // See how many points on the left are correct.  When they are correct,
    // increase the left counter accordingly.  When we encounter one that isn't
    // correct, stop.  We will switch it later.
    while ((data(splitDimension, left) < splitVal) && (left <= right))
      left++;

    // Now see how many points on the right are correct.  When they are correct,
    // decrease the right counter accordingly.  When we encounter one that isn't
    // correct, stop.  We will switch it with the wrong point we found in the
    // previous loop.
    while ((data(splitDimension, right) >= splitVal) && (left <= right))
      right--;
  }

  return left;
}

template<typename BoundType, typename MatType>
size_t MeanSplit<BoundType, MatType>::
    PerformSplit(MatType& data,
                 const size_t begin,
                 const size_t count,
                 const size_t splitDimension,
                 const double splitVal,
                 std::vector<size_t>& oldFromNew)
{
  // This method modifies the input dataset.  We loop both from the left and
  // right sides of the points contained in this node.  The points less than
  // splitVal should be on the left side of the matrix, and the points greater
  // than splitVal should be on the right side of the matrix.
  size_t left = begin;
  size_t right = begin + count - 1;

  // First half-iteration of the loop is out here because the termination
  // condition is in the middle.
  while ((data(splitDimension, left) < splitVal) && (left <= right))
    left++;
  while ((data(splitDimension, right) >= splitVal) && (left <= right))
    right--;

  while (left <= right)
  {
    // Swap columns.
    data.swap_cols(left, right);

    // Update the indices for what we changed.
    size_t t = oldFromNew[left];
    oldFromNew[left] = oldFromNew[right];
    oldFromNew[right] = t;  

    // See how many points on the left are correct.  When they are correct,
    // increase the left counter accordingly.  When we encounter one that isn't
    // correct, stop.  We will switch it later.
    while ((data(splitDimension, left) < splitVal) && (left <= right))
      left++;

    // Now see how many points on the right are correct.  When they are correct,
    // decrease the right counter accordingly.  When we encounter one that isn't
    // correct, stop.  We will switch it with the wrong point we found in the
    // previous loop.
    while ((data(splitDimension, right) >= splitVal) && (left <= right))
      right--;
  }

  return left;
}

}; // namespace tree
}; // namespace mlpack

#endif
