/**
 * @file decision_stump.hpp
 * @author Udit Saxena
 *
 * Definition of decision stumps.
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
#ifndef __MLPACK_METHODS_DECISION_STUMP_DECISION_STUMP_HPP
#define __MLPACK_METHODS_DECISION_STUMP_DECISION_STUMP_HPP

#include <mlpack/core.hpp>

namespace mlpack {
namespace decision_stump {

/**
 * This class implements a decision stump. It constructs a single level
 * decision tree, i.e., a decision stump. It uses entropy to decide splitting
 * ranges.
 *
 * The stump is parameterized by a splitting attribute (the dimension on which
 * points are split), a vector of bin split values, and a vector of labels for
 * each bin.  Bin i is specified by the range [split[i], split[i + 1]).  The
 * last bin has range up to \infty (split[i + 1] does not exist in that case).
 * Points that are below the first bin will take the label of the first bin.
 *
 * @tparam MatType Type of matrix that is being used (sparse or dense).
 */
template <typename MatType = arma::mat>
class DecisionStump
{
 public:
  /**
   * Constructor. Train on the provided data. Generate a decision stump from
   * data.
   *
   * @param data Input, training data.
   * @param labels Labels of training data.
   * @param classes Number of distinct classes in labels.
   * @param inpBucketSize Minimum size of bucket when splitting.
   */
  DecisionStump(const MatType& data,
                const arma::Row<size_t>& labels,
                const size_t classes,
                size_t inpBucketSize);

  /**
   * Classification function. After training, classify test, and put the
   * predicted classes in predictedLabels.
   *
   * @param test Testing data or data to classify.
   * @param predictedLabels Vector to store the predicted classes after
   *     classifying test data.
   */
  void Classify(const MatType& test, arma::Row<size_t>& predictedLabels);

  /**
   *
   *
   *
   *
   */
  DecisionStump(const DecisionStump<>& ds);

  /**
   *
   *
   *
   *
   *
   *
  ModifyData(MatType& data, const arma::Row<double>& D);
  */
  
  //! Access the splitting attribute.
  int SplitAttribute() const { return splitAttribute; }
  //! Modify the splitting attribute (be careful!).
  int& SplitAttribute() { return splitAttribute; }

  //! Access the splitting values.
  const arma::vec& Split() const { return split; }
  //! Modify the splitting values (be careful!).
  arma::vec& Split() { return split; }

  //! Access the labels for each split bin.
  const arma::Col<size_t> BinLabels() const { return binLabels; }
  //! Modify the labels for each split bin (be careful!).
  arma::Col<size_t>& BinLabels() { return binLabels; }

 private:
  //! Stores the number of classes.
  size_t numClass;

  //! Stores the value of the attribute on which to split.
  int splitAttribute;

  //! Size of bucket while determining splitting criterion.
  size_t bucketSize;

  //! Stores the splitting values after training.
  arma::vec split;

  //! Stores the labels for each splitting bin.
  arma::Col<size_t> binLabels;

  /**
   * Sets up attribute as if it were splitting on it and finds entropy when
   * splitting on attribute.
   *
   * @param attribute A row from the training data, which might be a
   *     candidate for the splitting attribute.
   */
  double SetupSplitAttribute(const arma::rowvec& attribute,
                             const arma::Row<size_t>& labels);

  /**
   * After having decided the attribute on which to split, train on that
   * attribute.
   *
   * @param attribute attribute is the attribute decided by the constructor
   *      on which we now train the decision stump.
   */
  template <typename rType> void TrainOnAtt(const arma::rowvec& attribute,
                                            const arma::Row<size_t>& labels);

  /**
   * After the "split" matrix has been set up, merge ranges with identical class
   * labels.
   */
  void MergeRanges();

  /**
   * Count the most frequently occurring element in subCols.
   *
   * @param subCols The vector in which to find the most frequently
   *     occurring element.
   */
  template <typename rType> rType CountMostFreq(const arma::Row<rType>& subCols);

  /**
   * Returns 1 if all the values of featureRow are not same.
   *
   * @param featureRow The attribute which is checked for identical values.
   */
  template <typename rType> int IsDistinct(const arma::Row<rType>& featureRow);

  /**
   * Calculate the entropy of the given attribute.
   *
   * @param attribute The attribute of which we calculate the entropy.
   * @param labels Corresponding labels of the attribute.
   */
  template <typename AttType, typename LabelType>
  double CalculateEntropy(arma::subview_row<LabelType> labels);
};

}; // namespace decision_stump
}; // namespace mlpack

#include "decision_stump_impl.hpp"

#endif
