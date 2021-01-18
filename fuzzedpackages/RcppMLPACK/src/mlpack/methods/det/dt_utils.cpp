/**
 * @file dt_utils.cpp
 * @author Parikshit Ram (pram@cc.gatech.edu)
 *
 * This file implements functions to perform different tasks with the Density
 * Tree class.
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
#include "dt_utils.hpp"

using namespace mlpack;
using namespace det;

void mlpack::det::PrintLeafMembership(DTree* dtree,
                                      const arma::mat& data,
                                      const arma::Mat<size_t>& labels,
                                      const size_t numClasses,
                                      const std::string leafClassMembershipFile)
{
  // Tag the leaves with numbers.
  int numLeaves = dtree->TagTree();

  arma::Mat<size_t> table(numLeaves, (numClasses + 1));
  table.zeros();

  for (size_t i = 0; i < data.n_cols; i++)
  {
    const arma::vec testPoint = data.unsafe_col(i);
    const int leafTag = dtree->FindBucket(testPoint);
    const size_t label = labels[i];
    table(leafTag, label) += 1;
  }

  if (leafClassMembershipFile == "")
  {
    Rcpp::Rcout << "Leaf membership; row represents leaf id, column represents "
        << "class id; value represents number of points in leaf in class."
        << std::endl << table;
  }
  else
  {
    // Create a stream for the file.
    std::ofstream outfile(leafClassMembershipFile.c_str());
    if (outfile.good())
    {
      outfile << table;
      Rcpp::Rcout << "Leaf membership printed to '" << leafClassMembershipFile
          << "'." << std::endl;
    }
    else
    {
      Rcpp::Rcout << "Can't open '" << leafClassMembershipFile << "' to write "
          << "leaf membership to." << std::endl;
    }
    outfile.close();
  }

  return;
}


void mlpack::det::PrintVariableImportance(const DTree* dtree,
                                          const std::string viFile)
{
  arma::vec imps;
  dtree->ComputeVariableImportance(imps);

  double max = 0.0;
  for (size_t i = 0; i < imps.n_elem; ++i)
    if (imps[i] > max)
      max = imps[i];

  Rcpp::Rcout << "Maximum variable importance: " << max << "." << std::endl;

  if (viFile == "")
  {
    Rcpp::Rcout << "Variable importance: " << std::endl << imps.t() << std::endl;
  }
  else
  {
    std::ofstream outfile(viFile.c_str());
    if (outfile.good())
    {
      outfile << imps;
      Rcpp::Rcout << "Variable importance printed to '" << viFile << "'."
          << std::endl;
    }
    else
    {
      Rcpp::Rcout << "Can't open '" << viFile << "' to write variable importance "
          << "to." << std::endl;
    }
    outfile.close();
  }
}


// This function trains the optimal decision tree using the given number of
// folds.
DTree* mlpack::det::Trainer(arma::mat& dataset,
                            const size_t folds,
                            const bool useVolumeReg,
                            const size_t maxLeafSize,
                            const size_t minLeafSize,
                            const std::string unprunedTreeOutput)
{
  // Initialize the tree.
  DTree* dtree = new DTree(dataset);

  // Prepare to grow the tree...
  arma::Col<size_t> oldFromNew(dataset.n_cols);
  for (size_t i = 0; i < oldFromNew.n_elem; i++)
    oldFromNew[i] = i;

  // Save the dataset since it would be modified while growing the tree.
  arma::mat newDataset(dataset);

  // Growing the tree
  double oldAlpha = 0.0;
  double alpha = dtree->Grow(newDataset, oldFromNew, useVolumeReg, maxLeafSize,
      minLeafSize);

  Rcpp::Rcout << dtree->SubtreeLeaves() << " leaf nodes in the tree using full "
      << "dataset; minimum alpha: " << alpha << "." << std::endl;

  // Compute densities for the training points in the full tree, if we were
  // asked for this.
  if (unprunedTreeOutput != "")
  {
    std::ofstream outfile(unprunedTreeOutput.c_str());
    if (outfile.good())
    {
      for (size_t i = 0; i < dataset.n_cols; ++i)
      {
        arma::vec testPoint = dataset.unsafe_col(i);
        outfile << dtree->ComputeValue(testPoint) << std::endl;
      }
    }
    else
    {
      Rcpp::Rcout << "Can't open '" << unprunedTreeOutput << "' to write computed"
          << " densities to." << std::endl;
    }

    outfile.close();
  }

  // Sequentially prune and save the alpha values and the values of c_t^2 * r_t.
  std::vector<std::pair<double, double> > prunedSequence;
  while (dtree->SubtreeLeaves() > 1)
  {
    std::pair<double, double> treeSeq(oldAlpha,
        dtree->SubtreeLeavesLogNegError());
    prunedSequence.push_back(treeSeq);
    oldAlpha = alpha;
    alpha = dtree->PruneAndUpdate(oldAlpha, dataset.n_cols, useVolumeReg);

    // Some sanity checks.
    //Log::Assert((alpha < std::numeric_limits<double>::max()) ||
    //    (dtree->SubtreeLeaves() == 1));
    //Log::Assert(alpha > oldAlpha);
    //Log::Assert(dtree->SubtreeLeavesLogNegError() < treeSeq.second);
  }

  std::pair<double, double> treeSeq(oldAlpha,
      dtree->SubtreeLeavesLogNegError());
  prunedSequence.push_back(treeSeq);

  Rcpp::Rcout << prunedSequence.size() << " trees in the sequence; maximum alpha:"
      << " " << oldAlpha << "." << std::endl;

  delete dtree;

  arma::mat cvData(dataset);
  size_t testSize = dataset.n_cols / folds;

  std::vector<double> regularizationConstants;
  regularizationConstants.resize(prunedSequence.size(), 0);

  // Go through each fold.
  for (size_t fold = 0; fold < folds; fold++)
  {
    // Break up data into train and test sets.
    size_t start = fold * testSize;
    size_t end = std::min((fold + 1) * testSize, (size_t) cvData.n_cols);

    arma::mat test = cvData.cols(start, end - 1);
    arma::mat train(cvData.n_rows, cvData.n_cols - test.n_cols);

    if (start == 0 && end < cvData.n_cols)
    {
      train.cols(0, train.n_cols - 1) = cvData.cols(end, cvData.n_cols - 1);
    }
    else if (start > 0 && end == cvData.n_cols)
    {
      train.cols(0, train.n_cols - 1) = cvData.cols(0, start - 1);
    }
    else
    {
      train.cols(0, start - 1) = cvData.cols(0, start - 1);
      train.cols(start, train.n_cols - 1) = cvData.cols(end, cvData.n_cols - 1);
    }

    // Initialize the tree.
    DTree* cvDTree = new DTree(train);

    // Getting ready to grow the tree...
    arma::Col<size_t> cvOldFromNew(train.n_cols);
    for (size_t i = 0; i < cvOldFromNew.n_elem; i++)
      cvOldFromNew[i] = i;

    // Grow the tree.
    oldAlpha = 0.0;
    alpha = cvDTree->Grow(train, cvOldFromNew, useVolumeReg, maxLeafSize,
        minLeafSize);

    // Sequentially prune with all the values of available alphas and adding
    // values for test values.  Don't enter this loop if there are less than two
    // trees in the pruned sequence.
    for (size_t i = 0;
         i < ((prunedSequence.size() < 2) ? 0 : prunedSequence.size() - 2); ++i)
    {
      // Compute test values for this state of the tree.
      double cvVal = 0.0;
      for (size_t j = 0; j < test.n_cols; j++)
      {
        arma::vec testPoint = test.unsafe_col(j);
        cvVal += cvDTree->ComputeValue(testPoint);
      }

      // Update the cv regularization constant.
      regularizationConstants[i] += 2.0 * cvVal / (double) dataset.n_cols;

      // Determine the new alpha value and prune accordingly.
      oldAlpha = 0.5 * (prunedSequence[i + 1].first +
          prunedSequence[i + 2].first);
      alpha = cvDTree->PruneAndUpdate(oldAlpha, train.n_cols, useVolumeReg);
    }

    // Compute test values for this state of the tree.
    double cvVal = 0.0;
    for (size_t i = 0; i < test.n_cols; ++i)
    {
      arma::vec testPoint = test.unsafe_col(i);
      cvVal += cvDTree->ComputeValue(testPoint);
    }

    if (prunedSequence.size() > 2)
      regularizationConstants[prunedSequence.size() - 2] += 2.0 * cvVal /
          (double) dataset.n_cols;

    test.reset();
    delete cvDTree;
  }

  double optimalAlpha = -1.0;
  long double cvBestError = -std::numeric_limits<long double>::max();

  for (size_t i = 0; i < prunedSequence.size() - 1; ++i)
  {
    // We can no longer work in the log-space for this because we have no
    // guarantee the quantity will be positive.
    long double thisError = -std::exp((long double) prunedSequence[i].second) +
        (long double) regularizationConstants[i];

    if (thisError > cvBestError)
    {
      cvBestError = thisError;
      optimalAlpha = prunedSequence[i].first;
    }
  }

  Rcpp::Rcout << "Optimal alpha: " << optimalAlpha << "." << std::endl;

  // Initialize the tree.
  DTree* dtreeOpt = new DTree(dataset);

  // Getting ready to grow the tree...
  for (size_t i = 0; i < oldFromNew.n_elem; i++)
    oldFromNew[i] = i;

  // Save the dataset since it would be modified while growing the tree.
  newDataset = dataset;

  // Grow the tree.
  oldAlpha = -DBL_MAX;
  alpha = dtreeOpt->Grow(newDataset, oldFromNew, useVolumeReg, maxLeafSize,
      minLeafSize);

  // Prune with optimal alpha.
  while ((oldAlpha < optimalAlpha) && (dtreeOpt->SubtreeLeaves() > 1))
  {
    oldAlpha = alpha;
    alpha = dtreeOpt->PruneAndUpdate(oldAlpha, newDataset.n_cols, useVolumeReg);

    // Some sanity checks.
    //Log::Assert((alpha < std::numeric_limits<double>::max()) ||
    //    (dtreeOpt->SubtreeLeaves() == 1));
    //Log::Assert(alpha > oldAlpha);
  }

  Rcpp::Rcout << dtreeOpt->SubtreeLeaves() << " leaf nodes in the optimally "
      << "pruned tree; optimal alpha: " << oldAlpha << "." << std::endl;

  return dtreeOpt;
}
