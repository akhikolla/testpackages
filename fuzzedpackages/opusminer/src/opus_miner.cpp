/* Open source implementation of the OPUS Miner algorithm which applies OPUS search for Filtered Top-k Association Discovery of Self-Sufficient Itemsets
** Copyright (C) 2012 Geoffrey I Webb
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
** Please report any bugs to Geoff Webb <geoff.webb@monash.edu>
*/

#include <Rcpp.h>

#include <cmath>
#include <limits>

#include "opus_miner.h"
#include "globals.h"
#include "find_itemsets.h"
#include "print_itemsets.h"
#include "filter_itemsets.h"
#include "utils.h"
#include "fisher.h"

std::priority_queue<itemsetRec> itemsets;

void init() {
  alpha = std::vector<double>();                 // globals.cpp
  tids = std::vector<tidset>();                  // globals.cpp
  itemsets = std::priority_queue<itemsetRec>();  // opus_miner.cpp
  itemNames = std::vector<std::string>();        // globals.cpp
  minValue = -std::numeric_limits<float>::max(); // find_itemsets.cpp
  TIDCount = std::map<itemset, int>();           // find_itemsets.cpp
}

Rcpp::GenericVector
#ifdef _WIN32
  __cdecl
#endif
opus(Rcpp::GenericVector tidList, int numItems, int numTrans, Rcpp::NumericVector k_, Rcpp::LogicalVector args) {

  init();
  tids = Rcpp::as< std::vector< tidset > >(tidList);
  noOfItems = numItems;
  noOfTransactions = numTrans;

  std::vector<itemsetRec> is;

  Rcpp::GenericVector output;

  k = Rcpp::as<unsigned int>(k_);

  printClosures = args[0];
  filter = args[1];
  searchByLift = args[2];
  correctionForMultCompare = args[3];
  redundancyTests = args[4];

  try {

    Rcpp::Rcout << "Finding itemsets ("
                << noOfTransactions << " transactions, "
                << noOfItems << " items)...\n\n";

    find_itemsets();

    Rcpp::Rcout << "\n\n";

    // extract the itemsets from the priority queue
    while (!itemsets.empty()) {
      is.push_back(itemsets.top());
      itemsets.pop();
    }

    if (filter) {
      Rcpp::Rcout << "Filtering itemsets...\n\n";
      filter_itemsets(is);
    }

    output = get_itemsets(is);
  }
  catch (const std::bad_alloc&) {
    Rcpp::Rcout << "Error: Out of memory.\n";
  }
  catch (...) {
    Rcpp::Rcout << "Error: Unhandled exception.\n";
  }

  return output;
}

// [[Rcpp::export(.opus_cpp)]]
Rcpp::GenericVector opus_cpp(Rcpp::GenericVector tidList, int numItems, int numTrans, Rcpp::NumericVector k_, Rcpp::LogicalVector args) {
  return opus(tidList, numItems, numTrans, k_, args);
}
