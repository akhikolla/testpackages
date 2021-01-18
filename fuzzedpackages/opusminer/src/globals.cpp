/* globals.cpp - a module of OPUS Miner providing global variable declarations.
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
*/

#include <cmath>
#include "globals.h"

unsigned int k = 100;  // the number of associations to return
bool filter = true;   // if true perform a filter for self-sufficiency
bool correctionForMultCompare = true; // if true we should correct alpha for the size of the search space
int noOfTransactions = 0;
itemID noOfItems = 0;
std::vector<tidset> tids;
std::vector<double> alpha;
std::vector<std::string> itemNames;

bool searchByLift = false;
bool redundancyTests = true;
bool printClosures = false;

void expandAlpha(const unsigned int depth) {
  if (alpha.empty()) {
    // alpha[0[ and [1] are not used.
    alpha.push_back(1.0);
    alpha.push_back(1.0);
    if (depth <= 1) return;
  }

  if (static_cast<int>(depth) > noOfItems) alpha.push_back(0.0);
  else if (static_cast<int>(depth) == noOfItems) alpha.push_back(alpha[depth-1]); // at deepest level so might as well use as much of the rest of the probability mass as possible
  else {
    unsigned int i;
    for (i = alpha.size(); i <= depth; i++) {
      alpha.push_back(std::min((std::pow(0.5, static_cast<int>(depth-1)) / std::exp(log_combin(noOfItems, depth))) * 0.05, alpha[depth-1]));
    }
  }
}
