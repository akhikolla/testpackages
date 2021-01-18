/* globals.h - header file for the globals.cpp module of OPUS Miner.
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


#ifndef __GLOBALS__
#define __GLOBALS__

#include <set>
#include <vector>
#include <queue>
#include <string>
#include <cmath>
#include "itemset.h"
#include "tidset.h"
#include "fisher.h"

typedef double p_value;

extern int noOfTransactions;
extern itemID noOfItems;
extern std::vector<tidset> tids;
extern std::priority_queue<itemsetRec> itemsets;
extern unsigned int k; // the maximum number of itemsets to find
extern bool filter;   // if true perform a filter for self-sufficiency
extern bool correctionForMultCompare; // if true we should correct alpha for the size of the search space
extern std::vector<double> alpha;
extern std::vector<std::string> itemNames;

extern void expandAlpha(const unsigned int depth);

inline double getAlpha(const unsigned int depth) {
  if (!correctionForMultCompare) {
    return 0.05;
  }

  if (depth >= alpha.size()) {
    expandAlpha(depth);
  }

  return alpha[depth];
}

extern bool searchByLift;
extern bool redundancyTests;
extern bool printClosures;

#endif // __GLOBALS__
