/* utils.h - header file for the utils.cpp module of OPUS Miner.
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

#include <set>
#include "globals.h"
#include "fisher.h"

#include <set>
#include "globals.h"

// true iff s1 is a subset of s2
// assumes that sets are stored in ascending order
template <class Type>
inline bool subset(Type &s1, Type &s2)
{ typename Type::const_iterator it1 = s1.begin();
  typename Type::const_iterator it2 = s2.begin();

  while (it1 != s1.end()) {
    if (it2 == s2.end()) return false;
    if (*it1 < *it2) return false;
    if (*it1 == *it2) {
      it1++;
    }
    it2++;
  }
  return true;
}

// get the tidset for an itemset
inline void gettids(const itemset &is, tidset &t) {

  itemset::const_iterator it = is.begin();

  if (is.size() == 1) {
    t = tids[*it];
  }
  else {
    const itemID item1 = *it++;
    const itemID item2 = *it++;
    intersection(t, tids[item1], tids[item2]);

    while (it != is.end()) {
      dintersection(t, tids[*it++]);
    }
  }
}

inline float countToSup(const int count) {
  return count / static_cast<float>(noOfTransactions);
}

extern int getNum(const char *str);

inline double itemSup(const itemID item) {
  return countToSup(tids[item].size());
}

// return the result of a fFisher exact test for an itemset i with support count count relative to support counts count1 and count2 for two subsets s1 and s2 that form a partition of i
inline double fisher(const int count, const int count1, const int count2) {
  return fisherTest(noOfTransactions-count1-count2+count, count1-count, count2-count, count);
}
