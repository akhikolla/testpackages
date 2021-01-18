/* find_itemsets.cpp - a module of OPUS Miner providing find_itemsets, a function to find productive non-redundant itemsets.
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

#include <Rcpp.h>

#include <limits>
#include <iterator>
#include <map>
#include <algorithm>

#include "globals.h"
#include "find_itemsets.h"
#include "utils.h"
#include "fisher.h"
#include "itemQClass.h"

// the minimum leverage of an itemset in the top-k so far
// any itemset whose leverage does not exceed this value cannot enter the top-k
float minValue = -std::numeric_limits<float>::max();

// for each itemset explored for which supersets might be in the best k, keep the count
std::map<itemset, int> TIDCount;

#ifdef SIXTEENBIT
// special map for pairs in order to save space
std::map<int, int> TIDPairCount;
#endif

// access function for TIDCount
inline bool getTIDCount(itemset is, int &count) {
  if (is.size() == 1) {
    count = tids[*is.begin()].size();
    return true;
  }
  #ifdef SIXTEENBIT
  else if (is.size() == 2) {
    std::map<int, int>::const_iterator it = TIDPairCount.find((*is.begin()<<16) + *is.rbegin());

    if (it == TIDPairCount.end()) {
      count = 0;
      return false;
    }
    else {
      count = it->second;
      return true;
    }
  }
  #endif
  else {
    std::map<itemset, int>::const_iterator it = TIDCount.find(is);

    if (it == TIDCount.end()) {
      count = 0;
      return false;
    }
    else {
      count = it->second;
      return true;
    }
  }
}

// check whether the count was stored. If not, we know that we have already determined that no superset can be in the best k
inline bool checkTIDCount(itemset is) {
  if (is.size() == 1) {
    // all single itemset counts are stored
    return true;
  }
  #ifdef SIXTEENBIT
  else if (is.size()==2) {
    std::map<int, int>::const_iterator it = TIDPairCount.find((*is.begin()<<16) + *is.rbegin());

    if (it == TIDPairCount.end()) return false;
    else return true;
  }
  #endif
  else {
    std::map<itemset, int>::const_iterator it = TIDCount.find(is);

    if (it == TIDCount.end()) return false;
    else return true;
  }
}

// array of element values used by itemgt
float *sortval;

// for sorting an array of items on sortval
int itemgt(const void *i1, const void *i2) {
  if (sortval[*static_cast<const itemID*>(i1)] > sortval[*static_cast<const itemID*>(i2)]) return -1;
  else return 1;
}

// for sorting an array of items on sortval
int itemlt(const void *i1, const void *i2) {
  if (sortval[*static_cast<const itemID*>(i1)] < sortval[*static_cast<const itemID*>(i2)]) return -1;
  else return 1;
}

void checkImmediateSubsets(itemset &is, const int isCnt, bool &redundant, bool &apriori) {
  itemset subset = is;
  itemset::const_iterator it;

  redundant = false;
  apriori = false;

  for (it = is.begin(); it != is.end(); it++) {
    int subsetCnt;

    subset.erase(*it);

    if (!getTIDCount(subset, subsetCnt)) {
      redundant = true;
      apriori = true;
      return;
    }

    if (redundancyTests && subsetCnt == isCnt) {
      redundant = true;
    }

    subset.insert(*it);
  }

  return;
}

// calculates leverage, p, whether the itemset is is redundant and whether it is possible to determine that all supersets of is will be redundant
// return true iff is is not redundant, val > minValue and p <= alpha
bool checkSubsetsX(itemset &sofar, itemset &remaining, const itemID limit, const int cnt, const double new_sup, float &val, double &p, const double alpha) {
  int sofarCnt;
  int remainingCnt;

  if (!getTIDCount(sofar, sofarCnt) || !getTIDCount(remaining, remainingCnt)) {
    return false;
  }

  // do test for sofar against remaining
  const float this_val = searchByLift ? new_sup / (countToSup(remainingCnt) * countToSup(sofarCnt))
                                      : new_sup - countToSup(remainingCnt) * countToSup(sofarCnt);

  if (this_val < val) {
    val = this_val;
    if (this_val <= minValue) return false;
  }

  const double this_p = fisher(cnt, sofarCnt, remainingCnt);

  if (this_p > p) {
    p = this_p;
    if (p > alpha) {
      return false;
    }
  }


  if (remaining.size() > 1) {
    itemset new_remaining(remaining);

    itemset::const_iterator it;

    for (it = remaining.begin(); it != remaining.end() && *it < limit; it++) {
      sofar.insert(*it);
      new_remaining.erase(*it);

      if (!checkSubsetsX(sofar, new_remaining, *it, cnt, new_sup, val, p, alpha)) {
        return false;
      }

      sofar.erase(*it);
      new_remaining.insert(*it);
    }
  }

  return p <= alpha  && val > minValue;
}

// calculates leverage, p, whether is is redundant and whether it is possible to determine that all supersets of is will be redundant
// return true iff is is not redundant, val > minValue and p <= alpha
bool checkSubsets(itemID item, itemset &is, const int cnt, const double new_sup, const int parentCnt, const double parentSup, float &val, double &p, const double alpha) {

  // do test for new item against the rest
  const int itemCnt = tids[item].size();

  val = searchByLift ? new_sup / (parentSup * itemSup(item))
                     : new_sup - parentSup * itemSup(item);

  if (val <= minValue) return false;

  p = fisher(cnt, itemCnt, parentCnt);

  if (p > alpha) return false;

  if (is.size() > 2) {
    itemset sofar;
    itemset remaining(is);

    sofar.insert(item);
    remaining.erase(item);

    itemset::const_iterator it;

    for (it = is.begin(); it != is.end(); it++) {
      if (*it != item) {
        sofar.insert(*it);
        remaining.erase(*it);

        if (!checkSubsetsX(sofar, remaining, *it, cnt, new_sup, val, p, alpha)) {
          return false;
        }

        sofar.erase(*it);
        remaining.insert(*it);
      }
    }
  }

  return p <= alpha  && val > minValue;
}

// insert is into the collection of k best itemsets
inline void insert_itemset(itemsetRec &is) {
  if (itemsets.size() >= k) {
    itemsets.pop();
  }
  itemsets.push(is);
  if (itemsets.size() == k) {
    const float newMin = itemsets.top().value;
    if (newMin > minValue) {
      minValue = newMin;
    }
  }
}

// perform OPUS search for specialisations of is (which covers cover) using the candidates in queue q
// maxItemSup is the maximum of the supports of all individual items in is
void opus(itemsetRec &is, tidset &cover, itemQClass &q, const int maxItemCount) {
  unsigned int i;
  const float parentSup = countToSup(cover.size());
  const int depth = is.size()+1;

  tidset newCover;
  itemQClass newQ;

  for (i = 0; i < q.size(); i++) {
    const itemID item = q[i].item;
    int count;

    // determine the number of TIDs that the new itemset covers
    intersection(newCover, cover, tids[item]);
    count = newCover.size();

    const int newMaxItemCount = std::max(maxItemCount, static_cast<int>(tids[item].size()));
    const float new_sup = countToSup(count);

    // this is a lower bound on the p value that may be obtained for this itemset or any superset
    const p_value lb_p = fisher(count, newMaxItemCount, count);

    // calculate an upper bound on the value that can be obtained by this itemset or any superset

    const float ubval = searchByLift ? ((count == 0) ? 0.0 : (1.0 / countToSup(maxItemCount)))
                                      : new_sup - new_sup * countToSup(maxItemCount);

    // performing OPUS pruning - if this test fails, the item will not be included in any superset of is
    if (lb_p <= getAlpha(depth) && ubval > minValue) {
      // only continue if there is any possibility of this itemset or its supersets entering the list of best itemsets
      float val;
      double p;
      bool redundant;
      bool apriori;

      is.insert(item);

      checkImmediateSubsets(is, count, redundant, apriori);

      if (!apriori) {
        if (checkSubsets(item, is, count, new_sup, cover.size(), parentSup, val, p, getAlpha(depth))) {
          is.count = count;
          is.value = val;
          is.p = p;
          insert_itemset(is);
        }

        // performing OPUS pruning - if this test fails, the item will not be included in any superset of is
        if (!redundant) {
	  #ifdef SIXTEENBIT
	  if (is.size()==2)
	    TIDPairCount[(*is.begin()<<16) + *is.rbegin()] = count;
	  else
	  #endif
          TIDCount[is] = count;

          if (!newQ.empty()) {
            // there are only more nodes to expand if there is a queue of items to consider expanding it with
            opus(is, newCover, newQ, newMaxItemCount);
          }

          newQ.insert(ubval, item);
        }
      }

      is.erase(item);
    }
  }
}

void find_itemsets() {
  itemQClass q; // a queue of items, to be sorted on an upper bound on value
  int i;
  unsigned int j;

  // initalise q - the queue of items ordered on an upper bound on value
  for (i = 0; i < noOfItems; i++) {
    const int c = tids[i].size();

    const float sup = countToSup(c);
    const float ubVal = searchByLift ? 1.0 / sup
                                     : sup - sup * sup;

    // make sure that the support is high enough for it to be possible to create a significant itemset
    if (fisher(c, c, c) <= getAlpha(2)) {
      q.append(ubVal, i); // it is faster to sort the q once it is full rather than doing an insertion sort
    }
  }

  q.sort();

  itemQClass newq;  // this is the queue of items that will be available for the item currently being explored

  newq.insert(q[0].ubVal, q[0].item);

  float prevMinVal = minValue; // remember the current minValue, and output an update if it improves in this iteration of the loop

  itemsetRec is;

  // we are stepping through all associations of i with j<i, so the first value of i that will have effect is 1
  for (j = 1; j < q.size() && q[j].ubVal > minValue; j++) {
    const itemID item = q[j].item;

    is.clear();
    is.insert(item);

    opus(is, tids[item], newq, tids[item].size());

    newq.append(q[j].ubVal, item);

    if (prevMinVal < minValue) {
      Rcpp::Rcout << "<" << minValue << ">";
      prevMinVal = minValue;
    }
    else Rcpp::Rcout << ".";
  }

}
