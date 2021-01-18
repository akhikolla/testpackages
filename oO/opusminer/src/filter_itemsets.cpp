/* filter_itemsets.cpp - a module of OPUS Miner providing filter_itemsets, a function to filter itemsets that are not indpendently productive.
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

#include <algorithm>

#include "filter_itemsets.h"
#include "globals.h"
#include "utils.h"

bool sizegt(itemset i1,  itemset i2) {
  return i1.size() > i2.size();
}

// check all combinations of intersections of partitions of tidsavail moved to either tidsleft or tidsright
bool checkSS2(std::vector<tidset> &uniqueTids, const int no, tidset &tidsleft, tidset &tidsright, const int availabletids,const int count, double alpha) {
  if (no == 0) {
    if (fisherTest(availabletids-tidsleft.size()-tidsright.size()+count, tidsleft.size()-count, tidsright.size()-count, count) > alpha) {
      return false;
    }
    else return true;
  }

   // first try with the tidset committed to the left then try with it committed to the right

  tidset newtids;
  intersection(newtids, uniqueTids[no-1], tidsleft);

  if (!checkSS2(uniqueTids, no-1, newtids, tidsright, availabletids, count, alpha)) {
    return false;
  }

  intersection(newtids, uniqueTids[no-1], tidsright);

  if (!checkSS2(uniqueTids, no-1, tidsleft, newtids, availabletids, count, alpha)) {
    return false;
  }

  return true;
}

// check whether itemset is is self sufficient given that it has supersets that cover the TIDs in supsettids
bool checkSS(itemset &is, tidset &supsettids) {
  bool result = true;

  // find for each item in is the TIDs that it covers that are not in supsettids
  std::vector<tidset> uniqueTids;
  uniqueTids.resize(is.size());


  int i;
  itemset::const_iterator it;

  for (it = is.begin(), i = 0; it != is.end(); it++, i++) {
    uniqueTids[i].resize(tids[*it].size());
    TID *ut_end = std::set_difference(tids[*it].begin(), tids[*it].end(), supsettids.begin(), supsettids.end(), &uniqueTids[i][0]);
    uniqueTids[i].resize(ut_end-&uniqueTids[i][0]);
    if (uniqueTids[i].size() == 0) {
      // there cannot be a significant association from adding this tidset
      result = false;
      break;
    }
  }

  if (result) {
    // set up a process that will check whether uniqueCov.size() is significantly greater than can be predicted by assuming independence between any partition of is
    tidset uniqueCov(uniqueTids[0]); // this is the TIDs covered by is that are not in supsettids
    int is_size = (int)is.size();
    // calculate uniqueCov
    for (i = 1; i < is_size; i++) {
      dintersection(uniqueCov, uniqueTids[i]);
    }

    tidset tidsright(uniqueTids[uniqueTids.size()-1]);  // this is the cover of the items committed to the right - initialise it to the last unique TID

    // start with the last item committed to the right, then successively commit eeach item first to the left then to the right
    for (i = uniqueTids.size()-2; i >= 0; i--) {
      result = checkSS2(uniqueTids, i, uniqueTids[i], tidsright, noOfTransactions-supsettids.size(), uniqueCov.size(), getAlpha(is.size()));

      if (result == false) return false;

      if (i > 0) {
        dintersection(tidsright, uniqueTids[i]);
      }
    }
  }

  return result;
}

// check whether itemsets can be explained by their supersets
void filter_itemsets(std::vector<itemsetRec> &is) {
  if (!is.empty()) {
    // Sort the itemsets so that the largest are first.
    // This way we only need to look at the itemsets before one that we are processing to find superset.
    // Also, we will determine whether a superset is self sufficient before trying to determine whether its subsets are
    std::sort(is.begin(), is.end(), sizegt);

    std::vector<itemsetRec>::iterator subset_it;

    itemset supitems;   // the additional items n the supersets of the current itemset
    tidset supsettids;  // the tids covered by the supitems
    tidset thissupsettids;  // the tids covered by the supitems

    for (subset_it = is.begin()+1; subset_it != is.end(); subset_it++) {
      // get the TIDs that are covered by the current itemset's supersets
      supsettids.clear();

      std::vector<itemsetRec>::iterator supset_it;
      for (supset_it = is.begin(); supset_it != subset_it; supset_it++) {
        if (supset_it->self_sufficient) {
          supitems.clear();
          if (subset(*subset_it, *supset_it)) {
            itemsetRec::const_iterator it;

            for (it = supset_it->begin(); it != supset_it->end(); it++) {
              if (subset_it->find(*it) == subset_it->end()) {
                supitems.insert(*it);
              }
            }

            if (!supitems.empty()) {
              gettids(supitems, thissupsettids);

              if (supsettids.empty()) {
                supsettids = thissupsettids;
              }
              else {
                dunion(supsettids, thissupsettids);
              }
            }
          }
        }
      }

      if (!supsettids.empty() && !checkSS(*subset_it, supsettids)) {
        // only call chechSS if one or more supersets were found (and hence TIDs retrieved
        subset_it->self_sufficient = false;
      }
    }
  }
}
