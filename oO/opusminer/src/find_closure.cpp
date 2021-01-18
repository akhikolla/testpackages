/* find_closure.cpp - a module of OPUS Miner providing find_closure, a function to find the closure of an itemset.
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

#include "find_closure.h"
#include "utils.h"

void find_closure(const itemset &is, itemset &closure) {
  tidset thistids;

  closure = is;

  gettids(is, thistids);

  itemID item;

  for (item = 0; item < noOfItems; item++) {
    if (tids[item].size() >= thistids.size() && is.find(item) == is.end() && count_intersection(thistids, tids[item]) == thistids.size()) {
      closure.insert(item);
    }
  }
}

