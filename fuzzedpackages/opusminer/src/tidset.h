/* tidset.h - header file for the tidset.cpp module of OPUS Miner.
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

#pragma once

#include <vector>
#include <algorithm>

#include "opus_miner.h"

typedef std::vector<TID> tidset;

inline unsigned int count_intersection(tidset &s1, tidset &s2) {
  // count the size of the intersection
  // relies on the sets both being stored in ascending order

  if (s1.size() == 0 || s2.size() == 0) {
    return 0;
  }

  tidset::const_iterator it1 = s1.begin();
  TID v1 = *it1;
  const tidset::const_iterator end1 = s1.end();
  tidset::const_iterator it2 = s2.begin();
  TID v2 = *it2;
  const tidset::const_iterator end2 = s2.end();

  int count = 0;

  while (true) {
    if (v1 == v2) {
      count++;
      it1++;
      if (it1 == end1) break;
      v1 = *it1;
      it2++;
      if (it2 == end2) break;
      v2 = *it2;
    }
    else if (v1 < v2) {
      it1++;
      if (it1 == end1) break;
      v1 = *it1;
    }
    else {
      it2++;
      if (it2 == end2) break;
      v2 = *it2;
    }
  }

  return count;
}


// find the intersection of two tidsets
// relies on the sets both being stored in ascending order
inline void intersection(tidset &result, tidset &s1, tidset &s2) {
  result.clear();
  result.reserve(std::min(s1.size(), s2.size()));

  if (s1.size() == 0 || s2.size() == 0) {
    return;
  }

  tidset::const_iterator it1 = s1.begin();
  TID v1 = *it1;
  const tidset::const_iterator end1 = s1.end();
  tidset::const_iterator it2 = s2.begin();
  TID v2 = *it2;
  const tidset::const_iterator end2 = s2.end();

  while (true) {
    if (v1 == v2) {
      result.push_back(v1);
      it1++;
      if (it1 == end1) break;
      v1 = *it1;
      it2++;
      if (it2 == end2) break;
      v2 = *it2;
    }
    else if (v1 < v2) {
      it1++;
      if (it1 == end1) break;
      v1 = *it1;
    }
    else {
      it2++;
      if (it2 == end2) break;
      v2 = *it2;
    }
  }
}

// destructively update s1 to its intersection with s2
inline void dintersection(tidset &s1, tidset &s2) {
  if (s1.size() == 0) {
    return;
  }

  if (s2.size() == 0) {
    s1.clear();
    return;
  }

  unsigned int from = 0;
  unsigned int to = 0;
  TID v1 = s1[0];
  const tidset::const_iterator end1 = s1.end();
  tidset::const_iterator it2 = s2.begin();
  TID v2 = *it2;
  const tidset::const_iterator end2 = s2.end();

  while (true) {
    if (v1 == v2) {
      s1[to++] = s1[from++];
      if (from == s1.size()) break;
      v1 = s1[from];
      it2++;
      if (it2 == end2) break;
      v2 = *it2;
    }
    else if (v1 < v2) {
      from++;
      if (from == s1.size()) break;
      v1 = s1[from];
    }
    else {
      it2++;
      if (it2 == end2) break;
      v2 = *it2;
    }
  }

  s1.resize(to);
}

// destructively update s1 to its union with s2
inline void dunion(tidset &s1, tidset &s2) {
  tidset result;

  tidset::const_iterator it1 = s1.begin();
  tidset::const_iterator it2 = s2.begin();

  while (true) {
    if (it1 == s1.end()) {
      while (it2 != s2.end()) {
        result.push_back(*it2);
        it2++;
      }
      break;
    }
    else if (it2 == s2.end()) {
      while (it1 != s1.end()) {
        result.push_back(*it1);
        it1++;
      }
      break;
    }
    else if (*it1 == *it2) {
      result.push_back(*it1);
      it1++;
      it2++;
    }
    else if (*it1 < *it2) {
      result.push_back(*it1);
      it1++;
    }
    else {
      result.push_back(*it2);
      it2++;
    }
  }

  s1 = result;
}
