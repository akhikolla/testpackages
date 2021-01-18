/* itemQClass.h - header file for the itemQClass.cpp module of OPUS Miner.
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

// This is a queue of items that is sorted in descending order on 
#include <vector>
#include <algorithm>

#include "opus_miner.h"

class itemQElem {
public:
  float ubVal;
  itemID item;
};

extern bool iqeGreater(itemQElem iqe1, itemQElem iqe2);

class itemQClass :
  public std::vector<itemQElem>
{
public:
  itemQClass(void);
  ~itemQClass(void);

  void insert(float ubVal, itemID item);
  inline void append(float ubVal, itemID item) {
    const int initialSize = size();

    resize(initialSize+1);

    at(initialSize).ubVal = ubVal;
    at(initialSize).item = item;
  }

  inline void sort() {
    std::sort(begin(), end(), iqeGreater);
  }
};
