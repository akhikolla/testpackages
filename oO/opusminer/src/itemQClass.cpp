/* itemQClass.cpp - a module of OPUS Miner providing the methods for the itemQClass class.
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

#include "itemQClass.h"
#include "globals.h"

bool iqeGreater(itemQElem iqe1, itemQElem iqe2) {
  return iqe1.ubVal > iqe2.ubVal;
}

itemQClass::itemQClass(void)
{
}

itemQClass::~itemQClass(void)
{
}

void itemQClass::insert(float ubVal, itemID item) {
  const int initialSize = size();

  resize(initialSize+1);

  if (initialSize == 0) {
    at(0).ubVal = ubVal;
    at(0).item = item;
  }
  else {
    int first = 0;
    int last = initialSize - 1;

    while (first < last) {
      const int mid = first + (last - first) / 2;
      if (ubVal <= at(mid).ubVal) {
        first = mid + 1;
      }
      else {
        last = mid;
      }
    }

    if (at(first).ubVal >= ubVal) {
      // this should only happen if all items in the queue have lower value than the new item
      first++;
    }

    for (last = initialSize; last > first; last--) {
      at(last) = at(last-1);
    }

    at(first).ubVal = ubVal;
    at(first).item = item;
  }
}
