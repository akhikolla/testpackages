/* itemset.h - header file for the itemset.cpp module of OPUS Miner.
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

#include <set>

#include "opus_miner.h"

class itemset :
  public std::set<itemID>
{
public:
  itemset(void);
  ~itemset(void);

};

class itemsetRec :
  public itemset
{
public:
  itemsetRec(void);
  ~itemsetRec(void);

  // used for sorting itemsets
  const bool operator <(const itemsetRec& pI) const
  {
    return (value >  pI.value);
  }

  int count;
  float value;
  double p;
  bool self_sufficient;
};
