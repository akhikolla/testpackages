// enums.h: Imprecise Classification Trees
//
// Copyright (C) 2018  Paul Fink
//
// This file is part of imptree.
//
// imptree is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// imptree is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with imptree.  If not, see <https://www.gnu.org/licenses/>.

#ifndef RCPP_IMPTREE_ENUMS_H
#define RCPP_IMPTREE_ENUMS_H

#include "translation.h"
#include <Rcpp.h>

enum class EntropyCorrection {
  no = 0,
  strobl,
  abellan,
};

enum class SplitMetric {
  entropyMax = 0,
  entropyRange,
}; 

enum class IpType {
  idm = 0,
  npi,
  npiapprox,
};

namespace IpTypeLookup {
  inline const char* toString(IpType v) {
    switch (v) {
      case IpType::idm:   return "IDM";
      case IpType::npi:   return "NPI";
      case IpType::npiapprox: return "NPIapprox";
    }
    throw Rcpp::exception(_("Only 'IDM', 'NPI' and 'NPIapprox' are supported as IpType"));
  }
}

enum class Dominance {
  interval = 0,
  maximality
};

#endif /*RCPP_IMPTREE_ENUMS_H*/
