#pragma once

#include <RcppCommon.h>
#include <boost/uuid/uuid.hpp>

namespace Rcpp {
  using boost::uuids::uuid;
  template<> SEXP wrap(const uuid&);
  template<> SEXP wrap(const std::vector<uuid>&);
} // Rcpp
