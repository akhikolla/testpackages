#include "wrap.h"

#include <Rcpp.h>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

namespace Rcpp {
  using boost::uuids::uuid;

  template<>
  SEXP wrap(const std::vector<uuid>& x) {
    StringVector res = no_init(x.size());
    std::transform(x.begin(), x.end(), res.begin(), [](const uuid& u) { return to_string(u); });
    return wrap(res);
  }
} // Rcpp
