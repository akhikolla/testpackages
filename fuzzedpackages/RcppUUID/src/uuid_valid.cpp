#include "wrap.h"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/string_generator.hpp>
#include <Rcpp.h>

using namespace Rcpp;
using boost::uuids::string_generator;
using boost::uuids::uuid;

inline bool try_uuid(const char* s) {
  try {
    string_generator gen;
    uuid id = gen(s);
    return id.is_nil() || id.version() != uuid::version_unknown;
  } catch(...) {
    return false;
  }
}

//' @title Validate UUIDs
//'
//' @description
//' Function validates uuids.
//'
//' @param x Character vector.
//' @return Logical vector.
//'
//' @export
//'
//' @examples
//' # validate UUIDs
//' uuid_validate(uuid_generate_random(2))
//' uuid_validate(uuid_generate_nil(2))
//' uuid_validate(uuid_generate_name(c("one", "two")))
//' uuid_validate(c("a", ""))
//'
// [[Rcpp::export(rng=false)]]
LogicalVector uuid_validate(StringVector x) {
  LogicalVector res = no_init(x.size());
  std::transform(x.begin(), x.end(), res.begin(), try_uuid);
  return res;
}
