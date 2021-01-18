#include "wrap.h"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/nil_generator.hpp>
#include <Rcpp.h>

using namespace Rcpp;
using boost::uuids::nil_generator;
using boost::uuids::uuid;

//' @title Generates Nil UUIDs
//'
//' @description
//' Function generates nil uuids.
//'
//' @param n Number of generated UUIDs.
//' @return Character vector with UUIDs.
//'
//' @references
//' <https://www.boost.org/doc/libs/1_72_0/libs/uuid/doc/uuid.html#Nil%20Generator>
//'
//' @export
//'
//' @examples
//' # generate nil UUIDs
//' uuid_generate_nil(2)
//'
// [[Rcpp::export(rng=false)]]
StringVector uuid_generate_nil(size_t n = 1)  {
  std::vector<uuid> res(n);
  std::generate(res.begin(), res.end(), nil_generator());
  return wrap(res);
}
