#include "wrap.h"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/name_generator.hpp>
#include <Rcpp.h>

using namespace Rcpp;
using boost::uuids::name_generator_sha1;
using boost::uuids::ns::x500dn;
using boost::uuids::uuid;

//' @title Generate UUIDs Version 5
//'
//' @description
//' Function generates name-based uuid is derived from content in a namespace.
//' A uuid with identical content shall yield the same uuid.
//' Hashing algorithm is SHA1. Namespace is X.500 DN.
//'
//' @param x Character vector.
//' @return Character vector with UUIDs.
//'
//' @note
//' This function generates valid uuids for the `NA` and empty strings.
//'
//' @export
//'
//' @references
//' <https://www.boost.org/doc/libs/1_72_0/libs/uuid/doc/uuid.html#Name%20Generator>
//'
//' @export
//'
//' @examples
//' # generate name UUIDs
//' uuid_generate_name(c("one", "two"))
//'
// [[Rcpp::export(rng=false)]]
StringVector uuid_generate_name(StringVector x)  {
  std::vector<uuid> res(x.size());
  name_generator_sha1 gen(x500dn());
  std::transform(x.begin(), x.end(), res.begin(), gen);
  return wrap(res);
}
