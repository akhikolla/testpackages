#include "utils.h"

using namespace Rcpp;

namespace armspp {

DottedPair listToDottedPair(List input, int index) {
  DottedPair output(R_NilValue);
  Nullable<CharacterVector> names = input.names();
  for (int i = input.size() - 1; i >= 0; --i) {
    RObject value;
    if (::Rf_isVector(input[i]) && !::Rf_isMatrix(input[i])) {
      GenericVector vectorValue = static_cast<GenericVector>(input[i]);
      value = vectorValue[index % vectorValue.size()];
    } else {
      value = input[i];
    }

    if (names.isNotNull() && names.as()[i] != "") {
      output.push_front(Named(
        as<std::string>(names.as()[i]),
        value
      ));
    } else {
      output.push_front(value);
    }
  }

  return output;
}

};
