#include <Rcpp.h>

using namespace Rcpp;

// helper function for breaking apart macro data lists per sub geography
// param macro_data a list of macro datasets from \code{pull_synth_data}
// param n_geo The number of sub geographies; equivalent to nrow(macro_data[[1]])
// [[Rcpp::export]]
List disaggregate_mdCPP(const List macro_data) {
  R_xlen_t n_list = macro_data.size(), n_geo = as<NumericMatrix>(macro_data[0]).nrow();
  List macro_data2(n_geo);

  for (R_xlen_t nr = 0; nr < n_geo; nr++) {
    List temp (n_list);
    for (R_xlen_t ix = 0; ix < n_list; ix++) {
      NumericMatrix mtemp = as<NumericMatrix>(macro_data[ix]);
      temp[ix] = mtemp.row(nr);
    }
    macro_data2[nr] = temp;
  }

  return macro_data2;
}

