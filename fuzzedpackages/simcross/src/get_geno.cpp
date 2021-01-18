#include <Rcpp.h>
using namespace Rcpp;

#include "get_geno.h"

// [[Rcpp::export(".get_geno")]]
IntegerMatrix get_geno(const List xodat, const double position)
{
    int n_ind = xodat.size();
    IntegerMatrix result(n_ind, 2);

    for(int i=0; i<n_ind; i++) {
        List ind = xodat[i];
        List mat = ind[0];
        List pat = ind[1];

        // convert mat and pat chr to genotype matrices
        result[i] = get_one_geno(mat, position);
        result[i+n_ind] = get_one_geno(pat, position);
    }

    return result;
}


int get_one_geno(const List chr, const double position)
{
    IntegerVector alleles = chr[0];
    NumericVector locations = chr[1];

    int n_loc = locations.size();

    // find closest position to right and take allele there
    for(int j=0; j<n_loc; j++) {
        if(locations[j] >= position)
            return alleles[j];
    }

    // shouldn't get here
    throw std::range_error("No location to right of position");
    return -1; // really won't get here
}
