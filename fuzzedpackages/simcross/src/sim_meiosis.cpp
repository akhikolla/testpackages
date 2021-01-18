#include <Rcpp.h>
using namespace Rcpp;

#include "random.h"
#include "sim_meiosis.h"

// [[Rcpp::export(".sim_crossovers")]]
NumericVector sim_crossovers(const double L, const int m, const double p,
                             const bool obligate_chiasma, const double Lstar)
{
    if(m==0) { // no-interference model is a lot easier
        int n_xo;

        if(obligate_chiasma) {
            // rejection sampling to get at least one chiasma
            while((n_xo = R::rpois(Lstar/50.0)) == 0);

            n_xo = R::rbinom((double)n_xo, 0.5);
        }
        else {
            n_xo = R::rpois(L/100.0);
        }

        NumericVector tmp(0);
        if(n_xo > 0) tmp = runif(n_xo, 0.0, L);
        return tmp.sort();
    }

    int n_points, first, n_nichi, n_ichi;

    double lambda1 = Lstar/50.0 * (double)(m+1) * (1.0 - p);
    double lambda2 = Lstar/50.0 * p;

    while(1) {
        // chiasma and intermediate points
        n_points = R::rpois(lambda1);

        // which point is the first chiasma?
        first = random_int(0, m);
        if(first > n_points) n_ichi = 0;
        else n_ichi = n_points/(m+1) + (int)(first < (n_points % (m+1)));

        // no. chiasma from no interference process
        if(p > 0) n_nichi = R::rpois(lambda2);
        else n_nichi = 0;

        if(!obligate_chiasma || n_ichi + n_nichi > 0) break;
    }

    // locations of chiasmata and intermediate points for process w/ interference
    NumericVector point_locations(0);
    if(n_points > 0) point_locations = runif(n_points, 0.0, L);
    point_locations.sort();

    // move every (m+1)st point back to front
    int n_chi=0;
    for(int j=first; j < n_points; j += (m+1), n_chi++)
        point_locations[n_chi] = point_locations[j];

    // chiasma locations from non-interference process
    NumericVector nichi_locations(0);
    if(n_nichi > 0) nichi_locations = runif(n_nichi, 0.0, L);

    // combine interference and no interference chiasma locations
    NumericVector chi_locations(n_chi + n_nichi);
    std::copy(point_locations.begin(), point_locations.begin()+n_chi, chi_locations.begin());
    std::copy(nichi_locations.begin(), nichi_locations.end(), chi_locations.begin()+n_chi);
    chi_locations.sort();

    // thin by 1/2
    int n_xo=0;
    for(int i=0; i<n_chi+n_nichi; i++) {
        if(R::unif_rand() < 0.5) { // flip coin -> chiasma
            chi_locations[n_xo] = chi_locations[i];
            n_xo++;
        }
    }

    NumericVector xo_locations(n_xo);
    std::copy(chi_locations.begin(), chi_locations.begin()+n_xo, xo_locations.begin());

    return xo_locations;
}



// [[Rcpp::export(".sim_meiosis")]]
List sim_meiosis(const List parent, const int m, const double p,
                 const bool obligate_chiasma, const double Lstar)
{
    const double tol=1e-12; // for comparison of chr lengths in parents

    List mat, pat;
    mat = parent[0];
    pat = parent[1];

    IntegerVector matalle = mat[0];
    NumericVector matloc  = mat[1];

    IntegerVector patalle = pat[0];
    NumericVector patloc  = pat[1];

    double L = max(matloc);
    if(fabs(L - max(patloc)) > tol)
        throw std::range_error("parent's two chromosomes are not the same length");

    // simulate crossover locations; add -1 to the beginning
    NumericVector tmp = sim_crossovers(L, m, p, obligate_chiasma, Lstar);
    NumericVector product(tmp.size() + 1);
    product[0] = -1.0;
    std::copy(tmp.begin(), tmp.end(), product.begin()+1);

    int cur_allele = random_int(0, 1); // first allele (0 or 1)

    int biggest_length = product.size() + matloc.size() + patloc.size();
    NumericVector loc(biggest_length);
    IntegerVector alle(biggest_length);

    int curpos = 0;
    if(product.size()==1) {
        if(cur_allele==0) return mat;
        else return pat;
    }
    else {
        int i, j;
        for(i=1; i<product.size(); i++) {

            if(cur_allele==0) { // mat chr
                for(j=0; j<matloc.size(); j++) {
                    if(matloc[j] >= product[i-1] && matloc[j] < product[i]) {
                        loc[curpos] = matloc[j];
                        alle[curpos] = matalle[j];
                        curpos++;
                    }
                    else if(matloc[j] > product[i]) break;
                }
                loc[curpos] = product[i];
                alle[curpos] = matalle[j];
                curpos++;
            }
            else { // pat chr
                for(j=0; j<patloc.size(); j++) {
                    if(patloc[j] >= product[i-1] && patloc[j] < product[i]) {
                        loc[curpos] = patloc[j];
                        alle[curpos] = patalle[j];
                        curpos++;
                    }
                    else if(patloc[j] > product[i]) break;
                }
                loc[curpos] = product[i];
                alle[curpos] = patalle[j];
                curpos++;
            }

            cur_allele = 1 - cur_allele;

        }

        double lastxo = max(product);

        if(cur_allele==0) { // mat chr
            for(j=0; j<matloc.size(); j++) {
                if(matloc[j] > lastxo) {
                    loc[curpos] = matloc[j];
                    alle[curpos] = matalle[j];
                    curpos++;
                }
            }
        }
        else { // pat chr
            for(j=0; j<patloc.size(); j++) {
                if(patloc[j] > lastxo) {
                    loc[curpos] = patloc[j];
                    alle[curpos] = patalle[j];
                    curpos++;
                }
            }
        }
    }

    if(curpos > 1) { // clean up repeated alleles

        NumericVector loc_clean(curpos);
        IntegerVector alle_clean(curpos);

        loc_clean[0] = loc[0];
        alle_clean[0] = alle[0];
        int lastpos=0;

        for(int i=1; i<curpos; i++) {
            if(alle_clean[lastpos] == alle[i]) {
                loc_clean[lastpos] = loc[i];
            }
            else {
                lastpos++;
                loc_clean[lastpos] = loc[i];
                alle_clean[lastpos] = alle[i];
            }
        }
        curpos = lastpos+1;
        loc = loc_clean;
        alle = alle_clean;
    }

    // copy over to short vectors
    NumericVector loc_result(curpos);
    IntegerVector alle_result(curpos);
    std::copy(loc.begin(), loc.begin()+curpos, loc_result.begin());
    std::copy(alle.begin(), alle.begin()+curpos, alle_result.begin());

    return List::create(Named("alleles")= alle_result, Named("locations")=loc_result);
}
