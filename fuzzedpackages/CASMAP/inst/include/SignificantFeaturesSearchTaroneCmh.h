/*
 * SignificantFeaturesSearchTaroneCmh.h
 *
 *  Created on: 8 May 2017
 *      Author: mikolajr
 */

#ifndef SIGNIFICANTINTERVALSEARCHTARONECMH_H_
#define SIGNIFICANTINTERVALSEARCHTARONECMH_H_

#include <array>
#include <vector>

/* LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */

#include "SignificantFeaturesSearchWithCovariates.h"
#include "Exception.h"
#include "types.h"

namespace SignificantPattern
{

class SignificantFeaturesSearchTaroneCmh : public SignificantFeaturesSearchWithCovariates
{
private:
    // super class pattern for code independence of changes in inheritance
    typedef SignificantFeaturesSearchWithCovariates super;

    void execute_constructor_taronecmh();
    void execute_destructor_taronecmh();

    // flag for marking locally whether initial arrays should be re-computed
    // due to the change of the covariates
    // Note: flag has to be local to class, since sub-, or super- class has to
    //       be able to mark itself that the file change has been (locally)
    //       accounted for
    long long lastCovariatesInitTime = 0;
    inline bool haveCovariatesChanged() const {
        return ( getCovariates().getLastInitTime() - lastCovariatesInitTime ) > 0;
    };
    inline void markCovariatesSeen() {
        lastCovariatesInitTime = getCovariates().getLastInitTime();
    }

protected:
    // Number of tables
    unsigned short K;
    // Number of observations per table
    std::vector<longint> Nt;
    // Number of observations in positive class per table
    std::vector<longint> nt;
    // Number of observations per table
    std::vector<longint> cum_Nt; //Cumulative sum of Nt
    // Now some precomputed quantities to save time when computing the CMH test
    // statistic
    std::vector<longint> Nt_nt; // Number of observations in zero class per table (Ni-ni for each of the K tables)
    std::vector<longint> hypercorner_bnd; // max(ni,Ni-ni) for each of the K tables
    inline virtual longint hypercorner_bnd_k(unsigned short k) {
        return (nt[k] > Nt_nt[k]) ? nt[k] : Nt_nt[k];
    }
    inline virtual bool notprunable_k(longint x_k, unsigned short k) {
        // If for any of the K tables, its margin x is smaller than the maximum of n
        // and N-n, then we cannot prune the interval (we are not in the "top-right"
        // hypercorner)
        return (x_k < hypercorner_bnd[k]);
    }
    inline virtual longint dim_margin_k(longint x_k, unsigned short k) {
        return (Nt[k]-x_k);
    }
    std::vector<double> gammat; // ni/Ni for each of the K tables
    std::vector<double> gammabint; // (ni/Ni)*(1-ni/Ni) for each of the K tables
    // And some others to save allocation time when computing the maximum CMH
    // test statistic in the "top right" hypercorner
    std::vector<double> f_vals, g_vals, betas;
    std::vector<unsigned short> idx_betas_sorted;


    static constexpr unsigned short NGRID=500; //Number of samples in the grid of tentative corrected significance thresholds, not counting 1 as it is a trivial threshold
    static constexpr double LOG10_MIN_PVAL=-30.0; //Minimum tentative corrected significance threshold = 10^{LOG10_MIN_PVAL}
    static constexpr double log10_p_step=-LOG10_MIN_PVAL/NGRID; // Step size in the grid

    // Grid of logarithmically spaced corrected significance thresholds. The sequence is ordered
    // from larger thresholds to smaller thresholds.
    std::array<double, NGRID+1> pgrid;
    // Index in the grid of the current tentative corrected significance threshold
    unsigned short idx_pgrid;

    // A LxK-dimensional matrix storing the frequency of each interval in each
    // table as they are processed
    longint **freq_par_cov;
    // A (NGRID+1)-dimensional vector such that freq_cnt_cmh[j] = #intervals
    // with maximum attainable CMH statistic in the j-th bucket
    longint *freq_cnt_cmh;

    /* -------------------------------- INITIALISATION AND TERMINATION METHODS ----------------------------------------- */
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    virtual void algorithm_init() override;
    virtual void algorithm_end() override;

    void freq_init() override;
    void freq_clear() override;
    void freq_constructor() override;
    void freq_destructor() override;

    /* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */
    double compute_pval(longint a, longint* x);

    double compute_score(longint a, longint* x);

    double compute_odds_ratio(longint* at, longint* x);

    double score_to_pval(double score);

    unsigned short bucket_idx(double pval);

    double compute_minpval(longint *x);

    void idx_betas_sort(unsigned short j);
    /*
     * Given the margins of the K tables and the minimum attainable CMH p-value,
     * checks whether the interval and all other intervals containing it can be
     * pruned from the search space
     * */
    double compute_envelope_minpval(longint *x);
    inline bool isprunable_freqcov(longint *x) {
        return isprunable(compute_envelope_minpval(x));
    };
    inline bool istestable_freqcov(longint *x) {
        return istestable(compute_minpval(x));
    };

    void decrease_threshold() override;

public:
    SignificantFeaturesSearchTaroneCmh();
    virtual ~SignificantFeaturesSearchTaroneCmh();
};




} /* namespace SignificantPattern */

#endif /* SIGNIFICANTINTERVALSEARCHTARONECMH_H_ */
