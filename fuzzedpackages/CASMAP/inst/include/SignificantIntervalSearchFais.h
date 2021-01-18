/*
 * SignificantIntervalSearchFais.h
 *
 *  Created on: 31 Mar 2017
 *      Author: mikolajr
 */

#ifndef LIBSIGINTERVALSEARCH_SIGNIFICANTINTERVALSEARCHFAIS_H_
#define LIBSIGINTERVALSEARCH_SIGNIFICANTINTERVALSEARCHFAIS_H_

#include "SignificantIntervalSearch.h"
#include "FeatureSet.h"
#include "Summary.h"
#include "types.h"

namespace SignificantPattern {

class SignificantIntervalSearchFais: public SignificantIntervalSearch {
private:
    // super class pattern for code independence of changes in inheritance
    typedef SignificantIntervalSearch super;

    void execute_constructor_fais();
    void execute_destructor_fais();

    //IntervalSetWithFreq pValsTestableInts;
    //IntervalSetWithFreq pValsSigInts;

    IntervalSetWithOddsRatio pValsTestableInts;
    IntervalSetWithOddsRatio pValsSigInts;

    SummaryFais summary;

protected:
    inline virtual SummaryFais& getSummaryRef() override { return summary; }
    inline virtual void summary_constructor() override { summary = SummaryFais(); }

    /// @todo
    /// All table pointers which have fixed size, known once input files have
    /// been read (i.e. N or L-dependent), can be re-factored in the
    /// execute-dependent class that uses stack instead of heap mem for them
    /// (`eltype arr[n]` instead of `eltype *arr`); alternatively: use
    /// std::vector (`std::vector<eltype> arr(n)`) just to remove the
    /// error-prone manual mem management
    /**
     * A L-dimensional vector storing the frequency of each interval as they are
     * processed.
     */
    longint *freq_par;
    /**
     * A (N+1)-dimensional vector such that freq_cnt[j] = #intervals with
     * x_{i}=j processed so far.
     */
    longint *freq_cnt;

    /**
     * Array with all values of minimum attainable P-value in [0,N]
     * pre-computed.
     */
    double *psi;

    /**
     * @name Region thresholds
     *
     * Sigma_k = [sl1,sl2] U [su1,su2]
     *
     * We always have `su1==N-sl2` and `su2==N-sl1`, but keeping each variable
     * separate reduces the number of computations at the expense of a tiny
     * amount of memory.
     */
    ///@{
    longint sl1, sl2, su1, su2;
    /**
     * Flag variable to keep track of the last change done to region Sigma_k:
     * if `flag==1`, the last change was `sl1++` (shrink on extremes of the W);
     * if `flag==0`, the last change was `sl2--` (shrink on center of the W).
     */
    int flag;
    ///@}

    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    /**
     * @copydoc super::algorithm_init()
     *
     * Additionally init minimum attainable p-values and region thresholds.
     */
    virtual void algorithm_init() override;
    virtual void algorithm_end() override;

    void freq_init() override;
    void freq_clear() override;
    void freq_constructor() override;
    void freq_destructor() override;

    inline void saveSignificantInterval(double pval, double score, double odds_ratio, longint tau, longint l, longint a) override {
        pValsSigInts.addFeature(tau, tau+l, a, score, odds_ratio, pval);
    }
    inline void saveTestableInterval(double pval, double score, double odds_ratio, longint tau, longint l, longint a) override {
        pValsTestableInts.addFeature(tau, tau+l, a, score, odds_ratio, pval);
    }

    /**
     * Allocate memory and initialise (pre-compute) minimum attainable p-values
     * psi(x) for all x in [0,N] and store them in array psi.
     *
     * Called in
     * algorithm_init().
     */
    void psi_init();
    /**
     * Precompute minimum attainable p-values #psi(x) for all x in [0,N].
     *
     * Called in psi_init().
     */
    virtual void psi_clear() = 0;
    /**
     * Cleanup memory after minimum attainable p-values psi(x) for all x in
     * [0,N] and zero-initialise the pointer.
     *
     * Called in execute_destructor_int().
     */
    void psi_destructor();
    /**
     * Zero-initialise pointers to the minimum attainable p-values table.
     *
     * Called in
     * execute_constructor_fais()
     * and in psi_destructor().
     */
    void psi_constructor();
    ///@}

    /**
     * @copydoc super::istestable()
     *
     * Check minimum attainable p-value #psi for frequency #freq_par of a given
     * interval `tau`.
     *
     * @param tau Interval from #testable_queue
     * @return `true` if p-value for `freq` is under the current threshold
     */
    inline bool istestable_int(longint tau) override {
        return istestable(psi[freq_par[tau]]);
    }
    /**
     * @copydoc super::isprunable()
     *
     * Directly check interval's frequency #freq_par against the current region
     * threshold #su2.
     *
     * @param tau Interval from #testable_queue
     * @return `true` if interval can be pruned
     */
    inline bool isprunable_int(longint tau) override {
        return (freq_par[tau] > su2);
    };

    inline double compute_interval_pval(longint a, longint tau) override {
        return compute_pval(a,freq_par[tau]);
    };

    inline double compute_interval_score(longint a, longint tau) override {
        return compute_score(a,freq_par[tau]);
    };

    inline double compute_interval_odds_ratio(longint a, longint tau) override {
        return compute_odds_ratio(a,freq_par[tau]);
    };

    /**
     * Evaluate Fisher's test on a table with margins `x`, #n and #N, and cell
     * count `a`.
     *
     * The p-value is defined as a two-tailed p-value which adds up the
     * probabilities of all tables less or equally likely to occur than the
     * one we observed
     *
     * @param a Cell count intended as a number of observation in positive class
     *          in the layer above
     * @param x Table margin intended as frequency of an interval
     * @return p-value
     */
    virtual double compute_pval(longint a, longint x) = 0;

    virtual double compute_score(longint a, longint x) = 0;

    virtual double compute_odds_ratio(longint a, longint x) = 0;

    virtual double score_to_pval(double score) = 0;

    /**
     * @copydoc super::decrease_threshold()
     *
     * Main operations that need to be performed are:
     * 1) Figure out whether we have to shrink "the W" on the left side or the right side, that is, if the current region
     *    is Sigma_{k} = [sl1,sl2] U [N-sl2,N-sl1], we need to figure out if Sigma_{k+1} is of the form
     *    Sigma_{k+1} = [sl1+1,sl2] U [N-sl2,N-sl1-1] (shrink left side) or Sigma_{k+1} = [sl1,sl2-1] U [N-sl2+1,N-sl1-1]
     *    (shrink right side). This is done with help of a binary flag that remembers which of the two types of region
     *    change happened the last time the threshold was decreased.
     * 2) Update variables sl1, sl2 and delta accordingly
     * 3) Recompute the number of testable items by removing those being excluded from the testable region due to the
     *    threshold change
     */
    virtual void decrease_threshold() override;

    void process_first_layer_pvalues() override;
    void process_intervals_pvalues() override;

    void process_first_layer_threshold() override;
    void process_intervals_threshold() override;

public:
    SignificantIntervalSearchFais();
    virtual ~SignificantIntervalSearchFais();

    inline IntervalSetWithOddsRatio const& getPValsTestableInts() const override {
        return pValsTestableInts;
    }
    inline IntervalSetWithOddsRatio const& getPValsSigInts() const override {
        return pValsSigInts;
    }

    inline virtual SummaryFais const& getSummary() const override { return summary; }
};

} /* namespace SignificantPattern */

#endif /* LIBSIGINTERVALSEARCH_SIGNIFICANTINTERVALSEARCHFAIS_H_ */
