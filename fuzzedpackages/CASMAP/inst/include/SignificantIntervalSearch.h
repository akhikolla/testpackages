/**
 * @file    SignificantIntervalSearch.h
 * @author  mabaker
 * @date    Aug 18, 2016
 *
 * Single class header file.
 *
 */

#ifndef SIGNIFICANTINTERVALSEARCH_H_
#define SIGNIFICANTINTERVALSEARCH_H_

#include "SignificantFeaturesSearch.h"
#include "FeatureSet.h"
#include "FilterIntervals.h"
#include "Summary.h"


namespace SignificantPattern
{

/**
 * An abstract base class for interval search methods.
 */
class SignificantIntervalSearch: virtual public SignificantFeaturesSearch
{
private:
    /**
     * Super class pattern for code independence of changes in inheritance
     */
    typedef SignificantFeaturesSearch super;

    void execute_constructor_int();
    void execute_destructor_int();

    /**
     * Summary object encapsulates execution summary including timing, number of
     * intervals tested and found significant, corrected threshold etc.
     */
    SummaryInt summary;

protected:
    inline virtual SummaryInt& getSummaryRef() override { return summary; }
    inline virtual void summary_constructor() override { summary = SummaryInt(); }

    /// @todo
    /// Simplify to selecting a subset of the IntervalSet
    /// instance (+ mv filtering method to this class:
    /// IntervalSet::cpp_filterIntervalsFromMemory())
    /**
     * Filtered (clustered) intervals found.
     *
     * Collected at the end in
     * process_significant_features()
     */
    FilterIntervals filter;
    /// @todo
    /// Unnecessary copy of an already existing data (1-to-1
    /// SignificantPattern::FeatureSet, but only with p-vals).
    /**
     * Significant intervals found.
     *
     * Collected at the end in
     * process_significant_features()
     */
    SignificantIntervals intervals;

    /**
     * Temporary LxN data matrix, for storing the values of nodes in the layer
     * above.
     */
    Genotype genotype_par;

    /**
     * Auxiliary variable to keep track of the current layer
     */
    longint last_tau;

    /// @todo
    /// Use std::queue
    /**
     * @name Queue of testable intervals in the layer below
     *
     * Manually-managed cyclic array.
     */
    ///@{
    longint *testable_queue;
    longint testable_queue_front;
    longint testable_queue_length;
    ///@}

    /**
     * @name INITIALISATION AND TERMINATION METHODS
     */
    ///@{
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    /**
     * @copydoc super::algorithm_init()
     *
     * Additionally, initialise queue of testable intervals and temporary data
     * layer.
     */
    virtual void algorithm_init() override;
    virtual void algorithm_end() override;

    /**
     * Allocate memory and initialise testable intervals queue.
     *
     * Called in algorithm_init().
     */
    virtual void testable_queue_init();
    /**
     * Clear (re-set) testable intervals queue.
     *
     * Called in
     * compute_corrected_significance_threshold(),
     * in find_significant_features()
     * and in testable_queue_init().
     */
    virtual void testable_queue_clear();
    /**
     * Zero-initialise pointer to the testable intervals queue.
     *
     * Called in execute_constructor_int() and in testable_queue_destructor().
     */
    virtual void testable_queue_constructor();
    /**
     * Cleanup the testable queue memory and zero-initialise the pointer.
     *
     * Called in execute_destructor_sisint().
     */
    virtual void testable_queue_destructor();
    ///@}

    /**
     * @name METHODS TO FIND THE CORRECTED SIGNIFICANCE THRESHOLD AND THE INTERVALS
     */
    ///@{

    virtual void compute_corrected_significance_threshold() override;
    virtual void find_significant_features() override;
    virtual void process_significant_features() override;

    /**
     * Process p-values of the first layer of intervals (of length 1) filling
     * initially up the queue of testable intervals.
     *
     * Called in find_significant_features().
     */
    virtual void process_first_layer_pvalues() = 0;
    /**
     * Process p-values of queued intervals, adding to the queue new layers of
     * testable intervals, until the queue is empty due to the pruning.
     *
     * Called last in
     * find_significant_features().
     */
    virtual void process_intervals_pvalues() = 0;

    /**
     * Process current p-value threshold for the first layer of intervals (of
     * length 1) filling initially up the queue of testable intervals.
     *
     * Called in compute_corrected_significance_threshold().
     */
    virtual void process_first_layer_threshold() = 0;
    /**
     * Process current p-value threshold for queued intervals, adding to the
     * queue new layers of testable intervals, until the queue is empty due to
     * the pruning.
     *
     * Called last in
     * compute_corrected_significance_threshold().
     */
    virtual void process_intervals_threshold() = 0;

    /**
     * Compute p-value for a testable interval `tau`, given current cell count
     * `a` for that interval.
     *
     * @param a Cell count intended as a number of observation in positive class
     *          in the layer above
     * @param tau testable interval from #testable_queue
     * @return p-value
     */
    virtual double compute_interval_pval(longint a, longint tau) = 0;

    /**
     * Compute test statistic for a testable interval `tau`, given current cell count
     * `a` for that interval.
     *
     * @param a Cell count intended as a number of observation in positive class
     *          in the layer above
     * @param tau testable interval from #testable_queue
     * @return test statistic
     */
    virtual double compute_interval_score(longint a, longint tau) = 0;

    /**
     * Compute odds ratio for a testable interval `tau`, given current cell count
     * `a` for that interval.
     *
     * @param a Cell count intended as a number of observation in positive class
     *          in the layer above
     * @param tau testable interval from #testable_queue
     * @return The odds ratio for the contingency table(s)
     */
    virtual double compute_interval_odds_ratio(longint a, longint tau) = 0;

    /**
     * Test if an interval is testable under the current significance threshold.
     *
     * @param tau interval from #testable_queue
     * @return `true` if interval can be pruned
     */
    virtual bool istestable_int(longint tau) = 0;
    /**
     * Test if an interval can be pruned from a search.
     *
     * @param tau interval from #testable_queue
     * @return `true` if interval can be pruned
     */
    virtual bool isprunable_int(longint tau) = 0;

    /**
     * Test interval and save it in memory if significant or if
     * #outputForTestableInts.
     *
     * @param threshold Significance threshold
     * @param pval P-value
     * @param score Test statistic
     * @param odds_ratio Odds ratio for the contingency table(s)
     * @param tau Starting position of the interval
     * @param l Length of the interval
     * @param a Cell count
     * @return
     */
    virtual bool testAndSaveInterval(double threshold, double score, double odds_ratio, double pval, longint tau, longint l, longint a);
    /**
     * Save in memory a significant interval found during search.
     *
     * @param pval P-value
     * @param score Test statistic
     * @param odds_ratio Odds ratio for the contingency table(s)
     * @param tau Starting position of the interval
     * @param l Length of the interval
     * @param a Cell count
     */
    virtual void saveSignificantInterval(double pval, double score, double odds_ratio, longint tau, longint l, longint a) = 0;
    /**
     * Save in memory an interval that was tested during search.
     *
     * @param pval P-value
     * @param score Test statistic
     * @param odds_ratio Odds ratio for the contingency table(s)
     * @param tau Starting position of the interval
     * @param l Length of the interval
     * @param a Cell count
     */
    virtual void saveTestableInterval(double pval, double score, double odds_ratio, longint tau, longint l, longint a) = 0;
    ///@}

public:
    SignificantIntervalSearch();
    virtual ~SignificantIntervalSearch();

    /// @todo
    /// Deprecate in favor of a IntervalSet object
    /**
     * Filtered intervals getter.
     *
     * @copydoc #filter
     *
     * @return Reference to #filter
     */
    virtual inline FilterIntervals& getFilteredIntervals() { return filter; }
    /// @todo
    /// Deprecate in favor of getPValsSigInts()
    /**
     * Significant intervals getter.
     *
     * @copydoc #intervals
     *
     * @return Reference to #intervals
     */
    virtual inline SignificantIntervals& getSignificantIntervals() { return intervals; }
    /**
     * Get all in-memory testable intervals (with their p-values etc).
     *
     * @return Object encapsulating multiple equal-length vectors.
     */
    virtual IntervalSetWithOddsRatio const& getPValsTestableInts() const = 0;
    /**
     * Get all in-memory significant intervals (with their p-values etc).
     *
     * @return Object encapsulating multiple equal-length vectors.
     */
    virtual IntervalSetWithOddsRatio const& getPValsSigInts() const = 0;

    inline virtual SummaryInt const& getSummary() const override { return summary; }
};




} /* namespace SignificantPattern */

#endif /* SIGNIFICANTINTERVALSEARCH_H_ */
