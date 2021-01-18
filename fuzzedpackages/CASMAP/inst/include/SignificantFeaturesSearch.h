/**
 * @file    SignificantFeaturesSearch.h
 * @author  mikolajr
 * @date    Jun 27, 2017
 *
 * Single class header file.
 *
 */

#ifndef SIGNIFICANTFEATURESSEARCH_H_
#define SIGNIFICANTFEATURESSEARCH_H_

#include <string>

#include "Genotype.h"
#include "Phenotype.h"
#include "Profiler.h"
#include "Summary.h"
#include "types.h"

namespace SignificantPattern
{

/**
 * An abstract base class for all search methods.
 */
class SignificantFeaturesSearch
{
private:
    /**
     * Flag to avoid multiple algorithm initializations, possible due to this
     * class being a tip in a "deadly diamond of death".
     *
     * Checked and set: `true` in algorithm_init(), `false` in algorithm_end()
     * and in execute_constructor_feat().
     */
    bool algorithm_initialised = false;

    /**
     * Zero-initialisation of variables and objects related to a single search
     * execution.
     *
     * Zero-initialisation means no memory allocation and no values
     * initialisation. Instead, per execution objects are deleted if they have
     * not been zero/NULL before, and are set as such.
     *
     * Called last in both execute_constructor() and execute_destructor_feat().
     */
    void execute_constructor_feat();
    /**
     * Cleanup of variables and objects related to a single search execution.
     *
     * Called last in execute_destructor().
     */
    void execute_destructor_feat();

    /**
     * Actions taken at the beginning of a single search execution including
     * a full cleanup and zero-initialisation, as well as algorithm
     * initialisation.
     *
     * Called at the beginning of
     * execute().
     * Invokes execute_destructor(),
     * followed by
     * execute_constructor() and
     * algorithm_init()
     *
     * @param alpha Significance threshold
     * @param L_max Sequence prefix length to take into account (0 = whole seq)
     */
    void execute_init(double alpha, longint L_max);
    /**
     * Clean-up and summarize output specific to a single search execution.
     * Actions taken at the end of a single search execution call including
     * algorithm-specific cleanup, significant features post-processing and
     * summary setup.
     *
     * Called at the end of execute(); invokes first algorithm_end(), then
     * process_significant_features().
     */
    void execute_end();

protected:
    /**
     * Writable summary getter.
     *
     * @return Reference to a summary object.
     */
    virtual Summary& getSummaryRef() = 0;
    /**
     * Zero-initialize summary object.
     *
     * Called last in execute_constructor_feat().
     *
     * @return Reference to a summary object.
     */
    virtual void summary_constructor() = 0;

    Phenotype phenotype;
    Genotype genotype;

    /**
     * Data getter.
     *
     * @return Data samples encapsulating object (a binary genotype)
     */
    inline Genotype const &getGenotype() const {
        return genotype;
    }
    /**
     * Labels getter.
     *
     * @return Sample labels encapsulating object (a binary phenotype)
     */
    inline Phenotype const &getPhenotype() const {
        return phenotype;
    }

    /**
     * Flag indicating whether to collect and output all testable intervals.
     */
    bool outputForTestableInts;
    /**
     * Profiler object times program execution or calculates memory usage on
     * demand.
     */
    Profiler profiler;

    /// @todo
    /// Fix an unlikely but possible BUG in WY methods, where
    /// #N is implicitly casted to
    /// `int` when passed to
    /// SignificantPattern::SignificantIntervalSearchWy#randperm (and
    /// SignificantPattern::SignificantIntervalSearchWy#randint); unlikely,
    /// because `int` would give maximum 1GB
    /// line in the ETH data file, so it won't happen for ETH data files smaller
    /// than L GB (L=#lines; NB: this limit can be raised 4 times with `unsigned
    /// int`); see: http://www.cplusplus.com/reference/climits/ .
    /**
     * Number of observations.
     */
    longint N;
    /**
     * Midpoint of interval [0,N], floor(N/2).
     */
    longint N_over_2;
    /**
     * Number of observations in the positive class (with label 1).
     */
    longint n;
    /**
     * Sequence length.
     */
    longint L;
    /**
     * Maximal sequence length.
     */
    longint L_max;

    /**
     * Current interval length.
     */
    longint l;
    /**
     * Number of testable intervals.
     */
    longint m;
    /**
     * Target FWER
     */
    double alpha;

    /**
     * Final corrected significance threshold
     */
    double delta_opt;

    /**
     * Current p-value (significance threshold)
     */
    double delta;

    /**
     * Logarithm of 1/binom(N,n). This terms appears for every evaluation of the
     * hypergeometric PDF, so it makes sense to precompute it
     */
    double log_inv_binom_N_n;

    /**
     * @name Profiling variables
     */
    ///@{
    longint n_featuresets_processed;
    longint n_pvalues_computed;
    longint n_significant_featuresets;
    ///@}

    ///@todo
    /// Use [`enum`](http://en.cppreference.com/w/cpp/language/enum) instead of
    /// a `bool` flag for format specification (to be able to support more than
    /// two formats).
    /**
     * Read data (genome) and labels (phenotype) files in the ETH or PLINK
     * format.
     *
     * Keep previous files if any of the files cannot be read or is invalid.
     *
     * @param xfilename Data file path
     * @param yfilename Labels file path
     * @param plinkFormat Flag indicating PLINK (or ETH) file format
     */
    virtual void readFiles(const std::string& xfilename, const std::string& yfilename, bool plinkFormat, const std::string& encoding);
    /**
     * Add file name extension of the PLINK data file.
     *
     * @param basefilename Path and name of the PLINK file w/o extension
     * @return Path and file name with extension
     */
    std::string getPlinkDataFilename(const std::string& basefilename) const;
    /**
     * Add file name extension of the PLINK labels file.
     *
     * @param basefilename Path and name of the PLINK file w/o extension
     * @return Path and file name with extension
     */
    std::string getPlinkLabelsFilename(const std::string& basefilename) const;
    /**
     * Read labels file to a temporary buffer.
     *
     * @param yfilename Labels file path
     * @param plinkFormat Flag indicating PLINK (or ETH) file format
     * @return Temporary labels file object (phenotype)
     */
    Phenotype readLabelsFileToBuffer(const std::string& yfilename, bool plinkFormat);
    /**
     * Read data file.
     *
     * @param xfilename Data file path
     * @param plinkFormat Flag indicating PLINK (or ETH) file format
     * @param phenotype_buf Labels file object corresponding to the data file,
     *                      to check for size consistency
     */
    void readDataFile(const std::string& xfilename, bool plinkFormat, Phenotype& phenotype_buf, const std::string& encoding);

    /**
     * @name INITIALISATION AND TERMINATION METHODS
     */
    ///@{
    /// @todo
    /// Consider re-factoring single-execute vars & methods into own
    /// classes; For instance, using templates, declaration of the search obj
    /// could look like this:
    /// ```
    ///     SignificantPattern::SignificantIntervalSearch<MethodInterval{Chi,Exact},ErrorFwer{Taron,Wy}> search;
    /// ```
    /**
     * Allocate and zero-initialise execute-specific memory.
     *
     * Execute-specific variables idempotent pseudo-constructor; called second
     * in execute_init(), and in the default constructor.
     */
    virtual void execute_constructor();
    /**
     * Free execute-specific memory.
     *
     * Execute-specific variables idempotent pseudo-destructor; called first in
     * execute_init(), on error during execution, and in the default destructor.
     */
    virtual void execute_destructor();

    // Note: freq_* declared here just for the common docs
    /**
     * Allocate memory and initialise frequency counters for intervals.
     *
     * Called in algorithm_init().
     */
    virtual void freq_init() = 0;
    /**
     * Clear (re-set) frequency counters.
     *
     * Called in
     * find_significant_features()
     * and in freq_init().
     */
    virtual void freq_clear() = 0;
    /**
     * Zero-initialise pointers to the frequency counters.
     *
     * Called in
     * execute_constructor_feat()
     * and in freq_destructor().
     */
    virtual void freq_constructor() = 0;
    /**
     * Cleanup the frequency counters memory and zero-initialise the pointers.
     *
     * Called in algorithm_end() and in
     * execute_destructor_feat().
     */
    virtual void freq_destructor() = 0;

    /**
     * Hooks for pre-execute algorithm-specific initialisations.
     *
     * Called last in execute_init().
     */
    inline virtual void algorithm_init() {
        if (!algorithm_initialised) {
            summary_constructor();
        }
        algorithm_initialised = true;
    };
    /**
     * Hooks for post-execute algorithm-specific cleanups.
     *
     * Called first in execute_end().
     */
    inline virtual void algorithm_end() {
        algorithm_initialised = false;
    };
    ///@}

    /**
     * @name METHODS TO FIND THE CORRECTED SIGNIFICANCE THRESHOLD AND THE INTERVALS
     */
    ///@{
    /**
     * Test if a feature set with given minimal p-value is testable under the
     * current significance threshold.
     *
     * @param minpval Minimum p-value of the feature set
     * @return `true` if minpval is under the current p-value threshold
     */
    inline bool istestable(double minpval) { return minpval <= delta; }
    /**
     * Test if a feature set with given minimal p-value can be pruned from
     * a search.
     *
     * @param minpval Minimum p-value of the feature set
     * @return `true` if minpval is over the current p-value threshold
     */
    inline bool isprunable(double minpval) { return !istestable(minpval); }

    /**
     * Wrapper function that encapsulates the functionality required to find
     * the corrected significance threshold.
     *
     * Called first in execute().
     */
    virtual void compute_corrected_significance_threshold() = 0;
    /**
     * Compute a final corrected significance threshold #delta_opt.
     */
    void compute_final_corrected_significance_threshold() { delta_opt = alpha/m; };

    /**
     * Wrapper function that encapsulates the functionality required to find
     * significant intervals.
     *
     * Called second in execute().
     */
    virtual void find_significant_features() = 0;

    /**
     * Decrease the current p-value threshold by one level/step.
     */
    virtual void decrease_threshold() = 0;
    /**
     *  Update the current minimum p-value threshold until FWER condition is
     *  satisfied again.
     */
    inline virtual void update_threshold() {
        while ((m * delta) > alpha) decrease_threshold();
    }

    /**
     * Post-processing of significant features found in the search (setting up
     * output).
     *
     * Called second in execute_end().
     */
    virtual void process_significant_features() = 0;
    ///@}

public:
    SignificantFeaturesSearch();
    virtual ~SignificantFeaturesSearch();

    /**
     * Get number of observations (samples).
     *
     * @return Number of observations
     */
    inline longint getNumObservations() const {
        // Important: from phenotype, not genotype, because it's read first
        return phenotype.getNumObservations();
    }
    /**
     * Get number of observations (samples) in the positive class (labelled 1).
     *
     * @return Number of positive observations
     */
    longint getNumPositiveObservations() const;
    /**
     * Get number of (binary) features considered for each sample.
     *
     * @return Number of features
     */
    inline longint getNumFeatures() const { return genotype.getNumFeatures(); }

    /**
     * Read data and labels files in the ETH file format into the memory.
     *
     * @param xfilename Data file path
     * @param yfilename Labels file path
     */
    virtual void readETHFiles(const std::string& xfilename, const std::string& yfilename, const std::string& encoding);
    /**
     * Write data and labels from memory to files in the ETH file format.
     *
     * @param xfilename Data file path
     * @param yfilename Labels file path
     */
    void writeETHFiles(const std::string& xfilename, const std::string& yfilename);
    /**
     * Read data and labels files in the ETH file format into the memory.
     *
     * @param basefilename Same path and name for both data and labels files in
     *                     PLINK file format, w/o file extension
     */
    virtual void readPlinkFiles(const std::string& basefilename, const std::string& encoding);
    /**
     * Execute search for significant features.
     *
     * @param alpha Significance threshold
     * @param L_max Sequence prefix length to take into account (0 means take
     *              the whole sequence)
     */
    void execute(double alpha, longint L_max);
    /// @todo
    /// Change return value to `Profiler const&` and add protected
    /// `getProfilerRef()`; update calls respectively.
    /**
     * Profiler getter.
     *
     * @copydoc #profiler
     *
     * @return Reference to #profiler.
     */
    inline Profiler& getProfiler() { return profiler; }

   /**
     * Read-only summary getter.
     *
     * @return Constant reference to a summary object.
     */
    virtual Summary const& getSummary() const = 0;
};




} /* namespace SignificantPattern */

#endif /* SIGNIFICANTFEATURESSEARCH_H_ */
