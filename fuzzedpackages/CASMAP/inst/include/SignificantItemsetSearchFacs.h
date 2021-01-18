/*
 * SignificantItemsetSearchFacs.h
 *
 *  Created on: 8 May 2017
 *      Author: mikolajr
 */

#ifndef SIGNIFICANTITEMSETSEARCHFACS_H_
#define SIGNIFICANTITEMSETSEARCHFACS_H_

#include "SignificantItemsetSearch.h"
#include "SignificantFeaturesSearchTaroneCmh.h"
#include "chi2.h"

#include <vector>

#include "Exception.h"
#include "Summary.h"
#include "TransactionsData.h"

namespace SignificantPattern
{

class SignificantItemsetSearchFacs : public SignificantItemsetSearch,
                                     public SignificantFeaturesSearchTaroneCmh
{
  private:
    /**
     * super class pattern for code independence of changes in inheritance
     */
    typedef SignificantFeaturesSearchTaroneCmh super_cov;
    typedef SignificantItemsetSearch super_feat;

    void execute_constructor_facs();
    void execute_destructor_facs();

    // Database object in which the final transaction database is stored. This
    // is a vertical representation, with perfect extensions of the root removed
    TransactionsData db;
    // Mapping of internal item labeling (consecutive ints) to original item
    // label
    std::vector<longint> item_label_map;
    // Vector of integers where the internal item labels of all items which are
    // perfect extensions of the root are stored
    std::vector<longint> pexs;

    /**
     * List holding the minimum attainable p-values
     */
    std::vector<double> minpvals;
    /**
     * List holding the corresponding values of the minimum attainable p-value
     * lower envelope
     */
    std::vector<double> lower_envelope_minpvals;

    /**
     * min(ni,Ni-ni) for each of the K tables
     *
     * @param k
     * @return
     */
    inline longint hypercorner_bnd_k(unsigned short k) override {
        return (nt[k] <= Nt_nt[k]) ? nt[k] : Nt_nt[k];
    }
    /**
     * If for any of the k tables its margin x_k is larger than the minimum of
     * n_k and N_k-n_k, then we cannot prune the itemset (we are not in the
     * "bottom-left" hypercorner).
     *
     * @param x_k k-th table margin
     * @param k table number
     * @return `bool` stating if itemset is not prunable due to k-th table
     */
    inline bool notprunable_k(longint x_k, unsigned short k) override {
        return (x_k > hypercorner_bnd[k]);
    }
    inline longint dim_margin_k(longint x_k, unsigned short k) override {
        return (x_k - 0);
    }


    /**
     * @name Global variables of the depth-first recursive search.
     */
    ///@{
    /**
     * Number of enumerated itemsets
     */
    longint n_enumerated;
    /**
     * Number of closed enumerated itemsets
     */
    longint n_enumerated_closed;
    /**
     * Output of recursion
     */
    longint r_out;
    ///@}

    SummaryFacs summary;

protected:
    inline virtual SummaryFacs& getSummaryRef() override { return summary; }
    inline virtual void summary_constructor() override { summary = SummaryFacs(); }

    /* -------------------------------- INITIALISATION AND TERMINATION METHODS ----------------------------------------- */
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    virtual void algorithm_init() override;
    virtual void algorithm_end() override;

    /* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */
    virtual void buildItemset(const std::vector<longint> &x_t,
                              const std::vector<longint> &iset,
                              const std::vector<longint> &pexs,
                              /* out */ std::vector<longint> &itemset) override;

    /**
     * @copydoc super::compute_corrected_significance_threshold()
     *
     * **Attention**: currently, the method is greedy and nothing is done here
     * (the threshold is not corrected before the actual search is executed).
     *
     */
    void compute_corrected_significance_threshold() override {};
    void find_significant_features() override;

    /**
     * `if_compute_pvals==true`, p-values for testable itemsets are computed on
     * the fly, using a greedy choice for the significance threshold.
     * `false` in the orig. code if file for significant itemset (intervals) is
     * not given
     */
    static constexpr bool if_compute_pvals = true;
    /**
     * List of tentative decreasing significance thresholds.
     */
    std::vector<double> tentative_sig_ths;

    void decrease_threshold() override;

    /**
     * Depth-first recursive search of itemset lattice.
     *
     * @param db
     * @param minpvals
     * @param lower_envelope_minpvals
     * @param iset
     * @param pexs
     * @param elim
     * @return Output of recursion
     */
    longint depth(TransactionsData &db, std::vector<double> &minpvals,
                  std::vector<double> &lower_envelope_minpvals,
                  std::vector<longint> &iset, std::vector<longint> &pexs,
                  std::vector<std::vector<bool> *> &elim);

    /**
     * @copydoc super::process_significant_features()
     *
     * Currently does nothing, as in this method tested and found significant
     * patterns are not post-processed (e.g. clustered), and can be accessed
     * directly via `getPValsTestableIsets()` and `getPValsSigIsets()` methods
     * from super.
     */
    virtual void process_significant_features() override;

    inline double score_to_pval(double score){
        return Chi2_sf(score, 1);
    };

  public:
    SignificantItemsetSearchFacs();
    virtual ~SignificantItemsetSearchFacs();

    inline virtual SummaryFacs const& getSummary() const override { return summary; }
};




} /* namespace SignificantPattern */

#endif /* SIGNIFICANTITEMSETSEARCHFACS_H_ */
