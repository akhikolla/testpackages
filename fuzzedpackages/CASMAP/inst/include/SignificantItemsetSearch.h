/*
 * SignificantItemsetSearch.h
 *
 *  Created on: 15 May 2017
 *      Author: mikolajr
 */

#ifndef SIGNIFICANTITEMSETSEARCH_H_
#define SIGNIFICANTITEMSETSEARCH_H_

#include "SignificantFeaturesSearch.h"

// #include <algorithm> // std::max
#include <vector>

#include "Exception.h"
#include "FeatureSet.h"
#include "Summary.h"
#include "TransactionsData.h"

namespace SignificantPattern
{

/**
 * An abstract base class for item set search methods.
 */
class SignificantItemsetSearch : virtual public SignificantFeaturesSearch
{
private:
    /**
     * super class pattern for code independence of changes in inheritance
     */
    typedef SignificantFeaturesSearch super;

    ItemsetSetWithOddsRatio pValsTestableIsets;
    inline void setPValsTestableIsets(ItemsetSetWithOddsRatio iset) {
        pValsTestableIsets = std::move(iset);
    }
    ItemsetSetWithOddsRatio pValsSigIsets;
    inline void setPValsSigIsets(ItemsetSetWithOddsRatio iset) {
        pValsSigIsets = std::move(iset);
    }

    SummaryIset summary;

    void execute_constructor_iset();
    void execute_destructor_iset();

protected:
    inline virtual SummaryIset& getSummaryRef() override { return summary; }
    inline virtual void summary_constructor() override { summary = SummaryIset(); }

    /* -------------------------------- INITIALISATION AND TERMINATION METHODS ----------------------------------------- */
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    /* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */
    virtual void buildItemset(const std::vector<longint> &x_t,
                              const std::vector<longint> &iset,
                              const std::vector<longint> &pexs,
                              /* out */ std::vector<longint> &itemset) = 0;
    bool testAndSaveItemset(double threshold, double score, double odds_ratio, double pval,
                            const std::vector<longint> &x_t, longint a,
                            const std::vector<longint> &iset,
                            const std::vector<longint> &pexs);
    inline void saveSignificantItemset(double pval, double score, double odds_ratio,
                                       const std::vector<longint> &itemset,
                                       longint a) {
        pValsSigIsets.addFeature(itemset, a, score, odds_ratio, pval);
    }
    inline void saveTestableItemset(double pval, double score, double odds_ratio,
                                    const std::vector<longint> &itemset,
                                    longint a) {
        pValsTestableIsets.addFeature(itemset, a, score, odds_ratio, pval);
    }

    // Post-processing of significant features found in the search
    virtual void process_significant_features() override; // called in execute_end()

  public:
    SignificantItemsetSearch();
    virtual ~SignificantItemsetSearch();

    /**
     * Get all in-memory testable itemsets (with their p-values etc).
     *
     * @return Object encapsulating multiple equal-length std::vectors.
     */
    inline ItemsetSetWithOddsRatio const& getPValsTestableIsets() const {
        return pValsTestableIsets;
    }
    /**
     * Get all in-memory significant itemsets (with their p-values etc).
     *
     * @return Object encapsulating multiple equal-length std::vectors.
     */
    inline ItemsetSetWithOddsRatio const& getPValsSigIsets() const {
        return pValsSigIsets;
    }

    inline virtual SummaryIset const& getSummary() const override { return summary; }
};

} /* namespace SignificantPattern */

#endif /* SIGNIFICANTITEMSETSEARCH_H_ */
