/*
 * SignificantIntervalSearchFastCmh.h
 *
 *  Created on: 2 Mar 2017
 *      Author: mikolajr
 */

#ifndef SIGNIFICANTINTERVALSEARCHFASTCMH_H_
#define SIGNIFICANTINTERVALSEARCHFASTCMH_H_

#include "SignificantIntervalSearch.h"
#include "SignificantFeaturesSearchTaroneCmh.h"
#include "chi2.h"

namespace SignificantPattern
{

class SignificantIntervalSearchFastCmh : public SignificantIntervalSearch,
                                         public SignificantFeaturesSearchTaroneCmh
{
private:
    // super class pattern for code independence of changes in inheritance
    typedef SignificantFeaturesSearchTaroneCmh super_cov;
    typedef SignificantIntervalSearch super_feat;

    IntervalSetWithOddsRatio pValsTestableInts;
    IntervalSetWithOddsRatio pValsSigInts;


    void execute_constructor_fastcmh();
    void execute_destructor_fastcmh();

protected:
    //OPT: Record frequencies for each class freq_par_cov[tau] (define and use
    //     IntervalSetWithFreqInClasses)
    inline void saveSignificantInterval(double pval, double score, double odds_ratio, longint tau, longint l, longint a) override {
        pValsSigInts.addFeature(tau, tau+l, a, score, odds_ratio, pval);
    }
    inline void saveTestableInterval(double pval, double score, double odds_ratio, longint tau, longint l, longint a) override {
        pValsTestableInts.addFeature(tau, tau+l, a, score, odds_ratio, pval);
    }

    /* -------------------------------- INITIALISATION AND TERMINATION METHODS ----------------------------------------- */
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    inline virtual void algorithm_init() override {
        super_cov::algorithm_init();
        super_feat::algorithm_init();
    };
    inline virtual void algorithm_end() override {
        super_feat::algorithm_end();
        super_cov::algorithm_end();
    }

    /* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */
    inline double compute_interval_pval(longint a, longint tau) override {
        return compute_pval(a, freq_par_cov[tau]);
    };

    inline double compute_interval_score(longint a, longint tau) override {
        return compute_score(a, freq_par_cov[tau]);
    };

    inline double compute_interval_odds_ratio(std::vector<longint> &at, longint tau){
        return compute_odds_ratio(at.data(), freq_par_cov[tau]);
    };

    inline double compute_interval_odds_ratio(longint a, longint tau) override {
        return -1;  //Dummy function, never should be called. Use version above instead
    };

    inline double score_to_pval(double score){
        return Chi2_sf(score, 1);
    };

    inline bool istestable_int(longint tau) override {
        return istestable_freqcov(freq_par_cov[tau]);
    }
    inline bool isprunable_int(longint tau) override {
        return isprunable_freqcov(freq_par_cov[tau]);
    };

    void process_first_layer_pvalues() override;
    void process_intervals_pvalues() override;

    void process_first_layer_threshold() override;
    void process_intervals_threshold() override;

public:
    using SignificantFeaturesSearchTaroneCmh::score_to_pval;  //TODO: Placement
    SignificantIntervalSearchFastCmh();
    virtual ~SignificantIntervalSearchFastCmh();

    inline IntervalSetWithOddsRatio const& getPValsTestableInts() const override {
        return pValsTestableInts;
    }
    inline IntervalSetWithOddsRatio const& getPValsSigInts() const override {
        return pValsSigInts;
    }
};




} /* namespace SignificantPattern */

#endif /* SIGNIFICANTINTERVALSEARCHFASTCMH_H_ */
