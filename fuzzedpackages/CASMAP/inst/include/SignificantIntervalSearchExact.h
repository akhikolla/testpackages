/*
 * SignificantIntervalSearchExact.h
 *
 *  Created on: Aug 18, 2016
 *      Author: mabaker
 */

#ifndef SIGNIFICANTINTERVALSEARCHEXACT_H_
#define SIGNIFICANTINTERVALSEARCHEXACT_H_

#include <string>
#include <iostream>
#include <fstream>

/* LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */
#include <time.h> //Already included in original LCM source code above
//#include <sys/time.h>
//#include <sys/resource.h>
/* END OF LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */

#include "SignificantIntervalSearchFais.h"
#include "Summary.h"
#include "chi2.h"


namespace SignificantPattern
{

class SignificantIntervalSearchExact : public SignificantIntervalSearchFais
{
private:
    // super class pattern for code independence of changes in inheritance
    typedef SignificantIntervalSearchFais super;

    double *loggamma;

    /* -------------------------------- INITIALISATION AND TERMINATION METHODS ----------------------------------------- */
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    void execute_constructor_exact();
    void execute_destructor_exact();

    virtual void algorithm_init() override;
    // virtual void algorithm_end() override;

    void loggamma_init();
    void loggamma_constructor();
    /**
     * Precompute values of log(x!) storing them in the array #loggamma
     */
    void loggamma_clear();
    void loggamma_destructor();

    virtual void psi_clear() override;


    /* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */
    /**
     * Exact test
     *
     * @copydoc super::compute_pval()
     */
    double compute_pval(longint a, longint x) override;

    double compute_score(longint a, longint x) override;

    double compute_odds_ratio(longint a, longint x) override;

    double score_to_pval(double score) override;

public:
    SignificantIntervalSearchExact();
    virtual ~SignificantIntervalSearchExact();

};




} /* namespace SignificantPattern */

#endif /* SIGNIFICANTINTERVALSEARCHEXACT_H_ */
