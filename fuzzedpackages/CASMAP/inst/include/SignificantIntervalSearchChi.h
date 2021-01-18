/*
 * SignificantIntervalSearchExact.h
 *
 *  Created on: Aug 18, 2016
 *      Author: mabaker
 */

#ifndef SIGNIFICANTINTERVALSEARCHCHI_H_
#define SIGNIFICANTINTERVALSEARCHCHI_H_

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

namespace SignificantPattern
{

class SignificantIntervalSearchChi : public SignificantIntervalSearchFais
{
private:
    // super class pattern for code independence of changes in inheritance
    typedef SignificantIntervalSearchFais super;

    double class_ratio, class_ratio_bin;

    /* -------------------------------- INITIALISATION AND TERMINATION METHODS ----------------------------------------- */
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    void execute_constructor_chi();
    void execute_destructor_chi();

    virtual void algorithm_init() override;
    virtual void algorithm_end() override;

    virtual void psi_clear() override;


    /* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */
    /**
     * Chi2 approximation
     *
     * @copydoc super::compute_pval()
     */
    double compute_pval(longint a, longint x) override;

    double compute_score(longint a, longint x) override;

    double compute_odds_ratio(longint a, longint x) override;

    virtual double score_to_pval(double score) override;

public:
    SignificantIntervalSearchChi();
    virtual ~SignificantIntervalSearchChi();

};




} /* namespace SignificantPattern */

#endif /* SIGNIFICANTINTERVALSEARCHCHI_H_ */
