/*
 * SignificantFeaturesSearchWithCovariates.h
 *
 *  Created on: 23 Mar 2017
 *      Author: mikolajr
 */

#ifndef SIGNIFICANTINTERVALSEARCHWITHCOVARIATES_H_
#define SIGNIFICANTINTERVALSEARCHWITHCOVARIATES_H_

#include "SignificantFeaturesSearch.h"
#include "FeatureSet.h"
#include "types.h"
#include "Phenotype.h"

namespace SignificantPattern {

class SignificantFeaturesSearchWithCovariates: virtual public SignificantFeaturesSearch {
private:
    // super class pattern for code independence of changes in inheritance
    typedef SignificantFeaturesSearch super;

    Phenotype covariates;

    Phenotype readCovariatesFileToBuffer(
            const std::string& covfilename, bool plinkFormat,
            const Phenotype& labels);

    void execute_constructor_cov();
    void execute_destructor_cov();

protected:
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    void readFilesWithCovariates(
            const std::string& xfilename, const std::string& yfilename, bool plinkFormat,
            const std::string& covfilename, bool covPlinkFormat = false, const std::string& encoding = "dominant");

    inline Phenotype const &getCovariates() const {
        return covariates;
    }
    inline virtual void setCovariates(const Phenotype& covariates_buf) {
        covariates = covariates_buf;
    }
    /**
     * Initialise covariates with one value for all observations.
     */
    inline virtual void initCovariates() {
        covariates.initialiseMatrix(getNumObservations());
    }

    /**
     * @copydoc super::algorithm_init()
     *
     * Additionally, if not initialised, first call initCovariates().
     */
    inline virtual void algorithm_init() override {
        if (!getCovariates().isInitialised()) {
            profiler.markStartInitialisation();
            initCovariates();
            profiler.markEndInitialisation();
        }
        super::algorithm_init();
    }

public:
    SignificantFeaturesSearchWithCovariates();
    virtual ~SignificantFeaturesSearchWithCovariates();

    inline unsigned short getNumCovariates() const {
            return getCovariates().getNumClasses();
    }
    inline std::vector<longint> getNumObservationsInClasses() const {
        return getCovariates().getNumObservationsInClasses();
    }
    std::vector<longint> getNumPositiveObservationsInClasses() const;

    /**
     * @copydoc super::readETHFiles()
     *
     * Additionally, after successful read, call initCovariates().
     */
    inline void readETHFiles(const std::string& xfilename, const std::string& yfilename, const std::string& encoding) override {
        super::readETHFiles(xfilename, yfilename, encoding);
        initCovariates();
    };
    /**
     * @copydoc super::readPlinkFiles()
     *
     * Additionally, after successful read, call initCovariates().
     */
    inline void readPlinkFiles(const std::string& basefilename, const std::string& encoding) override {
        super::readPlinkFiles(basefilename, encoding);
        initCovariates();
    };
    void readETHFilesWithCovariates(
            const std::string& xfilename, const std::string& yfilename,
            const std::string& covfilename, bool covPlinkFormat = false, const std::string& encoding = "dominant" );
    void readPlinkFilesWithCovariates(
            const std::string& basefilename,
            const std::string& covfilename, bool covPlinkFormat = true, const std::string& encoding = "dominant");
    void readCovariatesFile(
            const std::string& covfilename, bool plinkFormat = false);
    void writeETHFilesWithCovariates(
            const std::string& xfilename, const std::string& yfilename,
            const std::string& covfilename);
    void writeCovariatesFile(
            const std::string& covfilename);
};

} /* namespace SignificantPattern */

#endif /* SIGNIFICANTINTERVALSEARCHWITHCOVARIATES_H_ */
