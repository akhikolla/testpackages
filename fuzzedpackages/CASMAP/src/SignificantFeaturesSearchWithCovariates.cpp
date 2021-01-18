/*
 * SignificantFeaturesSearchWithCovariates.cpp
 *
 *  Created on: 23 Mar 2017
 *      Author: mikolajr
 */

#include "SignificantFeaturesSearchWithCovariates.h"

#include <numeric> // accumulate
#include <cassert>

#include "Exception.h"



namespace SignificantPattern {

SignificantFeaturesSearchWithCovariates::SignificantFeaturesSearchWithCovariates()
    : SignificantFeaturesSearch() // explicitly construct a virtual base class
{
    setCovariates(Phenotype()); // initialised with a single covariate value (= no covariates)
}
SignificantFeaturesSearchWithCovariates::~SignificantFeaturesSearchWithCovariates()
{
}

void SignificantFeaturesSearchWithCovariates::execute_constructor() {
    super::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchWithCovariates::execute_constructor()\n");
    #endif
    execute_constructor_cov();
}
void SignificantFeaturesSearchWithCovariates::execute_destructor() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchWithCovariates::execute_destructor()\n");
    #endif
    execute_destructor_cov();
    super::execute_destructor();
}
void SignificantFeaturesSearchWithCovariates::execute_constructor_cov() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchWithCovariates::execute_constructor_cov()\n");
    #endif
}
void SignificantFeaturesSearchWithCovariates::execute_destructor_cov() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchWithCovariates::execute_destructor_cov()\n");
    #endif
}



std::vector<longint> SignificantFeaturesSearchWithCovariates::getNumPositiveObservationsInClasses() const
{
    // OPT cache results and use the change mark to re-calculate on-demand
    unsigned char *labelsPtr = phenotype.getVectorPtr();
    unsigned char *classesPtr = getCovariates().getVectorPtr();
    std::vector<longint> positiveLablelsInClassCounters(getCovariates().getNumClasses(), 0);
    for (longint i = 0; i < getNumObservations(); ++i) {
        if(labelsPtr[i]) positiveLablelsInClassCounters[classesPtr[i]]++;
    }
    return positiveLablelsInClassCounters;
}



Phenotype SignificantFeaturesSearchWithCovariates::readCovariatesFileToBuffer(
        const std::string& covfilename, bool plinkFormat, const Phenotype& labels)
{
    longint N_expected = labels.getNumObservations();
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchWithCovariates::readCovariatesFileToBuffer('%s'), N=%ld\n", covfilename.c_str(), N_expected);
    #endif
    Phenotype covariates_new;
    profiler.markStartFileIO();
    if (plinkFormat)
        covariates_new.readPlinkCovFile(covfilename, labels);
    else
        covariates_new.readETHFile(covfilename, N_expected);
    profiler.markEndFileIO();
    //size sanity checks
    assert(N_expected==covariates_new.getNumObservations());
    std::vector<longint> v = covariates_new.getNumObservationsInClasses();
    assert(N_expected==std::accumulate(v.begin(), v.end(), 0));
    return covariates_new;
}
void SignificantFeaturesSearchWithCovariates::readCovariatesFile(
        const std::string& covfilename, bool plinkFormat)
{
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchWithCovariates::readCovariatesFile('%s'), N=%ld\n", covfilename.c_str(), getNumObservations());
    #endif
    const Phenotype& labels = getPhenotype();
    if (labels.getNumObservations() <= 0)
        throw Exception("Unknow number of observations - read labels first");
    setCovariates(readCovariatesFileToBuffer(covfilename, plinkFormat, labels));
}



void SignificantFeaturesSearchWithCovariates::readFilesWithCovariates(
        const std::string& xfilename, const std::string& yfilename, bool plinkFormat,
        const std::string& covfilename, bool covPlinkFormat, const std::string& encoding)
{
    Phenotype phenotype_buf = readLabelsFileToBuffer(yfilename, plinkFormat);
    Phenotype covariates_buf = readCovariatesFileToBuffer(covfilename, covPlinkFormat, phenotype_buf);
    readDataFile(xfilename, plinkFormat, phenotype_buf, encoding);
    // everything went smooth, memoize the buffers
    phenotype = phenotype_buf;
    setCovariates(covariates_buf);
}
void SignificantFeaturesSearchWithCovariates::readETHFilesWithCovariates(
        const std::string& xfilename, const std::string& yfilename,
        const std::string& covfilename, bool covPlinkFormat, const std::string& encoding)
{
    readFilesWithCovariates(xfilename, yfilename, false, covfilename, covPlinkFormat, encoding);
}
void SignificantFeaturesSearchWithCovariates::readPlinkFilesWithCovariates(
        const std::string& basefilename,
        const std::string& covfilename, bool covPlinkFormat, const std::string& encoding )
{
    readFilesWithCovariates(getPlinkDataFilename(basefilename), getPlinkLabelsFilename(basefilename), true,
            covfilename, covPlinkFormat, encoding);
}

void SignificantFeaturesSearchWithCovariates::writeETHFilesWithCovariates(
        const std::string& xfilename, const std::string& yfilename,
        const std::string& covfilename)
{
    super::writeETHFiles(xfilename, yfilename);
    writeCovariatesFile(covfilename);
}

void SignificantFeaturesSearchWithCovariates::writeCovariatesFile(
        const std::string& covfilename)
{
    profiler.markStartFileIO();
    getCovariates().writeETHFile(covfilename);
    profiler.markEndFileIO();
}

} /* namespace SignificantPattern */
