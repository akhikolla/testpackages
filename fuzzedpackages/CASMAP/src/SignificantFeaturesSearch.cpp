/**
 * @file    SignificantFeaturesSearch.cpp
 * @author  mikolajr
 * @date    Jun 28, 2017
 *
 * Single class source file.
 *
 */
#include "SignificantFeaturesSearch.h"

#include <stdlib.h>
#include <iostream>
#include <sstream> //stringstream
#include <math.h>
#include <numeric> // accumulate
#include <iterator> // next

#include "Exception.h"



using namespace std;

namespace SignificantPattern
{

SignificantFeaturesSearch::SignificantFeaturesSearch()
    : phenotype(),
      genotype(),
      profiler()
{
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearch()\n");
    #endif
    execute_constructor_feat();
}
SignificantFeaturesSearch::~SignificantFeaturesSearch() {
    #ifdef DEBUG
    fprintf(stderr, "~SignificantFeaturesSearch()\n");
    #endif
    execute_destructor_feat();
}

void SignificantFeaturesSearch::execute_constructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearch::execute_constructor()\n");
    #endif
    execute_constructor_feat();
}
void SignificantFeaturesSearch::execute_constructor_feat(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearch::execute_constructor_feat()\n");
    #endif
    outputForTestableInts = false;

    N=0; N_over_2=0; n=0; L=0; L_max=0;

    l=0; m=0; alpha=0;

    delta=0; delta_opt=0;

    log_inv_binom_N_n=0;

    n_featuresets_processed=0; n_pvalues_computed=0; n_significant_featuresets=0;

}
void SignificantFeaturesSearch::execute_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearch::execute_destructor()\n");
    #endif
    execute_destructor_feat();
}
void SignificantFeaturesSearch::execute_destructor_feat(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearch::execute_destructor_feat()\n");
    #endif
    algorithm_initialised = false;
}



void SignificantFeaturesSearch::execute_init(double alpha, longint L_max){
    if (!(phenotype.isInitialised() && genotype.isInitialised())) {
        throw Exception("Genotype and phenotype files have to be read first.");
    }

    // Cleanup after previous executes and re-init execute vars
    execute_destructor();
    execute_constructor();

    // Compute total number of observations and number of observations in
    // minority class
    N = getNumObservations();
    n = getNumPositiveObservations();

    // Store core constants
    N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)

    // Compute number of features
    L = getNumFeatures();

    // Set target FWER and nr of features cut-off
    this->alpha = alpha;
    this->L_max = L_max;

    // Algorithm-specific initialisation
    algorithm_init();
}

void SignificantFeaturesSearch::execute_end(){
    algorithm_end();

    // set up output features
    process_significant_features();

    // set up summary output object
    Summary& summary = getSummaryRef();
    summary.setN(N);
    summary.setn(n);
    summary.setL(L);
    summary.setm(m);
    summary.setL_max(L_max);
    summary.setNumSignificantFeatures(n_significant_featuresets);
    summary.setNumFeaturesProcessed(n_featuresets_processed);
    summary.setDelta(delta);
    summary.setDelta_opt(delta_opt);
    summary.setAlpha(alpha);
}

void SignificantFeaturesSearch::execute(double alpha, longint L_max)
{

    try {
        // Get time when program started to run
        profiler.reset();
        profiler.markStartExecution();

        // INITIALISATION
        profiler.markStartInitialisation();
        execute_init(alpha, L_max);
        profiler.markEndInitialisation();

        // MAIN FUNCTIONALITY
        profiler.markStartSignificanceThresholdCompute();
        compute_corrected_significance_threshold();
        profiler.markEndSignificanceThresholdCompute();
        profiler.markStartSignificantIntervalsCompute();
        find_significant_features();
        profiler.markEndSignificantIntervalsCompute();

        // POST-PROCESSING AND CLEANUP
        profiler.markStartPostprocessingAndCleanup();
        execute_end();
        profiler.markEndPostprocessingAndCleanup();

        // Get time when program finished
        profiler.markEndExecution();

        // Record peak memory usage
        profiler.markPeakMemory();
    } catch (const std::exception& e) {
        #ifdef DEBUG
        fprintf(stderr, "SignificantFeaturesSearch::execute() exception:\n\t%s\n", e.what());
        #endif
        execute_destructor();
        throw;
    }
}



longint SignificantFeaturesSearch::getNumPositiveObservations() const {
    std::vector<longint> v = phenotype.getNumObservationsInClasses();
    return std::accumulate(std::next(v.begin()), v.end(), 0);
}



Phenotype SignificantFeaturesSearch::readLabelsFileToBuffer(
        const std::string& yfilename, bool plinkFormat)
{
    Phenotype phenotype_buf = Phenotype();
    profiler.markStartFileIO();
    if (plinkFormat) phenotype_buf.readPlinkFamFile(yfilename);
    else phenotype_buf.readETHFile(yfilename);
    profiler.markEndFileIO();
    // TODO init as Phenotype(2) and throw exeception on file check, otherwise,
    //      if exepcted nr of classes not give, read nr on check
    if (phenotype_buf.getNumClasses() > 2) {
        std::stringstream errmsgstream;
        errmsgstream << "Too many phenotypes (" << phenotype_buf.getNumClasses()
                << ") in the labels file";
        throw Exception(errmsgstream.str());
    }
    return phenotype_buf;
}
void SignificantFeaturesSearch::readDataFile(
        const std::string& xfilename, bool plinkFormat, Phenotype& phenotype_buf, const std::string& encoding)
{
    profiler.markStartFileIO();
    if (plinkFormat) genotype.readPlinkRawFile(xfilename, phenotype_buf);
    else genotype.readETHFile(xfilename, phenotype_buf.getNumObservations(), encoding);
    profiler.markEndFileIO();
}

void SignificantFeaturesSearch::readFiles(const std::string& xfilename, const std::string& yfilename, bool plinkFormat, const std::string& encoding){
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearch::readFiles(): BEGIN, current data size: LxN=%ldx%ld\n", genotype.getNumFeatures(), genotype.getNumObservations());
    #endif
    // read{ETH,PlinkFam,PlinkRaw}File functions run check before allocating
    // memory, hence, to keep both files we need local buffer only for the
    // labels file (Y), which has to be read first. Label files are small, so
    // mem overhead is no issue.
    Phenotype phenotype_buf = readLabelsFileToBuffer(yfilename, plinkFormat);
    readDataFile(xfilename, plinkFormat, phenotype_buf, encoding);
    // no errors, genotype is already set, set phenotype from the buffer
    phenotype = phenotype_buf;
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearch::readFiles(): END, current data size: LxN=%ldx%ld\n", genotype.getNumFeatures(), genotype.getNumObservations());
    #endif
}

void SignificantFeaturesSearch::readETHFiles(const std::string& xfilename, const std::string& yfilename, const std::string& encoding)
{
    readFiles(xfilename, yfilename, false, encoding);
}

std::string SignificantFeaturesSearch::getPlinkDataFilename(const std::string& basefilename) const {
    return basefilename+".raw";
}
std::string SignificantFeaturesSearch::getPlinkLabelsFilename(const std::string& basefilename) const {
    return basefilename+".fam";
}
void SignificantFeaturesSearch::readPlinkFiles(const std::string& basefilename, const std::string& encoding)
{
    readFiles(getPlinkDataFilename(basefilename), getPlinkLabelsFilename(basefilename), true, encoding);
}

void SignificantFeaturesSearch::writeETHFiles(const std::string& xfilename, const std::string& yfilename) {
    profiler.markStartFileIO();
    phenotype.writeETHFile(yfilename);
    genotype.writeETHFile(xfilename);
    profiler.markEndFileIO();
}

} /* namespace SignificantPattern */

