#include "types.h"

#include "FeatureSet.h"
#include "Summary.h"

#include "SignificantFeaturesSearch.h"
#include "SignificantItemsetSearch.h"
#include "SignificantItemsetSearchFacs.h"
#include "SignificantIntervalSearch.h"
#include "SignificantIntervalSearchFais.h"
#include "SignificantIntervalSearchExact.h"
#include "SignificantIntervalSearchChi.h"
#include "SignificantFeaturesSearchWithCovariates.h"
#include "SignificantIntervalSearchFastCmh.h"

#include <Rcpp.h>
using namespace Rcpp;

//TODO: use Rcpp modules to create S4/RefClasses directly; cf.
//      http://dirk.eddelbuettel.com/code/rcpp/Rcpp-modules.pdf

SEXP _get_intervals(std::vector<SignificantPattern::Interval>& intervals) {
    size_t size = intervals.size();

    IntegerVector starts(size);
    IntegerVector ends(size);
    DoubleVector scores(size);
    DoubleVector odds_ratios(size);
    DoubleVector pvalues(size);

    for (size_t i = 0; i < size; i++) {
        starts[i] = intervals[i].getStart();
        ends[i] = intervals[i].getEnd();
        scores[i] = intervals[i].getScore();
        odds_ratios[i] = intervals[i].getOddsRatio();
        pvalues[i] = intervals[i].getPvalue();
    }

    return DataFrame::create(Named("start") = starts, Named("end") = ends,
                             Named("score") = scores, Named("odds_ratio") = odds_ratios, Named("pvalue") = pvalues);
}

SEXP _get_itemsets(const SignificantPattern::ItemsetSetWithOddsRatio itemsets) {

    size_t size = itemsets.getLength();
    const std::vector<double> scores  = itemsets.getScoreVector();
    const std::vector<double> odds_ratios  = itemsets.getOddsRatioVector();
    const std::vector<double> pvals  = itemsets.getPValueVector();
    const std::vector< std::vector<longint> > itemsets_vec = itemsets.getItemsetsVector();
    IntegerVector rownamesVec(size);
    // NumericVector pvalsVec(size);
    // List itemsetsList(size);
    // size_t iset_size;
    for (size_t i = 0; i < size; ++i) {
        // pvalsVec[i] = pvals[i];
        // std::vector<longint> itemset_vec = itemsets_vec[i];
        // iset_size = itemset_vec.size();
        // IntegerVector itemsetVec(iset_size);
        // // fprintf(stdout, "\tiset_size=%zu\n", iset_size);
        // for(int j = 0; j < iset_size; ++j)
        //     itemsetVec[j] = itemset_vec[j];
        // itemsetsList[i] = itemsetVec;
        rownamesVec[i] = i;
    }

    // Hack: create data frame from a list by a class atrr change, because
    // as.data.frame (or data.frame) R functions try to convert items in list
    // elements to separate columns
    // List ret = List::create(Named("itemsets") = itemsetsList,
    //                         Named("pvalue") =  pvalsVec);
    List ret =  List::create(Named("itemsets") = itemsets_vec,
                             Named("score") =  scores,
                             Named("odds_ratio") =  odds_ratios,
                             Named("pvalue") =  pvals);
    ret.attr("class") = "data.frame";
    ret.attr("row.names") = rownamesVec;
    return ret;
}


//' Internal function
//' @keywords internal
// [[Rcpp::export]]
SEXP lib_new_search_e() {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchExact> ptr( new SignificantPattern::SignificantIntervalSearchExact(), true );
    return ptr;
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
SEXP lib_new_search_chi() {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchChi> ptr( new SignificantPattern::SignificantIntervalSearchChi(), true );
    return ptr;
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
SEXP lib_new_search_fastcmh() {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchFastCmh> ptr( new SignificantPattern::SignificantIntervalSearchFastCmh(), true );
    return ptr;
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
SEXP lib_new_search_facs() {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearchFacs> ptr( new SignificantPattern::SignificantItemsetSearchFacs(), true );
    return ptr;
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_delete_search_e(SEXP inst) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchExact> ptr(inst);
    ptr.release();
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_delete_search_chi(SEXP inst) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchChi> ptr(inst);
    ptr.release();
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_delete_search_fastcmh(SEXP inst) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchFastCmh> ptr(inst);
    ptr.release();
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_delete_search_facs(SEXP inst) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearchFacs> ptr(inst);
    ptr.release();
}

// Because of multiple inheritance, a first non-virtual class that actually
// constructs the virtual public class SignificantFeaturesSearch has to be used
// for calling SignificantFeaturesSearch non-virtual methods; hence, duplicates.
// Moreover, for overriden methods with covariant (return) types, we have to
// declare a separate function for each different subtype, as for summaries, or
// results.

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_read_eth_files(SEXP inst, std::string x_filename, std::string y_filename, std::string encoding) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
    ptr->readETHFiles(x_filename, y_filename, encoding);
}

//TODO: read_eth_files_{int,iset}

// // [[Rcpp::export]]
// void read_eth_files_int(SEXP inst, std::string x_filename, std::string y_filename) {
//     Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
//     ptr->readETHFiles(x_filename, y_filename);
// }

// // [[Rcpp::export]]
// void read_eth_files_iset(SEXP inst, std::string x_filename, std::string y_filename) {
//     Rcpp::XPtr<SignificantPattern::SignificantItemsetSearch> ptr(inst);
//     ptr->readETHFiles(x_filename, y_filename);
// }

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_read_eth_files_with_cov_fastcmh(SEXP inst, std::string x_filename, std::string y_filename, std::string covfilename, std::string encoding) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchFastCmh> ptr(inst);
    ptr->readETHFilesWithCovariates(x_filename, y_filename, covfilename, false, encoding);
}


//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_read_eth_files_with_cov_facs(SEXP inst, std::string x_filename, std::string y_filename, std::string covfilename, std::string encoding) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearchFacs> ptr(inst);
    ptr->readETHFilesWithCovariates(x_filename, y_filename, covfilename, false, encoding);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_read_plink_files(SEXP inst, std::string base_filename, std::string encoding) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
    ptr->readPlinkFiles(base_filename, encoding);
}

//TODO: read_plink_files_{int,iset}

// // [[Rcpp::export]]
// void read_plink_files_int(SEXP inst, std::string base_filename) {
//     Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
//     ptr->readPlinkFiles(base_filename);
// }

// // [[Rcpp::export]]
// void read_plink_files_iset(SEXP inst, std::string base_filename) {
//     Rcpp::XPtr<SignificantPattern::SignificantItemsetSearch> ptr(inst);
//     ptr->readPlinkFiles(base_filename);
// }

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_read_plink_files_with_cov_fastcmh(SEXP inst, std::string base_filename, std::string covfilename, std::string encoding) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchFastCmh> ptr(inst);
    ptr->readPlinkFilesWithCovariates(base_filename, covfilename, true, encoding);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_read_plink_files_with_cov_facs(SEXP inst, std::string base_filename, std::string covfilename, std::string encoding) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearchFacs> ptr(inst);
    ptr->readPlinkFilesWithCovariates(base_filename, covfilename, true, encoding);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_read_covariates_file_fastcmh(SEXP inst, std::string cov_filename) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchFastCmh> ptr(inst);
    ptr->readCovariatesFile(cov_filename);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_read_covariates_file_facs(SEXP inst, std::string cov_filename) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearchFacs> ptr(inst);
    ptr->readCovariatesFile(cov_filename);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_write_eth_files_iset(SEXP inst, std::string x_filename, std::string y_filename) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearch> ptr(inst);
    ptr->writeETHFiles(x_filename, y_filename);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_write_eth_files_int(SEXP inst, std::string x_filename, std::string y_filename) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
    ptr->writeETHFiles(x_filename, y_filename);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_write_eth_files_with_cov_fastcmh(SEXP inst, std::string x_filename, std::string y_filename, std::string covfilename) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchFastCmh> ptr(inst);
    ptr->writeETHFilesWithCovariates(x_filename, y_filename, covfilename);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_write_eth_files_with_cov_facs(SEXP inst, std::string x_filename, std::string y_filename, std::string covfilename) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearchFacs> ptr(inst);
    ptr->writeETHFilesWithCovariates(x_filename, y_filename, covfilename);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_execute_iset(SEXP inst, double alpha, longint l_max) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearch> ptr(inst);
    ptr->execute(alpha, l_max);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_execute_int(SEXP inst, double alpha, longint l_max) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
    ptr->execute(alpha, l_max);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_summary_write_to_file_fais(SEXP inst, std::string output_file) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchFais> ptr(inst);
    ptr->getSummary().writeToFile(output_file);
}


//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_summary_write_to_file_fastcmh(SEXP inst, std::string output_file) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchFastCmh> ptr(inst);
    ptr->getSummary().writeToFile(output_file);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_summary_write_to_file_facs(SEXP inst, std::string output_file) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearchFacs> ptr(inst);
    ptr->getSummary().writeToFile(output_file);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_profiler_write_to_file(SEXP inst, std::string output_file) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
    ptr->getProfiler().writeToFile(output_file);
}
//TODO: profiler_write_to_file_{int,iset}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_filter_intervals_write_to_file(SEXP inst, std::string output_file) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
    ptr->getFilteredIntervals().writeToFile(output_file);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_pvals_testable_ints_write_to_file(SEXP inst, std::string output_file) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
    ptr->getPValsTestableInts().writeToFile(output_file);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_pvals_significant_ints_write_to_file(SEXP inst, std::string output_file) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
    ptr->getPValsSigInts().writeToFile(output_file);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_pvals_testable_isets_write_to_file(SEXP inst, std::string output_file) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearch> ptr(inst);
    ptr->getPValsTestableIsets().writeToFile(output_file);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
void lib_pvals_significant_isets_write_to_file(SEXP inst, std::string output_file) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearch> ptr(inst);
    ptr->getPValsSigIsets().writeToFile(output_file);
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
SEXP lib_get_significant_intervals(SEXP inst) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
    return _get_intervals(ptr->getSignificantIntervals().getSigInts());
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
SEXP lib_get_filtered_intervals(SEXP inst) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
    return _get_intervals(ptr->getFilteredIntervals().getSigInts());
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
SEXP lib_get_significant_itemsets(SEXP inst) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearch> ptr(inst);
    return _get_itemsets(ptr->getPValsSigIsets());
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
SEXP lib_get_result_fais(SEXP inst) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearchFais> ptr(inst);
    SignificantPattern::SummaryFais sum = ptr->getSummary();

    SEXP start = IntegerVector::create(sum.getSl1(), sum.getSu1());
    SEXP end = IntegerVector::create(sum.getSl2(), sum.getSu2());
    SEXP region = List::create(Named("start") = start,
                               Named("end") = end);

    return List::create(Named("n.int.processed") = sum.getNumFeaturesProcessed(),
                        Named("n.int.testable") = sum.getm(),
                        Named("testability.region") = region,
                        Named("testability.threshold") = sum.getDelta(),
                        Named("target.fwer") = sum.getAlpha(),
                        Named("corrected.significance.threshold") = sum.getDelta_opt());
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
SEXP lib_get_result_int(SEXP inst) {
    Rcpp::XPtr<SignificantPattern::SignificantIntervalSearch> ptr(inst);
    SignificantPattern::SummaryInt sum = ptr->getSummary();

    return List::create(Named("n.int.processed") = sum.getNumFeaturesProcessed(),
                        Named("n.int.testable") = sum.getm(),
                        Named("testability.threshold") = sum.getDelta(),
                        Named("target.fwer") = sum.getAlpha(),
                        Named("corrected.significance.threshold") = sum.getDelta_opt());
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
SEXP lib_get_result_iset(SEXP inst) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearch> ptr(inst);
    SignificantPattern::SummaryIset sum = ptr->getSummary();

    return List::create(Named("n.iset.processed") = sum.getNumFeaturesProcessed(),
                        Named("n.iset.testable") = sum.getm(),
                        Named("testability.threshold") = sum.getDelta(),
                        Named("target.fwer") = sum.getAlpha(),
                        Named("corrected.significance.threshold") = sum.getDelta_opt());
}

//' Internal function
//' @keywords internal
// [[Rcpp::export]]
SEXP lib_get_result_facs(SEXP inst) {
    Rcpp::XPtr<SignificantPattern::SignificantItemsetSearchFacs> ptr(inst);
    SignificantPattern::SummaryFacs sum = ptr->getSummary();

    return List::create(Named("n.iset.processed") = sum.getNumFeaturesProcessed(),
                        Named("n.iset.closed.processed") = sum.getNumItemsetsClosedProcessed(),
                        Named("n.iset.testable") = sum.getm(),
                        Named("testability.threshold") = sum.getDelta(),
                        Named("target.fwer") = sum.getAlpha(),
                        Named("corrected.significance.threshold") = sum.getDelta_opt());
}
