/*
 * SignificantItemsetSearchFacs.cpp
 *
 *  Created on: 8 May 2017
 *      Author: mikolajr
 */

#include "SignificantItemsetSearchFacs.h"

#include "algorithm" // min_element
#include <numeric>
// #include "ostream" // cout

namespace SignificantPattern
{

SignificantItemsetSearchFacs::SignificantItemsetSearchFacs()
    : SignificantFeaturesSearch(), // explicitly construct a virtual base class
      SignificantItemsetSearch(), SignificantFeaturesSearchTaroneCmh()
{
    #ifdef DEBUG
    fprintf(stderr, "SignificantItemsetSearchFacs()\n");
    #endif
    execute_constructor_facs();
}
SignificantItemsetSearchFacs::~SignificantItemsetSearchFacs() {
    #ifdef DEBUG
    fprintf(stderr, "~SignificantItemsetSearchFacs()\n");
    #endif
    execute_destructor_facs();
}

void SignificantItemsetSearchFacs::execute_constructor() {
    super_cov::execute_constructor();
    super_feat::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantItemsetSearchFacs::execute_constructor()\n");
    #endif
    execute_constructor_facs();
}
void SignificantItemsetSearchFacs::execute_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantItemsetSearchFacs::execute_destructor()\n");
    #endif
    execute_destructor_facs();
    super_feat::execute_destructor();
    super_cov::execute_destructor();
}

void SignificantItemsetSearchFacs::execute_constructor_facs() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantItemsetSearchFacs::execute_constructor_facs()\n");
    #endif
    n_enumerated = 0; n_enumerated_closed = 0; r_out = 0;

    db.clear();
}
void SignificantItemsetSearchFacs::execute_destructor_facs(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantItemsetSearchFacs::execute_destructor_facs()\n");
    #endif

    minpvals.clear();
    lower_envelope_minpvals.clear();

    tentative_sig_ths.clear();
}



void SignificantItemsetSearchFacs::algorithm_init() {
    super_cov::algorithm_init();
    super_feat::algorithm_init();

    db.initFromData(getGenotype(), getCovariates(), item_label_map, pexs);

    for (auto support: db.supports) {
        minpvals.push_back(compute_minpval(support.data()));
        lower_envelope_minpvals.push_back(
            compute_envelope_minpval(support.data()));
    }

    tentative_sig_ths.clear();
    if (if_compute_pvals) tentative_sig_ths.push_back(alpha);


    #ifdef DEBUG
    // std::cerr << "supports:" << std::endl;
    // for (auto support: db.supports) {
    //     for (auto item: support) std::cerr << item << ",";
    //     std::cerr << std::endl;
    // }
    // std::cerr << "minpvals:" << std::endl;
    // for (auto minpval: minpvals) std::cerr << minpval << std::endl;
    // std::cerr << "lower_envelope_minpvals:" << std::endl;
    // for (auto lempv: lower_envelope_minpvals) std::cerr << lempv << std::endl;
    fprintf(stderr, "SignificantItemsetSearchFacs::algorithm_init(), tentative_sig_ths=%g\n", tentative_sig_ths.back());
    #endif

}
void SignificantItemsetSearchFacs::algorithm_end(){
    super_feat::algorithm_end();
    super_cov::algorithm_end();

    SummaryFacs& summary = getSummaryRef();
    summary.setNumItemsetsClosedProcessed(n_enumerated_closed);
}





/* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */

void SignificantItemsetSearchFacs::buildItemset(const std::vector<longint> &x_t,
                                                const std::vector<longint> &iset,
                                                const std::vector<longint> &pexs,
                                                /* out */ std::vector<longint> &itemset)
{
    // Note: originally was writting to "std::cout", not std::cout
    //       and for testable items wrote only pval (output function empty)

    // for (auto x : x_t)
    //     std::cout << x << ' ';
    // for (size_t i = 0; i < iset.size() - 1; ++i)
    //     std::cout << item_label_map[iset[i]] << ' ';
    // if (pexs.size() > 0) {
    //     std::cout << item_label_map[iset[iset.size() - 1]] << ' ';
    //     std::cout << "pexs: ";
    //     for (size_t i = 0; i < pexs.size() - 1; ++i)
    //         std::cout << item_label_map[pexs[i]] << ' ';
    //     std::cout << item_label_map[pexs[pexs.size() - 1]] << std::endl;
    // } else
    //     std::cout << item_label_map[iset[iset.size() - 1]] << std::endl;

    itemset.clear();
    //itemset.reserve(x_t.size()+iset.size()+pexs.size());
    //itemset.insert(itemset.end(), x_t.begin(), x_t.end());
    itemset.reserve(iset.size()+pexs.size());
    for (auto item : iset) itemset.push_back(item_label_map[item]);
    for (auto pex : pexs) itemset.push_back(item_label_map[pex]);
}

void SignificantItemsetSearchFacs::decrease_threshold() {
        super_cov::decrease_threshold();
        if (if_compute_pvals) tentative_sig_ths.push_back(alpha/m);
}

void SignificantItemsetSearchFacs::find_significant_features() {
	n_enumerated = 0;  // Initialize number of enumerated itemsets to zero
	n_enumerated_closed = 0;  // Same for number of enumerated closed itemsets
	std::vector<longint> iset;
	std::vector<std::vector<bool> *> elim;

    // Clear the current layer frequency counters
    freq_clear();

	r_out = depth(db, minpvals, lower_envelope_minpvals, iset, pexs, elim);

    // Compute final corrected significance threshold only after depth-first
    // search which finds number of testable intervals
    compute_final_corrected_significance_threshold();

    // Significant itemsets are filtered again against the final corrected
    // threshold in process_significant_features()

    n_featuresets_processed = n_enumerated;
}

longint SignificantItemsetSearchFacs::depth(
    TransactionsData &db, std::vector<double> &minpvals,
    std::vector<double> &lower_envelope_minpvals, std::vector<longint> &iset,
    std::vector<longint> &pexs, std::vector<std::vector<bool> *> &elim) {

    // Declare variables for loops
    longint i, k, kk;
    // Declare variables to hold supports
    size_t current_support, max_support = 0;
    // Variable to check if itemsets are closed
    bool closed;
    // Variable to hold the value of the minimum attainable p-value and its
    // corresponding lower envelope
    double minpval, lower_envelope_minpval;
    // Number of items in the current (conditional) transaction database
    longint db_size = db.items.size();
    // Number of perfect extensions found during the processing of a new itemset
    longint n_proj_pexs;
    // Output of recursion
    size_t r_out;
    // Variables to compute itemset pvalue (only used if_compute_pvals==true)
    longint a;
    double score, odds_ratio, pval;

    unsigned char *labels = getPhenotype().getVectorPtr();
    unsigned char *cats = getCovariates().getVectorPtr();
    unsigned short n_cat = getCovariates().getNumClasses();

    std::vector<longint> at(n_cat);

    // Insert all transaction lists which must be used during closure check
    for (kk = (db_size - 1); kk > 0; --kk)
        elim.push_back(&(db.transactions_bool[kk]));

    // Start processing all possible child nodes
    for (k = 0; k < db_size; ++k) {
        current_support = db.transactions[k].size();
        // Between the time the item {db.items[k]} was added to the current
        // conditional database during the processing of the parent itemset, and
        // now, the prunability threshold might have changed, due to the earlier
        // depth-first processing of its siblings. Therefore, it's worth
        // checking again if the itemset can be pruned.
        if (not istestable(minpvals[k])) {
            // If the itemset is not testable, it might be prunable. In this
            // case, we must compute the minpval lower envelope if it was not
            // yet computed
            if (lower_envelope_minpvals[k] == -2)
                lower_envelope_minpvals[k] =
                    compute_envelope_minpval(db.supports[k].data());
            // Based on the value of the minpval lower envelope, check if
            // itemset is prunable
            if (isprunable(lower_envelope_minpvals[k])) {
                if (k < (db_size - 1))
                    elim.pop_back();
                continue;
            }
        }
        n_enumerated++; // Increase count of enumerated nodes
        max_support =
            (current_support > max_support)
                ? current_support
                : max_support; // max_support=max(current_support,max_support)

        // Check closure
        closed = true;
        for (kk = (elim.size() - 1); kk >= 0; --kk) {
            closed = false;
            for (auto trans_idx : db.transactions[k]) {
                if (not(*elim[kk])[trans_idx]) {
                    closed = true;
                    break;
                }
            }
            if (not closed)
                break;
        }
        if (not closed) {
            if (k < (db_size - 1))
                elim.pop_back();
            continue;
        }

        // Build new conditional transaction database
        TransactionsData proj_db;
        std::vector<double> proj_minpvals;
        std::vector<double> proj_lower_envelope_minpvals;
        n_proj_pexs = 0;
        for (kk = 0; kk < k; ++kk) {
            proj_db.transactions.push_back(std::vector<longint>());
            std::vector<longint> &U = proj_db.transactions.back();
            proj_db.transactions_bool.push_back(
                std::vector<bool>(db.transactions_bool[k]));
            std::vector<bool> &U_hash = proj_db.transactions_bool.back();
            // Compute intersection
            for (auto trans_idx : db.transactions[k]) {
                if (db.transactions_bool[kk][trans_idx])
                    U.push_back(trans_idx);
                else
                    U_hash[trans_idx] = false;
            }
            // Check if item is a perfect extension
            if (U.size() == current_support) {
                pexs.push_back(db.items[kk]);
                n_proj_pexs++;
                proj_db.transactions.pop_back();
                proj_db.transactions_bool.pop_back();
            }
            // If not, we might want to add the item to the next conditional
            // transaction database
            else {
                proj_db.supports.push_back(std::vector<longint>(n_cat));
                std::vector<longint> &proj_supports = proj_db.supports.back();
                // Compute supports for each category
                for (auto trans_idx : U)
                    proj_supports[cats[trans_idx]]++;
                minpval = compute_minpval(proj_supports.data());
                // If the itemset is testable, it cannot be prunable
                if (istestable(minpval)) {
                    proj_db.items.push_back(db.items[kk]);
                    proj_minpvals.push_back(minpval);
                    // Avoid computing the minpval lower envelope. Set its value
                    // to -2 (an arbitrary value outside its normal range), to
                    // indicate that the computation was skipped
                    proj_lower_envelope_minpvals.push_back(-2);
                }
                // If it is not testable, it could be prunable or not. We need
                // to compute the minpval lower envelope to answer that
                else {
                    lower_envelope_minpval =
                        compute_envelope_minpval(proj_supports.data());
                    if (not isprunable(lower_envelope_minpval)) {
                        proj_db.items.push_back(db.items[kk]);
                        proj_minpvals.push_back(minpval);
                        proj_lower_envelope_minpvals.push_back(
                            lower_envelope_minpval);
                    } else {
                        proj_db.transactions.pop_back();
                        proj_db.transactions_bool.pop_back();
                        proj_db.supports.pop_back();
                    }
                }
            }
        }
        // Invoke next step of depth-first search recursively
        iset.push_back(db.items[k]);
        r_out = (proj_db.items.size() > 0)
                    ? depth(proj_db, proj_minpvals,
                            proj_lower_envelope_minpvals, iset, pexs, elim)
                    : 0;

        #ifdef DEBUG
        fprintf(stderr, "depth(), r_out=%ld, current_support=%ld, k=%ld, minpvals[k]=%g, istestable(...)=%d\n",
            r_out, current_support, k, minpvals[k], istestable(minpvals[k]));
        #endif

        // If the itemset has no perfect extension and is testable, process it
        if ((r_out < current_support) and istestable(minpvals[k])) {
            // Increase counter of testable itemsets, and upgrade histogram of
            // supports
            freq_cnt_cmh[bucket_idx(minpvals[k])]++;
            m++;
            update_threshold();

            // If we have to compute pvalues for testable, and the pattern is
            // still testable, we assume that it will remain testable. This
            // won't be true for some patterns, leading to a few unnecessary
            // pvalues being computed. However, the extra runtime will in
            // practice be still smaller than that of executing the entire
            // enumeration process a second time, once the final testability
            // threshold is known
            if (if_compute_pvals and istestable(minpvals[k])) {
                // Compute cell-count
                //a = 0;
                //for (auto trans_idx : db.transactions[k])
                //    if (labels[trans_idx]) a++;
                std::fill(at.begin(), at.end(), 0);
                for (auto trans_idx : db.transactions[k])
                    if (labels[trans_idx]) at[cats[trans_idx]]++;
                a = std::accumulate(at.begin(), at.end(), 0);

                score = compute_score(a, db.supports[k].data());
                odds_ratio = compute_odds_ratio(at.data(), db.supports[k].data());
                pval = score_to_pval(score);
                //pval = compute_pval(a, db.supports[k].data());
                // TODO: This step is greedy. In pathological cases, some
                // significant itemsets could be lost...
                testAndSaveItemset(tentative_sig_ths.back(), score, odds_ratio, pval,
                                   db.supports[k], a, iset, pexs);
            }
        }
        if (r_out < current_support)
            n_enumerated_closed++; // Increase count of enumerated closed
                                   // itemsets

        // Undo changes to variables shared across recursion iterations (iset,
        // pexs and elim)
        for (i = 0; i < n_proj_pexs; ++i)
            pexs.pop_back();
        iset.pop_back();
        if (k < (db_size - 1))
            elim.pop_back();
    }
    return max_support;
}


/// @todo
/// Add list of tentative significance thresholds, or maybe only the too strict
/// ones to the results?
void SignificantItemsetSearchFacs::process_significant_features() {
    super_feat::process_significant_features();

    // // Output list of tentative significance thresholds for debugging purposes
    // std::cout << "TENTATIVE SIGNIFICANCE THRESHOLDS:" << std::endl;
    // for (int i = 0; i < tentative_sig_ths.size() - 1; ++i)
    //     std::cout << tentative_sig_ths[i] << "\t";
    // std::cout << tentative_sig_ths[tentative_sig_ths.size() - 1] << std::endl;
    // // Finally, check if some of the tentative significance thresholds were too
    // // strict, and thus patterns might have been missed
    // double min_tentative_sig_th =
    //     *min_element(tentative_sig_ths.begin(), tentative_sig_ths.end());
    // if (min_tentative_sig_th < delta_opt)
    //     std::cout
    //         << "WARNING!: The greedy p-value evaluation approach used "
    //            "threshold "
    //         << min_tentative_sig_th
    //         << " during the procedure. Some significant patterns might be lost."
    //         << std::endl;
    // else
    //     std::cout << "Greedy p-value evaluation approach was successful. Have "
    //                  "a nice day :)"
    //               << std::endl;
}

} /* namespace SignificantPattern */
