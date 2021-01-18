/*
 * TransactionsData.cpp
 *
 *  Created on: 9 May 2017
 *      Author: mikolajr
 */

#include <cassert>
#include <numeric> // iota
#include <set>
#include <vector>
#include<algorithm> //sort

#include "TransactionsData.h"
#include "types.h"

namespace SignificantPattern {

TransactionsData::TransactionsData() {}
TransactionsData::~TransactionsData() {}

void TransactionsData::clear() {
    transactions.clear();
    transactions_bool.clear();
    items.clear();
    supports.clear();
}

void TransactionsData::initFromData(const Genotype &data,
                                    const Phenotype &covariates,
                                    std::vector<longint> &item_label_map,
                                    std::vector<longint> &pexs) {

    // 1. Convert ETH internal binary data format to transactions format
    // (transposed; sparse array instead of boolean array)
    unsigned char **X_tr = data.getMatrixPtr();
    longint L = data.getNumFeatures(), N = data.getNumObservations();

    // List of transactions
    std::vector<std::vector<longint>> T_list(N);
    // Set of all items in original transaction database
    // indices of rows containing at least one "1" (OR over column vectors)
    std::set<longint> items_set;

    for (longint j = 0; j < L; j++) {
        for (longint i = 0; i < N; i++) {
            if (X_tr[j][i])
                T_list[i].push_back(j); // Attention: in some examples transactions indices started at 1
            items_set.insert(T_list[i].begin(), T_list[i].end());
        }
    }

    // 2. Build temporary transaction database

    // Total number of items in transaction database
    size_t n_items = items_set.size();

    // Temporal transaction database (in vertical representation)
    TransactionsData tmp_db;
    tmp_db.transactions.resize(n_items);

    // Relabel items as consecutive integers, keeping track of the mapping
    // back to their original labels.
    for (auto item : items_set)
        item_label_map.push_back(item);

    // sum_of_transaction_sizes[i]=sum of sizes of all transactions containing
    // item i
    std::vector<longint> sum_of_transaction_sizes(n_items);
    longint n_trans = 0, tmp_idx;
    for (auto trans : T_list) {
        // Add transaction index to transaction lists of all items contained in
        // transaction. For that, we need to map the original label index to the
        // consecutive indexing using internally, using find() and distance()
        // for that purpose
        for (auto item : trans) {
            tmp_idx = std::distance(
                item_label_map.begin(),
                std::find(item_label_map.begin(), item_label_map.end(), item));
            tmp_db.transactions[tmp_idx].push_back(n_trans);
            sum_of_transaction_sizes[tmp_idx] += trans.size();
        }
        // Increment transaction index
        n_trans++;
    }
    assert(N == n_trans);

    // Sort items in decreasing order of cumulative size of all transactions
    // they belong to (rationale, items contained in very large transactions
    // will lead to denser conditional transaction databases, and should be
    // processed in the "shallow" part of the enumeration tree for extra
    // performance)
    std::vector<longint> idx_sort(n_items);
    std::iota(idx_sort.begin(), idx_sort.end(), 0);
    std::sort(idx_sort.begin(), idx_sort.end(), [&](size_t i, size_t j) {
        return sum_of_transaction_sizes[i] > sum_of_transaction_sizes[j];
    });

    // 3. Create final transaction database, ordered according to idx_sort above
    // and storing perfect extensions of the root separately

    // Category (integer between 0 and n_cat-1) of each sample
    unsigned char *cats = covariates.getVectorPtr();
    // std::vector<unsigned short> cats(N);
    // unsigned char *cats_char = covariates.getVectorPtr();
    // std::copy(cats_char, cats_char + N, &cats[0]);

    // Number of different categories for the covariate
    short n_cat = covariates.getNumClasses();

    // Clear this before adding new elements
    this->clear();
    longint current_support; // temporal variable to store support of each item
    for (auto i_sorted : idx_sort) {
        current_support = tmp_db.transactions[i_sorted].size();
        if (current_support == N)
            pexs.push_back(i_sorted);
        else {
            this->transactions.push_back(tmp_db.transactions[i_sorted]);
            this->transactions_bool.push_back(std::vector<bool>(N));
            for (auto trans_idx : this->transactions.back())
                this->transactions_bool.back()[trans_idx] = true;
            this->items.push_back(i_sorted);
            // Compute supports for each category
            this->supports.push_back(std::vector<longint>(n_cat));
            for (auto trans_idx : this->transactions.back())
                this->supports.back()[cats[trans_idx]]++;
        }
    }

}

} /* namespace SignificantPattern */
