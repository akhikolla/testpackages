/*
 * TransactionsData.h
 *
 *  Created on: 9 May 2017
 *      Author: mikolajr
 */

#ifndef TRANSACTIONSDATA_H_
#define TRANSACTIONSDATA_H_

#include <vector>

#include "Genotype.h"
#include "Phenotype.h"
#include "types.h"

namespace SignificantPattern
{

class TransactionsData
{
public:
    // List of lists holding the transactions in which each item appears
    std::vector< std::vector<longint> > transactions;
    // List of boolean vectors, such that transactions_bool[i][j]==true iff
    // j is contained in transactions[i]. This redundancy will in fact help
    // make the intersection of transaction lists more efficiently
    std::vector< std::vector<bool> > transactions_bool;
    // List holding the corresponding item IDs
    std::vector<longint> items;
    // List holding the corresponding supports
    std::vector< std::vector<longint> > supports;

    void clear();
    void initFromData(const Genotype &data, const Phenotype &covariates,
                      std::vector<longint> &item_label_map,
                      std::vector<longint> &pexs);

    TransactionsData();
    virtual ~TransactionsData();
};




} /* namespace SignificantPattern */

#endif /* TRANSACTIONSDATA_H_ */
