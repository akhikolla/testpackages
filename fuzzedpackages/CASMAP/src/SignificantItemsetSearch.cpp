/*
 * SignificantItemsetSearch.cpp
 *
 *  Created on: 15 May 2017
 *      Author: mikolajr
 */

#include "SignificantItemsetSearch.h"

// #include <ostream> // std::cout


using namespace std;


namespace SignificantPattern
{

SignificantItemsetSearch::SignificantItemsetSearch()
    : SignificantFeaturesSearch() // explicitly construct a virtual base class
{
    #ifdef DEBUG
    fprintf(stderr, "SignificantItemsetSearch()\n");
    #endif
    execute_constructor_iset();
}
SignificantItemsetSearch::~SignificantItemsetSearch() {
    #ifdef DEBUG
    fprintf(stderr, "~SignificantItemsetSearch()\n");
    #endif
    execute_destructor_iset();
}

void SignificantItemsetSearch::execute_constructor() {
    super::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantItemsetSearch::execute_constructor()\n");
    #endif
    execute_constructor_iset();
}
void SignificantItemsetSearch::execute_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantItemsetSearch::execute_destructor()\n");
    #endif
    execute_destructor_iset();
    super::execute_destructor();
}

void SignificantItemsetSearch::execute_constructor_iset() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantItemsetSearch::execute_constructor_iset()\n");
    #endif
    setPValsTestableIsets(ItemsetSetWithOddsRatio());
    setPValsSigIsets(ItemsetSetWithOddsRatio());
}
void SignificantItemsetSearch::execute_destructor_iset(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantItemsetSearch::execute_destructor_iset()\n");
    #endif
}



bool SignificantItemsetSearch::testAndSaveItemset(
    double threshold, double score, double odds_ratio, double pval, const vector<longint> &x_t, longint a,
    const vector<longint> &iset, const vector<longint> &pexs)
{
    std::vector<longint> itemset; bool itemsetBuilt = false;
    if (outputForTestableInts) {

        // std::cout << pval << ',' << a << ',';

        buildItemset(x_t, iset, pexs, itemset); itemsetBuilt = true;
        saveTestableItemset(pval, score, odds_ratio, itemset, a);
    }
    bool isSignificant = (pval <= threshold);
    if(isSignificant)
    {

        // std::cout << pval << ',' << a << ',';

        if (!itemsetBuilt) buildItemset(x_t, iset, pexs, itemset);
        saveSignificantItemset(pval, score, odds_ratio, itemset, a);
        n_significant_featuresets++;
    }
    return isSignificant;
}

void SignificantItemsetSearch::process_significant_features() {

    // Filter sig itemsets against the final threshold
    ItemsetSetWithOddsRatio finalPValsSigIsets;

    const ItemsetSetWithOddsRatio pValsSigIsets = getPValsSigIsets();
    std::vector< std::vector<longint> > itemsets = pValsSigIsets.getItemsetsVector();
    std::vector<double> score = pValsSigIsets.getScoreVector();
    std::vector<double> odds_ratio = pValsSigIsets.getOddsRatioVector();
    std::vector<double> pvals = pValsSigIsets.getPValueVector();
    std::vector<longint> alphas = pValsSigIsets.getAlphaVector();

    for(size_t i = 0; i < pValsSigIsets.getLength(); i++)
        if (pvals[i] <= delta_opt)
            finalPValsSigIsets.addFeature(itemsets[i], alphas[i], score[i], odds_ratio[i], pvals[i]);
    setPValsSigIsets(finalPValsSigIsets);
    this->n_significant_featuresets = getPValsSigIsets().getLength();
}


} /* namespace SignificantPattern */
