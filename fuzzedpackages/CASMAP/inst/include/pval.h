/*
 * pval.h
 *
 *  Created on: 27 Mar 2017
 *      Author: mikolajr
 */

#ifndef PVAL_H
#define PVAL_H

#include "types.h"

/**
 * Precompute minimum attainable p-values psi(x) for all x in [0,N].
 */
void fisher_minpvals(longint N, longint n, longint N_over_2,
                     /*out*/ double *psi);

/**
 * Evaluate Fisher's exact test on a table with margins x, n and N and cell count a.
 * The p-value is defined as a two-tailed p-value which adds up the probabilities
 * of all tables less or equally likely to occur than the one we observed
 */
double fisher_pval(longint a, longint x, longint N, longint n, double *loggamma, double log_inv_binom_N_n);

/**
 * Precompute minimum attainable p-values psi(x) for all x in [0,N] using Chi2
 * CDF.
 */
void chi2_minpvals(longint N, longint n, longint N_over_2, double class_ratio,
                   double class_ratio_bin,
                   /*out*/ double *psi);

/**
 * Evaluate Fisher's Chi2-approximated P-value for a setting like in fisher_pval().
 */
double chi2_pval(longint a, longint x, longint N, longint n, double class_ratio_bin);

/**
 * Evaluate Fisher's Chi2-approximated test statistic for a setting like in fisher_pval().
 */
double chi2_score(longint a, longint x, longint N, longint n, double class_ratio_bin);

/**
 * Compute odds ratio for a single contingency table.
 */
double odds_ratio(longint a, longint x, longint N, longint n);

/**
 * Computes the minimum attainable CMH p-value depending on the margins x, n and
 * N for the K tables
 */
double cmh_minpval(longint *x, unsigned short K, longint* Nt, longint* nt,
        //precomputed quantities to save time when computing the CMH test statistic
        longint* Nt_nt, double* gammat, double* gammabint);

/**
 * Computes the CMH p-value as a function of the margins x, n and N and the cell
 * counts a for the K tables
 */
double cmh_pval(longint a, longint *x, unsigned short K, longint* Nt,
        double* gammat, double* gammabint);

/**
 * Computes the CMH test statistic as a function of the margins x, n and N and the cell
 * counts a for the K tables
 */
double cmh_score(longint a, longint *x, unsigned short K, longint* Nt,
                double* gammat, double* gammabint);

/**
 * Compute odds ratio for K contingency tables.
 */
double cmh_odds_ratio(longint* at, longint* x, unsigned short K, longint* Nt, longint* nt);

#endif


