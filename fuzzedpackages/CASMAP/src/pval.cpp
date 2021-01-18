/*
 * pval.h
 *
 *  Created on: 27 Mar 2017
 *      Author: mikolajr
 */

#include "pval.h"

#include <math.h>

#include <limits>

#include "chi2.h" // Chi2_sf
#include "double_comp.h" // doubleeq

//#define RATIO_TH -23
#define LOGPVAL_EQ_RELTOL 1E-12 // exp(-27.6) \approx 1E-12

using namespace std;

void chi2_minpvals(longint N, longint n, longint N_over_2, double class_ratio,
                   double class_ratio_bin,
                   /*out*/ double *psi) {
    double num, den;
    longint x;

    // First, compute the left side of "the W", i.e. the range [0,n]
    psi[0] = 1;
    // In this range, the recursion psi(x)=$\psi(x-1)$*[(n-(x-1))/(N-(x-1))] can
    // be seen to be correct
    for(x=1; x<=n; x++) {
        num = x*(1-class_ratio); num *= num;
        den = x*(1-((double)x)/N)*class_ratio_bin;
        psi[x] = Chi2_sf(num/den,1);
    }

    // Now, compute the minimum attainable p-values in the range [N-N_over_2,N]
    for(x=(n+1); x<=N_over_2; x++) {
        num = n*(1-((double)x)/N); num *= num;
        den = x*(1-((double)x)/N)*class_ratio_bin;
        psi[x] = Chi2_sf(num/den,1);
    }

    // Finally, since psi(x) is symmetric around N_over_2, complete the right half by symmetry
    for(x=(N_over_2+1); x<=N; x++) psi[x] = psi[N-x];
}

double chi2_pval(longint a, longint x, longint N, longint n, double class_ratio_bin)
{
    double aux, num, den;
    aux = ((double)x)/N;
    num = a-n*aux; 
    //num = pow(num,2); //this is simply num squared
    //need to replace pow, because pow(longint, int)
    num = num*num;
    den = x*(1-aux)*class_ratio_bin;
    return Chi2_sf(num/den,1);
}

double chi2_score(longint a, longint x, longint N, longint n, double class_ratio_bin)
{
    double aux, num, den;
    aux = ((double)x)/N;
    num = a-n*aux; 
    //num = pow(num,2); //this is simply num squared
    //need to replace pow, because pow(longint, int)
    num = num*num;
    den = x*(1-aux)*class_ratio_bin;
    return num/den;
}

double odds_ratio(longint a, longint x, longint N, longint n)
{
    double num, den;
    num = a * ((N + a) - (x + n));
    den = (n - a) * (x - a);
    if (den==0) return std::numeric_limits<double>::infinity();  //TODO: Test
    else return num/den;
}


void fisher_minpvals(longint N, longint n, longint N_over_2,
                     /*out*/ double *psi) {
    double xi1;
    longint x, x_init;

    // First, compute the left side of "the W", i.e. the range [0,n]
    psi[0] = 1;
    // In this range, the recursion psi(x)=$\psi(x-1)$*[(n-(x-1))/(N-(x-1))] can
    // be seen to be correct
    for(x=1; x<=n; x++) psi[x] = (((double)(n-(x-1)))/(N-(x-1)))*psi[x-1];

    // Now, start computing xi1 in the range [N-N_over_2,N] using another
    // recursion, this time starting in N.
    // Note that we don't need to store all values, since this will be used only
    // to initialise psi[N_over_2]
    x_init = N-N_over_2;
    xi1 = 1;
    // In this range, the recursion
    // $\xi_{1}(x)$=$\xi_{1}(x+1)$*[((x-1)-n)/(x+1)] can be seen to be correct
    for(x=(N-1); x>=x_init; x--) xi1 = (((double)((x+1)-n))/(x+1))*xi1;

    // Now, use the value of $\xi_{1}(N-N_over_2)$=xi1[0] to get
    // $\psi(N_over_2)$=psi[N_over_2] using the same recursion if N is odd, or
    // simply copy the value of xi1[0] since
    // $\xi_{1}(N-N_over_2)=\psi(N_over_2)$ if N is even
    if (N % 2) psi[N_over_2] = (((double)(x_init-n))/x_init)*xi1;
    else psi[N_over_2] = xi1;

    // Now, using psi[N_over_2] compute the right side of "the W", i.e. the
    // range [n+1,N_over_2] using the same recursion as for $\xi_{1}$
    for(x=(N_over_2-1); x > n; x--) psi[x] = (((double)((x+1)-n))/(x+1))*psi[x+1];

    // Finally, since psi(x) is symmetric around N_over_2, complete the right
    // half by symmetry
    for(x=x_init; x<=N; x++) psi[x] = psi[N-x];
}

double fisher_pval(longint a, longint x, longint N, longint n, double *loggamma, double log_inv_binom_N_n)
{
    longint a_min, a_max;
    double log_p_left, log_p_right;
    double p_left, p_right, pval;
    //double p_diff, p_ratio_left, p_ratio_right;
    double pre_comp_xterms;

    // Compute the contribution of all terms depending on x but not on a
    pre_comp_xterms = loggamma[x] + loggamma[N-x];
    a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
    a_max = (x > n) ? n : x;//min(x,n)


    // The algorithm to compute the p-value proceeds as follows. We inspect simultaneously probability values on the left and right tails of the
    // hypergeometric distribution, "accepting" each time the one that is smaller. When that happens, we move the index in the appropriate direction,
    // that is, a_min++ if we "accept" a value on the left and a_max-- if we "accept" a value on the right. As soon as the "accepted" value is located
    // in index a, we know that we have already explored all values of the hypergeometric probability mass whose probabilities are smaller or equal
    // than the probability of a. The only tricky case occurs when the probability values on both sides happen to be identical. The way to handle
    // that case is by "accepting" both values simultaneously.
    pval = 0; //Accumulate probabilities in this variable
    while(a_min<a_max){
            log_p_left = pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_min] + loggamma[n-a_min] + loggamma[x-a_min] + loggamma[(N-n)-(x-a_min)]);
            log_p_right = pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_max] + loggamma[n-a_max] + loggamma[x-a_max] + loggamma[(N-n)-(x-a_max)]);
            p_left = exp(log_p_left);
            p_right = exp(log_p_right);


            // Note: in Felipe's p-values equality test implementation
            //
            //    p_diff = p_left - p_right;
            //    p_diff = (p_diff > 0) ? p_diff : -p_diff;
            //    p_diff = log(p_diff);
            //    p_ratio_left = p_diff - log_p_left;
            //    p_ratio_right = p_diff - log_p_right;
            //
            //    if( (p_ratio_left <= RATIO_TH) || (p_ratio_right <= RATIO_TH) ) {
            //        ...
            //
            //  we have:
            //
            //    min(p_ratio_left, p_ratio_right) == log(1-exp( -fabs(log_p_right - log_p_left) ))
            //
            // so it simplifies to a more numerically stable (`1-exp(a-b)`
            // instead of `exp(a)-exp(b)`), and more concise:
            //
            //     if ( log( 1-exp(-fabs(log_p_right - log_p_left)) ) <= RATIO_TH ) {
            //         ...
            //
            // BEWARE: if actually p_left == p_right
            //         (equivalently log_p_right == log_p_left), then above one
            //         would get log(0), which likely will raise an error (dep.
            //         on C lib impl.)
            //
            // Instead, use plain relative error of exponents instead (since
            // exp() is strictly increasing monotone function, and exponents
            // here are not small negative numbers, so no worries reg. use of
            // `/` operator):
            if (doubleeq(log_p_left, log_p_right, LOGPVAL_EQ_RELTOL)) {
                    pval += (p_left+p_right);
                    if((a==a_min) || (a==a_max)) return pval;
                    a_min++; a_max--;
            }
            else if(p_left < p_right){
                    pval += p_left;
                    if(a==a_min) return pval;
                    a_min++;
            }
            else{
                    pval += p_right;
                    if(a==a_max) return pval;
                    a_max--;
            }
    }
    // If we get to this part of the code, it means is the mode of the distribution and, therefore, its associated p-value is 1
    return 1;
}




double cmh_minpval(longint *x, unsigned short K, longint* Nt, longint* nt,
        //precomputed quantities to save time when computing the CMH test statistic
        longint* Nt_nt, double* gammat, double* gammabint) {
    double left_tail_num = 0, right_tail_num = 0, den = 0;
    double aux1, aux2;
    for(unsigned short k=0; k<K; k++){
        aux1 = x[k]-Nt_nt[k]; aux2 = x[k]*gammat[k];
        left_tail_num += ((aux1 > 0) ? aux1 : 0) - aux2;
        right_tail_num += ((x[k] > nt[k]) ? nt[k] : x[k]) - aux2;
        den += x[k]*(1-((double)x[k])/Nt[k])*gammabint[k];
    }
    left_tail_num *= left_tail_num; right_tail_num *= right_tail_num;
    if(den==0) return 1;
    else return Chi2_sf(((left_tail_num > right_tail_num) ? left_tail_num : right_tail_num)/den,1);
}

double cmh_pval(longint a, longint *x, unsigned short K, longint* Nt,
        double* gammat, double* gammabint){
    long long k;
    double num = a, den = 0;
    for(k=0; k<K; k++){
        num -= x[k]*gammat[k];
        den += x[k]*(1-((double)x[k])/Nt[k])*gammabint[k];
    }
    num *= num;
    if(den==0) return 1;
    else return Chi2_sf(num/den,1);
}

double cmh_score(longint a, longint *x, unsigned short K, longint* Nt,
                double* gammat, double* gammabint){
    long long k;
    double num = a, den = 0;
    for(k=0; k<K; k++){
        num -= x[k]*gammat[k];
        den += x[k]*(1-((double)x[k])/Nt[k])*gammabint[k];
    }
    num *= num;
    if(den==0) return 0;
    else return num/den;
}

double cmh_odds_ratio(longint* at, longint* x, unsigned short K, longint* Nt, longint* nt)
{
    longint k;
    double num = 0, den = 0;
    for(k=0; k<K; k++){
        if(Nt[k]==0) continue;
        num += at[k] * ((Nt[k] + at[k]) - (x[k] + nt[k]));
        den += (nt[k] - at[k]) * (x[k] - at[k]);
    }
    if (den==0) return std::numeric_limits<double>::infinity();  //TODO: Test
    else return num/den;
}
