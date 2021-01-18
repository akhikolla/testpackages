#include "functions.h"

//function to calculate variance-covariance matrix from posterior samples
void calcPost(int i, int npars, 
                  arma::vec *tempmn, arma::mat *meanmat, 
                  arma::mat *meanmat1, NumericMatrix posterior, 
                  arma::mat *propcov)
{
    int j, k, m;

    //calculates variance-covariance matrix
    //first: means
    for(j = 0; j < npars; j++) {
        (*tempmn)[j] = 0;
        for(k = 0; k <= i; k++) {
            (*tempmn)[j] += posterior(k, j);
        }
        (*tempmn)[j] = (*tempmn)[j] / ((double) i + 1);
    }
    //matrix of product of means
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            (*meanmat)(j, k) = (*tempmn)[j] * (*tempmn)[k];
        }
    }
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            (*propcov)(j, k) = 0.0;
        }
    }
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            for(m = 0; m <= i; m++) {
                (*propcov)(j, k) += posterior(m, j) * posterior(m, k);
            }
            (*propcov)(j, k) -= (i + 1) * (*meanmat)(j, k);
            (*propcov)(j, k) = (*propcov)(j, k) / ((double) i);
        }
    }
//        Rprintf("Initial cov:\n");
//        for(k = 0; k < npars; k++)
//        {
//            for(j = 0; j < npars; j++) Rprintf("%f ", (*propcov)(j, k));
//            Rprintf("\n");
//        }
    return;
}

//function to update variance-covariance matrix based on current samples
void adaptUpdate(int i, int npars, 
                  arma::vec *tempmn, arma::mat *meanmat, 
                  arma::mat *meanmat1, NumericVector posterior, 
                  arma::mat *propcov)
{
    int j, k;

    //recursively update mean and covariance matrix
    for(j = 0; j < npars; j++) {
        (*tempmn)[j] = ((*tempmn)[j] * i + posterior[j]) / ((double) i + 1);
    }
    //new matrix of product of means
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            (*meanmat1)(j, k) = (*tempmn)[j] * (*tempmn)[k];
        }
    }
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            (*propcov)(j, k) = (((double) (i - 1)) / ((double) i)) * (*propcov)(j, k) + (1.0 / ((double) i)) * (i * (*meanmat)(j, k) - (i + 1) * (*meanmat1)(j, k) + posterior[j] * posterior[k]);
        }
    }
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            (*meanmat)(j, k) = (*meanmat1)(j, k);
        }
    }
    return;
}

