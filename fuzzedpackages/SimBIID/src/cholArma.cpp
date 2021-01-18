#include "functions.h"

// function to compute Cholesky decomposition and to deal with
// cases where it fails
arma::mat cholArma(arma::mat sigma, double *scale)
{
    bool success = false;
    int j = 0;
    arma::mat cholSigma (sigma.n_rows, sigma.n_cols);
    while(success == false && j < 20) {
        success = arma::chol(cholSigma, sigma);
        if(success == false) {
            sigma += arma::eye(sigma.n_rows, sigma.n_cols) * 1e-6;
        }
        j++;
    }
//    if(j > 1) {
//        Rprintf("chol %d\n", j);
//    }
    // return scale
    *scale = 1e-6 * (j - 1);
    if(success == false) {
        // // write to file as check
        // FILE *filemcmc;
        // char outfile [64] = "output_chol.txt";
        // //append current iteration to runinfo.txt
        // filemcmc = fopen(outfile, "w");
        // for(j = 0; j < sigma.n_rows; j++) {
        //     for(k = 0; k < sigma.n_cols; k++) {
        //         fprintf(filemcmc, "%f ", sigma(j, k));
        //     }
        //     fprintf(filemcmc, "\n");
        // }
        // fclose(filemcmc);
        stop("Error in Cholesky decomposition");
    }
    return cholSigma;
}

