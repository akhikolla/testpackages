

#include <RcppArmadillo.h>
extern "C" {
#include "declarations.h"
}
#include "customsdp.h"


//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::dvec csdpArma(
              int n_p,
              int nconstraints_p,
              int nblocks_p,
              const arma::ivec& blocktypes_p,
              const arma::ivec& blocksizes_p,
              const Rcpp::List& C_p,
              const Rcpp::List& A_p,
              const arma::dvec& b_p,
              const arma::cube& car,
              const int printlevel = 0)
{

    struct blockmatrix C;
    struct constraintmatrix *constraints;
    double *b;
    double pobj, dobj;
    int status;
    arma::dvec out(car.n_rows);


    /*
     * setup C
     */
    C = blkmatrix_R2csdpArma(C_p);

    /*
     * setup constraints
     */
    constraints = constraints_R2csdpArma(A_p);

    /*
     * Allocate storage for RHS
     */
    b = double_vector_R2csdpArma(nconstraints_p,b_p);

    /*
     * Solve the problem
     */
    status = custom_sdpCpp(n_p,nconstraints_p,C,b,constraints,0.0,&pobj,&dobj, car, out, printlevel);


    // free_prob(n_p,nconstraints_p,C,b,constraints,X,y,Z);
    /*
     * freeing up the memory, copied from free_prob
     */
    free(b);
    free_mat(C);

    int i;
    struct sparseblock *ptr;
    struct sparseblock *oldptr;
    if (constraints != NULL)
    {
        for (i=1; i<=nconstraints_p; i++)
        {
            /*
             * Get rid of constraint i.
             */

            ptr=constraints[i].blocks;
            while (ptr != NULL)
            {
                free(ptr->entries);
                free(ptr->iindices);
                free(ptr->jindices);
                oldptr=ptr;
                ptr=ptr->next;
                free(oldptr);
            };
        };
        /*
         * Finally, free the constraints array.
         */

        free(constraints);
    };

    return out;
}
