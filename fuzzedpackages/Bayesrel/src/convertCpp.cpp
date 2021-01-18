//
//
//
//  Created by Julius Pfadt on 26.08.20.
//

#include <stdio.h>
#include <RcppArmadillo.h>
extern "C" {
#include "blockmat.h"
}


//[[Rcpp::depends(RcppArmadillo)]]

arma::ivec int_vector_csdp2RArma(int n, int *y)
{
    return arma::ivec(y, n+1);
}

arma::dvec double_vector_csdp2RArma(int n, double *y)
{
    return arma::dvec(y, n+1);
}

int * int_vector_R2csdpArma(int n, const arma::ivec& y)
{
  // return y.memptr(); // <- this should also work
  int *ret;
  int i;
  ret = (int *) malloc((n+1) * sizeof(int));
  if (ret == NULL)
    return NULL;
  for (i=1; i<=n; i++)
    ret[i] = y[i];
  return ret;
}



double * double_vector_R2csdpArma(int n, const arma::dvec& y)
{
  // return y.memptr(); // <- this should also work
  double *ret;
  int i;
  ret = (double *) malloc((n+1) * sizeof(double));
  if (ret == NULL)
    return NULL;
  for (i=1; i<=n; i++)
    ret[i] = y[i];
  return ret;
}


/*
 the input of this function is a kind of list but of class csdpBlkMat
 dont know if a Rcpp type list can also work
*/
struct blockmatrix blkmatrix_R2csdpArma(const Rcpp::List& X)
{
    struct blockmatrix ret;
    int blksize, blktype, allocsize, j, k;
    Rcpp::List cur_block;

    int nblocks = X["nblocks"];
    Rcpp::List blocks = X["blocks"];
    ret.nblocks = nblocks;
    ret.blocks = (struct blockrec *) malloc((nblocks + 1) * sizeof(struct blockrec));
    for (j=1; j<=nblocks; j++) {
      cur_block = blocks[j-1];
      blksize = cur_block["blocksize"];
      ret.blocks[j].blocksize = blksize;
      blktype = cur_block["blockcategory"];
      ret.blocks[j].blockcategory = (blktype == 1) ? MATRIX : DIAG;
      if (blktype == 1) {
        allocsize = blksize*blksize;
        ret.blocks[j].data.mat = (double *) malloc(allocsize * sizeof(double));
        arma::dvec dblvec = cur_block["data"];
        for (k=0; k<allocsize; k++)
            ret.blocks[j].data.mat[k] = dblvec(k);
      }
      else {
        ret.blocks[j].data.vec = double_vector_R2csdpArma(blksize, cur_block["data"]);
      }
    }
    return ret;
}


Rcpp::List blkmatrix_csdp2RArma(const blockmatrix& X)
{
    Rcpp::List ret;
    arma::dvec data;
    Rcpp::List blocks;
    Rcpp::List tmp = Rcpp::List::create(Rcpp::_["blocksize"], Rcpp::_["blockcategory"], Rcpp::_["data"]);

    int j,k, allocsize, blocksize, blockcategory;

    int nblocks = X.nblocks;
    ret.push_back(nblocks, "n.blocks");

    for (j=1; j<=X.nblocks; j++) {
        blocksize = X.blocks[j].blocksize;
        blockcategory = (X.blocks[j].blockcategory == MATRIX) ? 1 : 2;


        if (X.blocks[j].blockcategory == MATRIX) {
            allocsize = X.blocks[j].blocksize * X.blocks[j].blocksize;
            data.set_size(allocsize);
            for (k=0; k<allocsize; k++)
                data(k) = X.blocks[j].data.mat[k];
        } else {
            data = double_vector_csdp2RArma(X.blocks[j].blocksize, X.blocks[j].data.vec);
        }
        tmp["blocksize"] = blocksize;
        tmp["blockcategory"] = blockcategory;
        tmp["data"] = data;
        blocks.insert(j-1, tmp);

    }
    ret.push_back(blocks, "blocks");
    return ret;
}


struct constraintmatrix *constraints_R2csdpArma(const Rcpp::List& A)
{
    struct constraintmatrix *constraints;
    struct sparseblock *blockptr;

    int nblocks;
    Rcpp::List Ai, Aij;

    int i,j;

    int nconstraints = A.length();
    constraints = (struct constraintmatrix *) malloc((nconstraints + 1) * sizeof(struct constraintmatrix));

    for (i=1; i<=nconstraints; i++) {
        Ai = A[i-1];
        /*
         * Terminate block linked list with NULL
         */
        constraints[i].blocks = NULL;
        nblocks = Ai.length();

        for (j=nblocks; j>=1; j--) {
            Aij = Ai[j-1];
            /*
             * Allocate block data structure
             */
            blockptr = (struct sparseblock *) malloc(sizeof(struct sparseblock));
            /*
             * Initialize block data structure
             */
            blockptr->blocknum=Aij["blocknum"];
            blockptr->blocksize=Aij["blocksize"];
            blockptr->constraintnum=Aij["constraintnum"];
            blockptr->next=NULL;
            blockptr->nextbyblock=NULL;
            blockptr->numentries=Aij["numentries"];

            /*
             * Enter data
             */
            blockptr->iindices = int_vector_R2csdpArma(blockptr->numentries, Aij["iindices"]);

            blockptr->jindices = int_vector_R2csdpArma(blockptr->numentries, Aij["jindices"]);

            blockptr->entries = double_vector_R2csdpArma(blockptr->numentries, Aij["entries"]);

            /*
             * Insert block into linked list
             */
            blockptr->next=constraints[i].blocks;
            constraints[i].blocks=blockptr;
        }
    }
    return constraints;
}


// void printBlockMat(const blockmatrix& C) {
//     Rcpp::Rcout << "nblocks: " << C.nblocks << std::endl;
//     int j;
//     for (j=1; j<=C.nblocks; j++) {
//         Rcpp::Rcout << "blocksize: " << C.blocks[j].blocksize << std::endl;
//         Rcpp::Rcout << "blockcategory: " << C.blocks[j].blockcategory << std::endl;
//         Rcpp::Rcout << "blockdatarec.vec: " << *C.blocks[j].data.vec << std::endl;
//         Rcpp::Rcout << "blockdatarec.mat: " << *C.blocks[j].data.mat << std::endl;
//     }
//     return;
// }
//
//
// void printConstMat (const constraintmatrix& S, int k) {
//     int i;
//     for (i=1; i<=k; i++) {
//         Rcpp::Rcout << "stuff: " << S.blocks[i].blocksize << std::endl;
//     }
//     return;
// }
