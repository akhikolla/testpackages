#include <Rcpp.h>
using namespace Rcpp;
/*================================================================================================*/

/* Node status */
#define NODE_TERMINAL -1

IntegerVector predictClassTreeOOB (NumericMatrix x, int nsample, int mdim,
                                IntegerVector lDaughter, IntegerVector rDaughter,
                                IntegerVector nodestatus, NumericVector xbestsplit,
                                IntegerVector nodeclass,IntegerVector bestvar,
                                int treeSize,IntegerVector cat, int maxcat);

NumericVector predictRegTreeOOB (NumericMatrix x, int nsample, int mdim,
                                 IntegerVector lDaughter, IntegerVector rDaughter,
                                 IntegerVector nodestatus, NumericVector split,
                                 NumericVector nodepred, IntegerVector splitVar,
                                 int treeSize,IntegerVector cat, int maxcat);


void permuteOOB(int m, NumericMatrix x, IntegerVector inv, int nsample, int mdim);
/*================================================================================================*/
/*                                       Regression                                               */
/*================================================================================================*/

// [[Rcpp::export]]
List Rcpp_compVarImpReg (NumericMatrix x, NumericVector y,
                		     int nsample, int mdim,int nTree, int nPerm,
                         IntegerMatrix lDaughter, IntegerMatrix rDaughter,
                         IntegerMatrix nodestatus, NumericMatrix split,
                         NumericMatrix nodepred,  IntegerMatrix splitVar,
                         IntegerMatrix inm, IntegerVector ndbigtree,
                         IntegerVector cat, int maxcat) {


  double  ooberr, ooberrperm,delta,resOOB,r;
  int  nOOB;
  NumericVector ytr(nsample), xtmp(nsample),
                errimp(mdim), impSD(mdim);

  IntegerVector inv(nTree);

  GetRNGstate();
  /*************************************
   * Start the loop over trees.
   *************************************/
  for (int j = 0; j < nTree; j++) {

            /***********************************************************************************/
            /* predict the OOB data with the current tree */
          	/* ytr is the prediction on OOB data by the current tree */
            ytr=predictRegTreeOOB( x, nsample,  mdim, lDaughter(_,j), rDaughter(_,j),
                                   nodestatus(_,j),  split(_,j), nodepred(_,j),  splitVar(_,j),
                                   ndbigtree(j), cat, maxcat);

        		ooberr = 0.0;
        		nOOB = 0; /* nOOB is the number of OOB samples for this tree */
            inv=inm(_,j); /* which samples are “in-bag” in which trees  */

        		for (int n = 0; n < nsample; n++) {
            			if (inv[n] == 0) {
                      nOOB++;
              				resOOB = ytr[n] - y[n];
                      ooberr += resOOB * resOOB;
            			}

        		}

            /***********************************************************************************/
            /* permutation Variable importance */

		        for (int mr = 0; mr < mdim; mr++) {

                  /* make a copy of the m-th variable into xtmp */
                  xtmp = x(mr,_);
                  ooberrperm = 0.0;

                  for (int k = 0; k < nPerm; ++k) {
                      permuteOOB(mr, x, inv, nsample, mdim);
                      ytr=predictRegTreeOOB( x, nsample,  mdim, lDaughter(_,j), rDaughter(_,j),
                                        nodestatus(_,j),  split(_,j), nodepred(_,j),  splitVar(_,j),
                                        ndbigtree(j), cat, maxcat);

                      for (int n = 0; n < nsample; ++n) {
                          if (inv[n] == 0) {
                              r = ytr[n] - y[n];
                              ooberrperm += r * r;

                          }
                      }
                  }
                  delta = (ooberrperm / nPerm - ooberr) / nOOB;
                  errimp[mr] += delta;
                  impSD[mr] += delta * delta;
                  /* copy original data back */
                   x(mr,_)= xtmp;

          }// end of permutation importance


  } //end of tree iterations

  PutRNGstate();

  for (int m = 0; m < mdim; ++m) {
  		errimp[m] = errimp[m] / nTree;
		  impSD[m] = sqrt( ((impSD[m] / nTree) - (errimp[m] * errimp[m])) / nTree );
  }

  return List ::create( Named ("importance") = errimp,
                        Named ("importanceSD") = impSD );


}// end of compVarImp

/*================================================================================================*/
/*                                     Classification                                             */
/*================================================================================================*/

// [[Rcpp::export]]
List  Rcpp_compVarImpCL (NumericMatrix x, IntegerVector y,
                         int nsample, int mdim,
                         int nTree,int nclass,
                         IntegerMatrix lDaughter,
                         IntegerMatrix rDaughter,
                         IntegerMatrix nodestatus,
                         NumericMatrix xbestsplit,
                         IntegerMatrix nodeclass,
                         IntegerMatrix bestvar,
                         IntegerMatrix inbag,
                         IntegerVector ndbigtree,
                         IntegerVector cat,
                         int maxcat) {


  int noutall,nrightall,nrightimpall;
  double av=0.0, delta=0.0;
  IntegerVector yts(nsample),yvs(nsample), jin(nsample);

  NumericVector  xtmp(nsample);
  NumericMatrix  imprt( mdim,(nclass+1)),impsd( mdim,(nclass+1));


  GetRNGstate();
  /*************************************
  * Start the loop over trees.
  *************************************/
  for (int j = 0; j < nTree; j++) {
    /***********************************************************************************/
    /* predict the OOB data with the current tree */
    /* yts is the prediction on OOB data by the current tree */
    yts = predictClassTreeOOB( x, nsample, mdim,
                            lDaughter(_,j),rDaughter(_,j),
                            nodestatus(_,j),xbestsplit(_,j),
                            nodeclass(_,j),bestvar(_,j),
                            ndbigtree(j),cat,maxcat);

    /***********************************************************************************/
    /* count number of OOB cases in the current iteration.
    nout[n] is the number of OOB cases for the n-th class.
    noutall is the number of OOB cases overall. */

    jin=inbag(_,j); /* which samples are “in-bag” in which trees  */
    IntegerVector nout(nclass);
    noutall=0;
    for (int n = 0; n < nsample; n++) {
      if (jin[n] == 0) {
        nout[y[n] - 1]++;
        noutall++;
      }
    }//end of nsample iterations


    /***********************************************************************************/
    /* permutation Variable importance */

    /* Count the number of correct prediction by the current tree
    among the OOB samples, by class. */
    nrightall = 0;
    IntegerVector nright(nclass);
    for (int n = 0; n < nsample; n++) {
      /* out-of-bag and predicted correctly: */
      if (jin[n] == 0 && yts[n] == y[n]) {
        nright[y[n] - 1]++;
        nrightall++;
      }
    }

    for (int m = 0; m < mdim; m++) {
      nrightimpall = 0;
      IntegerVector nrightimp(nclass);
      /* make a copy of the m-th variable into xtmp */
      xtmp = x(m,_);
      /* Permute the m-th variable. */
      permuteOOB(m, x, jin, nsample, mdim);
      /* Predict the modified data using the current tree. */
      yvs = predictClassTreeOOB( x, nsample, mdim,
                              lDaughter(_,j),rDaughter(_,j),
                              nodestatus(_,j),xbestsplit(_,j),
                              nodeclass(_,j),bestvar(_,j),
                              ndbigtree(j),cat,maxcat);
      /* copy original data back */
      x(m,_)= xtmp;

      /* Count how often correct predictions are made with
      the modified data. */
      for (int n = 0; n < nsample; n++) {
        if (jin[n] == 0 ){
          if (yvs[n] == y[n])  {
            nrightimp[y[n] - 1]++;
            nrightimpall++;
          }
        }
      }
      /* Accumulate decrease in proportions of correct
      predictions. */
      /* class-specific measures first: */
      for (int n = 0; n < nclass; n++) {
        if (nout[n] > 0) {
          delta = ((double) (nright[n] - nrightimp[n])) / nout[n];
          imprt(m, n) += delta;
          impsd(m, n)+= delta * delta;
        }
      }
      /* overall measure, across all classes: */
      if (noutall > 0) {
        delta = ((double)(nrightall - nrightimpall)) / noutall;
        imprt(m,nclass) += delta;
        impsd(m,nclass) += delta * delta;
      }

    }//end of mdim iterations

  }//end of tree iterations


  PutRNGstate();

  for (int m = 0; m < mdim; m++) {
    /* class-specific measures */
    for (int k = 0; k < nclass; ++k) {
      av = imprt(m,k) / nTree;
      impsd(m,k) = sqrt(((impsd(m,k)/ nTree) - av*av) / nTree);
      imprt(m,k) = av;
    }
    /* overall measures */
    av = imprt(m,nclass)/ nTree;
    impsd(m,nclass) = sqrt(((impsd(m,nclass)/ nTree) - av*av) / nTree);
    imprt(m,nclass) = av;
  }
  return List ::create( Named ("importance") =imprt,
                        Named ("importanceSD") =impsd );
}
/*================================================================================================*/
/**************************************************************************************************/
/*================================================================================================*/
NumericVector predictRegTreeOOB(NumericMatrix x, int nsample, int mdim,
                                IntegerVector lDaughter, IntegerVector rDaughter,
                                IntegerVector nodestatus, NumericVector split,
                                NumericVector nodepred, IntegerVector splitVar,
                                int treeSize,IntegerVector cat, int maxcat) {

  int  k, m;
  double dpack;
  NumericVector ypred(nsample);
  IntegerVector cbestsplit(maxcat * treeSize);
  /* decode the categorical splits */
  if (maxcat > 1) {
    for (int i = 0; i < treeSize; i++) {
      if (nodestatus[i] != NODE_TERMINAL && cat[splitVar[i] - 1] > 1) {
          dpack = split[i];
          /* unpack `npack' into bits */
          /* unpack(dpack, maxcat, cbestsplit + i * maxcat); */
          for (int j = 0; j < cat[splitVar[i] - 1]; ++j) {
            cbestsplit[j + i*maxcat] = ((unsigned long) dpack & 1) ? 1 : 0;
            dpack = dpack / 2.0 ;
            /* cbestsplit[j + i*maxcat] = npack & 1; */
          }
      }
    }
  }


  for (int i = 0; i < nsample; ++i) {
    k = 0;
    while (nodestatus[k] != NODE_TERMINAL) { /* go down the tree */
      m = splitVar[k] - 1;
      if (cat[m] == 1) {
         /* Split by a numerical predictor */
          k = (x(m , i) <= split[k]) ?lDaughter[k] - 1 : rDaughter[k] - 1;
      }else {
          /* Split by a categorical predictor */
          k = cbestsplit[(int) x(m , i) - 1 + k * maxcat] ? lDaughter[k] - 1 : rDaughter[k] - 1;
      }
    }
    /* terminal node: assign prediction and move on to next */
    ypred[i] = nodepred[k];
  }
  return ypred;

}
/*================================================================================================*/
IntegerVector predictClassTreeOOB( NumericMatrix x, int nsample, int mdim,
                                IntegerVector lDaughter, IntegerVector rDaughter,
                                IntegerVector nodestatus, NumericVector xbestsplit,
                                IntegerVector nodeclass,IntegerVector bestvar,
                                int treeSize,IntegerVector cat, int maxcat) {
  int m, k;
  IntegerVector yp(nsample),cbestsplit(maxcat * treeSize);
  double dpack;

  /* decode the categorical splits */
  if (maxcat > 1) {
    for (int i = 0; i < treeSize; i++) {
      if (nodestatus[i] != NODE_TERMINAL) {
        if (cat[bestvar[i] - 1] > 1) {
          dpack = xbestsplit[i];
          /* unpack `dpack' into bits */
          /* unpack(dpack, maxcat, cbestsplit + i * maxcat); */
          for (int j = 0; j < cat[bestvar[i] - 1]; j++) {
            cbestsplit[j + i*maxcat] = ((unsigned long) dpack & 1) ? 1 : 0;
            dpack = dpack / 2;
          }
        }
      }
    }
  }

  for (int i = 0; i < nsample; i++) {
    k = 0;
    while (nodestatus[k] != (NODE_TERMINAL)) {
      m = bestvar[k] - 1;
      if (cat[m] == 1) {
        /* Split by a numerical predictor */
        k = (x(m , i)  <= xbestsplit[k]) ? lDaughter[k] - 1 : rDaughter[k]- 1;
      } else {
        /* Split by a categorical predictor */
        k = cbestsplit[(int) x(m,i) - 1 + k * maxcat] ? lDaughter[k] - 1 : rDaughter[k]- 1;
      }
    }
    /* Terminal node: assign class label */
    yp[i] = nodeclass[k];
  }
  return yp;

}


/*================================================================================================*/

void permuteOOB(int m, NumericMatrix x, IntegerVector inv, int nsample, int mdim) {
  /* Permute the OOB part of a variable in x.
  * Argument:
  *   m: the variable to be permuted
  *   x: the data matrix (variables in rows)
  *   inv: vector indicating which case is OOB
  *   nsample: number of cases in the data
  *   mdim: number of variables in the data
  */
  double  tmp;
  int i, last, k, nOOB = 0;
  NumericVector tp(nsample);

  for (i = 0; i < nsample; i++) {
    /* make a copy of the OOB part of the data into tp (for permuting) */
    if (inv[i] == 0) {
      tp[nOOB] = x(m, i);
      nOOB++;
    }
  }
  /* Permute tp */
  last = nOOB;
  for (i = 0; i < nOOB; i++) {
    k = (int) last * unif_rand();
    tmp = tp[last - 1];
    tp[last - 1] = tp[k];
    tp[k] = tmp;
    last--;
  }

  /* Copy the permuted OOB data back into x. */
  nOOB = 0;
  for (i = 0; i < nsample; i++) {
    if (inv[i] == 0) {
      x(m, i)= tp[nOOB];
      nOOB++;
    }
  }

}
