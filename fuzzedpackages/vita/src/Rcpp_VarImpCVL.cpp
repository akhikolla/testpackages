#include <Rcpp.h>
using namespace Rcpp;

/*================================================================================================*/

/* Node status */
#define NODE_TERMINAL -1

NumericVector predictRegTreeSL (NumericMatrix x, int nsample, int mdim,
                                IntegerVector lDaughter, IntegerVector rDaughter,
                                IntegerVector nodestatus, NumericVector split,
                                NumericVector nodepred, IntegerVector splitVar,
                                int treeSize, IntegerVector cat, int maxcat);

IntegerVector predictClassTreeSL (NumericMatrix x, int nsample, int mdim,
                                  IntegerVector lDaughter, IntegerVector rDaughter,
                                  IntegerVector nodestatus, NumericVector xbestsplit,
                                  IntegerVector nodeclass,IntegerVector bestvar,
                                  int treeSize,IntegerVector cat, int maxcat);

NumericMatrix permuteSL (int m, NumericMatrix x,  int nsample, int mdim);
/*================================================================================================*/
/*                                       Regression                                               */
/*================================================================================================*/

// [[Rcpp::export]]
List Rcpp_VarImpCVLReg (NumericMatrix x_l, NumericVector y_l,
              		      int nsample, int mdim,int nTree, int nPerm,
                        IntegerMatrix lDaughter, IntegerMatrix rDaughter,
                        IntegerMatrix nodestatus, NumericMatrix split,
                        NumericMatrix nodepred,  IntegerMatrix splitVar,
                        IntegerVector ndbigtree, IntegerVector cat, int maxcat) {

/*************************************************************************
* Compute fold-specific variabele iportance for the l-th data set
* Argument:
*   x_l,y_l ... observations from the l-th data set
*   mdim ... number of variables in data set
*   nsample ...  number of cases in l-th data set
*
*   nPerm ... Number of times the l-th data are permuted per tree for
*           assessing  fold specific variable importance.
*   cat ... integer vector of number of classes in the predictor;
*   maxcat ... maximum of cat
*
*   nTree ... number of trees in the l-th Forest
*
*   lDaughter,rDaughter, nodestatus,
*   split, nodepred, splitVar,ndbigtree ... l-th Forest based on observations that
*                                           are not part of the l-th set.
*************************************************************************/


  double  res_sl, slerr, slerrperm,delta,r;

  NumericVector ytr(nsample), xtmp(nsample), errimp(mdim);


  GetRNGstate();
  /*************************************
   * Start the loop over trees.
   *************************************/
  for (int j = 0; j < nTree; j++) {

            /***********************************************************************************/
            /* predict the l-th data set with the current tree */
          	/* ytr is the prediction on l-th data set by the current tree */
            ytr=predictRegTreeSL(x_l, nsample,  mdim, lDaughter(_,j), rDaughter(_,j),
                                nodestatus(_,j),  split(_,j), nodepred(_,j),  splitVar(_,j),
                                ndbigtree(j), cat, maxcat);
            /***********************************************************************************/


        		slerr = 0.0;
            res_sl = 0.0;

        		for (int n = 0; n < nsample; n++) {
              				res_sl = ytr[n] - y_l[n];
                      slerr += res_sl * res_sl;
        		}

            /***********************************************************************************/
            /*  fold specific permutation Variable importance */

		        for (int mr = 0; mr < mdim; mr++) {

                  /* make a copy of the m-th variable into xtmp */
                  xtmp =x_l(mr,_);

                  slerrperm = 0.0;

                  for (int k = 0; k < nPerm; ++k) {
                      permuteSL(mr,x_l, nsample, mdim);

                      ytr=predictRegTreeSL(x_l, nsample,  mdim, lDaughter(_,j), rDaughter(_,j),
                                        nodestatus(_,j),  split(_,j), nodepred(_,j),  splitVar(_,j),
                                        ndbigtree(j), cat, maxcat);

                      for (int n = 0; n < nsample; ++n) {
                              r = ytr[n] -y_l[n];
                              slerrperm += r * r;
                      }
                  }
                  delta = (slerrperm / nPerm - slerr) / nsample;
                  errimp[mr] += delta;
                  /* copy original data back */
                  x_l(mr,_)= xtmp;

          }// end of permutation importance


  } //end of tree iterations

  PutRNGstate();

  for (int m = 0; m < mdim; ++m) {
  		errimp[m] = errimp[m] / nTree;
  }

  return List ::create( Named ("fold_importance") =errimp);


}// end of Rcpp_VarImpCVLReg


/*================================================================================================*/
/*                                     Classification                                             */
/*================================================================================================*/

// [[Rcpp::export]]
List  Rcpp_VarImpCVLCL (NumericMatrix x_l, IntegerVector y_l,
                        int nsample, int mdim, int nTree, int nclass,
                        IntegerMatrix lDaughter, IntegerMatrix rDaughter,
                        IntegerMatrix nodestatus, NumericMatrix xbestsplit,
                        IntegerMatrix nodeclass, IntegerMatrix bestvar,
                        IntegerVector ndbigtree, IntegerVector cat, int maxcat) {

/*************************************************************************
 * Compute fold-specific variabele iportance for the l-th data set
 * Argument:
 *   x_l,y_l ... observations from the l-th data set
 *   mdim ... number of variables in data set
 *   nsample ...  number of cases in l-th data set
 *
 *   nPerm ... Number of times the l-th data are permuted per tree for
 *           assessing  fold specific variable importance.
 *   cat ... integer vector of number of classes in the predictor;
 *   maxcat ... maximum of cat
 *
 *   nTree ... number of trees in the l-th Forest
 *
 *   lDaughter,rDaughter, nodestatus,
 *   xbestsplit, nodepred, bestvar,ndbigtree ... l-th Forest based on observations that
 *                                               are not part of the l-th set.
 *************************************************************************/

  int nrightall, nrightimpall;
  double delta=0.0;
  IntegerVector yts(nsample), yvs(nsample);
  NumericVector xtmp(nsample), imprt( mdim);


  GetRNGstate();
  /*************************************
  * Start the loop over trees.
  *************************************/
  for (int j = 0; j < nTree; j++) {
    /***********************************************************************************/
    /* predict the OOB data with the current tree */
    /* yts is the prediction  on l-th data set by the current tree */
    yts = predictClassTreeSL (x_l, nsample, mdim,
                              lDaughter(_,j),rDaughter(_,j),
                              nodestatus(_,j),xbestsplit(_,j),
                              nodeclass(_,j),bestvar(_,j),
                              ndbigtree(j),cat,maxcat);

    /***********************************************************************************/
    /* permutation Variable importance */

    /* Count the number of correct prediction by the current tree
    among the l-th data set. */
    nrightall = 0;

    for (int n = 0; n < nsample; n++) {
      /*l-th data set and predicted correctly: */
      if (yts[n] == y_l[n]) {
        nrightall++;
      }
    }

    for (int m = 0; m < mdim; m++) {
      nrightimpall = 0;
      /* make a copy of the m-th variable into xtmp */
      xtmp = x_l(m,_);
      /* Permute the m-th variable. */
      permuteSL(m, x_l, nsample, mdim);
      /* Predict the modified data using the current tree. */
      yvs = predictClassTreeSL (x_l, nsample, mdim,
                                lDaughter(_,j),rDaughter(_,j),
                                nodestatus(_,j),xbestsplit(_,j),
                                nodeclass(_,j),bestvar(_,j),
                                ndbigtree(j),cat,maxcat);
      /* copy original data back */
      x_l(m,_)= xtmp;

      /* Count how often correct predictions are made with
      the modified data. */
      for (int n = 0; n < nsample; n++) {
          if (yvs[n] == y_l[n])  {
            nrightimpall++;
          }
      }
      /* Accumulate decrease in proportions of correct
      predictions. */
      /* overall measure, across all classes: */
      delta = ((double)(nrightall - nrightimpall)) / nsample;
      imprt(m) += delta;

    }//end of mdim iterations

  }//end of tree iterations

  PutRNGstate();

  for (int m = 0; m < mdim; m++) {
    /* overall measures */
    imprt(m) =  imprt(m)/ nTree;
  }
  return List ::create( Named ("fold_importance") =imprt);

}
/*================================================================================================*/
/**************************************************************************************************/
/*================================================================================================*/

NumericVector predictRegTreeSL(NumericMatrix x, int nsample, int mdim,
                               IntegerVector lDaughter, IntegerVector rDaughter,
                               IntegerVector nodestatus, NumericVector split,
                               NumericVector nodepred, IntegerVector splitVar,
                               int treeSize, IntegerVector cat, int maxcat) {
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
IntegerVector predictClassTreeSL( NumericMatrix x, int nsample, int mdim,
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

NumericMatrix permuteSL(int m, NumericMatrix x,  int nsample, int mdim) {
  /*************************************************************************
   * Permute the l-th data set part of a variable in x.
   * Argument:
   *   m: the variable to be permuted
   *   x: the data matrix (variables in rows)
   *   nsample: number of cases in the data
   *   mdim: number of variables in the data
   *************************************************************************/

  double  tmp;
  int i, last, k;
  NumericVector tp(nsample);


  /* make a copy of the l-th data into tp (for permuting) */

  tp = x(m,_);

  /* Permute tp */
  last = nsample;
  for (i = 0; i < nsample; ++i) {
    k = (int) last * unif_rand();
    tmp = tp[last - 1];
    tp[last - 1] = tp[k];
    tp[k] = tmp;
    last--;
  }

  /* Copy the permuted   l-th data back into x. */
  x(m,_)= tp;
  return x;
}
/*================================================================================================*/
