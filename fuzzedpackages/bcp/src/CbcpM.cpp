/*
 bcp: an R package for performing a Bayesian analysis
 of change point problems.

 Copyright (C) 2013 Susan Wang and John W. Emerson

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, a copy is available at
 http://www.r-project.org/Licenses/

 -------------------
 FILE: CbcpM.cpp  */

/*  LIBRARIES  */

#include "bcp.h"


/* BEGIN BCP-M stuff */

// b is the number of blocks if position i is not a change point
double getprob(double B0, double B1, double W0, double W1,
                int b, Params& params) {
  double xmax1 = (B1*params.w[0]/W1)/(1+(B1*params.w[0]/W1));
  double xmax0 = (B0*params.w[0]/W0)/(1+(B0*params.w[0]/W0));
  double ratio, p;
  ratio = params.priors[b-1];
  // Rprintf("B0:%.10f, B1:%.10f, W0:%.10f, W1:%.10f, b:%d\n", B0, B1, W0, W1, b);
  if (b >= params.nn - 4 / params.kk) {
    p = 0;
  } else if (B0 == 0) {
    ratio *= exp(Rf_lbeta((double) (params.kk * (b+1) + 1) / 2,
                          (double) ((params.nn2 - b-1) * params.kk - 2) / 2))*
              Rf_pbeta(xmax1,
                         (double) (params.kk * (b+1) + 1) / 2,
                         (double) ((params.nn2 - b-1) * params.kk - 2) / 2, 1, 0);
    ratio *= exp((params.nn2*params.kk-1) * log(W0)/2 + log((b*params.kk+1)/2)
              - 0.5*((b*params.kk+1)*log(params.w[0])+
                  ((b+1)*params.kk+1) * log(B1) +
                  ((params.nn2-b-1)*params.kk-2)*log(W1)
              ));
    p = ratio/(1+ratio);
  } else {
    ratio *= exp(0.5*(
             ((params.nn2 - b) * params.kk - 2) * log(W0/W1) +
             (params.kk * b + 1) * log(B0/B1) +
             params.kk*log(W1/B1)
                   ));
    ratio *= exp(Rf_lbeta((double) (params.kk * (b+1) + 1) / 2,
                          (double) ((params.nn2 - b-1) * params.kk - 2) / 2))*
             Rf_pbeta(xmax1,
                      (double) (params.kk * (b+1) + 1) / 2,
                      (double) ((params.nn2 - b-1) * params.kk - 2) / 2, 1, 0);
    ratio /= (exp(Rf_lbeta((double) (params.kk * b + 1) / 2,
                           (double) ((params.nn2 - b) * params.kk - 2) / 2))*
              Rf_pbeta(xmax0,
                      (double) (params.kk * b + 1) / 2,
                      (double) ((params.nn2 - b) * params.kk - 2) / 2, 1, 0));
    p = ratio/(1+ratio);
  }
  // Rprintf("p:%0.2f\n", p);
  return p;
}

MCMCStepSeq pass(MCMCStepSeq &step, HelperVariables &helpers, Params &params)
{

  int i, j, cp;

  int bsize1, bsize2, bsize3;
  double tmp;
  DoubleVec bmean1(params.kk, 0);
  DoubleVec bmean2(params.kk, 0);
  DoubleVec bmean3(params.kk, 0);
  double bZ1, bZ2, bZ3;

  IntVec bvals(2);
  DoubleVec Wvals(2);
  DoubleVec Bvals(2);

  DoubleVec bmeanlast(params.kk);
  double bZlast = 0;

  // this is used to reference the current block in the MCMCStep we came in with
  int prevblock = 0;

  // this is used to reference the current MCMCStep we build from scratch
  MCMCStepSeq stepnew(step);
  int currblock = 0;

  // some other variables to denote stuff in current block and previous block
  // Note that "last" refers to the immediately left-adjacent block
  // whereas "prev" refers to the same block in the variable step
  double thisblockZ = step.bZ[0];
  int thisbend = step.bend[0];

  double lastblockZ = 0;
  int lastbend = -1; // this is simply a notational convenience

  // start the loop
  for (i = 0; i < params.nn - 1; i++) {
    if (i == step.bend[prevblock]) {
      // we're at an old change point, so we need to refresh "this" to be the
      // immediately following block
      lastblockZ = thisblockZ;
      prevblock++;

      thisbend = step.bend[prevblock];
      thisblockZ = step.bZ[prevblock];
    }
    /****
     * consider merging blocks if currently a change point
     */
    bvals[0] = stepnew.b;
    if (step.rho[i] == 0) { // not a change point at the moment
      Bvals[0] = stepnew.B;
      Wvals[0] = stepnew.W;
    } else { // it is a change point, so let's try removing the change point
      bvals[0]--;
      tmp = thisblockZ + lastblockZ;
      bZ3 = 0;
      if (lastbend > -1) {
        bsize3 = helpers.cumksize[thisbend] - helpers.cumksize[lastbend];
        for (j = 0; j < params.kk; j++) {
          bmean3[j] = (helpers.cumymat[j][thisbend] - helpers.cumymat[j][lastbend]) / bsize3;
          bZ3 += pow(bmean3[j], 2) * bsize3;
        }
      } else {
        bsize3 = helpers.cumksize[thisbend];
        for (j = 0; j < params.kk; j++) {
          bmean3[j] = helpers.cumymat[j][thisbend] / bsize3;
          bZ3 += pow(bmean3[j], 2) * bsize3;
        }
      }
      if (params.kk == 1 && bvals[0] == 1) Bvals[0] = 0; // force this to avoid rounding errs
      else
        Bvals[0] = stepnew.B - tmp + bZ3;
      Wvals[0] = stepnew.W + tmp - bZ3;
    }


    /****
     * consider breaking blocks if not a change point
     */
    bvals[1] = stepnew.b;
    if (step.rho[i] == 1) {
      Bvals[1] = stepnew.B;
      Wvals[1] = stepnew.W;
    } else {
      bZ1 = 0;
      bZ2 = 0;
      bvals[1]++;
      bsize2 = helpers.cumksize[thisbend] - helpers.cumksize[i];

      if (lastbend > -1)
        bsize1 = helpers.cumksize[i] - helpers.cumksize[lastbend];
      else
        bsize1 = helpers.cumksize[i];
      tmp = thisblockZ;
      for (j = 0; j < params.kk; j++) {
        bmean2[j] = (helpers.cumymat[j][thisbend] - helpers.cumymat[j][i]) / bsize2;
        if (lastbend > -1) {
          bmean1[j] = (helpers.cumymat[j][i] - helpers.cumymat[j][lastbend]) / bsize1;
        } else {
          bmean1[j] = helpers.cumymat[j][i] / bsize1;
        }
        bZ1 += pow(bmean1[j], 2) * bsize1;
        bZ2 += pow(bmean2[j], 2) * bsize2;
      }
      Bvals[1] = stepnew.B - tmp + bZ1 + bZ2;
      Wvals[1] = stepnew.W + tmp - bZ1 - bZ2;
    }

    // if (i == 4122) return(stepnew);
    double p = getprob(Bvals[0], Bvals[1], Wvals[0], Wvals[1], bvals[0], params);
    // do the sampling and then updates
    double myrand = Rf_runif(0.0, 1.0);

    if (myrand < p) {
      cp = 1;
    } else {
      cp = 0;
    }
    // Rprintf("i:%d  p=%0.4f, myrand=%0.2f, cp=%d\n", i, p, myrand, cp);

    stepnew.B = Bvals[cp];
    stepnew.W = Wvals[cp];
    stepnew.b = bvals[cp];

    if (cp != step.rho[i]) { // we modified the change point status
      if (cp == 0) {
        // removed a change point
        // update last block's stuff since the last block is now further back
        thisblockZ = bZ3;
        if (currblock > 0) {
          lastbend = stepnew.bend[currblock - 1];
          lastblockZ = stepnew.bZ[currblock - 1];
        } else {
          lastblockZ = 0;
          lastbend = -1; // this is simply a notational convenience
        }
      } else { // added a change point
        thisblockZ = bZ2;
        lastblockZ = bZ1;
      }
    }
    stepnew.rho.push_back(cp);

    if (stepnew.rho[i] == 1) {
      if (step.rho[i] == 1) { // never calculated these quantities yet; do it now
        lastblockZ = 0;

        if (lastbend > -1)
          bsize1 = helpers.cumksize[i] - helpers.cumksize[lastbend];
        else
          bsize1 = helpers.cumksize[i];
        for (j = 0; j < params.kk; j++) {
          if (lastbend > -1) {
            bmean1[j] = (helpers.cumymat[j][i] - helpers.cumymat[j][lastbend]) / bsize1;
          } else {
            bmean1[j] = helpers.cumymat[j][i] / bsize1;
          }
          lastblockZ += pow(bmean1[j], 2) * bsize1;
        }
      }
      // we've added a change point, so we want to record some stuff
      stepnew.bsize.push_back(bsize1);
      stepnew.bend.push_back(i);
      stepnew.bmean.push_back(bmean1);
      stepnew.bZ.push_back(lastblockZ);
      currblock++;
      lastbend = i;

    }
  }
  // done with a full pass, now let's add info on the final block
  if (lastbend > -1)
    stepnew.bsize.push_back(params.nn2 - helpers.cumksize[lastbend]);
  else
    stepnew.bsize.push_back(params.nn2);
  for (j = 0; j < params.kk; j++) {
    if (lastbend > -1) {
      bmeanlast[j] = (helpers.cumymat[j][params.nn - 1] - helpers.cumymat[j][lastbend]) / stepnew.bsize[currblock];
    } else {
      bmeanlast[j] = helpers.cumymat[j][params.nn - 1] / params.nn2;
    }
    bZlast += pow(bmeanlast[j], 2) * stepnew.bsize[currblock];
  }
  stepnew.bmean.push_back(bmeanlast);
  stepnew.bZ.push_back(bZlast);
  stepnew.bend.push_back(params.nn - 1);
  stepnew.rho.push_back(1);
  return stepnew;
}

// [[Rcpp::export]]
SEXP rcpp_bcpM(SEXP pdata, SEXP pid, SEXP pmcmcreturn, SEXP pburnin, SEXP pmcmc,
                         SEXP pa, SEXP pw)
{

  NumericMatrix data(pdata);
  int mcmcreturn = INTEGER_DATA(pmcmcreturn)[0];
  int burnin = INTEGER_DATA(pburnin)[0];
  int mcmc = INTEGER_DATA(pmcmc)[0];

  // INITIALIZATION OF LOCAL VARIABLES
  int i, j, m, k;
  double wstar, xmax;

  // INITIALIZATION OF OTHER OBJECTS
  HelperVariables helpers(data, pid);
  Params params(pw, helpers.cumksize.size(), data.nrow(), pa, false, false,
                0, 0, data.ncol());
  //params.print();
  //helpers.print();
  int MM = burnin + mcmc;

  //helpers.print();
  //params.print();

  MCMCStepSeq step(helpers, params);

  int MM2, nn2;
  if (mcmcreturn == 0) {
    MM2 = 1;
    nn2 = 1;
  } else {
    nn2 = params.nn;
    MM2 = MM;
  }
  // Things to be returned to R:
  NumericMatrix pmean(params.nn, params.kk);
  NumericMatrix ss(params.nn, params.kk);
  NumericMatrix pvar(params.nn, params.kk);
  NumericVector pchange(params.nn);
  NumericVector blocks(burnin + mcmc);
  NumericMatrix rhos(nn2, MM2);
  // NumericVector liks(MM2);
  NumericMatrix results(nn2*MM2,params.kk);

  double tmpMean;

  // Rprintf("starting\n");
  GetRNGstate(); // Consider Dirk's comment on this.
  // step.print();
  for (i = 0; i < params.nn; i++) {
    pchange[i] = 0;
    for (j = 0; j < params.kk; j++) {
      pmean(i, j) = 0;
    }
  }
  for (m = 0; m < MM; m++) {
    // Rprintf("Step %d -- ", m);
    step = pass(step, helpers, params);
    // Rprintf("blocks:%d, B:%0.2f\n", step.b, step.B);
    blocks[m] = step.b;
    if (m >= burnin || mcmcreturn == 1) {
      // compute posteriors
      if (step.B == 0) {
        wstar = params.w[0] * (step.b*params.kk + 1) / (step.b * params.kk +3);
      } else {

        xmax = step.B * params.w[0] / step.W / (1 + step.B * params.w[0] / step.W);
        // Rprintf("xmax:%0.2f\n", xmax);
        // wstar = log(step.W) - log(step.B)
        //   + Rf_lbeta((double) (step.b* params.kk + 3) / 2, (double) ((params.nn2 - step.b)*params.kk - 4) / 2)
        //   + Rf_pbeta(xmax, (double) (step.b*params.kk + 3) / 2, (double) ((params.nn2  - step.b)*params.kk - 4) / 2, 1, 1)
        //   - Rf_lbeta((double) (step.b*params.kk + 1) / 2, (double) ((params.nn2  - step.b)*params.kk - 2) / 2)
        //   - Rf_pbeta(xmax, (double) (step.b * params.kk+ 1) / 2, (double) ((params.nn2  - step.b)*params.kk - 2) / 2, 1, 1);
        // wstar = exp(wstar);
        wstar = (step.W/step.B)*
          Rf_beta((double) (step.b* params.kk + 3) / 2, (double) ((params.nn2 - step.b)*params.kk - 4) / 2) *
          Rf_pbeta(xmax, (double) (step.b*params.kk + 3) / 2, (double) ((params.nn2  - step.b)*params.kk - 4) / 2, 1, 0) /
          Rf_beta((double) (step.b*params.kk + 1) / 2, (double) ((params.nn2  - step.b)*params.kk - 2) / 2) /
          Rf_pbeta(xmax, (double) (step.b * params.kk+ 1) / 2, (double) ((params.nn2  - step.b)*params.kk - 2) / 2, 1, 0);
        // Rprintf("wstar:%0.2f\n", wstar);

      }
      // for posterior estimate of overall noise variance
      // if (m >= burnin)
        // pvar += (step.W + wstar*step.B)/(params.nn2 * params.kk-3);
      k = 0;
      for (j = 0; j < params.nn; j++) {
        // Rprintf("j:%d out of %d (%d, %d)  | ", j, params.nn, pchange.size(), step.rho.size());
        // Rprintf("pchange[%d]: %0.2f, step.rho:%d\n", j, pchange[j], step.rho[j]);
        if (m >= burnin)
          pchange[j] += (double) step.rho[j];
        for (i = 0; i < params.kk; i++) {
          tmpMean = step.bmean[k][i] * (1 - wstar) + helpers.ybar * wstar;
          // Rprintf("i:%d -- tmpMean:%0.2f, wstar:%0.2f, bmean:%0.2f, ybar:%0.2f\n",
                  // i, tmpMean, wstar, step.bmean[k][i], helpers.ybar);
          if (m >= burnin) {
            pmean(j, i) += tmpMean;
            ss(j, i) += tmpMean * tmpMean;
            // Rprintf("pmean:%0.2f, ss:%0.2f\n", pmean(j,i), ss(j,i));
          }
          if (mcmcreturn == 1)
            results(m*params.nn+j, i) = tmpMean;
        }

        if (mcmcreturn == 1)
          rhos(j, m) = step.rho[j];
        if (step.rho[j] == 1) k++;
      }
    }
  }
  // Rprintf("post processing\n");
  // step.print();
  // post processing
  for (j = 0; j < params.nn; j++) {
    pchange[j] /= mcmc;
    for (i = 0; i < params.kk; i++) {
      pmean(j, i) /= mcmc;
      pvar(j, i) = (ss(j, i) / mcmc - pmean(j,i)*pmean(j,i))*(mcmc/(mcmc-1));
    }
  }
  // Rprintf("ending\n");

  PutRNGstate();

  List z;
  z["posterior.mean"] = pmean;
  z["posterior.var"] = pvar;
  z["posterior.prob"] = pchange;
  z["blocks"] = blocks;
  z["mcmc.rhos"] = rhos;
  z["mcmc.means"] = results;
  // z["lik"] = liks;
  return z;

} /* END MAIN  */

/* END BCP-M STUFF*/
