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
FILE: Cbcp.cpp  */

#include "bcp.h"

int sampleFromLikelihoods(DoubleVec &likvals, double maxlik) {
  int i;
  int k = likvals.size();
  DoubleVec cumlikvals(k);
  cumlikvals[0] = exp(likvals[0] - maxlik);
  for (i = 1; i < k; i++) {
    cumlikvals[i] = cumlikvals[i - 1] + exp(likvals[i] - maxlik);
  }
  double p = Rf_runif(0.0, 1.0);
  for (i = 0; i < k; i++) {
    if (p < cumlikvals[i] / cumlikvals[k - 1])
      return i;
  }
  return -1; // this won't be triggered anyway
}

DoubleVec matrixCalcs(HelperVariables& helpers,
                      Params& params,
                      DoubleVec &w,
                      int start, int end) {
  int n;
  if (start > 0) {
    n = helpers.cumksize[end]-helpers.cumksize[start-1];
  } else
    n = helpers.cumksize[end];

  mat Winv = zeros(params.kk, params.kk);
  mat Pmat = eye(n, n) - ones(n,n)/n;
  if (start > 0)
    start = helpers.cumksize[start-1];
  end = helpers.cumksize[end]-1;
  mat Xtilde = Pmat*helpers.X(span(start, end), span(1,helpers.X.n_cols-1));
  mat XXW = Xtilde.t()*Xtilde;
  bool ok = FALSE;
  while(!ok) {
    for (int i = 0; i < params.kk; i++) {
      if (XXW(i,i) < 1e-12) {
        Xtilde = Pmat*(helpers.X(span(start, end), span(1,helpers.X.n_cols-1))+ mvrnormArma(n, params));
        XXW = Xtilde.t()*Xtilde;
        break;
      }
      Winv(i,i) = XXW(i,i)*w[i+1]/(1-w[i+1]);
      if (i == params.kk-1) ok = TRUE;
    }
  }
  XXW = XXW + Winv;
  mat sumxy = Xtilde.t()*helpers.Y.subvec(start,end);
  DoubleVec ret(2); // (Z, logC)
  ret[0] = as_scalar(sumxy.t()*XXW.i()*sumxy);
  double detval, detsign;
  mat tmp = XXW*Winv.i();
  log_det(detval, detsign, tmp);
  ret[1] = -0.5*detval;
  return ret;
}

// printtimp defaults to 0
MCMCStepSeq pass(MCMCStepSeq &step, HelperVariables &helpers,
                 Params &params, int printtmp=0) {

  int i, cp, tau, s;

  int bsize1, bsize2, bsize3;
  double bmean1, bmean2, bmean3, maxlik;
  double bK1=0, bK2=0, bK3=0, bK4=0;
  double tmp;

  DoubleVec likvals(6);
  DoubleVec logCvals(6);
  DoubleVec Kvals(6);
  DoubleVec Qvals(6);
  IntVec bvals(2);
  DoubleVec Wvals(2);
  DoubleVec Bvals(2);
  double Qval0, logCval0, Kval0;

  DoubleVec matOutput(2);
  DoubleVec matOutput1(2);
  DoubleVec matOutput2(2);
  int prevblock = 0;

  // this is used to reference the current MCMCStep we build from scratch
  MCMCStepSeq stepnew(step);
  int currblock = 0;

  // some other variables to denote stuff in current block and previous block
  // Note that "last" refers to the immediately left-adjacent block
  // whereas "prev" refers to the same block in the variable step
  double thisblockZ = step.bZ[0];
  int thisbend = step.bend[0];
  // double thisbetasqv = step.bbetasqv[0];
  int thisbtau = step.btau[0];
  double thisblogC = step.blogC[0];
  double thisbK = step.bK[0];
  double thisbQ = step.bQ[0];

  double lastblockZ = 0;
  int lastbend = -1; // this is simply a notational convenience
  // double lastbetasqv = 0;
  int lastbtau = 0;
  double lastblogC = 0;
  double lastbK = 0;
  double lastbQ = 0;

  // start the loop
  for (i = 0; i < params.nn - 1; i++) {//
    // Rprintf("i:%d\n", i);

    maxlik = -DBL_MAX;
    if (i == step.bend[prevblock]) {
      // we're at an old change point, so we need to refresh "this" to be the
      // immediately following block
      lastblockZ = thisblockZ;
      lastbtau = thisbtau;
      lastblogC = thisblogC;
      lastbK = thisbK;
      lastbQ = thisbQ;
      prevblock++;

      thisbend = step.bend[prevblock];
      thisblockZ = step.bZ[prevblock];
      thisbtau = step.btau[prevblock];
      thisblogC = step.blogC[prevblock];
      thisbK = step.bK[prevblock];
      thisbQ = step.bQ[prevblock];
    }

    // set the defaults
    if (step.rho[i] == 1) {
      Qval0 = stepnew.Q - thisbQ - lastbQ;
      logCval0 = stepnew.logC - thisblogC - lastblogC;
      Kval0 = stepnew.K - thisbK - lastbK;
    } else {
      Qval0 = stepnew.Q - thisbQ;
      logCval0 = stepnew.logC - thisblogC;
      Kval0 = stepnew.K - thisbK;
    }

    /****
    * consider if cp = 0 (not a change point)
    */
    bvals[0] = stepnew.b - 1 * (step.rho[i] == 1);
    if (lastbend > -1) {
      bsize3 = helpers.cumksize[thisbend] - helpers.cumksize[lastbend];
      bmean3 = (helpers.cumy[thisbend] - helpers.cumy[lastbend]) / bsize3;
    } else {
      bsize3 = helpers.cumksize[thisbend];
      bmean3 = helpers.cumy[thisbend] / bsize3;
    }
    tmp = 0;
    if (step.rho[i] == 1)
      tmp = thisblockZ + lastblockZ - pow(bmean3, 2) * bsize3;
    if (bvals[0] == 1)
      Bvals[0] = 0;
    else
      Bvals[0] = stepnew.B - tmp;
    Wvals[0] = stepnew.W + tmp;

    // now consider the change point type (tau)
    for (tau = 0; tau <= 1; tau++) {
      if (tau == 1 && bsize3 < params.nreg) {
        likvals[tau] = -DBL_MAX;
        break;
      }
      Qvals[tau] = Qval0;
      logCvals[tau] = logCval0;
      if (tau == 1) {
        matOutput = matrixCalcs(helpers, params, stepnew.w, lastbend+1, thisbend);
        Qvals[tau] += matOutput[0];
        logCvals[tau] += matOutput[1];
      }

      Kvals[tau] = Kval0 + logKcalc(bsize3, tau, params);
      likvals[tau] = likelihood(Bvals[0], Wvals[0], bvals[0], params, logCvals[tau],
                                Qvals[tau], Kvals[tau]);
      if (likvals[tau] > maxlik)
        maxlik = likvals[tau];
    }

    /****
    * consider if cp = 1 (make a change point)
    */
    bvals[1] = stepnew.b + (step.rho[i] != 1);
    bsize2 = helpers.cumksize[thisbend] - helpers.cumksize[i];
    if (lastbend > -1)
      bsize1 = helpers.cumksize[i] - helpers.cumksize[lastbend];
    else
      bsize1 = helpers.cumksize[i];
    tmp = 0;
    bmean2 = (helpers.cumy[thisbend] - helpers.cumy[i]) / bsize2;
    if (lastbend > -1) {
      bmean1 = (helpers.cumy[i] - helpers.cumy[lastbend]) / bsize1;
    } else {
      bmean1 = helpers.cumy[i] / bsize1;
    }
    if (step.rho[i] == 0)
      tmp = thisblockZ - pow(bmean1, 2) * bsize1 - pow(bmean2, 2) * bsize2;

    // Rprintf("thisbZ:%0.2f, bmean1:%0.2f, bsize1:%d, bmean2:%0.2f, bsize2:%d, tmp:%0.4f\n",
    //   thisblockZ, bmean1, bsize1,
    //   bmean2, bsize2, tmp);
    Bvals[1] = stepnew.B - tmp;
    Wvals[1] = stepnew.W + tmp;
    // now consider the change point type for 2 blocks [tau = (left, right)]
    // tau = 2: (0, 0), tau = 3: (0, 1), tau = 4: (1, 0), tau = 5: (1, 1)
    bK1 = logKcalc(bsize1, 0, params);
    if (bsize1 >= params.nreg) {
      bK2 = logKcalc(bsize1, 1, params);
      matOutput1 = matrixCalcs(helpers, params, stepnew.w, lastbend+1, i);
    }
    bK3 = logKcalc(bsize2, 0, params);
    if (bsize2 >= params.nreg) {
      matOutput2 = matrixCalcs(helpers, params, stepnew.w, i+1, thisbend);
      bK4 = logKcalc(bsize2,1, params);
    }
    // Rprintf("bK1:%0.2f, bK2:%0.2f bK3:%0.2f bK4:%0.2f\n", bK1, bK2, bK3,bK4);
    for (tau = 2; tau < 6; tau++) {
      if ((bsize1 < params.nreg && tau > 3) || (bsize2 < params.nreg &&
          (tau == 3 || tau == 5))) {
        likvals[tau] = -DBL_MAX;
        continue;
      }

      Kvals[tau] = Kval0 + bK1*(tau < 4) + bK2*(tau >= 4)+
        bK3*(tau == 2 || tau == 4)+
        bK4*(tau == 3 || tau == 5);
      Qvals[tau] = Qval0;
      logCvals[tau] = logCval0;
      if (tau >= 4) {
        Qvals[tau] += matOutput1[0];
        logCvals[tau] += matOutput1[1];
      }
      if (tau == 3 || tau == 5) {
        Qvals[tau] += matOutput2[0];
        logCvals[tau] += matOutput2[1];
      }
      likvals[tau] = likelihood(Bvals[1], Wvals[1], bvals[1], params,
                                logCvals[tau], Qvals[tau],
                                                    Kvals[tau]);
      if (likvals[tau] > maxlik)
        maxlik = likvals[tau];
    }
    // do the sampling and then updates

    s = sampleFromLikelihoods(likvals, maxlik);
    if (s < 2)
      cp = 0;
    else
      cp = 1;
    // Rprintf("i:%d s:%d likvals:%0.2f Qvals:%0.2f K:%0.2f C:%0.2f\n", i, s, likvals[s], Qvals[s], Kvals[s], logCvals[s]);
    if (cp != step.rho[i]) { // we modified the change point status
      if (cp == 0) {
        // removed a change point
        // update last block's stuff since the last block is now further back
        thisblockZ = pow(bmean3, 2) * bsize3;
        thisbK += Kvals[s] - stepnew.K+lastbK;

        if (currblock > 0) {
          lastbend = stepnew.bend[currblock - 1];
          lastblockZ = stepnew.bZ[currblock - 1];
          // lastbtau = stepnew.btau[currblock - 1];
          lastblogC = stepnew.blogC[currblock - 1];
          lastbQ = stepnew.bQ[currblock - 1];
          lastbK = stepnew.bK[currblock - 1];
        } else {
          lastblockZ = 0;
          lastbend = -1; // this is simply a notational convenience
          lastbQ = 0;
          // lastbtau = 0;
          lastblogC = 0;
          lastbK = 0;
        }
      } else { // added a change point
        thisblockZ = pow(bmean2, 2) * bsize2;
        lastblockZ = pow(bmean1, 2) * bsize1;

      }
    }
    stepnew.rho.push_back(cp);

    if (stepnew.rho[i] == 0) {
      if(step.rho[i]==0) thisbK += Kvals[s] - stepnew.K;
      if (s == 1) {
        thisbQ = matOutput[0];
        thisblogC = matOutput[1];
      } else {
        thisbQ = 0;
        thisblogC = 0;
      }
      thisbtau = s;


    } else if (stepnew.rho[i] == 1) {
      // we've added a change point, so we want to record some stuff
      if (s < 4) {
        lastbtau = 0;
        lastbQ = 0;
        lastblogC = 0;
      } else {
        lastbtau = 1;
        lastbQ = matOutput1[0];
        lastblogC = matOutput1[1];
      }

      if (s == 2 || s == 4) {
        thisbtau = 0;
        thisbQ = 0;
        thisblogC = 0;
      } else {
        thisbtau = 1;
        thisbQ = matOutput2[0];
        thisblogC = matOutput2[1];
      }

      lastbK = bK1*(s == 2 || s == 3) + bK2*(s >= 4);
      thisbK = bK3*(s == 2 || s == 4) + bK4*(s == 3 || s == 5);

      stepnew.bsize.push_back(bsize1);
      stepnew.bend.push_back(i);
      // stepnew.bmean.push_back(bmean1);
      stepnew.bZ.push_back(lastblockZ);
      stepnew.btau.push_back(lastbtau);
      stepnew.blogC.push_back(lastblogC);
      stepnew.bQ.push_back(lastbQ);
      stepnew.bK.push_back(lastbK);
      currblock++;
      lastbend = i;
    }
    stepnew.lik = likvals[s];
    stepnew.B = Bvals[cp];
    stepnew.W = Wvals[cp];
    stepnew.b = bvals[cp];
    stepnew.logC = logCvals[s];
    stepnew.Q = Qvals[s];
    stepnew.K = Kvals[s];
  }
  // done with a full pass, now let's add info on the final block
  if (lastbend < 0)
    stepnew.bsize.push_back(params.nn2);
  else
    stepnew.bsize.push_back(params.nn2 - helpers.cumksize[lastbend]);

  stepnew.bQ.push_back(thisbQ);
  stepnew.bend.push_back(params.nn - 1);
  stepnew.btau.push_back(thisbtau);
  stepnew.blogC.push_back(thisblogC);
  stepnew.bZ.push_back(thisblockZ);
  stepnew.bK.push_back(thisbK);

  // sample w2
  double lik;
  for (i = 1; i <= params.kk; i++) {
    DoubleVec w2cand = stepnew.w;
    w2cand[i] = stepnew.w[i] + Rf_runif( - 0.05 * params.w[i], 0.05 * params.w[i]);
    if (w2cand[i] < 0 || w2cand[i] > params.w[i]) continue;
    // Rprintf("w2cand:%0.2f\n", w2cand[i]);
    DoubleVec logCvec(stepnew.b);
    DoubleVec Qvec(stepnew.b);
    double Q0 = 0;
    double logC0 = 0;
    lastbend = -1;

    for (s = 0; s < stepnew.b; s++) {
      thisbend = stepnew.bend[s];
      // Rprintf("s:%d\n",s);
      if (stepnew.btau[s] == 1) {
        matOutput = matrixCalcs(helpers, params, w2cand, lastbend+1, thisbend);
        logCvec[s] = matOutput[1];
        Qvec[s] = matOutput[0];
        logC0 += logCvec[s];
        Q0 += Qvec[s];
      } else {
        logCvec[s] = 0;
        Qvec[s] = 0;
      }

      lastbend = thisbend;
    }
    lik = likelihood(stepnew.B, stepnew.W, stepnew.b, params, logC0,
                     Q0, stepnew.K);
    // Rprintf("lik:%0.9f\n", lik);
    double p = exp(lik - stepnew.lik);

    p = p / (1 + p);
    double u = Rf_runif(0.0, 1.0);
    // Rprintf("p:%0.8f u: %0.8f\n", p, u);
    if (u < p) {
      stepnew.w[i] = w2cand[i];
      stepnew.lik = lik;
      stepnew.blogC = logCvec;
      stepnew.bQ = Qvec;
      stepnew.logC = logC0;
      stepnew.Q = Q0;
    }
  }
  // Rprintf("w:%0.10f\n", stepnew.w[1]);
  stepnew.rho.push_back(1);
  return stepnew;
}

// [[Rcpp::export]]
SEXP rcpp_bcpR(SEXP py, SEXP px, SEXP pgrpinds, SEXP pid, SEXP pmcmcreturn, SEXP pburnin, SEXP pmcmc, SEXP pa,
               SEXP pw, SEXP pba, SEXP pnreg) {
  // INITIALIZATION OF LOCAL VARIABLES
  int i, j, m, start, end, start2, end2, resultStart, resultEnd;
  double Wtilde, wstar, bmean, xmax, tmpAlpha;

  // INITIALIZATION OF OTHER OBJECTS
  HelperVariables helpers(py, px, pid);
  mat grpInds = as<mat>(pgrpinds);
  // helpers.print();
  Params params(pw, helpers.cumksize.size(), helpers.Y.n_rows, pa, true, false,
                NUMERIC_DATA(pba)[0], INTEGER_DATA(pnreg)[0]);
  MCMCStepSeq step(helpers, params);
  mat Winv = zeros(params.kk, params.kk);
  int mcmcreturn = INTEGER_DATA(pmcmcreturn)[0];
  int burnin = INTEGER_DATA(pburnin)[0];
  int mcmc = INTEGER_DATA(pmcmc)[0];
  int MM = burnin + mcmc;

  int nn2 = params.nn;
  int MM2 = burnin + mcmc;
  if (mcmcreturn == 0) {
    MM2 = 1;
    nn2 = 1;
  }

  // Things to be returned to R:
  // NumericVector pmean(params.nn);
  vec pmean = zeros<vec>(params.nn);
  vec pvar = zeros<vec>(params.nn);
  vec ss = zeros<vec>(params.nn);
  NumericVector pchange(params.nn);
  NumericVector blocks(burnin + mcmc);
  mat betaposts = zeros(params.nn, 1+params.kk); // not calculating variances
  NumericMatrix rhos(nn2, MM2);
  // store (fitted, intercept, slopes)
  mat results = zeros(MM2*nn2, 2+params.kk); // store conditional means
  // NumericVector liks(MM2);
  GetRNGstate(); // Consider Dirk's comment on this.
  // step.print();

  for (m = 0; m < MM; m++) {
    step = pass(step, helpers, params);
    blocks[m] = step.b;
    if (mcmcreturn == 1) {
      // liks[m] = step.lik;
      for (j = 0; j < params.nn; j++) {
        rhos(j, m) = step.rho[j];
      }
    }
    if (m >= burnin || mcmcreturn == 1) {
      // Rprintf("m:%d\n", m);
      // step.print();
      // compute posteriors
      Wtilde = step.W - step.Q;
      if (step.b > 1) {
        if (step.B <= 0)
          step.B = 0.00001;
        xmax = step.B * params.w[0] / Wtilde / (1 + step.B * params.w[0] / Wtilde);
        wstar = log(Wtilde) - log(step.B)
          + Rf_lbeta((double) (step.b + 3) / 2, (double) (params.nn2 - step.b - 4) / 2)
          + Rf_pbeta(xmax, (double) (step.b + 3) / 2, (double) (params.nn2 - step.b - 4) / 2, 1, 1)
          - Rf_lbeta((double) (step.b + 1) / 2, (double) (params.nn2 - step.b - 2) / 2)
          - Rf_pbeta(xmax, (double) (step.b + 1) / 2, (double) (params.nn2 - step.b - 2) / 2, 1, 1);
          wstar = exp(wstar);
      } else
        wstar = params.w[0]/2;
      // pvar += (Wtilde + wstar*step.B)/(params.nn2 - 3);

      start = 0;
      for (i = 0; i < step.b; i++) {
        end = step.bend[i];

        if (start > 0) {
          bmean = (helpers.cumy[end]-helpers.cumy[start-1])/step.bsize[i];
          start2 = helpers.cumksize[start-1];
        } else  {
          bmean = helpers.cumy[end]/step.bsize[i];
          start2 = 0;
        }
        end2 = helpers.cumksize[end]-1;
        mat grpIndMat = grpInds(span(start, end), span(start2, end2));

        tmpAlpha = bmean * (1 - wstar) + helpers.ybar * wstar;
        vec tmpMean = grpIndMat*(ones(step.bsize[i], 1)* tmpAlpha);
        vec fitted = tmpMean;
        if (mcmcreturn == 1) {
          resultStart = params.nn*m+start;
          resultEnd = params.nn*m+end;
          results(span(resultStart, resultEnd), 1) = tmpMean;
        }
        if (m >= burnin)
          betaposts(span(start, end), 0) += tmpMean;

        if (step.btau[i] == 1) {
          mat Pmat = eye(step.bsize[i], step.bsize[i])
          - ones(step.bsize[i],step.bsize[i])/step.bsize[i];

          mat Xtilde = Pmat*helpers.X(span(start2, end2), span(1, params.kk));
          mat XXW = Xtilde.t()*Xtilde;
          bool ok = FALSE;
          while(!ok) {
            for (int ii = 0; ii < params.kk; ii++) {
              if (XXW(ii,ii) < 1e-12) {
                Xtilde = Pmat*(helpers.X(span(start2, end2), span(1,params.kk))+
                  mvrnormArma(step.bsize[i], params));
                XXW = Xtilde.t()*Xtilde;
                break; // this breaks the loop and since !ok, checks for XXW=0 all over again
              }
              Winv(ii,ii) = XXW(ii,ii)*step.w[ii+1]/(1-step.w[ii+1]);
              if (ii == params.kk-1) ok = TRUE;
            }
          }
          XXW = XXW + Winv;
          mat sumxy = Xtilde.t()*helpers.Y.subvec(start2,end2);
          vec bhats = XXW.i()*sumxy;
          mat bhatMat = repmat(bhats.t(), end-start+1, 1);
          fitted += grpIndMat*Xtilde*bhats;
          if (m >= burnin) {
            betaposts(span(start, end), span(1, params.kk)) += bhatMat;
            betaposts(span(start, end), 0) -= bhatMat*helpers.X(span(start2, end2), span(1, params.kk)).t()*
              ones(step.bsize[i], 1)/step.bsize[i];
          }
          if (mcmcreturn == 1) {
            results(span(resultStart, resultEnd), 0) = fitted;
            results(span(resultStart, resultEnd), span(2, params.kk+1)) = bhatMat;
            results(span(resultStart, resultEnd), 1) -= bhatMat*helpers.X(span(start2, end2), span(1, params.kk)).t()*
              ones(step.bsize[i], 1)/step.bsize[i];
          }

        }
        if (m >= burnin) {
          pmean.subvec(start,end) += fitted;
          ss.subvec(start, end) += fitted % fitted;
        }

        start = end+1;

      }
      if (m >= burnin)
        for (j = 0; j < params.nn; j++) {
          pchange[j] += (double) step.rho[j];
        }
    }
  }

  // post processing
  for (j = 0; j < params.nn; j++) {
    pchange[j] = pchange[j] / (double) mcmc;
    pmean[j] = pmean[j] / mcmc;
    pvar[j] = (ss[j] / mcmc - pmean[j] * pmean[j])*(mcmc/(mcmc-1));
  }
  // pvar = pvar / mcmc;
  betaposts /= (double) mcmc;

  PutRNGstate();

  List z;
  z["posterior.mean"] = wrap(pmean);
  z["posterior.var"] = wrap(pvar);
  z["posterior.prob"] = pchange;
  z["blocks"] = blocks;
  z["mcmc.rhos"] = rhos;
  z["betaposts"] = wrap(betaposts);
  z["mcmc.means"] = wrap(results);
  // z["lik"] = liks;
  return z;
} /* END MAIN  */
