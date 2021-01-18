#ifndef BESTROUND_HPP
#define BESTROUND_HPP

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

#include "intCalib.h"
#include "Losses/Losses.h"

template <typename Type> colvec bestRound(const Type& A, const colvec& y, colvec& w,
                                          const mat& Bnds, const colvec& scale,
                                          const std::string lossType) {
    /* Optimal Rounding to Integer Algorithm
     * INPUT:
     *   A, matrix of comodities per each farm
     *   y, vector of targets
     *   w, vector of initial real weights to be rounded
     *   Rc, vector of indicators values for Restricted weights
     *   Bnds, matrix of target-bounds centered for zeros point targets
     *   scale, vector of rescaling factors
     *   lossType, string indicating which loss function to optimize (see "intCalib.hpp" for additional details)
     *
     * OUTPUT:
     *   w, vector of calibrated integer weights
     */

    int i, j;
    double mantissa, ffcm, ffcp;
    LogicalVector ismod(w.size(), true); // Is the i-th weight modifiable?
    colvec grad(w.size());
    colvec ee(y.size());
    colvec tmpEm(y.size());
    colvec tmpEp(y.size());
    colvec s, tau, lambda, iscale;
    umat ord;
#if _DEBUG
    double ffnew;
    int oobo;
    time_t now;
#endif // _DEBUG

    int lt = 0;
    if(lossType == "L1") lt = LOSS_L1;
    if(lossType == "aL1") lt = LOSS_aL1;
    if(lossType == "rL1") lt = LOSS_rL1;
    if(lossType == "LB1") lt = LOSS_LB1;
    if(lossType == "rB1") lt = LOSS_rB1;
    if(lossType == "L2") lt = LOSS_L2;
    if(lossType == "aL2") lt = LOSS_aL2;
    if(lossType == "rL2") lt = LOSS_rL2;
    if(lossType == "LB2") lt = LOSS_LB2;
    if(lossType == "rB2") lt = LOSS_rB2;

    switch (lt) {
    case LOSS_aL1: // bounds asymmetry factor
    case LOSS_aL2:
        tau = Bnds.col(1) / (Bnds.col(1) - Bnds.col(0));
        break;
    case LOSS_rL1: // scale factor for relative errors
    case LOSS_rL2:
        iscale = 1.0 / scale;
        break;
    default:
        break;
    }

    // Calculate the current loss and the gradient
    ee = y - A * w;
    switch (lt) {
    case LOSS_L1: // L1-norm minimization
        s = sign(ee);
        grad = L1::ffGrd(A, s);
        break;
    case LOSS_aL1: // Asymmetric L1-norm
        lambda = tau - (ee < 0.0);
        grad = aL1::ffGrd(A, lambda);
        break;
    case LOSS_rL1: // L1-norm with relative errors
        s = sign(ee);
        grad = rL1::ffGrd(A, s, iscale);
        break;
    case LOSS_LB1: // L1-norm outside the boundaries and zero within
        s = conv_to<colvec>::from(ee > Bnds.col(1)) - conv_to<colvec>::from(ee < Bnds.col(0));
        grad = LB1::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
        break;
    case LOSS_rB1: // relative L1-norm outside the boundaries and zero within
        s = conv_to<colvec>::from(ee > Bnds.col(1)) - conv_to<colvec>::from(ee < Bnds.col(0));
        grad = rB1::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
        break;
    case LOSS_L2: // L2-norm minimization
        grad = L2::ffGrd(A, ee);
        break;
    case LOSS_aL2: // Asymmetric L2-norm
        lambda = abs(tau - (ee < 0.0));
        grad = aL2::ffGrd(A, ee, lambda);
        break;
    case LOSS_rL2: // L2-norm with relative errors
        grad = rL2::ffGrd(A, ee, iscale);
        break;
    case LOSS_LB2: // L2-norm outside the boundaries and zero within
        grad = LB2::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
        break;
    case LOSS_rB2: // relative L2-norm outside the boundaries and zero within
        grad = rB2::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
        break;
    default:
        grad.zeros();
        break;
    }
    ord = stable_sort_index(abs(grad), "descend"); // sort the gradient

    for (j = 0; ; j = 0) {
        for (; (size_t) j < w.size(); j++) if(ismod[ord[j]]) break; // find roundable position
        if ((size_t) j >= w.size()) break;
        i = ord[j];
        ismod[i] = false;
        mantissa = w[i] - trunc(w[i]);
        if (mantissa == 0.0) continue;
        tmpEm = ee + mantissa * A.col(i);
        tmpEp = ee - (1.0 - mantissa) * A.col(i);
        switch (lt) {
        case LOSS_L1: // L1-norm minimization
            ffcm = L1::ff(tmpEm);
            ffcp = L1::ff(tmpEp);
            break;
        case LOSS_aL1: // Asymmetric L1-norm
            lambda = tau - (tmpEm < 0.0);
            ffcm = aL1::ff(tmpEm, lambda);
            lambda = tau - (tmpEp < 0.0);
            ffcp = aL1::ff(tmpEp, lambda);
            break;
        case LOSS_rL1: // L1-norm with relative errors
            ffcm = rL1::ff(tmpEm, iscale);
            ffcp = rL1::ff(tmpEp, iscale);
            break;
        case LOSS_LB1: // L1-norm outside the boundaries and zero within
            ffcm = LB1::ff(Bnds.col(0), Bnds.col(1), tmpEm);
            ffcp = LB1::ff(Bnds.col(0), Bnds.col(1), tmpEp);
            break;
        case LOSS_rB1: // relative L1-norm outside the boundaries and zero within
            ffcm = rB1::ff(Bnds.col(0), Bnds.col(1), tmpEm);
            ffcp = rB1::ff(Bnds.col(0), Bnds.col(1), tmpEp);
            break;
        case LOSS_L2: // L2-norm minimization
            ffcm = L2::ff(tmpEm);
            ffcp = L2::ff(tmpEp);
            break;
        case LOSS_aL2: // Asymmetric L2-norm
            lambda = tau - (tmpEm < 0.0);
            ffcm = aL2::ff(tmpEm, lambda);
            lambda = tau - (tmpEp < 0.0);
            ffcp = aL2::ff(tmpEp, lambda);
            break;
        case LOSS_rL2: // L2-norm with relative errors
            ffcm = rL2::ff(tmpEm, iscale);
            ffcp = rL2::ff(tmpEp, iscale);
            break;
        case LOSS_LB2: // L2-norm outside the boundaries and zero within
            ffcm = LB2::ff(Bnds.col(0), Bnds.col(1), tmpEm);
            ffcp = LB2::ff(Bnds.col(0), Bnds.col(1), tmpEp);
            break;
        case LOSS_rB2: // relative L2-norm outside the boundaries and zero within
            ffcm = rB2::ff(Bnds.col(0), Bnds.col(1), tmpEm);
            ffcp = rB2::ff(Bnds.col(0), Bnds.col(1), tmpEp);
            break;
        default:
            ffcm = 0.0;
            ffcp = 0.0;
            break;
        }
        if (ffcm < ffcp) {
            ee = tmpEm;
            w[i] = trunc(w[i]);
#if _DEBUG
            ffnew = ffcm;
#endif // _DEBUG
        }
        else {
            if (ffcm > ffcp) {
                ee = tmpEp;
                w[i] = ceil(w[i]);
#if _DEBUG
                ffnew = ffcp;
#endif // _DEBUG
            }
            else {
                if (mantissa < 0.5) {
                    ee = tmpEm;
                    w[i] = trunc(w[i]);
#if _DEBUG
                    ffnew = ffcm;
#endif // _DEBUG
                }
                else {
                    ee = tmpEp;
                    w[i] = ceil(w[i]);
#if _DEBUG
                    ffnew = ffcp;
#endif // _DEBUG
                }
            }
        }
        // Update the gradient
        switch (lt) {
        case LOSS_L1: // L1-norm minimization
            L1::updategrd(A, s, ee, grad, ord, j);
            s = sign(ee);
            break;
        case LOSS_aL1: // Asymmetric L1-norm
            aL1::updategrd(A, s, ee, grad, ord, j);
            lambda = tau - (ee < 0.0);
            break;
        case LOSS_rL1: // L1-norm with relative errors
            rL1::updategrd(A, iscale, s, ee, grad, ord, j);
            s = sign(ee);
            break;
        case LOSS_LB1: // L1-norm outside the boundaries and zero within
            LB1::updategrd(A, Bnds, s, ee, grad, ord, j);
            s = conv_to<colvec>::from(ee > Bnds.col(1)) - conv_to<colvec>::from(ee < Bnds.col(0));
            break;
        case LOSS_rB1: // relative L1-norm outside the boundaries and zero within
            rB1::updategrd(A, Bnds, s, ee, grad, ord, j);
//            Rprintf("MaxDiffGrad: %f, ", max(abs(grad - rB1::ffGrd(A, ee, Bnds.col(0), Bnds.col(1)))));
            s = conv_to<colvec>::from(ee > Bnds.col(1)) - conv_to<colvec>::from(ee < Bnds.col(0));
            break;
        case LOSS_L2: // L2-norm minimization
            grad = L2::ffGrd(A, ee);
            break;
        case LOSS_aL2: // Asymmetric L2-norm
            lambda = abs(tau - (ee < 0.0));
            grad = aL2::ffGrd(A, ee, lambda);
            break;
        case LOSS_rL2: // L2-norm with relative errors
            grad = rL2::ffGrd(A, ee, iscale);
            break;
        case LOSS_LB2: // L2-norm outside the boundaries and zero within
            grad = LB2::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
            break;
        case LOSS_rB2: // relative L2-norm outside the boundaries and zero within
            grad = rB2::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
            break;
        default:
            break;
        }
        if (!IS_L1BASED(lt)) {
            ord = stable_sort_index(abs(grad), "descend");
        }
#if _DEBUG
        now = std::time(0);
        for (i = 0, oobo = 0; (size_t) i < ee.size(); i++) {
            if (ee[i] < Bnds.col(0)[i] || ee[i] > Bnds.col(1)[i]) oobo++;
        }
#if _DEBUG != 2
        Rprintf("M: %d, F:%f, ", oobo, ffnew);
#else // IF _DEBUG == 2 THEN
        Rprintf("\033[1mM:\033[31m%d\033[39;0m, \033[1mF:\033[32m%f\033[39;0m, ", oobo, ffnew);
#endif // _DEBUG == 2
#if _POSCHK
        Rprintf("pos = %d # ", ord[j]);
#endif // _POSCHK
        Rprintf("%s", std::ctime(&now));
#endif // _DEBUG
    }
    return w;
}

#endif // BESTROUND_HPP
