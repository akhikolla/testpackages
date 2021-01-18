#ifndef INTCALIB_HPP
#define INTCALIB_HPP

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;


// Calculate the number of fitted targets outside the boundaries
int countOOB(const colvec& L, const colvec& U, const colvec& e) {
    int i, oob = 0;
    for (i = 0; (size_t) i < e.size(); i++) {
        if (L[i] > e[i] || U[i] < e[i]) oob++;
    }
    return oob;
}

#include "Losses/Losses.h"

template <typename Type> colvec IntProgCalib(const Type& A, const colvec& y, colvec& w, const colvec& dse,
                                             const vec& lower, const vec& upper, const mat& Bnds,
                                             const colvec& scale, const std::string lossType) {
    /* Integer Programming Calibration algorithm
     * INPUT:
     *   A, matrix of comodities per each farm
     *   y, vector of targets
     *   w, vector of initial weights
     *   dse, vector of dse weights as reference vector (used only for LASSO type objective)
     *   lower, vector of lower bounds for Restricted weights
     *   upper, vector of lower bounds for Restricted weights
     *   Bnds, matrix of target-bounds centered for zeros point targets
     *   scale, vector of rescaling factors
     *   lossType, string indicating which loss function to optimize:
     *     * "L1", summation of absolute errors
     *     * "aL1", asymmetric summation of absolute errors
     *     * "rL1", summation of absolute relative errors
     *     * "LB1", summation of absolute errors if outside the boundaries
     *     * "rB1", summation of absolute relative errors if outside the boundaries
     *     * "rbLasso1", summation of absolute relative errors if outside the boundaries + Lasso based on dse weights
     *     * "L2", summation of square errors
     *     * "aL2", asymmetric summation of square errors
     *     * "rL2", summation of square relative errors
     *     * "LB2", summation of square errors if outside the boundaries
     *     * "rB2", summation of square relative errors if outside the boundaries
     *     * "rbLasso2", summation of square relative errors if outside the boundaries + Lasso based on dse weights
     *
     * OUTPUT:
     *   w, vector of calibrated integer weights
     */
#if _DEBUG
#if _POSCHK
    int df, pos;
#endif // _POSCHK
    time_t now;
#endif // _DEBUG
    int i, j, ac = 0;
    colvec cc(w.size());
    colvec grad(w.size());
    colvec ee(y.size());
    colvec tmpE(y.size());
    double ffnew, ffcc;
    double minPen = 0.0, totPen = 0.0;
    int oob, oobo;
    colvec s(y.size());
    colvec sp(w.size());
    colvec pen(w.size()), tpen(w.size());
    umat ord, nord;
    colvec tau, lambda, iscale;

    int lt = 0;
    if(lossType == "L1") lt = LOSS_L1;
    if(lossType == "aL1") lt = LOSS_aL1;
    if(lossType == "rL1") lt = LOSS_rL1;
    if(lossType == "LB1") lt = LOSS_LB1;
    if(lossType == "rB1") lt = LOSS_rB1;
    if(lossType == "rbLasso1") lt = LOSS_RBLASSO1;
    if(lossType == "L2") lt = LOSS_L2;
    if(lossType == "aL2") lt = LOSS_aL2;
    if(lossType == "rL2") lt = LOSS_rL2;
    if(lossType == "LB2") lt = LOSS_LB2;
    if(lossType == "rB2") lt = LOSS_rB2;
    if(lossType == "rbLasso2") lt = LOSS_RBLASSO2;

    switch (lt) {
    case LOSS_aL1: // bounds asymmetry factor
    case LOSS_aL2:
        tau = Bnds.col(1) / (Bnds.col(1) - Bnds.col(0));
        break;
    case LOSS_rL1: // scale factor for relative errors
    case LOSS_rL2:
        iscale = 1.0 / scale;
        break;
    case LOSS_RBLASSO1:
    case LOSS_RBLASSO2:
        for (i = 0; (size_t) i < w.size(); i++) {
            if (w[i] < lower[i]) {
                w[i] = lower[i];
            }
            else {
                if (w[i] > upper[i]) {
                    w[i] = upper[i];
                }
            }
        }
        w = round(w);
        minPen = sum(abs(dse - w));
        break;
    default:
        break;
    }

    /* Calculate the errors, count the number
     * of missed targets, and create a temporary
     * copy of the vector of weights. */
    ee = y - A * w;
    oobo = countOOB(Bnds.col(0), Bnds.col(1), ee);
    cc = w;
    /* Calculate the currnet Loss */
    switch (lt) {
    case LOSS_L1: // L1-norm minimization
        ffnew = L1::ff(ee);
        s = sign(ee);
        break;
    case LOSS_aL1: // Asymmetric L1-norm
        s = sign(ee);
        lambda = tau - (ee < 0.0);
        ffnew = aL1::ff(ee, lambda);
        break;
    case LOSS_rL1: // L1-norm with relative errors
        ffnew = rL1::ff(ee, iscale);
        s = sign(ee);
        break;
    case LOSS_LB1: // L1-norm outside the boundaries and zero within
        ffnew = LB1::ff(Bnds.col(0), Bnds.col(1), ee);
        s = conv_to<colvec>::from(ee > Bnds.col(1)) - conv_to<colvec>::from(ee < Bnds.col(0));
        break;
    case LOSS_rB1: // relative L1-norm outside the boundaries and zero within
        ffnew = rB1::ff(Bnds.col(0), Bnds.col(1), ee);
        s = conv_to<colvec>::from((ee > Bnds.col(1)) - (ee < Bnds.col(0)));
        break;
    case LOSS_RBLASSO1: // relative L1-norm outside the boundaries and zero within + LASSO penalty based on DSE
        pen = dse - w;
        tpen = pen;
        totPen = L1::ff(pen) - minPen;
        ffnew = rB1::ff(Bnds.col(0), Bnds.col(1), ee) + totPen;
        s = conv_to<colvec>::from(ee > Bnds.col(1)) - conv_to<colvec>::from(ee < Bnds.col(0));
        sp = sign(pen);
        break;
    case LOSS_L2: // L2-norm minimization
        ffnew = L2::ff(ee);
        break;
    case LOSS_aL2: // Asymmetric L2-norm
        lambda = abs(tau - (ee < 0.0));
        ffnew = aL2::ff(ee, lambda);
        break;
    case LOSS_rL2: // L2-norm with relative errors
        ffnew = rL2::ff(ee, iscale);
        break;
    case LOSS_LB2: // L2-norm outside the boundaries and zero within
        ffnew = LB2::ff(Bnds.col(0), Bnds.col(1), ee);
        break;
    case LOSS_rB2: // relative L2-norm outside the boundaries and zero within
        ffnew = rB2::ff(Bnds.col(0), Bnds.col(1), ee);
        break;
    case LOSS_RBLASSO2: // relative L2-norm outside the boundaries and zero within + LASSO penalty based on DSE
        pen = dse - w;
        tpen = pen;
        totPen = L1::ff(pen) - minPen;
        ffnew = rB2::ff(Bnds.col(0), Bnds.col(1), ee) + totPen;
        s = conv_to<colvec>::from(ee > Bnds.col(1)) - conv_to<colvec>::from(ee < Bnds.col(0));
        sp = sign(pen);
        break;
    default:
        break;
    }

#if _DEBUG
                now = std::time(0);
                Rprintf("%s", std::ctime(&now));
#if _DEBUG != 2
                Rprintf("M:%d, F:%f, ", oobo, ffnew);
#else // IF _DEBUG == 2 THEN
                Rprintf("\033[1mM:\033[31m%d\033[39;0m, \033[1mF:\033[32m%f\033[39;0m\n", oobo, ffnew);
#endif // _DEBUG == 2
#endif // _DEBUG

    // Repeat the following until there are no more changes on "w"
    do {
        /* Calculate the gradient, and determine
         * the optimal ordered positions. */
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
        case LOSS_RBLASSO1: // relative L1-norm outside the boundaries and zero within + LASSO penalty based on DSE
            s = conv_to<colvec>::from((ee > Bnds.col(1)) - (ee < Bnds.col(0)));
            grad = rB1::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
            sp = sign(pen);
            grad -= sp;
            break;
        case LOSS_L2: // L2-norm minimization
            if (ac == 0) grad = L2::ffGrd(A, ee);
            break;
        case LOSS_aL2: // Asymmetric L2-norm
            lambda = abs(tau - (ee < 0.0));
            grad = aL2::ffGrd(A, ee, lambda);
            break;
        case LOSS_rL2: // L2-norm with relative errors
            if (ac == 0) grad = rL2::ffGrd(A, ee, iscale);
            break;
        case LOSS_LB2: // L2-norm outside the boundaries and zero within
            if (ac == 0) grad = LB2::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
            break;
        case LOSS_rB2: // relative L2-norm outside the boundaries and zero within
            if (ac == 0) grad = rB2::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
            break;
        case LOSS_RBLASSO2: // relative L2-norm outside the boundaries and zero within + LASSO penalty based on DSE
            if (ac == 0) {
                grad = rB2::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
                sp = sign(pen);
                grad -= sp;
            }
            break;
        default:
            grad.zeros();
            break;
        }
        ord = stable_sort_index(abs(grad), "descend");

        /* Reset the counter of any change in vector w */
        ac = 0;
#if _DEBUG
#if _POSCHK
        pos = 0;
#endif // _POSCHK
#endif // _DEBUG
        for (j = 0; (size_t) j < w.size(); j++) {
            R_CheckUserInterrupt();
            /* Consider the weights with the most effective reduction,
             * change them according to the gradient, and update
             * the values of the errors (which are store on a temporary vector). */
            i = ord[j];
            if (grad[i] == 0.0) continue;
            if (grad[i] > 0.0) {
                if (cc[i] < (lower[i] + 1.0)) continue; //{Rprintf("%d - w_i < 2\n", j); continue;}
                cc[i]--;
                tmpE = ee + A.col(i);
                tpen[i]++;
            }
            else {
                if (cc[i] > (upper[i] - 1.0)) continue; //{Rprintf("%d - w_i > 5\n", j); continue;}
                cc[i]++;
                tmpE = ee - A.col(i);
                tpen[i]--;
            }
            /* Counting the number of missed targets (Out Of Bounds),
             * and the new value of the loss function to optimize. */
            oob = countOOB(Bnds.col(0), Bnds.col(1), tmpE);
            switch (lt) {
            case LOSS_L1:
                ffcc = L1::ff(tmpE);
                break;
            case LOSS_aL1:
                lambda = tau - (tmpE < 0.0);
                ffcc = aL1::ff(tmpE, lambda);
                break;
            case LOSS_rL1:
                ffcc = rL1::ff(tmpE, iscale);
                break;
            case LOSS_LB1:
                ffcc = LB1::ff(Bnds.col(0), Bnds.col(1), tmpE);
                break;
            case LOSS_rB1:
                ffcc = rB1::ff(Bnds.col(0), Bnds.col(1), tmpE);
                break;
            case LOSS_RBLASSO1:
                ffcc = rB1::ff(Bnds.col(0), Bnds.col(1), tmpE) + totPen - fabs(pen[i]) + fabs(tpen[i]);
                break;
            case LOSS_L2:
                ffcc = L2::ff(tmpE);
                break;
            case LOSS_aL2:
                lambda = abs(tau - (tmpE < 0.0));
                ffcc = aL2::ff(tmpE, lambda);
                break;
            case LOSS_rL2:
                ffcc = rL2::ff(tmpE, iscale);
                break;
            case LOSS_LB2:
                ffcc = LB2::ff(Bnds.col(0), Bnds.col(1), tmpE);
                break;
            case LOSS_rB2:
                ffcc = rB2::ff(Bnds.col(0), Bnds.col(1), tmpE);
                break;
            case LOSS_RBLASSO2:
                ffcc = rB2::ff(Bnds.col(0), Bnds.col(1), tmpE) + totPen - fabs(pen[i]) + fabs(tpen[i]);
                break;
            default:
                break;
            }
            /* If both the loss function and the number of missed
             * targets are reduced, then the old point is changed,
             * the gradient is updated and the algorithm starts
             * from the beginning until convergence (based on the
             * number of changes). */
            if ((ffcc < ffnew) && OFF_WHEN_BND(oob <= oobo, lt)) {// -1 instead of lt ???
#if _DEBUG
                now = std::time(0);
#if _POSCHK
                df = j - pos;
#endif
#if _DEBUG != 2
                Rprintf("M:%d, F:%f, ", oob, ffcc);
                Rprintf("S:%s # ", grad[i] > 0.0 ? "-" : "+");
#if _POSCHK
                Rprintf("[after %d in %dG of %dW] # ", df, j, i);
#endif
#else // IF _DEBUG == 2 THEN
                Rprintf("\033[1mM:\033[31m%d\033[39;0m, \033[1mF:\033[32m%f\033[39;0m, ", oob, ffcc);
                Rprintf("S:\033[34m%s\033[39;0m # ", grad[i] > 0.0 ? "-" : "+");
#if _POSCHK
                Rprintf("[after %d in %dG of %dW] # ", df, j, i);
#endif
#endif // _DEBUG == 2
                Rprintf("%s", std::ctime(&now));
#endif // _DEBUG
                ac++; // Count how many changes of the weigths are considered
                ffnew = ffcc;
                w[i] = cc[i];
                ee = tmpE;
                switch (lt) {
                case LOSS_L1:
                    j = L1::updategrd(A, s, ee, grad, ord, j);
                    s = sign(ee);
                    break;
                case LOSS_aL1:
                    j = aL1::updategrd(A, s, ee, grad, ord, j);
                    s = sign(ee);
                    break;
                case LOSS_rL1:
                    j = rL1::updategrd(A, iscale, s, ee, grad, ord, j);
                    s = sign(ee);
                    break;
                case LOSS_LB1:
                    j = LB1::updategrd(A, Bnds, s, ee, grad, ord, j);
                    s = conv_to<colvec>::from(ee > Bnds.col(1)) - conv_to<colvec>::from(ee < Bnds.col(0));
                    break;
                case LOSS_rB1:
                    j = rB1::updategrd(A, Bnds, s, ee, grad, ord, j);
                    s = conv_to<colvec>::from(ee > Bnds.col(1)) - conv_to<colvec>::from(ee < Bnds.col(0));
                    break;
                case LOSS_RBLASSO1:
                    totPen += fabs(tpen[i]) - fabs(pen[i]); //update the total penalty
                    pen[i] = tpen[i]; // update the partial penalties before absolute value
                    j = rB1::updategrd(A, Bnds, s, ee, grad, ord, j);
                    if (sp[i] != d_sign(pen[i])) {
                        sp[i] = d_sign(pen[i]); //update sign for partial penalties
                        grad[i] -= sp[i]; //update gradient
                        j = -1;
                    }
                    grad -= sp;
                    break;
                case LOSS_L2:
                    grad = L2::ffGrd(A, ee);
                    break;
                case LOSS_aL2:
                    grad = aL2::ffGrd(A, ee, lambda);
                    break;
                case LOSS_rL2:
                    grad = rL2::ffGrd(A, ee, iscale);
                    break;
                case LOSS_LB2:
                    grad = LB2::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
                    break;
                case LOSS_rB2:
                    grad = rB2::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
                    break;
                case LOSS_RBLASSO2:
                    totPen += fabs(tpen[i]) - fabs(pen[i]); //update the total penalty
                    pen[i] = tpen[i]; // update the partial penalties before absolute value
                    grad = rB2::ffGrd(A, ee, Bnds.col(0), Bnds.col(1));
                    if (sp[i] != d_sign(pen[i])) {
                        sp[i] = d_sign(pen[i]); //update sign for partial penalties
                        grad[i] -= sp[i]; //update gradient
                    }
                    break;
                default:
                    break;
                }
                if (!IS_L1BASED(lt)) { // if the objective function is based on L2-norm (or L0-norm)
                    nord = stable_sort_index(abs(grad), "descend"); //sort the gradient
                    if (any(conv_to<vec>::from(nord != ord))) { // and calculate the position according to the changes
                        j = -1;
                        ord = nord;
                    }
                }
                if (oob < oobo) oobo = oob;
#if _DEBUG
#if _POSCHK
                if (j < 0) {
                    pos ^= pos;
                }
                else {
                    pos = j;
                }
#endif // _POSCHK
#endif // _DEBUG
            }
            else {
                cc[i] = w[i];
            }
        }
    } while(ac);
#if _DEBUG
    now = std::time(0);
    Rprintf("%s", std::ctime(&now));
#endif
    return w;
}

#endif // INTCALIB_HPP
