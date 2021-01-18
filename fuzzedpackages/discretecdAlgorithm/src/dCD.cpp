#include "dBasic.h"
#include "dCD.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
#include <math.h>
#define MATHLIB_STANDALONE
#define RMATH_DLL
#include "Rmath.h"

using namespace Eigen;
using namespace std;



// int i, int j, MatrixXi nzIndex(?,2), MatrixXd beta(di, rj), bool isBetaZero
// MatrixXd pm(nh, rj), MatrixXd logitm(nh, rj), MatrixXd ym(nh, rj), MatrixXd dm(nh, di),
void NewtonIterUpdate(const int& nh, const int& di, const int& rj, const int& nrows, const MatrixXd& dm, const MatrixXd& ym, const MatrixXi& nzIndex,
                     MatrixXd& logitm, MatrixXd& beta, bool& isBetaZero, const double& hval, const double& penalty, const double& qtol)
{

	int it1, rInd, cInd;
	double norm, mad;
	MatrixXd yMp(nh, rj), dif(di, rj), grad(di, rj), tmp(di, rj);
	VectorXd rowMax;

    // one Newton step
    grad = MatrixXd::Zero(di, rj);

    // compute the yMp matrix to prepare for gradient computation
    yMp = logitm;
    rowMax = yMp.rowwise().maxCoeff();
    yMp.colwise() -= rowMax;
    yMp = yMp.array().exp();
    for (it1 = 0; it1 < nh; ++it1)
    {
        yMp.row(it1) = ym.row(it1) - yMp.row(it1) / yMp.row(it1).sum();
    }

    // compute the gradient
    for (it1 = 0; it1 < nrows; ++it1)
    {
        rInd = nzIndex(it1, 0);
        cInd = nzIndex(it1, 1);
        grad.row(cInd) +=  yMp.row(rInd) * dm(rInd, cInd);
    }

    // update beta
    if(isBetaZero)
    {
        norm = grad.norm();
        if (norm <= penalty)
        {
            return;
        }
        else
        {
            isBetaZero = false;
            beta = dif = ((penalty / norm - 1.0) / hval) * grad;
            mad = dif.lpNorm<Infinity>();
        }
    }
    else
    {
        grad -= hval * beta;
        norm = grad.norm();
        if (norm <= penalty)
        {
            isBetaZero = true;
            dif = -beta;
            mad = dif.lpNorm<Infinity>();
            beta = MatrixXd::Zero(di, rj);
        }
        else
        {
            tmp = ((penalty / norm - 1.0) / hval) * grad;
            dif = tmp - beta;
            mad = dif.lpNorm<Infinity>();
            if (mad < qtol)	return;
            else	beta = tmp;
        }
    }

    // update the logit matrix
    for (it1 = 0; it1 < nrows; ++it1)
    {
        rInd = nzIndex(it1, 0);
        cInd = nzIndex(it1, 1);
        logitm.row(rInd) += dif.row(cInd) * dm(rInd, cInd);
    }

    // check convergence
    if (mad < qtol)	return;

}

void NewtonIterTmp(const int& nh, const int& di, const int& rj, const int& nrows, const MatrixXd& dm, const MatrixXd& ym, const MatrixXi& nzIndex,
                   const MatrixXd& logitm, const MatrixXd& beta, const bool& isBetaZero, const double& hval, const double& penalty,
                   const double& qtol, MatrixXd& logitm_Tmp, MatrixXd& beta_Tmp, bool& isBetaZero_Tmp, bool& needUpdate)
{

	int it1, rInd, cInd;
	double norm, mad;
	MatrixXd yMp(nh, rj), grad(di, rj), dif(di, rj), tmp(di, rj);
	VectorXd rowMax;


	/* the first iteration starts */

	grad = MatrixXd::Zero(di, rj);

	// compute the yMp matrix to prepare for gradient computation
	yMp = logitm;
	rowMax = yMp.rowwise().maxCoeff();
	yMp.colwise() -= rowMax;
	yMp = yMp.array().exp();
	for (it1 = 0; it1 < nh; ++it1)
	{
		yMp.row(it1) = ym.row(it1) - yMp.row(it1) / yMp.row(it1).sum();
	}

	// compute the gradient
	for (it1 = 0; it1 < nrows; ++it1)
	{
		rInd = nzIndex(it1, 0);
		cInd = nzIndex(it1, 1);
		grad.row(cInd) +=  yMp.row(rInd) * dm(rInd, cInd);
	}

	// update beta_Tmp
	if(isBetaZero)
	{
		norm = grad.norm();
		if (norm <= penalty)
		{
			needUpdate = false;
			return;
		}
		else
		{
			isBetaZero_Tmp = false;
			needUpdate = true;
			beta_Tmp = dif = ((penalty / norm - 1.0) / hval) * grad;
			mad = dif.lpNorm<Infinity>();
		}
	}
	else
	{
		grad -= hval * beta;
		norm = grad.norm();
		if (norm <= penalty)
		{
			isBetaZero_Tmp = true;
			needUpdate = true;
			dif = -beta;
			mad = dif.lpNorm<Infinity>();
			beta_Tmp = MatrixXd::Zero(di, rj);
		}
		else
		{
			tmp = ((penalty / norm - 1.0) / hval) * grad;
			dif = tmp - beta;
			mad = dif.lpNorm<Infinity>();
			if (mad < qtol)
			{
				needUpdate = false;
				return;
			}
			else
			{
				isBetaZero_Tmp = false;
				needUpdate = true;
				beta_Tmp = tmp;
			}
		}
	}

	// update the logitm_Tmp matrix
	logitm_Tmp = logitm;
	for (it1 = 0; it1 < nrows; ++it1)
	{
		rInd = nzIndex(it1, 0);
		cInd = nzIndex(it1, 1);
		logitm_Tmp.row(rInd) += dif.row(cInd) * dm(rInd, cInd);
	}

	// check convergence
	if (mad < qtol)	return;

}

void NewtonIterTmp_GD(const int& nh, const int& di, const int& rj, const int& nrows, const MatrixXd& dm, const MatrixXd& ym, const MatrixXi& nzIndex,
                   const MatrixXd& logitm, const MatrixXd& beta, const bool& isBetaZero, const double& hval, const double& penalty,
                   const double& qtol, MatrixXd& grad, MatrixXd& dif, MatrixXd& logitm_Tmp, MatrixXd& beta_Tmp, bool& isBetaZero_Tmp, bool& needUpdate)
{

	int it1, rInd, cInd;
	double norm, mad;
	MatrixXd yMp(nh, rj), tmp(di, rj);
	VectorXd rowMax;


	/* the first iteration starts */

	grad = MatrixXd::Zero(di, rj);

	// compute the yMp matrix to prepare for gradient computation
	yMp = logitm;
	rowMax = yMp.rowwise().maxCoeff();
	yMp.colwise() -= rowMax;
	yMp = yMp.array().exp();
	for (it1 = 0; it1 < nh; ++it1)
	{
		yMp.row(it1) = ym.row(it1) - yMp.row(it1) / yMp.row(it1).sum();
	}

	// compute the gradient
	for (it1 = 0; it1 < nrows; ++it1)
	{
		rInd = nzIndex(it1, 0);
		cInd = nzIndex(it1, 1);
		grad.row(cInd) +=  yMp.row(rInd) * dm(rInd, cInd);
	}

	// update beta_Tmp
	if(isBetaZero)
	{
		norm = grad.norm();
		if (norm <= penalty)
		{
			needUpdate = false;
			return;
		}
		else
		{
			isBetaZero_Tmp = false;
			needUpdate = true;
			beta_Tmp = dif = ((penalty / norm - 1.0) / hval) * grad;
			mad = dif.lpNorm<Infinity>();
		}
	}
	else
	{
		grad -= hval * beta;
		norm = grad.norm();
		if (norm <= penalty)
		{
			isBetaZero_Tmp = true;
			needUpdate = true;
			dif = -beta;
			mad = dif.lpNorm<Infinity>();
			beta_Tmp = MatrixXd::Zero(di, rj);
		}
		else
		{
			tmp = ((penalty / norm - 1.0) / hval) * grad;
			dif = tmp - beta;
			mad = dif.lpNorm<Infinity>();
			if (mad < qtol)
			{
				needUpdate = false;
				return;
			}
			else
			{
				isBetaZero_Tmp = false;
				needUpdate = true;
				beta_Tmp = tmp;
			}
		}
	}

	// update the logitm_Tmp matrix
	logitm_Tmp = logitm;
	for (it1 = 0; it1 < nrows; ++it1)
	{
		rInd = nzIndex(it1, 0);
		cInd = nzIndex(it1, 1);
		logitm_Tmp.row(rInd) += dif.row(cInd) * dm(rInd, cInd);
	}

	// check convergence
	if (mad < qtol)	return;

}

// Update the intercept, constrain beta_{j00} to zero
void InterceptUpdate(const int& nh, const MatrixXd& ym, MatrixXd& logitm, MatrixXd& beta, const double& hval, const double& qtol)
{

	int it1;
	double mad;
	MatrixXd yMp, dif;
	VectorXd rowMax;

    // compute the yMp matrix to prepare for updating the intercepts
    yMp = logitm;
    rowMax = yMp.rowwise().maxCoeff();
    yMp.colwise() -= rowMax;
    yMp = yMp.array().exp();
    for (it1 = 0; it1 < nh; ++it1)
    {
        yMp.row(it1) = ym.row(it1) - yMp.row(it1) / yMp.row(it1).sum();
    }

    // update the intercepts
    dif = yMp.colwise().sum() * (-1.0 / hval);
    dif(0) = 0.0;			// do not update beta_{j00}
    mad = dif.lpNorm<Infinity>();

    if (mad < qtol)	return;

    beta += dif;

    // update the logit matrix
    for (it1 = 0; it1 < nh; ++it1)
        logitm.row(it1) += dif;

}



void dmFetch(MatrixXd& dmt, MatrixXi& nzIndt, int& rCount, const MatrixXd& dm, const int& nrow, const int& ncol, const VectorXi& rIn, const VectorXi& cIn,
             const VectorXd& scaling)
{
	double value;
	rCount = -1;
	int colInd;
	for (int it1 = 0; it1 < ncol; ++it1)
	{
		colInd = cIn(it1);
		for (int it2 = 0; it2 < nrow; ++it2)
		{
			value = dm(rIn(it2), colInd);
			if (value != 0)
			{
				dmt(it2, it1) = value * scaling(it1);
				rCount += 1;
				nzIndt(rCount, 0) = it2;
				nzIndt(rCount, 1) = it1;
			}
			else
				dmt(it2, it1) = value;
		}
	}
	rCount += 1;
}



void firstDMFetch(MatrixXd& dmt, MatrixXi& nzIndt, int& rCount, const MatrixXd& dm, const int& nrow, const int& ncol, const VectorXi& rIn, const VectorXi& cIn,
				 VectorXd& scaling)
{
	double value;
	rCount = -1;
	int colInd;
	for (int it1 = 0; it1 < ncol; ++it1)
	{
		double nzc = 0.0;
		colInd = cIn(it1);
		for (int it2 = 0; it2 < nrow; ++it2)
		{
			dmt(it2, it1) = value = dm(rIn(it2), colInd);
			if (value != 0)
			{
				rCount += 1;
				nzIndt(rCount, 0) = it2;
				nzIndt(rCount, 1) = it1;
				nzc += 1.0;
			}
		}
		scaling(it1) = 1.0/sqrt(nzc);
		dmt.col(it1) *= scaling(it1);
	}
	rCount += 1;
}

// line search to update beta value.
void lineSearch(const MatrixXd& beta, MatrixXd& beta_tmp, const MatrixXd& logitm, MatrixXd& logitm_tmp,
                const VectorXi& ynzind, const double& penalty, const MatrixXi& nzIndex, const MatrixXd& dm,
                const MatrixXd& dif, const MatrixXd& grad, double& alpha, const double& delta,
                const double& sigma, const bool& isBetaZero, bool& isBetaZero_Tmp, const int& nh,
                const int& nrows, const int& di, const int& rj, bool& needUpdate)
{
    // definition
    VectorXi yInd;
    int l = 0; // l is the tracker of the number of iterations of line search
    int rInd, cInd;
    double pLog, St_old = 0.0, St_new = 0.0, big_delta = 0.0;
    MatrixXd yMp, logitm_T, beta_T;
    VectorXd rowMax, log_values;

    // calculate old score
    pLog = 0.0;
    yInd = ynzind;
    yMp = logitm;
    rowMax = yMp.rowwise().maxCoeff();
    yMp.colwise() -= rowMax;
    log_values = yMp.array().exp().rowwise().sum().log();
    pLog -= penalty * (beta).norm();
    for (int it2 = 0; it2 < nh; ++it2) {
        pLog -= log_values(it2) - yMp(it2, yInd(it2));
    }
    St_old = -pLog;

    // calculate big_delta
    big_delta = - (dif * grad.transpose()).rowwise().sum().sum() + penalty * (beta_tmp).norm() - penalty * (beta).norm();

    // loop to find alpha
    while (true) {
        logitm_T = logitm;

        // recalculate dif
        MatrixXd dif_T = alpha * dif;
        beta_T = beta + dif_T;

        // recalculate logitm
        for (int ii = 0; ii < nrows; ++ii) {
            rInd = nzIndex(ii, 0);
            cInd = nzIndex(ii, 1);
            logitm_T.row(rInd) += dif_T.row(cInd) * dm(rInd, cInd);
        }

        // calculate new score
        pLog = 0.0;
        yMp = logitm_T;
        rowMax = yMp.rowwise().maxCoeff();
        yMp.colwise() -= rowMax;
        log_values = yMp.array().exp().rowwise().sum().log();
        pLog -= penalty * (beta_T).norm();
        for (int it2 = 0; it2 < nh; ++it2) {
            pLog -= log_values(it2) -yMp(it2, yInd(it2));
        }
        St_new = -pLog;

        double error_tol = pow(10.0, -6);
        // decide if stopping criterion is met
        if ((St_new - St_old) <= (alpha * sigma * big_delta) || abs(alpha) < error_tol) {

            if (abs(alpha) >= error_tol) {
                if (!isBetaZero) {
                    // If beta was not zero, and line search suggest it should be zero.
                    if (beta_T.lpNorm<Infinity>() <= pow(10.0, -6)){
                        isBetaZero_Tmp = TRUE;
                        beta_T = MatrixXd::Zero(di, rj);
                        dif_T = beta_T - beta;
                        logitm_T = logitm;
                        for (int ii = 0; ii < nrows; ++ii) {
                            rInd = nzIndex(ii, 0);
                            cInd = nzIndex(ii, 1);
                            logitm_T.row(rInd) += dif_T.row(cInd) * dm(rInd, cInd);
                        }
                    }
                    else {
                        // if beta was not zero and line search agree
                        isBetaZero_Tmp = FALSE;
                    }
                }
                else {
                    // if beta was zero, but line search suggest to move.
                    isBetaZero_Tmp = FALSE;
                }
            }
            else if (abs(alpha) < error_tol) {
                // if line search suggest to stay where beta was.
                isBetaZero_Tmp = isBetaZero;
                beta_T = beta;
                logitm_T = logitm;
                needUpdate = FALSE;
            }
            beta_tmp = beta_T;
            logitm_tmp = logitm_T;
            break;
        }

        // for the next loop
        l = l + 1;
        alpha = alpha * delta;
    }
}

// linesearch added to beta update_by Jean
void innerLearning(const int& node, MatrixXi& G, const int& eor_nr, const MatrixXi& eor, vector<int>& active, double& MAD, const VectorXi& nobsVec,
                   const MatrixXi& ndfs, const VectorXMXd& dM, const VectorXMXd& yM, const VectorXVXi& yNZIndex, VectorXMXd& logitM, MatrixXMXd& betaM,
                   MatrixXb& IsBetaZeros, const MatrixXd& hvals, const MatrixXd& penalties, const double& qtol, const VectorXVXi& obsIndex,
                   const MatrixXVXi& levelIndex, const MatrixXVXd& scales, const int& maxRows, const int& maxCols)
{
    int rn, pn, rCount1, n1, d1, r1, it1;
    bool need_update1, isBetaZeroFlag1;
    MatrixXd dmt1(maxRows, maxCols),logitm_Tmp, logitm_Tmp1, beta_Tmp1, yMp;
    MatrixXi nzIndt1(maxRows * maxCols, 2);
    VectorXd rowMax, log_values;
    VectorXi yInd;
    MAD = 0.0;


    for (it1 = 0; it1 < eor_nr; ++it1){
        rn = eor(it1, 0)-1;
        pn = eor(it1, 1)-1;

        n1 = nobsVec(pn);
        d1 = ndfs(rn, pn);
        r1 = ndfs(pn, pn);

        dmFetch(dmt1, nzIndt1, rCount1, dM(rn), n1, d1, obsIndex(pn), levelIndex(rn, pn), scales(rn, pn));

        MatrixXd grad(d1, r1), dif(d1, r1);

        NewtonIterTmp_GD(n1, d1, r1, rCount1, dmt1, yM(pn), nzIndt1, logitM(pn), betaM(rn, pn), IsBetaZeros(rn, pn), hvals(rn, pn), penalties(rn, pn), qtol, grad, dif, logitm_Tmp1, beta_Tmp1, isBetaZeroFlag1, need_update1);

        if (need_update1) {

            // use line search to update beta and logit matrix.
            double alpha = 1.0;
            double delta = 0.5;
            double sigma = 0.1;

            lineSearch(betaM(rn, pn), beta_Tmp1, logitM(pn), logitm_Tmp1, yNZIndex(pn), penalties(rn, pn), nzIndt1, dmt1, dif, grad, alpha, delta, sigma, IsBetaZeros(rn, pn), isBetaZeroFlag1, n1, rCount1, d1, r1, need_update1);

            if (need_update1) {
                MAD = max(MAD, (beta_Tmp1-betaM(rn, pn)).lpNorm<Infinity>());
                betaM(rn, pn) = beta_Tmp1;
                logitM(pn) = logitm_Tmp1;
                IsBetaZeros(rn, pn) = isBetaZeroFlag1;

                auto beta_inter_Tmp = betaM(node, pn);
                InterceptUpdate(n1, yM(pn), logitM(pn), betaM(node, pn), hvals(node, pn), qtol);
                betaM(node, pn) = beta_inter_Tmp + alpha*(betaM(node, pn)-beta_inter_Tmp);
                MAD = max(MAD, (beta_inter_Tmp - betaM(node, pn)).lpNorm<Infinity>());

                if (IsBetaZeros(rn, pn)) {
                    G(rn, pn) = 0;
                }
                else {
                    G(rn, pn) = 1;
                }
            }
        }
    }
}

// CD algorithm, main part.
/*
 Note: certain restrictions apply to the data matrix.
 Let's define the jth problem as follows: using node j as the response and all other nodes as predictors (1 <= j <= p).
0. For any j, we enforce that the levels of node j are coded as 0, 1, 2, ..., . We also enforce that sample indices of the data matrix start from
   0 rather than 1.
1. For any j, the data for the jth problem (after excluding interventional samples of the jth node) should have enough samples such that every node
   has at least two levels, that is, at least one degrees of freedom.
2. We enforce using baseline (dummy variable) coding schemes.
3. ndfs: a p by p matrix. The jth column stores the number of levels for node j and the degrees of freedoms of all other nodes for the jth problem.
4. obsIndex:  a p by 1 col vector. Element j is a col vector of ordered observational sample indices for the jth problem.
5. levelIndex: a p by p matrix, where the (i, j)th element is a (col) vector of ordered levels of node i (excluding a baseline) that should be included in
   the design matrix for the jth problem.
6. yNZIndex: a p by 1 vector. Each element j is a vector whose hth element is the nonzero level for node j in sample h.
*/
void OneCDLoop(const int& node, MatrixXi& G, const int& eor_nr, const MatrixXi& eor, vector<int>& active, double& MAD, const VectorXi& nobsVec,
                          const MatrixXi& ndfs, const VectorXMXd& dM, const VectorXMXd& yM, const VectorXVXi& yNZIndex, VectorXMXd& logitM, MatrixXMXd& betaM,
                          MatrixXb& IsBetaZeros, const MatrixXd& hvals, const MatrixXd& penalties, const double& qtol, const VectorXVXi& obsIndex,
                          const MatrixXVXi& levelIndex, const MatrixXVXd& scales, const int& maxRows, const int& maxCols)
{
    /* betaM(node, j) and hvals(node, j) stores the intercept terms for the jth problem. */

    int rn, pn, rCount1, n1, d1, r1,
    it1, it2, rInd, cInd;
    double pLog_rn, pLog_pn, alpha, alpha1, alpha2;
    bool need_update1, isBetaZeroFlag1, need_update2, isBetaZeroFlag2;
    MatrixXd dmt1(maxRows, maxCols), logitm_Tmp, logitm_Tmp1, beta_Tmp1,
    logitm_Tmp2, beta_Tmp2, yMp;
    MatrixXi nzIndt1(maxRows * maxCols, 2);
    VectorXd rowMax, log_values;
    VectorXi yInd;
    MAD = 0.0;
    double delta = 0.5;
    double sigma = 0.1;

    for (it1 = 0; it1 < eor_nr; ++it1)
    {

        rn = eor(it1, 0) - 1;
        pn = eor(it1, 1) - 1;
        G(rn, pn) = G(pn, rn) = 0;

        if(Cycle(node, &G(0, 0), rn + 1, pn + 1))
        {
            n1 = nobsVec(rn);	d1 = ndfs(pn, rn);	r1 = ndfs(rn, rn);

            // create a local copy of the design matrix of node pn
            dmFetch(dmt1, nzIndt1, rCount1, dM(pn), n1, d1, obsIndex(rn), levelIndex(pn, rn), scales(pn, rn));

            // update beta
            MatrixXd grad(d1, r1), dif(d1, r1);
            NewtonIterTmp_GD(n1, d1, r1, rCount1, dmt1, yM(rn), nzIndt1, logitM(rn), betaM(pn, rn), IsBetaZeros(pn, rn), hvals(pn, rn), penalties(pn, rn), qtol, grad, dif, logitm_Tmp1, beta_Tmp1, isBetaZeroFlag1, need_update1);

            if (need_update1) {

                alpha = 1.0;
                // use line search to update beta and logit matrix.

                lineSearch(betaM(pn, rn), beta_Tmp1, logitM(rn), logitm_Tmp1, yNZIndex(rn), penalties(pn, rn), nzIndt1, dmt1, dif, grad, alpha, delta, sigma, IsBetaZeros(pn, rn), isBetaZeroFlag1, n1, rCount1, d1, r1, need_update1);

                if (need_update1) {
                    MAD = max(MAD, (beta_Tmp1-betaM(pn, rn)).lpNorm<Infinity>());
                    betaM(pn, rn) = beta_Tmp1;
                    logitM(rn) = logitm_Tmp1;
                    IsBetaZeros(pn, rn) = isBetaZeroFlag1;

                    //                    update intercept
                    auto beta_inter_Tmp = betaM(node, rn);
                    InterceptUpdate(n1, yM(rn), logitM(rn), betaM(node, rn), hvals(node, rn), qtol);
                    betaM(node, rn) = beta_inter_Tmp + alpha*(betaM(node, rn)-beta_inter_Tmp);
                    MAD = max(MAD, (beta_inter_Tmp - betaM(node, rn)).lpNorm<Infinity>());

                }
            }

            // update the graph and the active set
            if (!IsBetaZeros(pn, rn))
            {
                G(pn, rn) = 1;
                active.push_back(it1);
            }
        }

        else if (Cycle(node, &G(0, 0), pn + 1, rn + 1))
        {
            n1 = nobsVec(pn);	d1 = ndfs(rn, pn);	r1 = ndfs(pn, pn);

            // create a local copy of the design matrix of node rn
            dmFetch(dmt1, nzIndt1, rCount1, dM(rn), n1, d1, obsIndex(pn), levelIndex(rn, pn), scales(rn, pn));

            // update beta
            MatrixXd grad(d1, r1), dif(d1, r1);
            NewtonIterTmp_GD(n1, d1, r1, rCount1, dmt1, yM(pn), nzIndt1, logitM(pn), betaM(rn, pn), IsBetaZeros(rn, pn), hvals(rn, pn), penalties(rn, pn), qtol, grad, dif, logitm_Tmp1, beta_Tmp1, isBetaZeroFlag1, need_update1);

            if (need_update1) {

                alpha = 1.0;
                // use line search to update beta and logit matrix.

                lineSearch(betaM(rn, pn), beta_Tmp1, logitM(pn), logitm_Tmp1, yNZIndex(pn), penalties(rn, pn), nzIndt1, dmt1, dif, grad, alpha, delta, sigma, IsBetaZeros(rn, pn), isBetaZeroFlag1, n1, rCount1, d1, r1, need_update1);

                if (need_update1) {
                    MAD = max(MAD, (beta_Tmp1-betaM(rn, pn)).lpNorm<Infinity>());
                    betaM(rn, pn) = beta_Tmp1;
                    logitM(pn) = logitm_Tmp1;
                    IsBetaZeros(rn, pn) = isBetaZeroFlag1;

                    // update intercept
                    auto beta_inter_Tmp = betaM(node, pn);
                    InterceptUpdate(n1, yM(pn), logitM(pn), betaM(node, pn), hvals(node, pn), qtol);
                    betaM(node, pn) = beta_inter_Tmp + alpha*(betaM(node, pn)-beta_inter_Tmp);
                    MAD = max(MAD, (beta_inter_Tmp - betaM(node, pn)).lpNorm<Infinity>());

                }
            }

            // update the graph and the active set
            if (!IsBetaZeros(rn, pn))
            {
                G(rn, pn) = 1;
                active.push_back(it1);
            }
        }

        else
        {
            if (!IsBetaZeros(pn, rn))
            {
                // consider pn to rn
                n1 = nobsVec(rn);	d1 = ndfs(pn, rn);	r1 = ndfs(rn, rn);
                dmFetch(dmt1, nzIndt1, rCount1, dM(pn), n1, d1, obsIndex(rn), levelIndex(pn, rn), scales(pn, rn));
                MatrixXd grad1(d1, r1), dif1(d1, r1);
                NewtonIterTmp_GD(n1, d1, r1, rCount1, dmt1, yM(rn), nzIndt1, logitM(rn), betaM(pn, rn), IsBetaZeros(pn, rn), hvals(pn, rn), penalties(pn, rn), qtol, grad1, dif1, logitm_Tmp1, beta_Tmp1, isBetaZeroFlag1, need_update1);

                // line search
                if (need_update1) {

                    alpha1 = 1.0;
                    // use line search to update beta and logit matrix.

                    lineSearch(betaM(pn, rn), beta_Tmp1, logitM(rn), logitm_Tmp1, yNZIndex(rn), penalties(pn, rn), nzIndt1, dmt1, dif1, grad1, alpha1, delta, sigma, IsBetaZeros(pn, rn), isBetaZeroFlag1, n1, rCount1, d1, r1, need_update1);
                }

                if (!need_update1)
                {
                    // compute the rn^th penalized negative log-likelihood assuming that all beta for node pn are zero
                    logitm_Tmp = logitM(rn);
                    for (it2 = 0; it2 < rCount1; ++it2)
                    {
                        rInd = nzIndt1(it2, 0);
                        cInd = nzIndt1(it2, 1);
                        logitm_Tmp.row(rInd) -= betaM(pn, rn).row(cInd) * dmt1(rInd, cInd);
                    }
                    yInd = yNZIndex(rn);
                    yMp = logitm_Tmp;
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_rn = 0.0;
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_rn += log_values(it2) - yMp(it2, yInd(it2));
                    }

                    // compute the difference in the rn^th penalized negative log-likelihood
                    yMp = logitM(rn);
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_rn -= penalties(pn, rn) * (betaM(pn, rn)).norm();
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_rn -= log_values(it2) - yMp(it2, yInd(it2));
                    }
                }
                else if (!isBetaZeroFlag1)
                {
                    // compute the rn^th penalized negative log-likelihood assuming that all beta for node pn are zero
                    logitm_Tmp = logitM(rn);
                    for (it2 = 0; it2 < rCount1; ++it2)
                    {
                        rInd = nzIndt1(it2, 0);
                        cInd = nzIndt1(it2, 1);
                        logitm_Tmp.row(rInd) -= betaM(pn, rn).row(cInd) * dmt1(rInd, cInd);
                    }
                    yInd = yNZIndex(rn);
                    yMp = logitm_Tmp;
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_rn = 0.0;
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_rn += log_values(it2) - yMp(it2, yInd(it2));
                    }

                    // compute the difference in the rn^th penalized negative log-likelihood
                    yMp = logitm_Tmp1;
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_rn -= penalties(pn, rn) * (beta_Tmp1).norm();
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_rn -= log_values(it2) - yMp(it2, yInd(it2));
                    }
                }
                else
                {
                    logitm_Tmp = logitm_Tmp1;
                    pLog_rn = 0.0;
                }


                // consider rn to pn
                n1 = nobsVec(pn);	d1 = ndfs(rn, pn);	r1 = ndfs(pn, pn);
                dmFetch(dmt1, nzIndt1, rCount1, dM(rn), n1, d1, obsIndex(pn), levelIndex(rn, pn), scales(rn, pn));
                MatrixXd grad2(d1, r1), dif2(d1, r1);
                NewtonIterTmp_GD(n1, d1, r1, rCount1, dmt1, yM(pn), nzIndt1, logitM(pn), betaM(rn, pn), IsBetaZeros(rn, pn), hvals(rn, pn),
                                 penalties(rn, pn), qtol, grad2, dif2, logitm_Tmp2, beta_Tmp2, isBetaZeroFlag2, need_update2);

                // line search
                if (need_update2) {

                    alpha2 = 1.0;
                    // use line search to update beta and logit matrix.

                    lineSearch(betaM(rn, pn), beta_Tmp2, logitM(pn), logitm_Tmp2, yNZIndex(pn), penalties(rn, pn), nzIndt1, dmt1, dif2, grad2, alpha2, delta, sigma, IsBetaZeros(rn, pn), isBetaZeroFlag2, n1, rCount1, d1, r1, need_update2);
                }

                if (!need_update2)
                {
                    pLog_pn = 0.0;
                }
                else
                {
                    // compute the pn^th penalized negative log-likelihood assuming that all beta for node rn are zero
                    yInd = yNZIndex(pn);
                    yMp = logitM(pn);
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_pn = 0.0;
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_pn += log_values(it2) - yMp(it2, yInd(it2));
                    }

                    // compute the difference in the pn^th penalized negative log-likelihood
                    yMp = logitm_Tmp2;
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_pn -= penalties(rn, pn) * (beta_Tmp2).norm();
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_pn -= log_values(it2) - yMp(it2, yInd(it2));
                    }
                }

                if (pLog_pn <= pLog_rn)
                {
                    if (need_update1)
                    {
                        MAD = max(MAD, (beta_Tmp1 - betaM(pn, rn)).lpNorm<Infinity>());
                        logitM(rn) = logitm_Tmp1;
                        betaM(pn, rn) = beta_Tmp1;
                        IsBetaZeros(pn, rn) = isBetaZeroFlag1;

                        beta_Tmp1 = betaM(node, rn);
                        InterceptUpdate(nobsVec(rn), yM(rn), logitM(rn), betaM(node, rn), hvals(node, rn), qtol);
                        betaM(node, rn) = beta_Tmp1 + alpha1*(betaM(node, rn)-beta_Tmp1);
                        MAD = max(MAD, (beta_Tmp1 - betaM(node, rn)).lpNorm<Infinity>());
                    }
                    if (!IsBetaZeros(pn, rn))
                    {
                        G(pn, rn) = 1;
                        active.push_back(it1);
                    }
                }
                else
                {
                    MAD = max(MAD, betaM(pn, rn).lpNorm<Infinity>());
                    logitM(rn) = logitm_Tmp;
                    betaM(pn, rn) = MatrixXd::Zero(ndfs(pn, rn), ndfs(rn, rn));
                    IsBetaZeros(pn, rn) = true;

                    beta_Tmp1 = betaM(node, rn);
                    InterceptUpdate(nobsVec(rn), yM(rn), logitM(rn), betaM(node, rn), hvals(node, rn), qtol);
                    MAD = max(MAD, (beta_Tmp1 - betaM(node, rn)).lpNorm<Infinity>());

                    if (need_update2)
                    {
                        MAD = max(MAD, (beta_Tmp2 - betaM(rn, pn)).lpNorm<Infinity>());
                        logitM(pn) = logitm_Tmp2;
                        betaM(rn, pn) = beta_Tmp2;
                        IsBetaZeros(rn, pn) = isBetaZeroFlag2;

                        beta_Tmp2 = betaM(node, pn);
                        InterceptUpdate(nobsVec(pn), yM(pn), logitM(pn), betaM(node, pn), hvals(node, pn), qtol);
                        betaM(node, pn) = beta_Tmp2 + alpha2*(betaM(node, pn)-beta_Tmp2);
                        MAD = max(MAD, (beta_Tmp2 - betaM(node, pn)).lpNorm<Infinity>());
                    }
                    if (!IsBetaZeros(rn, pn))
                    {
                        G(rn, pn) = 1;
                        active.push_back(it1);
                    }
                }
            }
            else if (!IsBetaZeros(rn, pn))
            {
                // consider pn to rn
                n1 = nobsVec(rn);	d1 = ndfs(pn, rn);	r1 = ndfs(rn, rn);
                dmFetch(dmt1, nzIndt1, rCount1, dM(pn), n1, d1, obsIndex(rn), levelIndex(pn, rn), scales(pn, rn));
                MatrixXd grad1(d1, r1), dif1(d1, r1);
                NewtonIterTmp_GD(n1, d1, r1, rCount1, dmt1, yM(rn), nzIndt1, logitM(rn), betaM(pn, rn), IsBetaZeros(pn, rn), hvals(pn, rn),
                                 penalties(pn, rn), qtol, grad1, dif1, logitm_Tmp1, beta_Tmp1, isBetaZeroFlag1, need_update1);

                // line search
                if (need_update1) {

                    alpha1 = 1.0;
                    // use line search to update beta and logit matrix.

                    lineSearch(betaM(pn, rn), beta_Tmp1, logitM(rn), logitm_Tmp1, yNZIndex(rn), penalties(pn, rn), nzIndt1, dmt1, dif1, grad1, alpha1, delta, sigma, IsBetaZeros(pn, rn), isBetaZeroFlag1, n1, rCount1, d1, r1, need_update1);
                }


                if (!need_update1)
                {
                    pLog_rn = 0.0;
                }
                else
                {
                    // compute the rn^th penalized negative log-likelihood assuming that all beta for node pn are zero
                    yInd = yNZIndex(rn);
                    yMp = logitM(rn);
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_rn = 0.0;
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_rn += log_values(it2) - yMp(it2, yInd(it2));
                    }

                    // compute the difference in the rn^th penalized negative log-likelihood
                    yMp = logitm_Tmp1;
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_rn -= penalties(pn, rn) * (beta_Tmp1).norm();
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_rn -= log_values(it2) - yMp(it2, yInd(it2));
                    }
                }

                // consider rn to pn
                n1 = nobsVec(pn);	d1 = ndfs(rn, pn);	r1 = ndfs(pn, pn);
                dmFetch(dmt1, nzIndt1, rCount1, dM(rn), n1, d1, obsIndex(pn), levelIndex(rn, pn), scales(rn, pn));
                MatrixXd grad2(d1, r1), dif2(d1, r1);
                NewtonIterTmp_GD(n1, d1, r1, rCount1, dmt1, yM(pn), nzIndt1, logitM(pn), betaM(rn, pn), IsBetaZeros(rn, pn), hvals(rn, pn),
                                 penalties(rn, pn), qtol, grad2, dif2, logitm_Tmp2, beta_Tmp2, isBetaZeroFlag2, need_update2);

                // line search
                if (need_update2) {

                    alpha2 = 1.0;
                    // use line search to update beta and logit matrix.

                    lineSearch(betaM(rn, pn), beta_Tmp2, logitM(pn), logitm_Tmp2, yNZIndex(pn), penalties(rn, pn), nzIndt1, dmt1, dif2, grad2, alpha2, delta, sigma, IsBetaZeros(rn, pn), isBetaZeroFlag2, n1, rCount1, d1, r1, need_update2);
                }

                if (!need_update2)
                {
                    // compute the pn^th penalized negative log-likelihood assuming that all beta for node rn are zero
                    logitm_Tmp = logitM(pn);
                    for (it2 = 0; it2 < rCount1; ++it2)
                    {
                        rInd = nzIndt1(it2, 0);
                        cInd = nzIndt1(it2, 1);
                        logitm_Tmp.row(rInd) -= betaM(rn, pn).row(cInd) * dmt1(rInd, cInd);
                    }
                    yInd = yNZIndex(pn);
                    yMp = logitm_Tmp;
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_pn = 0.0;
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_pn += log_values(it2) - yMp(it2, yInd(it2));
                    }

                    // compute the difference in the pn^th penalized negative log-likelihood
                    yMp = logitM(pn);
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_pn -= penalties(rn, pn) * (betaM(rn, pn)).norm();
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_pn -= log_values(it2) - yMp(it2, yInd(it2));
                    }
                }
                else if (!isBetaZeroFlag2)
                {
                    // compute the pn^th penalized negative log-likelihood assuming that all beta for node rn are zero
                    logitm_Tmp = logitM(pn);
                    for (it2 = 0; it2 < rCount1; ++it2)
                    {
                        rInd = nzIndt1(it2, 0);
                        cInd = nzIndt1(it2, 1);
                        logitm_Tmp.row(rInd) -= betaM(rn, pn).row(cInd) * dmt1(rInd, cInd);
                    }
                    yInd = yNZIndex(pn);
                    yMp = logitm_Tmp;
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_pn = 0.0;
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_pn += log_values(it2) - yMp(it2, yInd(it2));
                    }

                    // compute the difference in the pn^th penalized negative log-likelihood
                    yMp = logitm_Tmp2;
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_pn -= penalties(rn, pn) * (beta_Tmp2).norm();
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_pn -= log_values(it2) - yMp(it2, yInd(it2));
                    }
                }
                else
                {
                    logitm_Tmp = logitm_Tmp2;
                    pLog_pn = 0.0;
                }

                if (pLog_rn <= pLog_pn)
                {
                    if (need_update2)
                    {
                        MAD = max(MAD, (beta_Tmp2 - betaM(rn, pn)).lpNorm<Infinity>());
                        logitM(pn) = logitm_Tmp2;
                        betaM(rn, pn) = beta_Tmp2;
                        IsBetaZeros(rn, pn) = isBetaZeroFlag2;

                        beta_Tmp2 = betaM(node, pn);
                        InterceptUpdate(nobsVec(pn), yM(pn), logitM(pn), betaM(node, pn), hvals(node, pn), qtol);
                        betaM(node, pn) = beta_Tmp2 + alpha2*(betaM(node, pn)-beta_Tmp2);
                        MAD = max(MAD, (beta_Tmp2 - betaM(node, pn)).lpNorm<Infinity>());
                    }
                    if (!IsBetaZeros(rn, pn))
                    {
                        G(rn, pn) = 1;
                        active.push_back(it1);
                    }
                }
                else
                {
                    MAD = max(MAD, betaM(rn, pn).lpNorm<Infinity>());
                    logitM(pn) = logitm_Tmp;
                    betaM(rn, pn) = MatrixXd::Zero(ndfs(rn, pn), ndfs(pn, pn));
                    IsBetaZeros(rn, pn) = true;

                    beta_Tmp2 = betaM(node, pn);
                    InterceptUpdate(nobsVec(pn), yM(pn), logitM(pn), betaM(node, pn), hvals(node, pn), qtol);
                    MAD = max(MAD, (beta_Tmp2 - betaM(node, pn)).lpNorm<Infinity>());

                    if (need_update1)
                    {
                        MAD = max(MAD, (beta_Tmp1 - betaM(pn, rn)).lpNorm<Infinity>());
                        logitM(rn) = logitm_Tmp1;
                        betaM(pn, rn) = beta_Tmp1;
                        IsBetaZeros(pn, rn) = isBetaZeroFlag1;

                        beta_Tmp1 = betaM(node, rn);
                        InterceptUpdate(nobsVec(rn), yM(rn), logitM(rn), betaM(node, rn), hvals(node, rn), qtol);
                        betaM(node, rn) = beta_Tmp1 + alpha1*(betaM(node, rn)-beta_Tmp1);
                        MAD = max(MAD, (beta_Tmp1 - betaM(node, rn)).lpNorm<Infinity>());
                    }
                    if (!IsBetaZeros(pn, rn))
                    {
                        G(pn, rn) = 1;
                        active.push_back(it1);
                    }
                }
            }
            else
            {
                // consider pn to rn
                n1 = nobsVec(rn);	d1 = ndfs(pn, rn);	r1 = ndfs(rn, rn);
                dmFetch(dmt1, nzIndt1, rCount1, dM(pn), n1, d1, obsIndex(rn), levelIndex(pn, rn), scales(pn, rn));
                MatrixXd grad1(d1, r1), dif1(d1, r1);
                NewtonIterTmp_GD(n1, d1, r1, rCount1, dmt1, yM(rn), nzIndt1, logitM(rn), betaM(pn, rn), IsBetaZeros(pn, rn), hvals(pn, rn),
                                 penalties(pn, rn), qtol, grad1, dif1, logitm_Tmp1, beta_Tmp1, isBetaZeroFlag1, need_update1);

                // line search
                if (need_update1) {

                    alpha1 = 1.0;
                    // use line search to update beta and logit matrix.

                    lineSearch(betaM(pn, rn), beta_Tmp1, logitM(rn), logitm_Tmp1, yNZIndex(rn), penalties(pn, rn), nzIndt1, dmt1, dif1, grad1, alpha1, delta, sigma, IsBetaZeros(pn, rn), isBetaZeroFlag1, n1, rCount1, d1, r1, need_update1);
                }

                if (!need_update1)
                {
                    pLog_rn = 0.0;
                }
                else
                {
                    // compute the rn^th penalized negative log-likelihood assuming that all beta for node pn are zero
                    yInd = yNZIndex(rn);
                    yMp = logitM(rn);
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_rn = 0.0;
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_rn += log_values(it2) - yMp(it2, yInd(it2));
                    }

                    // compute the difference in the rn^th penalized negative log-likelihood
                    yMp = logitm_Tmp1;
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_rn -= penalties(pn, rn) * (beta_Tmp1).norm();
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_rn -= log_values(it2) - yMp(it2, yInd(it2));
                    }
                }

                // consider rn to pn
                n1 = nobsVec(pn);	d1 = ndfs(rn, pn);	r1 = ndfs(pn, pn);
                dmFetch(dmt1, nzIndt1, rCount1, dM(rn), n1, d1, obsIndex(pn), levelIndex(rn, pn), scales(rn, pn));
                MatrixXd grad2(d1, r1), dif2(d1, r1);
                NewtonIterTmp_GD(n1, d1, r1, rCount1, dmt1, yM(pn), nzIndt1, logitM(pn), betaM(rn, pn), IsBetaZeros(rn, pn), hvals(rn, pn),
                                 penalties(rn, pn), qtol, grad2, dif2, logitm_Tmp2, beta_Tmp2, isBetaZeroFlag2, need_update2);

                // line search
                if (need_update2) {

                    alpha2 = 1.0;
                    // use line search to update beta and logit matrix.

                    lineSearch(betaM(rn, pn), beta_Tmp2, logitM(pn), logitm_Tmp2, yNZIndex(pn), penalties(rn, pn), nzIndt1, dmt1, dif2, grad2, alpha2, delta, sigma, IsBetaZeros(rn, pn), isBetaZeroFlag2, n1, rCount1, d1, r1, need_update2);
                }

                if (!need_update2)
                {
                    pLog_pn = 0.0;
                }
                else
                {
                    // compute the pn^th penalized negative log-likelihood assuming that all beta for node rn are zero
                    yInd = yNZIndex(pn);
                    yMp = logitM(pn);
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_pn = 0.0;
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_pn += log_values(it2) - yMp(it2, yInd(it2));
                    }

                    // compute the difference in the pn^th penalized negative log-likelihood
                    yMp = logitm_Tmp2;
                    rowMax = yMp.rowwise().maxCoeff();
                    yMp.colwise() -= rowMax;
                    log_values = yMp.array().exp().rowwise().sum().log();
                    pLog_pn -= penalties(rn, pn) * (beta_Tmp2).norm();
                    for (it2 = 0; it2 < n1; ++it2)
                    {
                        pLog_pn -= log_values(it2) - yMp(it2, yInd(it2));
                    }
                }

                if (pLog_pn <= pLog_rn)
                {
                    if (need_update1)
                    {
                        MAD = max(MAD, (beta_Tmp1 - betaM(pn, rn)).lpNorm<Infinity>());
                        logitM(rn) = logitm_Tmp1;
                        betaM(pn, rn) = beta_Tmp1;
                        IsBetaZeros(pn, rn) = isBetaZeroFlag1;

                        beta_Tmp1 = betaM(node, rn);
                        InterceptUpdate(nobsVec(rn), yM(rn), logitM(rn), betaM(node, rn), hvals(node, rn), qtol);
                        betaM(node, rn) = beta_Tmp1 + alpha1*(betaM(node, rn)-beta_Tmp1);
                        MAD = max(MAD, (beta_Tmp1 - betaM(node, rn)).lpNorm<Infinity>());
                    }
                    if (!IsBetaZeros(pn, rn))
                    {
                        G(pn, rn) = 1;
                        active.push_back(it1);
                    }
                }
                else
                {
                    if (need_update2)
                    {
                        MAD = max(MAD, (beta_Tmp2 - betaM(rn, pn)).lpNorm<Infinity>());
                        logitM(pn) = logitm_Tmp2;
                        betaM(rn, pn) = beta_Tmp2;
                        IsBetaZeros(rn, pn) = isBetaZeroFlag2;

                        beta_Tmp2 = betaM(node, pn);
                        InterceptUpdate(nobsVec(pn), yM(pn), logitM(pn), betaM(node, pn), hvals(node, pn), qtol);
                        betaM(node, pn) = beta_Tmp2 + alpha2*(betaM(node, pn)-beta_Tmp2);
                        MAD = max(MAD, (beta_Tmp2 - betaM(node, pn)).lpNorm<Infinity>());
                    }
                    if (!IsBetaZeros(rn, pn))
                    {
                        G(rn, pn) = 1;
                        active.push_back(it1);
                    }
                }
            }
        }
    }
}



void CDOnePoint(const int& node, MatrixXi& G, const int& eor_nr, const MatrixXi& eor, const double& eps, const VectorXi& nobsVec,
                const MatrixXi& ndfs, const VectorXMXd& dM, const VectorXMXd& yM, const VectorXVXi& yNZIndex, VectorXMXd& logitM, MatrixXMXd& betaM,
                MatrixXb& IsBetaZeros, const MatrixXd& hvals, const MatrixXd& penalties, const double& qtol, const VectorXVXi& obsIndex,
                const MatrixXVXi& levelIndex, const MatrixXVXd& scales, const int& maxRows, const int& maxCols)
{
    double MAD;
    unsigned int num;
    unsigned int times1, times2;
    times1 = 0;

    MatrixXi activeSet_old = eor; // Store the previous active set. And compare with the new active set, if the two active sets are the same, stop the algorithm.

    while(true)
    {
        vector<int> active;
        vector<int> numbers;

		OneCDLoop(node, G, eor_nr, eor, active, MAD, nobsVec, ndfs, dM, yM, yNZIndex, logitM, betaM,
				 IsBetaZeros, hvals, penalties, qtol, obsIndex, levelIndex, scales, maxRows, maxCols);

        times1++;

		num = active.size();
		MatrixXi activeSet(num, 2);

        int num_cnt = 0;

        for (int row = 0; row < node; row++) {
            for (int col = 0; col < node; col++) {
                if (G(row, col)) {
                    activeSet(num_cnt, 0) = row+1;
                    activeSet(num_cnt, 1) = col+1;
                    num_cnt++;
                }
            }
        }

        if (activeSet_old.size()/2 == num) {
            if ((activeSet - activeSet_old).norm() == 0) {
                // break if one more outer loop will not change the active set.
                // cout << "has " << times1 << "outer part" << endl;
                break;
            }
        }

        if(MAD < eps)	{
            // break if coefficient converges
            // cout << "has " << times1 << "outer part" << endl;
            break;
        }

        // inner loop
        // Algorithm 2
        times2 = 0;

		while(true)
		{
			vector<int> sactive;
			innerLearning(node, G, num, activeSet, sactive, MAD, nobsVec, ndfs, dM, yM, yNZIndex, logitM, betaM,
					  IsBetaZeros, hvals, penalties, qtol, obsIndex, levelIndex, scales, maxRows, maxCols);

            times2++;
			if(MAD < eps || times2 > 1000) {
                // throw out warning if inner loop reaches maximum number of iterations.
                if (times2 > 1000) {
                    Rprintf("\n");
                    Rprintf("the %d th inner iteration reaches maximum. See Adaptive Penalized Estimation of Directed Acyclic Graphs From Categorical Data (http://arxiv.org/abs/1403.2310). Chapter 3.2 Algorithm 2 for help. \n", times1);
                }
                break;
            }
		}

        // break if the outer loop reaches maximum number of iterations
        if(times1 > 20)	{
            Rprintf("\n");
            Rprintf("the outer iteration reaches maximum. See Adaptive Penalized Estimation of Directed Acyclic Graphs From Categorical Data (http://arxiv.org/abs/1403.2310). Chapter 3.2 Algorithm 1 for help. \n");
            break;
        }
    }
}

/*
Note: certain restrictions apply to the data matrix.
Let's define the jth problem as follows: using node j as the response and all other nodes as predictors (1 <= j <= p).
0. For any j, we enforce that the levels of node j are coded as 0, 1, 2, ..., . We also enforce that sample indices of the data matrix start from
   0 rather than 1.
1. For any j, the data for the jth problem (after excluding interventional samples of the jth node) should have enough samples such that every node
   has at least two levels, that is, at least one degrees of freedom.
2. We enforce using baseline (dummy variable) coding schemes.
3. ndfs: a p by p matrix. The jth column stores the number of levels for node j and the degrees of freedoms of all other nodes for the jth problem.
4. obsIndex:  a p by 1 col vector. Element j is a col vector of ordered observational sample indices for the jth problem.
5. levelIndex: a p by p matrix, where the (i, j)th element is a (col) vector of ordered levels of node i (excluding a baseline) that should be included in
   the design matrix for the jth problem.
6. yNZIndex: a p by 1 vector. Each element j is a vector whose hth element is the nonzero level for node j in sample h.
7. upperbound: a large positive value used to truncate the adaptive weights. A -1 value indicates that there is no truncation.
*/
void CDAlgo(int node, int dataSize, const MatrixXi& data, const VectorXi& nlevels, const VectorXVXi& obsIndex, const MatrixXVXi& levelIndex, int eor_nr, const MatrixXi& eor, int nlam, double eps, double convLb, double qtol, VectorXd& lambdaSeq, VectorXd& log_like, VectorXd& dur, MatrixXMXd& betaM, MatrixXMXd& betaN, MatrixXi& estimateG, MatrixXd& weights, double gamma, double upperbound, int threshold)
{
	int maxRows = dataSize, maxCols = nlevels.maxCoeff(); //nlevels record number of levels for each node.
	if(upperbound >= 0.0) // to calculate weights
	{
		MatrixXd umtr = MatrixXd::Constant(node, node, pow(upperbound, gamma)); // node by node matrix with elements all equal to upperbound to the power of gamma
		weights = weights.array().pow(gamma).min(umtr.array()); //
	}
	else if(upperbound == -1.0) // no truncation
		weights = weights.array().pow(gamma);
	else
		Rprintf("upperbound sould be positive!");

	// create the full design matrix dM for each node (including the baseline column which is assumed to be the last factor level)
	// initialize all relevant information
	VectorXi nobsVec(node);
	MatrixXi ndfs(node, node); // p by p matrix. The jth column stores the number of levels for node j and the degrees of freedoms of all other nodes for the jth problem.
	MatrixXd penalties(node, node);
	VectorXMXd dM(node), yM(node);
	VectorXVXi yNZIndex(node); //p by 1 vector. Each element j is a vector whose hth element is the nonzero level for node j in sample h
	MatrixXVXd scales(node, node);
	for (int it1 = 0; it1 < node; ++it1) // it1 -> i
	{
		int size;
		nobsVec(it1) = obsIndex(it1).size(); // ith element is the number of observations for the ith problem
		for (int it2 = 0; it2 < node; ++it2) // it2 -> j
		{
			ndfs(it2, it1) = size = levelIndex(it2, it1).size();
			scales(it2, it1) = VectorXd::Zero(size);
		}
		dM(it1) = MatrixXd::Zero(dataSize, nlevels(it1));// design matrix for each node, the ith element is the matrix that indicate which level the ith node have
		yM(it1) = MatrixXd::Zero(nobsVec(it1), ndfs(it1, it1)); // ith element is a matrix record that rows are for the observations of ith problem, colums are for levels of ith node.
		yNZIndex(it1) = VectorXi::Zero(nobsVec(it1));
		penalties.col(it1) = (ndfs.col(it1) * ndfs(it1, it1)).cast<double>().cwiseSqrt();
	}
	penalties = penalties.array() * weights.array();	// Adaptive Lasso weights
	for (int it1 = 0; it1 < dataSize; ++it1) // i to determine design matrix
	{
		for (int it2 = 0; it2 < node; ++it2) // j
			dM(it2)(it1, data(it1, it2)) = 1.0;
	}
	for (int it1 = 0; it1 < node; ++it1) // to
	{
		int nrows = nobsVec(it1), ncols = ndfs(it1, it1); // number of observations fro the ith problem, number of levels for the ith node
		double value;
		const VectorXi& oind = obsIndex(it1); //ordered observational sample indices for the ith problem
		const VectorXi& lind = levelIndex(it1, it1); // ordered levels of node i, that should be included in the design matrix fot the ith problem
		const MatrixXd& dm = dM(it1);
		for (int it2 = 0; it2 < nrows; ++it2)
		{
			for (int it3 = 0; it3 < ncols; ++it3)
			{
				value = dm(oind(it2), lind(it3));
				yM(it1)(it2, it3) = value;
				if (value != 0)
					yNZIndex(it1)(it2) = it3;
			}
		}
	}



	// compute the maximum lambda and sequence of lambdas
	// initialize betaM
	// double lambdaMax = 0.0;
	VectorXMXd logitM(node);
	MatrixXb IsBetaZeros(node, node);
	for (int j = 0; j < node; ++j)
	{
		MatrixXd yMp = yM(j);
		RowVectorXd colMean = yMp.colwise().mean();
		yMp.rowwise() -= colMean;

		int rCount, di, rj = ndfs(j, j), rInd, cInd;
		MatrixXd dmt(maxRows, maxCols);
		MatrixXi nzIndt(maxRows * maxCols, 2);
		for (int i = 0; i < node; ++i)
		{
			if(i != j)
			{
				di = ndfs(i, j);
				betaM(i, j) = MatrixXd::Zero(di, rj);
				IsBetaZeros(i, j) = true;
				firstDMFetch(dmt, nzIndt, rCount, dM(i), nobsVec(j), ndfs(i, j), obsIndex(j), levelIndex(i, j), scales(i, j));
				MatrixXd grad = MatrixXd::Zero(di, rj);
				for (int it1 = 0; it1 < rCount; ++it1)
				{
					rInd = nzIndt(it1, 0);
					cInd = nzIndt(it1, 1);
					grad.row(cInd) +=  yMp.row(rInd) * dmt(rInd, cInd);
				}
				// lambdaMax = max(lambdaMax, grad.norm()/penalties(i, j));
			}
		}
		/*  begin initializing intercepts  */
		RowVectorXd counts = yM(j).colwise().sum();
		if (counts(0) == 0.0) {
		  counts(0) = pow(10, -10);
		}
		counts /= counts(0);
		counts = counts.array().log();
		counts(0) = 0.0;
		betaM(node, j) = counts;
		logitM(j) = counts.replicate(nobsVec(j), 1);
		/*  end initializing intercepts  */
	}
    // generate lambda seq by CDAlgo.
//	lambdaSeq(0) = lambdaMax;
//	double ratio = pow(fmlam, 1.0 / (nlam - 1));
//	for (int it1 = 1; it1 < nlam; ++it1)
//		lambdaSeq(it1) = lambdaSeq(it1 - 1) * ratio;
    // edit penalty. Fei
//	penalties *= lambdaMax / ratio;

	MatrixXi G = MatrixXi::Zero(node, node);
	MatrixXd hvals(node + 1, node);
	for (int it1 = 0; it1 < nlam; ++it1)
	{
	  // Rcpp::Rcout << "\n the " << it1 << "th lambda \n";
		// update the Hessian
		/* we use the fact that the coding scheme is baseline coding (i.e., 0-1 dummy variables) to update the Hessian matrix */
		for (int j = 0; j < node; ++j)
		{
			MatrixXd yMp = logitM(j);
			VectorXd rowOP = yMp.rowwise().maxCoeff();
			yMp.colwise() -= rowOP;
			yMp = yMp.array().exp();
			rowOP = yMp.rowwise().sum();
			yMp = yMp.cwiseQuotient(rowOP.replicate(1, yMp.cols()));
			yMp = yMp.array() * (1.0 - yMp.array());

			hvals(node, j) = -max(yMp.colwise().sum().maxCoeff(), convLb);
			int di, nh = nobsVec(j), rj = ndfs(j, j), colInd;
			double value, scaleFactor;
			const VectorXi& rIn = obsIndex(j);
			for (int i = 0; i < node; ++i)
			{
				if (i != j)
				{
					// approximate the Hessian_{ji}
					di = ndfs(i, j);
					const MatrixXd& dm = dM(i);
					const VectorXi& cIn = levelIndex(i, j);
					const VectorXd& scaling = scales(i, j);
					MatrixXd hessianDiag = MatrixXd::Zero(di, rj);
					for (int it2 = 0; it2 < di; ++it2)
					{
						colInd = cIn(it2);
						scaleFactor = scaling(it2);
						for (int h = 0; h < nh; ++h)
						{
							value = dm(rIn(h), colInd);
							if (value != 0)
								hessianDiag.row(it2) += pow(value * scaleFactor, 2) * yMp.row(h);
						}
					}
					hvals(i, j) = -max(hessianDiag.maxCoeff(), convLb);
				}
			}
		}

//		penalties *= ratio; // Fei
        if (it1==0) {
            penalties *= lambdaSeq[0];
        }
        else {
            penalties *= lambdaSeq[it1]/lambdaSeq[it1-1];
        }

        // to keep record of time
        clock_t start, finish;
        double duration;

        start = clock();
		CDOnePoint(node, G, eor_nr, eor, eps, nobsVec, ndfs, dM, yM, yNZIndex, logitM, betaM,
				  IsBetaZeros, hvals, penalties, qtol, obsIndex, levelIndex, scales, maxRows, maxCols);
        finish = clock();

        //comput the log_likelihood for each lambda
        VectorXd pLog1(node);
        pLog1= VectorXd::Zero(node);
        for (int j1 = 0; j1 < node; j1++){
            VectorXi yInd1;
            MatrixXd yMp1;
            VectorXd rowMax1, log_values1;
            int n11;
            n11 = nobsVec(j1);
            yMp1 = logitM(j1);
            yInd1 = yNZIndex(j1);
            rowMax1 = yMp1.rowwise().maxCoeff();
            yMp1.colwise() -= rowMax1;
            log_values1 = (yMp1.array().exp().rowwise().sum()).log();
            for (int iti = 0; iti < n11 ; ++iti) {
                pLog1(j1) +=  yMp1(iti, yInd1(iti)) - log_values1(iti);
            }
        }
        log_like(it1) = 0.0;
        for (int j1 = 0; j1 < node; j1++) {
            log_like(it1) += pLog1(j1);
        }


        duration = (double)(finish - start) / CLOCKS_PER_SEC;

        dur(it1) = duration;

        int P = 0;

        for (int i = 0; i < node; i++) {
            for (int j = 0; j < node; j++) {
                if (G(i, j) == 1) {
                    P++;
                }
                estimateG((it1*node+i), j) = G(i, j);
            }
        }

        // to make sure sparsity of DAG
        if (P >= threshold*node) {
            break;
        }
	}
}

void maxLambda(int node, int dataSize, const MatrixXi& data, const VectorXi& nlevels, const VectorXVXi& obsIndex, const MatrixXVXi& levelIndex, MatrixXMXd& betaM, MatrixXd& weights, double& lambda, double gamma, double upperbound)
{
    int maxRows = dataSize, maxCols = nlevels.maxCoeff(); //nlevels record number of levels for each node.
    if(upperbound >= 0.0) // to calculate weights
    {
        MatrixXd umtr = MatrixXd::Constant(node, node, pow(upperbound, gamma)); // node by node matrix with elements all equal to upperbound to the power of gamma
        weights = weights.array().pow(gamma).min(umtr.array()); //
    }
    else if(upperbound == -1.0) // no truncation
        weights = weights.array().pow(gamma);
    else
        Rprintf("upperbound sould be positive!");

    // create the full design matrix dM for each node (including the baseline column which is assumed to be the last factor level)
    // initialize all relevant information
    VectorXi nobsVec(node);
    MatrixXi ndfs(node, node); // p by p matrix. The jth column stores the number of levels for node j and the degrees of freedoms of all other nodes for the jth problem.
    MatrixXd penalties(node, node);
    VectorXMXd dM(node), yM(node);
    VectorXVXi yNZIndex(node); //p by 1 vector. Each element j is a vector whose hth element is the nonzero level for node j in sample h
    MatrixXVXd scales(node, node);
    for (int it1 = 0; it1 < node; ++it1) // it1 -> i
    {
        int size;
        nobsVec(it1) = obsIndex(it1).size(); // ith element is the number of observations for the ith problem
        for (int it2 = 0; it2 < node; ++it2) // it2 -> j
        {
            ndfs(it2, it1) = size = levelIndex(it2, it1).size();
            scales(it2, it1) = VectorXd::Zero(size);
        }
        dM(it1) = MatrixXd::Zero(dataSize, nlevels(it1));// design matrix for each node, the ith element is the matrix that indicate which level the ith node have
        yM(it1) = MatrixXd::Zero(nobsVec(it1), ndfs(it1, it1)); // ith element is a matrix record that rows are for the observations of ith problem, colums are for levels of ith node.
        yNZIndex(it1) = VectorXi::Zero(nobsVec(it1));
        penalties.col(it1) = (ndfs.col(it1) * ndfs(it1, it1)).cast<double>().cwiseSqrt();
    }
    penalties = penalties.array() * weights.array();	// Adaptive Lasso weights
    for (int it1 = 0; it1 < dataSize; ++it1) // i to determine design matrix
    {
        for (int it2 = 0; it2 < node; ++it2) // j
            dM(it2)(it1, data(it1, it2)) = 1.0;
    }
    for (int it1 = 0; it1 < node; ++it1) // to
    {
        int nrows = nobsVec(it1), ncols = ndfs(it1, it1); // number of observations fro the ith problem, number of levels for the ith node
        double value;
        const VectorXi& oind = obsIndex(it1); //ordered observational sample indices for the ith problem
        const VectorXi& lind = levelIndex(it1, it1); // ordered levels of node i, that should be included in the design matrix fot the ith problem
        const MatrixXd& dm = dM(it1);
        for (int it2 = 0; it2 < nrows; ++it2)
        {
            for (int it3 = 0; it3 < ncols; ++it3)
            {
                value = dm(oind(it2), lind(it3));
                yM(it1)(it2, it3) = value;
                if (value != 0)
                    yNZIndex(it1)(it2) = it3;
            }
        }
    }



    // compute the maximum lambda and sequence of lambdas
    // initialize betaM
    double lambdaMax = 0.0;
    VectorXMXd logitM(node);
    MatrixXb IsBetaZeros(node, node);
    for (int j = 0; j < node; ++j)
    {
        MatrixXd yMp = yM(j);
        RowVectorXd colMean = yMp.colwise().mean();
        yMp.rowwise() -= colMean;

        int rCount, di, rj = ndfs(j, j), rInd, cInd;
        MatrixXd dmt(maxRows, maxCols);
        MatrixXi nzIndt(maxRows * maxCols, 2);
        for (int i = 0; i < node; ++i)
        {
            if(i != j)
            {
                di = ndfs(i, j);
                betaM(i, j) = MatrixXd::Zero(di, rj);
                IsBetaZeros(i, j) = true;
                firstDMFetch(dmt, nzIndt, rCount, dM(i), nobsVec(j), ndfs(i, j), obsIndex(j), levelIndex(i, j), scales(i, j));
                MatrixXd grad = MatrixXd::Zero(di, rj);
                for (int it1 = 0; it1 < rCount; ++it1)
                {
                    rInd = nzIndt(it1, 0);
                    cInd = nzIndt(it1, 1);
                    grad.row(cInd) +=  yMp.row(rInd) * dmt(rInd, cInd);
                }
                lambdaMax = max(lambdaMax, grad.norm()/penalties(i, j));
            }
        }
        /*  begin initializing intercepts  */
        RowVectorXd counts = yM(j).colwise().sum();
        if (counts(0) == 0.0) {
          counts(0) = pow(10, -10);
        }
        counts /= counts(0);
        counts = counts.array().log();
        counts(0) = 0.0;
        betaM(node, j) = counts;
        logitM(j) = counts.replicate(nobsVec(j), 1);
        /*  end initializing intercepts  */
    }

    lambda = lambdaMax;
}

