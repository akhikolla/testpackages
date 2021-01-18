/*
 * MvLrmaL1.cpp
 *
 *  Created on: Aug 25, 2015
 *      Author: Javon
 */

#include "MvLrmaL1.h"
#include <math.h>
//#include <assert.h>
#include "utils.h"

const double MvLrmaL1::THRESHOLD_OUTER = 0.00001;
const double MvLrmaL1::THRESHOLD_INNER = 0.00001;

MvLrmaL1::MvLrmaL1(const vector<mat> & datasets)
: MvBiClus(datasets) {
	m_lz = 0;
	m_lus.zeros();
	m_lvs.zeros();
	init();
}

MvLrmaL1::MvLrmaL1(const vector<mat>& datasets, float lz, const fvec& lus,
		const fvec& lvs)
: MvBiClus(datasets), m_lz(lz), m_lus(lus), m_lvs(lvs) {
	init();
}

void MvLrmaL1::init() {
	// initialize V, column cluster assignment
	m_pVecV = new vector<vec>(m_nView);
//	cout << m_pV->size() << endl;
	m_pFeatCluster = new vector<uvec>(m_nView);
	for (int i = 0; i < m_nView; i++) {
		// V
		vec v(m_nFeat[i]);
		v.zeros();
		m_pVecV->at(i) = v;

		// column cluster assignment
		uvec fc(m_nFeat[i]);
		fc.zeros();
		m_pFeatCluster->at(i) = fc;
	}

	// initialize Z
	m_vecZ.set_size(m_nSample);
	m_vecZ.zeros();

	// initialize U
	m_matU.set_size(m_nSample, m_nView);
	m_matU.zeros();

	// initialize row cluster assignment
	m_cluster.set_size(m_nSample);
	m_cluster.zeros();

	// initialize sigma
	m_sigmas.set_size(m_nView);
	m_sigmas.zeros();

	// initialize terminating criteria
	m_maxOuterIter = MvLrmaL1::MAX_OUTER_ITER;
	m_maxInnerIter = MvLrmaL1::MAX_INNER_ITER;
	m_thresholdOuter = MvLrmaL1::THRESHOLD_OUTER;
	m_thresholdInner = MvLrmaL1::THRESHOLD_INNER;
}

MvLrmaL1::~MvLrmaL1() {
	if (m_pVecV != NULL) {
		delete m_pVecV;
		m_pVecV = NULL;
	}
	if (m_pFeatCluster != NULL) {
		delete m_pFeatCluster;
		m_pFeatCluster = NULL;
	}
}

/** run clustering */
void MvLrmaL1::clustering() {
	// compute F-norm of all matrixes for initialization
	vec matFnormSquare(m_nView);
	for (int iView = 0; iView < m_nView; iView++) {
		matFnormSquare[iView] = pow(norm(m_datasets[iView], "fro"), 2);
	}
	sdecomp(0, m_lz * matFnormSquare[0] / sum(matFnormSquare), m_lvs[0]); // initialize with first view ssvd
	// initialize Z with U of first view
	for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
		if (m_matU.at(iSamp, 0) != 0) {
			m_vecZ[iSamp] = 1;
		}
		else {
			m_vecZ[iSamp] = 0;
		}
	}
	//cout << "initial Z: " << m_vecZ.t() << endl;
	// initialize rest views with Z
	for (uint8_t iView = 1; iView < m_nView; iView++) {
		m_matU.col(iView) = m_matU.col(0); // initialize U with that of first view
		solveUV(iView);

//		cout << "initial of view " << iView << ":" << endl;
//		cout << "U: " << m_matU.col(iView).t() << endl;
//		cout << "V: " << m_pVecV->at(iView).t() << endl;
	}
	// solve Z
	solveZ();
	if (m_debugLevel >= 1) {
		calcObj(m_obj); //  calculate the objective
		Rprintf("mvbc: Initial objective: %2.2e, err: %2.2e, penalty: %2.2e, nonzero in z: %d\n",
				m_obj.getOveralObj(), m_obj.getLoss(), m_obj.getPenalty(), countNoneZero(m_vecZ));
	}
//	cout << "second initial Z: " << m_vecZ.t();

	// start alternating
	vec utmp(m_nSample);
	uint32_t iIter = 0;
	for (; iIter < m_maxOuterIter; iIter++) {
		utmp = m_matU.col(0); // preserve previous V of first view

		// given Z, solve U, V
		for (int iView = 0; iView < m_nView; iView++) {
			solveU(iView, m_vecZ, m_lus[iView]); // solve U
			solveV(iView, m_vecZ, m_lvs[iView]); // solve V
//			solveUV(iView); // solve U, V to the end with fixed Z
		}

		// given U, V, solve Z
		solveZ();
//		cout << "Z: " << (*m_pZ).t() << endl;
		// calculate the different between current and last U, i.e., 2-norm of U - pre_U
		double diff = 0;
		for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
			diff += pow((utmp[iSamp] - m_matU.at(iSamp, 0)), 2);
		}
		diff = sqrt(diff);

		if (m_debugLevel >= 1) {
			calcObj(m_obj); // calculate the obj
			Rprintf("mvbc: Iter %d -- change in u1: %2.2e, nonzero in z: %d, objective: %2.2e, err: %2.2e, penalty: %2.2e\n",
					iIter, diff, countNoneZero(m_vecZ), m_obj.getOveralObj(), m_obj.getLoss(), m_obj.getPenalty());
		}

		if (diff < m_thresholdOuter) {
			// the threshold is reached, terminate the loop
			break;
		}
	}

	//assert(iIter <= m_maxOuterIter);
	if (iIter == m_maxOuterIter) {
		Rprintf("mvbc does not converge, try to increase the maximum iteration limit");
	}

	// produce clusters according to Z and V's
	// sample clusters
	for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
		if (m_vecZ[iSamp] == 0) {
			m_cluster[iSamp] = 0;
		}
		else {
			m_cluster[iSamp] = 1;
		}
	}
	// feature clusters
	for (int iView = 0; iView < m_nView; iView++) {
		for (uword iFeat = 0; iFeat < m_nFeat[iView]; iFeat++) {
			if (m_pVecV->at(iView).at(iFeat) != 0) {
				m_pFeatCluster->at(iView).at(iFeat) = 1;
			}
			else {
				m_pFeatCluster->at(iView).at(iFeat) = 0;
			}
		}
	}

	// calculate the objective, if havn't
	if (m_debugLevel < 1) {
		calcObj(m_obj);
	}
}

const uvec& MvLrmaL1::getCluster() {
	return m_cluster;
}

const vector<uvec>& MvLrmaL1::getFeatClus() {
	return *m_pFeatCluster;
}

/** calculate objective */
void MvLrmaL1::calcObj(Objective& obj) {
	// calculate the loss
	double loss = 0;
	for (int iView = 0; iView < m_nView; iView++) {
		for (uword iFeat = 0; iFeat < m_nFeat[iView]; iFeat++) {
			// armadillo stores matrix in column major, so we iterate column by column
			for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
				loss += pow((m_datasets[iView].at(iSamp, iFeat)
						- m_vecZ[iSamp] * m_matU.at(iSamp, iView) * m_pVecV->at(iView).at(iFeat)), 2);
			}
		}
	}
	// calculate penalty
	double penalty = m_lz * norm(m_vecZ, 1);
	for (int iView = 0; iView < m_nView; iView++) {
		penalty += m_lvs[iView] * norm(m_pVecV->at(iView), 1);
	}

	// set up output
	obj.setLoss(loss);
	obj.setPenalty(penalty);
	obj.setOveralObj(loss + penalty);
}

/** sparse low (one) rank approximation for single view */
void MvLrmaL1::sdecomp(uint8_t iView, double lu, double lv) {
	// initialize with svd
	mat U;
	vec s;
	mat V;
	svd_econ(U, s, V, m_datasets[iView]);
	double sigma2root = sqrt(s[0]);
	m_matU.col(iView) = U.col(0) * sigma2root;
	for (uint32_t iFeat = 0; iFeat < m_nFeat.at(iView); iFeat++) {
		m_pVecV->at(iView).at(iFeat) = V.at(iFeat, 0) * sigma2root;
	}
//	cout << "svd initial: " << endl;
//	cout << "U: " << m_matU.col(iView).t();
//	cout << "V: " << m_pVecV->at(iView).t();

	// start alternating
	vec Z(m_nSample);
	Z.ones(); // use Z with all ones
	vec utmp(m_nSample);
	for (uint32_t iIter = 0; iIter < m_maxInnerIter; iIter++) {
		utmp = m_matU.col(iView); // preserve previous V

		solveV(iView, Z, lv); // solve V
//		cout << "V: " << m_pVecV->at(iView).t();
		solveU(iView, Z, lu); // solve U
//		cout << "U: " << m_matU.col(iView).t();

		// calculate the different between current and last V, i.e., 2-norm of U - pre_U
		double diff = 0;
		for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
			diff += pow((utmp[iSamp] - m_matU.at(iSamp, iView)), 2);
		}
		diff = sqrt(diff);

		if (m_debugLevel >= 2) {
			Rprintf("    sdecomp: Iter %d -- change in u: %2.2e\n", iIter, diff);
		}

		if (diff < m_thresholdInner) {
			// the threshold is reached, terminate the loop
			break;
		}
	}

//	cout << "svd initial: " << endl;
//	cout << "U: " << m_matU.col(iView).t();
//	cout << "V: " << m_pVecV->at(iView).t();
}

/** solve U given Z and V (from class member) for the given view */
void MvLrmaL1::solveU(uint8_t iView, const vec& Z, double lu) {
	double v2sum = pow(norm(m_pVecV->at(iView), 2), 2); // sum of square of v

	if (v2sum == 0) {
		m_matU.col(iView).zeros(); // set all zeros
		return;
	}

	// compute U
	double beta = lu / 2;
	for (uint32_t i = 0; i < m_nSample; i++) {
		if (Z[i] == 0) {
			m_matU.at(i, iView) = 0; // if z at the position is 0, then set u 0
			continue;
		}

		// none 0 in Z
		// calculate alpha
		m_matU.at(i, iView) = (m_datasets[iView].row(i) * (m_pVecV->at(iView))).eval().at(0, 0);
		m_matU.at(i, iView) = m_matU.at(i, iView) / (Z[i] * v2sum);
		double tmpBeta = beta / (pow(Z[i], 2) * v2sum);
		// apply soft thresholding
		if (m_matU.at(i, iView) > tmpBeta) {
			m_matU.at(i, iView) -= tmpBeta;
		}
		else if (m_matU.at(i, iView) < -1 * tmpBeta){
			m_matU.at(i, iView) += tmpBeta;
		}
		else {
			m_matU.at(i, iView) = 0;
		}
	}
}

/** solve V given Z and U (from class member) for the given view */
void MvLrmaL1::solveV(uint8_t iView, const vec& Z, double lv) {
	// calculate element wise product of Z and U
	rowvec zu(m_nSample);
	double zu2sum = 0; // zu sum of square
	for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
		zu.at(iSamp) = Z[iSamp] * m_matU.at(iSamp, iView);
		zu2sum += pow(zu[iSamp], 2);
	}

	// calculate beta
	double beta = lv / (2 * zu2sum);
	// update V
	for (uword iFeat = 0; iFeat < m_nFeat.at(iView); iFeat++) {
		double alpha = (zu * m_datasets[iView].col(iFeat)).eval().at(0, 0) / zu2sum;
		// apply soft thresholding
		if (alpha > beta) {
			m_pVecV->at(iView).at(iFeat) = alpha - beta;
		}
		else if (alpha < -1 * beta){
			m_pVecV->at(iView).at(iFeat) = alpha + beta;
		}
		else {
			m_pVecV->at(iView).at(iFeat) = 0;
		}
	}
}

/** with given Z, optimize U and V to the end */
void MvLrmaL1::solveUV(uint8_t iView) {
	// solve V
	solveV(iView, m_vecZ, m_lvs.at(iView));
//	cout << "solveUV - V: " << m_pVecV->at(iView).t();
	// start alternating
	vec vtmp(m_nFeat.at(iView));
	for (uint32_t iIter = 0; iIter < m_maxInnerIter; iIter++) {
		vtmp = m_pVecV->at(iView); // preserve previous V
		solveU(iView, m_vecZ, m_lus[iView]); // solve U
//		cout << "solveUV - U: " << m_matU.col(iView).t();
		solveV(iView, m_vecZ, m_lvs[iView]); // solve V
//		cout << "solveUV - V: " << m_pVecV->at(iView).t();

		// calculate the different between current and last V, i.e., 2-norm of V - pre_V
		double diff = 0;
		for (uword iFeat = 0; iFeat < m_nFeat.at(iView); iFeat++) {
			diff += pow((vtmp[iFeat] - m_pVecV->at(iView).at(iFeat)), 2);
		}
		diff = sqrt(diff);

		if (m_debugLevel >= 2) {
			Rprintf("    solve_uv: Iter %d -- change in v: %2.2e\n", iIter, diff);
		}

		if (diff < m_thresholdInner) {
			// the threshold is reached, terminate the loop
			break;
		}
	}
}

/** solve Z, given U and V, which are retrieved from class member */
void MvLrmaL1::solveZ() {
	// calculate beta
	double beta = m_lz / 2;

	// update Z
	// compute alpha, and store the results in Z
	m_vecZ.zeros(); // initialize with all zeros
	vec e2sum(m_nSample); // sum of square of each row of matrix E, i.e., the concatenation of all UV
	e2sum.zeros(); // init with all zeros
	for (int iView = 0; iView < m_nView; iView++) {
		for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
			for (uword iFeat = 0; iFeat < m_nFeat.at(iView); iFeat++) {
				double e = m_matU.at(iSamp, iView) * m_pVecV->at(iView).at(iFeat);
				m_vecZ.at(iSamp) += e * m_datasets[iView].at(iSamp, iFeat);
				e2sum[iSamp] += pow(e, 2);
			}
		}
	}
	// we calculate Z via rules of soft thresholding
	for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
		if (e2sum[iSamp] == 0) {
			//assert(m_vecZ.at(iSamp) == 0); // m_vecZ.at(iSamp) is supposed to be 0
			continue;
		}
		m_vecZ.at(iSamp) /= e2sum(iSamp);
		double tmpBeta = beta / e2sum[iSamp];
		if (m_vecZ.at(iSamp) > tmpBeta) {
			m_vecZ.at(iSamp) -= tmpBeta;
		}
		else if (m_vecZ.at(iSamp) < -1 * tmpBeta) {
			m_vecZ.at(iSamp) += tmpBeta;
		}
		else {
			m_vecZ.at(iSamp) = 0;
		}
	}
}
