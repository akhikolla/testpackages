/*
 * MvSsvd.cpp
 *
 *  Created on: Aug 20, 2015
 *      Author: Javon
 */

#include "MvSsvd.h"
#include <math.h>
//#include <assert.h>

const double MvSsvd::THRESHOLD_OUTER = 0.00001;
const double MvSsvd::THRESHOLD_INNER = 0.00001;

MvSsvd::MvSsvd(const vector<mat> & datasets)
: MvBiClus(datasets) {
	m_lz = 0;
	m_pLvs = NULL;

	init();
}

MvSsvd::MvSsvd(const vector<mat> & datasets, float lz, const fvec & lvs)
: MvBiClus(datasets), m_lz(lz) {
	m_pLvs = new fvec(lvs);
	init();
}

MvSsvd::~MvSsvd() {
	if (m_pLvs != NULL) {
		delete m_pLvs;
		m_pLvs = NULL;
	}
	if (m_pZ != NULL) {
		delete m_pZ;
		m_pZ = NULL;
	}

	if (m_pV != NULL) {
		for (vector<vec *>::iterator it = m_pV->begin(); it != m_pV->end(); it++) {
			delete *it;
		}
		delete m_pV;
		m_pV = NULL;
	}

	if (m_pU != NULL) {
		delete m_pU;
		m_pU = NULL;
	}

	if (m_pSigma != NULL) {
		delete m_pSigma;
		m_pSigma = NULL;
	}
	if (m_pFeatCluster != NULL) {
		delete m_pFeatCluster;
		m_pFeatCluster = NULL;
	}
}

void MvSsvd::init() {
	// initialize V, column cluster assignment
	m_pV = new vector<vec *>(m_nView);
//	cout << m_pV->size() << endl;
	m_pFeatCluster = new vector<uvec>(m_nView);
	for (int i = 0; i < m_nView; i++) {
		// V
		vec * pVec = new vec(m_nFeat[i]);
		pVec->zeros();
		m_pV->at(i) = pVec;

		// column cluster assignment
		uvec fc(m_nFeat[i]);
		fc.zeros();
		m_pFeatCluster->at(i) = fc;
	}

//	cout << m_pV->size() << endl;

	// initialize Z
	m_pZ = new vec(m_nSample);
	m_pZ->zeros();

	// initialize U
	m_pU = new mat(m_nSample, m_nView);
	m_pU->zeros();

	// initialize row cluster assignment
	m_pCluster = new uvec(m_nSample);
	m_pCluster->zeros();

	// initialize sigma
	m_pSigma = new vec(m_nView);
	m_pSigma->zeros();

	// initialize terminating criteria
	m_maxOuterIter = MvSsvd::MAX_OUTER_ITER;
	m_maxInnerIter = MvSsvd::MAX_INNER_ITER;
	m_thresholdOuter = MvSsvd::THRESHOLD_OUTER;
	m_thresholdInner = MvSsvd::THRESHOLD_INNER;
}

/** run clustering*/
void MvSsvd::clustering() {
	// compute F-norm of all matrixes for initialization
	vec matFnormSquare(m_nView);
	for (int iView = 0; iView < m_nView; iView++) {
		matFnormSquare[iView] = pow(norm(m_datasets[iView], "fro"), 2);
	}
	ssvd(0, m_lz * matFnormSquare[0] / sum(matFnormSquare), m_pLvs->at(0)); // initialize with first view ssvd
	// initialize Z with U of first view
	for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
		if (m_pU->at(iSamp, 0) != 0) {
			m_pZ->at(iSamp) = 1;
		}
		else {
			m_pZ->at(iSamp) = 0;
		}
	}
//	cout << "initial Z: " << (*m_pZ).t() << endl;
	// initialize rest views with Z
	for (uint8_t iView = 1; iView < m_nView; iView++) {
		m_pU->col(iView) = m_pU->col(0); // initialize U with that of first view
		solveUV(iView);

//		cout << "initial of view " << iView << ":" << endl;
//		cout << "U: " << m_pU->col(iView).t() << endl;
//		cout << "V: " << (*m_pV->at(iView)).t() << endl;
//		cout << "sigma: " << m_pSigma->at(iView) << endl;
	}
	// solve Z
	solveZ();
	if (m_debugLevel >= 1) {
		calcObj(m_obj); //  calculate the objective
		Rprintf("mvbc: Initial objective: %2.2e, err: %2.2e, penalty: %2.2e, nonzero in z: %d\n",
				m_obj.getOveralObj(), m_obj.getLoss(), m_obj.getPenalty(), static_cast<int>(sum(*m_pZ)));
	}
//	cout << "second initial Z: " << (*m_pZ).t();
//	cout << "sigma: " << m_pSigma->at(0) << endl;
//	cout << "U1: " << m_pU->col(0).t();
//	cout << "sigma: " << m_pSigma->at(1) << endl;
//	cout << "U2: " << m_pU->col(1).t();

	// start alternating
	vec utmp(m_nSample);
	uint32_t iIter = 0;
	for (; iIter < m_maxOuterIter; iIter++) {
		utmp = m_pU->col(0); // preserve previous V of first view

		// given Z, solve U, V
		for (int iView = 0; iView < m_nView; iView++) {
			solveUV(iView);
		}

		// given U, V, solve Z
		solveZ();
//		cout << "Z: " << (*m_pZ).t() << endl;
		// calculate the different between current and last V, i.e., 2-norm of V - pre_V
		double diff = 0;
		for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
			diff += pow((utmp[iSamp] - m_pU->at(iSamp, 0)), 2);
		}
		diff = sqrt(diff);

		if (m_debugLevel >= 1) {
			calcObj(m_obj); // calculate the obj
			Rprintf("mvbc: Iter %d -- change in u1: %2.2e, nonzero in z: %d, objective: %2.2e, err: %2.2e, penalty: %2.2e\n",
					iIter, diff, static_cast<int>(sum(*m_pZ)), m_obj.getOveralObj(), m_obj.getLoss(), m_obj.getPenalty());
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
		if (m_pZ->at(iSamp) == 0) {
			m_pCluster->at(iSamp) = 0;
		}
		else {
			m_pCluster->at(iSamp) = 1;
		}
	}
	// feature clusters
	for (int iView = 0; iView < m_nView; iView++) {
		for (uword iFeat = 0; iFeat < m_nFeat[iView]; iFeat++) {
			if (m_pV->at(iView)->at(iFeat) != 0) {
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

/* retrieve cluster assignment after clustering has been performed */
const uvec & MvSsvd::getCluster() {
	return *m_pCluster;
}

const vector<uvec> & MvSsvd::getFeatClus() {
	return *m_pFeatCluster;
}

/** calculate objective value using all learned model parameters */
void MvSsvd::calcObj(Objective& obj) {
	// calculate the loss
	double loss = 0;
	for (int iView = 0; iView < m_nView; iView++) {
		for (uword iFeat = 0; iFeat < m_nFeat[iView]; iFeat++) {
			// armadillo stores matrix in column major, so we iterate column by column
			for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
				loss += pow((m_datasets[iView].at(iSamp, iFeat)
						- m_pSigma->at(iView) * m_pZ->at(iSamp) * m_pU->at(iSamp, iView) * m_pV->at(iView)->at(iFeat)), 2);
			}
		}
	}
	// calculate penalty
	double penalty = m_lz * norm(*m_pZ, 1);
	for (int iView = 0; iView < m_nView; iView++) {
		penalty += m_pLvs->at(iView) * m_pSigma->at(iView) * norm(*m_pV->at(iView), 1);
	}

	// set up output
	obj.setLoss(loss);
	obj.setPenalty(penalty);
	obj.setOveralObj(loss + penalty);
}

/** solve U given Z and V (from class member) for the given view */
void MvSsvd::solveU(uint8_t iView, const vec& Z, double lu) {
	m_pU->col(iView).zeros(); // initialize U by setting zero for all the components
	if (norm(*m_pV->at(iView), 1) == 0) { // test whether V contains all zeros, the testing may be optimized
		return;
	}

	// compute U
	double beta = lu / 2;
	for (uint32_t i = 0; i < m_nSample; i++) {
		if (Z[i] != 0) {
			// none 0 in Z, note here we ignore the factor of V's two norm, V is unit vector, so norm(V) = 1,
			// which won't affect the result
			m_pU->at(i, iView) = (m_datasets[iView].row(i) * (*m_pV->at(iView))).eval().at(0, 0);
			m_pU->at(i, iView) = m_pU->at(i, iView) / Z[i];

			// U is regularized, apply soft thresholding
			if (beta != 0) {
				double tmpBeta = beta / Z[i];
				if (m_pU->at(i, iView) > tmpBeta) {
					m_pU->at(i, iView) -= tmpBeta;
				}
				else if (m_pU->at(i, iView) < -1 * beta){
					m_pU->at(i, iView) += tmpBeta;
				}
				else {
					m_pU->at(i, iView) = 0;
				}
			}
		}
	}

	// calculate 2-norm of U to obtain sigma
	m_pSigma->at(iView) = norm(m_pU->col(iView), 2);
	if (m_pSigma->at(iView) != 0) {
		// normalize U
		for (uint32_t i = 0; i < m_nSample; i++) {
			m_pU->at(i, iView) /= m_pSigma->at(iView);
		}
	}
}

/** solve V given Z and U (from class member) for the given view */
void MvSsvd::solveV(uint8_t iView, const vec& Z, double lv) {
	// calculate element wise product of Z and U
	// we ignore the 2-norm of zu, as zu is unit vector
	shared_ptr<rowvec> pzu(new rowvec(m_nSample));
	for (uint32_t i = 0; i < m_nSample; i++) {
		pzu->at(i) = Z[i] * m_pU->at(i, iView);
	}
	// calculate beta
	double beta = lv / 2;
	// update V
	for (uword i = 0; i < m_nFeat.at(iView); i++) {
		double alpha = ((*pzu) * m_datasets[iView].col(i)).eval().at(0, 0);
		if (alpha > beta) {
			m_pV->at(iView)->at(i) = alpha - beta;
		}
		else if (alpha < -1 * beta){
			m_pV->at(iView)->at(i) = alpha + beta;
		}
		else {
			m_pV->at(iView)->at(i) = 0;
		}
	}
	// calculate 2-norm of V to obtain sigma
	m_pSigma->at(iView) = norm(*m_pV->at(iView), 2);
	if (m_pSigma->at(iView) != 0) {
		// normalize V
		for (uint32_t i = 0; i < m_nFeat.at(iView); i++) {
			m_pV->at(iView)->at(i) /= m_pSigma->at(iView);
		}
	}
}

/** solve U and V to the end given Z */
void MvSsvd::solveUV(uint8_t iView) {
	// solve V
	solveV(iView, *m_pZ, m_pLvs->at(iView));
	// start alternating
	vec vtmp(m_nFeat.at(iView));
	for (uint32_t iIter = 0; iIter < m_maxInnerIter; iIter++) {
		vtmp = *m_pV->at(iView); // preserve previous V
		solveU(iView, *m_pZ, 0); // solve U
		solveV(iView, *m_pZ, m_pLvs->at(iView)); // solve V

		// calculate the different between current and last V, i.e., 2-norm of V - pre_V
		double diff = 0;
		for (uword iFeat = 0; iFeat < m_nFeat.at(iView); iFeat++) {
			diff += pow((vtmp[iFeat] - m_pV->at(iView)->at(iFeat)), 2);
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

/** sparse svd of data matrix of given view,
 * the results are stored in class member U, V and sigma corresponding to given view */
void MvSsvd::ssvd(uint8_t iView, double lu, double lv) {
	// initialize with svd
	mat U;
	vec s;
	mat V;
	svd_econ(U, s, V, m_datasets[iView]);
//	svd(U, s, V, m_datasets[iView]);
	m_pU->col(iView) = U.col(0);
	for (uint32_t iFeat = 0; iFeat < m_nFeat.at(iView); iFeat++) {
		m_pV->at(iView)->at(iFeat) = V.at(iFeat, 0);
	}
	m_pSigma->at(iView) = s[0];

//	cout << "svd initial: " << endl;
//	cout << "U: " << m_pU->col(iView).t();
//	cout << "V: " << (*m_pV->at(iView)).t();
//	cout << "sigma: " << m_pSigma->at(iView) << endl;

	// start alternating
	vec Z(m_nSample);
	Z.ones(); // use Z with all ones
	vec utmp(m_nSample);
	for (uint32_t iIter = 0; iIter < m_maxInnerIter; iIter++) {
		utmp = m_pU->col(iView); // preserve previous V

		solveV(iView, Z, lv); // solve V
//		cout << "V: " << (*m_pV->at(iView)).t();
		solveU(iView, Z, lu); // solve U
//		cout << "U: " << m_pU->col(iView).t();

		// calculate the different between current and last V, i.e., 2-norm of V - pre_V
		double diff = 0;
		for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
			diff += pow((utmp[iSamp] - m_pU->at(iSamp, iView)), 2);
		}
		diff = sqrt(diff);

		if (m_debugLevel >= 2) {
			Rprintf("    ssvd: Iter %d -- change in u: %2.2e\n", iIter, diff);
		}

		if (diff < m_thresholdInner) {
			// the threshold is reached, terminate the loop
			break;
		}
	}

//	cout << "ssvd initial: " << endl;
//	cout << "U: " << m_pU->col(iView).t() << endl;
//	cout << "V: " << (*m_pV->at(iView)).t() << endl;
//	cout << "sigma: " << m_pSigma->at(iView) << endl;
}

/** solve Z, given U and V, which are retrieved from class member */
void MvSsvd::solveZ() {
	// calculate theta
	double theta = m_lz / 2;

	// ratio of gammas against the one of first view
	double* sigmaRatio = new double [m_nView];
	for (int i = 0; i < m_nView; i++) {
		sigmaRatio[i] = m_pSigma->at(i) / m_pSigma->at(0);
	}

	// update Z
	// compute gamma, and store the results in Z
	m_pZ->zeros(); // initialize with all zeros
	vec e2norm(m_nSample); // 2-norm of each row of matrix E, i.e., the concatenation of all UV
	e2norm.zeros(); // init with all zeros
	for (int iView = 0; iView < m_nView; iView++) {
		for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
			double coef = sigmaRatio[iView] * m_pU->at(iSamp, iView);
			for (uword iFeat = 0; iFeat < m_nFeat.at(iView); iFeat++) {
				double e = coef * m_pV->at(iView)->at(iFeat);
				// gamma[iSamp] = \sum_{iView, iFeat} (sigma[iView]/sigma[0] * U[iSamp, iView] * V[iView][iFeat] * M[iView][iSamp, iFeat])
				m_pZ->at(iSamp) += e * m_datasets[iView].at(iSamp, iFeat);
//				e2norm[iSamp] += pow(e, 2); // as V is a unit vector, so we can move this calculation out of this loop without affecting the result
			}
			e2norm[iSamp] += pow(coef, 2); // calculate e^2 here instead inside the most inner loop
		}
	}
	// we calculate \tilde{Z} via rules of soft thresholding
	for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
		if (e2norm[iSamp] == 0) {
			//assert(m_pZ->at(iSamp) == 0); // m_pZ->at(iSamp) is supposed to be 0
//			if (m_pZ->at(iSamp) != 0) {
//				cout << "error here in gamma calculation" << endl;
//				cout << "e2 = " << e2norm[iSamp] << ", " << "m_pZ->at(iSamp) = " << m_pZ->at(iSamp) << endl;
//			}
			continue;
		}
		m_pZ->at(iSamp) /= e2norm(iSamp);
		double tmpTheta = theta / e2norm[iSamp];
		if (m_pZ->at(iSamp) > tmpTheta) {
			m_pZ->at(iSamp) -= tmpTheta;
		}
		else if (m_pZ->at(iSamp) < -1 * tmpTheta) {
			m_pZ->at(iSamp) += tmpTheta;
		}
		else {
			m_pZ->at(iSamp) = 0;
		}
	}

	// update U and sigma, so when converting Z to indicator vector, the object remains unchanged
	for (int iView = 0; iView < m_nView; iView++) {
		for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
			m_pU->at(iSamp, iView) *= m_pZ->at(iSamp); // Z . U
		}
		double u2norm = norm(m_pU->col(iView), 2);
		// normalize U to unit vector
		for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
			m_pU->at(iSamp, iView) /= u2norm;
		}
		m_pSigma->at(iView) = sigmaRatio[iView] * u2norm; // update to sigma
	}
	delete [] sigmaRatio;
	// now it's ready to convert Z to indicator vector
	for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
		if (m_pZ->at(iSamp) != 0) {
			m_pZ->at(iSamp) = 1;
		}
	}
}
