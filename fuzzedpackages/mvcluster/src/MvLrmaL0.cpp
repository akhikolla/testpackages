/*
 * MvLrmaL0.cpp
 *
 *  Created on: Aug 27, 2015
 *      Author: Javon
 */

#include "MvLrmaL0.h"
//#include <assert.h>

const double MvLrmaL0::THRESHOLD = 0.00001;
const float MvLrmaL0::GAMMA = 1.2;

MvLrmaL0::MvLrmaL0(const vector<mat>& datasets)
: MvBiClus(datasets),
  m_debugLevel(0),
  m_maxIter(MvLrmaL0::MAX_ITER),
  m_threshold(MvLrmaL0::THRESHOLD),
  m_sz(0),
  m_gammaz(MvLrmaL0::GAMMA),
  m_gammau(MvLrmaL0::GAMMA),
  m_gammav(MvLrmaL0::GAMMA) {
	m_svs.set_size(datasets.size());
	m_svs.zeros();
	init();
	m_iSeedFeat = m_nFeat[0];
}

MvLrmaL0::MvLrmaL0(const vector<mat>& datasets, uint32_t sz, const uvec& sv)
: MvBiClus(datasets),
  m_debugLevel(0),
  m_maxIter(MvLrmaL0::MAX_ITER),
  m_threshold(MvLrmaL0::THRESHOLD),
  m_sz(sz),
  m_svs(sv),
  m_gammaz(MvLrmaL0::GAMMA),
  m_gammau(MvLrmaL0::GAMMA),
  m_gammav(MvLrmaL0::GAMMA) {
	init();
	m_iSeedFeat = m_nFeat[0]; // indicates seed feature is not specified
}

void MvLrmaL0::init() {
	// initialize V, column cluster assignment
	m_pVecV = new vector<vec>(m_nView);
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
}

MvLrmaL0::~MvLrmaL0() {
	if (m_pVecV != NULL) {
		delete m_pVecV;
		m_pVecV = NULL;
	}
	if (m_pFeatCluster != NULL) {
		delete m_pFeatCluster;
		m_pFeatCluster = NULL;
	}
}

void MvLrmaL0::clustering() {
	// initialization
	m_vecZ.ones(); // initialize Z with all ones
	m_matU.ones(); // initialize U with all ones
	// initialize V with all ones for all views
	for (uint8_t iView = 0; iView < m_nView; iView++) {
		m_pVecV->at(iView).ones();
	}

	if (m_iSeedFeat < m_nFeat[0]) {
		// seed feature is specified
		// initialize V with all zeros except the seed feature position
		m_pVecV->at(0).zeros();
		m_pVecV->at(0).at(m_iSeedFeat) = 1;
	}
	else {
		// seed feature is not specified
		// update V
		updateV(0);
//		cout << "V0: " << m_pVecV->at(0).t();
	}

	// update U of view 1
	updateU(0);
	// initialize Z according to none zero in  U of view 1
	vec tmp = m_matU.col(0);
//	cout << "U-tmp: " << tmp.t();
	map(tmp, m_sz);
	//cout << "none zero in U: " << countNoneZero(tmp) << endl;
	for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
		if (tmp[iSamp] == 0) {
			// if u is 0 then set corresponding z 0
			m_vecZ[iSamp]  = 0;
		}
	}
//	cout << "Z: " << m_vecZ.t();
	// initialize rest of the views
	for (uint8_t iView = 1; iView < m_nView; iView++) {
		updateV(iView); // update V
//		cout << "view" << iView << " - V: " << m_pVecV->at(iView).t();
		updateU(iView); // update U
//		cout << "view" << iView << " - U: " << m_matU.col(iView).t();
	}
	// update Z again
	updateZ();
//	cout << "Z: " << m_vecZ.t();

	if (m_debugLevel > 0) {
		m_obj = calcObj();
		Rprintf("initial -- objective: %2.2e\n", m_obj);
	}

	// start alternating
	uint32_t iter = 0;
	double preObj = 0;
	vec preZ(m_nSample);
	for (; iter < m_maxIter; iter++) {
		if (m_debugLevel > 0) {
			preObj = m_obj;
		}
		preZ = m_vecZ; // keep previous Z to help terminate the iteration
		for (uint8_t iView = 0; iView < m_nView; iView++) {
			updateV(iView); // update V
			updateU(iView); // update U
		}
		updateZ(); // update Z

		double diff = norm(preZ - m_vecZ, 2);
		if (m_debugLevel > 0) {
			m_obj = calcObj();
			Rprintf("Iter %d -- objective: %2.2e, change in obj: %2.2e, in z: %2.2e\n",
			            iter, m_obj, preObj - m_obj, diff);
		}

		if (diff < m_threshold) {
			// the difference is less than the threshold, end the iteration
			break;
		}
	}

	//assert(iter <= m_maxIter);

	if (iter == m_maxIter) {
		Rprintf("The algorithm does not converge, try to increase the maximum iteration limit");
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
		m_obj = calcObj();
	}
}

const uvec& MvLrmaL0::getCluster() {
	return m_cluster;
}

const vector<uvec>& MvLrmaL0::getFeatClus() {
	return *m_pFeatCluster;
}

/** calculate the objective */
double MvLrmaL0::calcObj() {
	double obj = 0;
	for (uint8_t iView = 0; iView < m_nView; iView++) {
		for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
			for (uint32_t iFeat = 0; iFeat < m_nFeat[iView]; iFeat++) {
				obj += pow(m_datasets[iView].at(iSamp, iFeat) - m_matU.at(iSamp, iView) * m_vecZ[iSamp] * m_pVecV->at(iView).at(iFeat), 2);
			}
		}
	}
	return obj;
}

/** update U */
void MvLrmaL0::updateU(uint8_t iView) {
	vec zu(m_nSample); // for the element wise production of z and u
	vec zz(m_nSample); // for element wise squre of z
	for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
		zu[iSamp] = m_vecZ[iSamp] * m_matU.at(iSamp, iView);
		zz[iSamp] = m_vecZ[iSamp] * m_vecZ[iSamp];
	}
	double v2sum = pow(norm(m_pVecV->at(iView), 2), 2);
	// calculate lipchitz constant
	double lip = v2sum * norm(zz, 2);
	if (lip == 0) {
		m_matU.col(iView).zeros(); // the constant is zero, so set U zero
		return;
	}

	// calculate gradient
	vec tmp1 = zu * v2sum;
	//cout << m_datasets[iView];
	//cout << m_pVecV->at(iView);
	vec tmp2 = m_datasets[iView] * m_pVecV->at(iView);
//	vec grad = zu * v2sum - m_datasets[iView] * m_pVecV->at(iView);
	vec grad = tmp1 -tmp2;
	for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
		grad[iSamp] *= m_vecZ[iSamp];
	}

	// update U
	m_matU.col(iView) -= grad / (m_gammau * lip);
}

void MvLrmaL0::updateV(uint8_t iView) {
	vec zu(m_nSample); // element wise product of z and u
//	cout << "m_matU.col(iView): " << m_matU.col(iView).t();
	for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
		zu.at(iSamp) = m_vecZ[iSamp] * m_matU.at(iSamp, iView);
	}
//	cout << "zu: " << zu.t();
	double lip = pow(norm(zu, 2), 2); // lipchiz constant
//	cout << "lip: " << lip << endl;
	if (lip == 0) {
		// the constant is zero, set V all zero
		m_pVecV->at(iView).zeros();
		return;
	}
	// calculate gradient
	vec grad = m_pVecV->at(iView) * lip - m_datasets[iView].t() * zu;
//	cout << "V(0)" << m_pVecV->at(iView) << ", lip: " << lip << endl;
//	cout << "M(:, 0): " << m_datasets[iView].col(0).t();
//	cout << "grad: " << grad.t();

	// update V, proximal step
	m_pVecV->at(iView) -= grad / (lip * m_gammav);
	// map to satisfy constraint
	map(m_pVecV->at(iView), m_svs[iView]);
}

void MvLrmaL0::updateZ() {
	vec grad(m_nSample); // gradient
	grad.zeros();
	double lip = 0; // lipschz constant
	vec liptv(m_nSample); // assist in calculating lipschz
	liptv.zeros();
	for (uint8_t iView = 0; iView < m_nView; iView++) {
//		vec xv = m_datasets[iView] * m_pVecV->at(iView); // X * V
		double v2sum = pow(norm(m_pVecV->at(iView), 2), 2);
//		cout << "view" << iView << " - v2sum: " << v2sum << endl;
		for (uint32_t iSamp = 0; iSamp < m_nSample; iSamp++) {
			double uu = m_matU.at(iSamp, iView) * m_matU.at(iSamp, iView); // square of components in U
//			grad[iSamp] += v2sum * (uu * m_vecZ[iSamp] - xv[iSamp] * m_matU[iSamp, iView]); // gradient
			grad[iSamp] += v2sum * uu * m_vecZ[iSamp] - (m_datasets[iView].row(iSamp) * m_pVecV->at(iView)).eval().at(0, 0) * m_matU.at(iSamp, iView); // gradient
			liptv[iSamp] += v2sum * uu;
		}
	}
	lip = norm(liptv, 2); // calculate lipschz

//	cout << "lip: " << lip << endl;
//	cout << "grad: " << grad.t();

	if (lip == 0) {
		// the constant is 0, set Z all zeros
		m_vecZ.zeros();
		return;
	}
	// update Z, proximal step
	m_vecZ -= grad / (lip * m_gammaz);
	// map to satisfy constraint
	map(m_vecZ, m_sz);
}

/** map u to its closest that has <= s none zeros */
void MvLrmaL0::map(vec& v, uint32_t s) {
	// calculate absolute of each component in v
	vector<double> av = conv_to<vector<double> >::from(abs(v));

//	cout << "av: ";
//	for (uint i = 0; i < av.size(); i++) {
//		cout<< av[i]<<" ";
//	}
//	cout << endl;

	// find the s-th largest of absolute of v
	double t = selectionByMedianOfMedians(av, s);
	// calculate the number of numbers that are larger than t
	uint32_t nLrg = 0;
	for (uint32_t i = 0; i < av.size(); i++) {
		if (av[i] > t) {
			nLrg++;
		}
	}
	uint32_t cnt = s - nLrg;
	// we need cnt number of numbers that are equal to t
	// set component that has an absolute value less than t to 0
	for (uint32_t i = 0; i < av.size(); i++) {
		if (av[i] < t) {
			v[i] = 0;
		}
		else if (av[i] == t) {
			if (cnt > 0) {
				cnt--;
			}
			else {
				v[i] = 0;
			}
		}
	}
}
