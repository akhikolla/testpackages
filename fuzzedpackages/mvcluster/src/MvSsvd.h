/*
 * MvSsvd.h
 *
 *  Created on: Aug 20, 2015
 *      Author: Javon
 */

#ifndef SRC_MVSSVD_H_
#define SRC_MVSSVD_H_

#include "MvBiClus.h"
#include "Objective.h"


class MvSsvd: public MvBiClus {
public:
	/** Constructor with only dataset as input */
	MvSsvd(const vector<mat> &);

	/** Constructor with inputs:
	 * 		datasets,
	 * 		lz - \lambda_z
	 * 		lvs - vector consists of \lambda_v's  */
	MvSsvd(const vector<mat> & datasets, float lz, const fvec & lvs);
	virtual ~MvSsvd();

	inline void setLz(float lz) {
		m_lz = lz;
	}

	inline float getLz() {
		return m_lz;
	}

	inline void setLvs(const fvec & lvs) {
		if (m_pLvs != NULL) {
			delete m_pLvs;
		}
		m_pLvs = new fvec(lvs);
	}

	inline const fvec & getLvs() {
		return *m_pLvs;
	}

	inline const vec & getZ() {
		return *m_pZ;
	}

	inline const mat & getU() {
		return *m_pU;
	}

	inline const vector<vec *> & getV() {
		return *m_pV;
	}

	inline void debugLevel(uint8_t debugLevel) {
		m_debugLevel = debugLevel;
	}

	inline void setDebugLevel(uint8_t debugLevel) {
		m_debugLevel = debugLevel;
	}

	inline void setMaxInnerIter(uint32_t maxInnerIter) {
		m_maxInnerIter = maxInnerIter;
	}

	inline void setMaxOuterIter(uint32_t maxOuterIter) {
		m_maxOuterIter = maxOuterIter;
	}

	inline void setThresholdInner(double thresholdInner) {
		m_thresholdInner = thresholdInner;
	}

	inline void setThresholdOuter(double thresholdOuter) {
		m_thresholdOuter = thresholdOuter;
	}

	/** run clustering*/
	virtual void clustering();

	/* retrieve cluster assignment after clustering has been performed */
	virtual const uvec & getCluster();

	/** retrieve feature cluster assignment after clustering is performed */
	virtual const vector<uvec> & getFeatClus();

	// constants
	static const uint32_t MAX_OUTER_ITER = 100000; // default maximum iterations of outer loop, i.e., alternating between Z and U,V
	static const double THRESHOLD_OUTER; // default terminating threshold of outer loop
	static const uint32_t MAX_INNER_ITER = 10000; // default maximum iterations of inner loop, i.e., alternating between U and V, with fixed Z
	static const double THRESHOLD_INNER; // default terminating threshold of inner loop

private:
	void init();
	void calcObj(Objective& obj);
	void ssvd(uint8_t iView, double lu, double lv);
	void solveU(uint8_t iView, const vec& Z, double lu);
	void solveV(uint8_t iView, const vec& Z, double lv);
	void solveUV(uint8_t iView);
	void solveZ();

	uint8_t m_debugLevel; // 0 - none-debugging mode, 1 - up to first level, 2 - up to second level

	/** preset parameter */
	vec * m_pSigma;
	uint32_t m_maxOuterIter; // maximum iterations of outer loop, i.e., alternating between Z and U,V
	double m_thresholdOuter; // terminating threshold of outer loop
	uint32_t m_maxInnerIter; // maximum iterations of inner loop, i.e., alternating between U and V, with fixed Z
	double m_thresholdInner; // terminating threshold of inner loop

	float m_lz;
	fvec * m_pLvs;
	vec * m_pZ;
	mat * m_pU; // n x v
	vector<vec *> * m_pV;
	Objective m_obj; // the value of objective function after optimization

	uvec * m_pCluster;
	vector<uvec> * m_pFeatCluster;
};

#endif /* SRC_MVSSVD_H_ */
