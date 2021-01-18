/*
 * MvLrmaL1.h
 *
 *  Created on: Aug 25, 2015
 *      Author: Javon
 */

#ifndef SRC_MVLRMAL1_H_
#define SRC_MVLRMAL1_H_

#include "MvBiClus.h"
#include "Objective.h"

class MvLrmaL1: public MvBiClus {
public:
	MvLrmaL1(const vector<mat> & datasets);

	MvLrmaL1(const vector<mat> & datasets, float lz, const fvec& lus, const fvec& lvs);

	virtual ~MvLrmaL1();

	/** run clustering*/
	virtual void clustering();

	/* retrieve cluster assignment after clustering has been performed */
	virtual const uvec & getCluster();

	/** retrieve feature cluster assignment after clustering is performed */
	virtual const vector<uvec> & getFeatClus();

	inline uint8_t getDebugLevel() const {
		return m_debugLevel;
	}

	inline void setDebugLevel(uint8_t debugLevel) {
		m_debugLevel = debugLevel;
	}

	inline const fvec& getLus() const {
		return m_lus;
	}

	inline void setLus(const fvec& lus) {
		m_lus = lus;
	}

	inline const fvec& getLvs() const {
		return m_lvs;
	}

	inline void setLvs(const fvec& lvs) {
		m_lvs = lvs;
	}

	inline float getLz() const {
		return m_lz;
	}

	inline void setLz(float lz) {
		m_lz = lz;
	}

	inline const mat& getMatU() const {
		return m_matU;
	}

	inline const vec& getVecZ() const {
		return m_vecZ;
	}

	inline uint32_t getMaxInnerIter() const {
		return m_maxInnerIter;
	}

	inline void setMaxInnerIter(uint32_t maxInnerIter) {
		m_maxInnerIter = maxInnerIter;
	}

	inline uint32_t getMaxOuterIter() const {
		return m_maxOuterIter;
	}

	inline void setMaxOuterIter(uint32_t maxOuterIter) {
		m_maxOuterIter = maxOuterIter;
	}

	inline const Objective& getObj() const {
		return m_obj;
	}

	inline void setObj(const Objective& obj) {
		m_obj = obj;
	}

	inline const vector<vec>& getV() const {
		return *m_pVecV;
	}

	inline const vec& getSigmas() const {
		return m_sigmas;
	}

	inline void setSigmas(const vec& sigmas) {
		m_sigmas = sigmas;
	}

	double getThresholdInner() const {
		return m_thresholdInner;
	}

	inline void setThresholdInner(double thresholdInner) {
		m_thresholdInner = thresholdInner;
	}

	inline double getThresholdOuter() const {
		return m_thresholdOuter;
	}

	inline void setThresholdOuter(double thresholdOuter) {
		m_thresholdOuter = thresholdOuter;
	}

	// constants
	static const uint32_t MAX_OUTER_ITER = 100000; // default maximum iterations of outer loop, i.e., alternating between Z and U,V
	static const double THRESHOLD_OUTER; // default terminating threshold of outer loop
	static const uint32_t MAX_INNER_ITER = 10000; // default maximum iterations of inner loop, i.e., alternating between U and V, with fixed Z
	static const double THRESHOLD_INNER; // default terminating threshold of inner loop

private:
	void init();
	void calcObj(Objective& obj);
	void sdecomp(uint8_t iView, double lu, double lv);
	void solveU(uint8_t iView, const vec& Z, double lu);
	void solveV(uint8_t iView, const vec& Z, double lv);
	void solveUV(uint8_t iView);
	void solveZ();

	uint8_t m_debugLevel; // 0 - none-debugging mode, 1 - up to first level, 2 - up to second level

	/** preset parameter */
	vec m_sigmas;
	uint32_t m_maxOuterIter; // maximum iterations of outer loop, i.e., alternating between Z and U,V
	double m_thresholdOuter; // terminating threshold of outer loop
	uint32_t m_maxInnerIter; // maximum iterations of inner loop, i.e., alternating between U and V, with fixed Z
	double m_thresholdInner; // terminating threshold of inner loop

	float m_lz;
	fvec m_lus;
	fvec m_lvs;
	vec m_vecZ;
	mat m_matU; // n x v
	vector<vec> * m_pVecV;
	Objective m_obj; // the value of objective function after optimization

	uvec m_cluster;
	vector<uvec> * m_pFeatCluster;
};

#endif /* SRC_MVLRMAL1_H_ */
