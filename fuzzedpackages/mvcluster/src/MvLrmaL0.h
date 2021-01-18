/*
 * MvLrmaL0.h
 *
 *  Created on: Aug 27, 2015
 *      Author: Javon
 */

#ifndef SRC_MVLRMAL0_H_
#define SRC_MVLRMAL0_H_

#include "MvBiClus.h"
#include "Objective.h"
#include "utils.h"

class MvLrmaL0: public MvBiClus {
public:
	MvLrmaL0(const vector<mat>&);

	/** sz -- maximum number of none zeros in Z
	 *  sv -- maximum number of none zeros in V for each view */
	MvLrmaL0(const vector<mat>& datasets, uint32_t sz, const uvec& sv);

	virtual ~MvLrmaL0();

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

	inline float getGammau() const {
		return m_gammau;
	}

	inline void setGammau(float gammau) {
		m_gammau = gammau;
	}

	inline float getGammav() const {
		return m_gammav;
	}

	inline void setGammav(float gammav) {
		m_gammav = gammav;
	}

	float getGammaz() const {
		return m_gammaz;
	}

	void setGammaz(float gammaz) {
		m_gammaz = gammaz;
	}

	inline const mat& getMatU() const {
		return m_matU;
	}

	inline uint32_t getMaxIter() const {
		return m_maxIter;
	}

	inline void setMaxIter(uint32_t maxIter) {
		m_maxIter = maxIter;
	}

	inline const double getObj() const {
		return m_obj;
	}

	inline const vector<vec>& getPVecV() const {
		return *m_pVecV;
	}

	inline const uvec& getSvs() const {
		return m_svs;
	}

	inline void setSvs(const uvec& svs) {
		m_svs = svs;
	}

	inline uint32_t getSz() const {
		return m_sz;
	}

	inline void setSz(uint32_t sz) {
		m_sz = sz;
	}

	inline double getThreshold() const {
		return m_threshold;
	}

	inline void setThreshold(double threshold) {
		m_threshold = threshold;
	}

	inline  const vec& getVecZ() const {
		return m_vecZ;
	}

	uint32_t getISeedFeat() const {
		return m_iSeedFeat;
	}

	void setISeedFeat(uint32_t iSeedFeat) {
		m_iSeedFeat = iSeedFeat;
	}

	// constants
	static const uint32_t MAX_ITER = 1000; // default maximum iterations
	static const double THRESHOLD; // default terminating threshold of outer loop
	static const float GAMMA;

private:
	void init();
	double calcObj();
	void updateU(uint8_t iView);
	void updateV(uint8_t iView);
	void updateZ();

	/** map u to its closest that has <= s none zeros */
	void map(vec& v, uint32_t s);

	uint8_t m_debugLevel; // 0 - none-debugging mode, 1 - up to first level

	/** preset parameter */
	uint32_t m_maxIter; // maximum iterations, i.e., alternating between Z and U,V
	double m_threshold; // terminating threshold
	uint32_t m_sz;
	uvec m_svs;
	uint32_t m_iSeedFeat; // the index of seed feature of view 1, V of view 1 is set to zero except for this feature during initialization
	float m_gammaz;
	float m_gammau;
	float m_gammav;

	vec m_vecZ;
	mat m_matU; // n x v
	vector<vec> * m_pVecV;
	double m_obj; // the value of objective function after optimization

	uvec m_cluster;
	vector<uvec> * m_pFeatCluster;
};

#endif /* SRC_MVLRMAL0_H_ */
