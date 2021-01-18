/*
 * MvClus.h
 *
 *  Created on: Aug 20, 2015
 *      Author: Javon
 */

#ifndef SRC_MVCLUS_H_
#define SRC_MVCLUS_H_

#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>
#include <memory>
#include <vector>
#include <stdint.h>
//#include <stdlib.h>

using namespace std;
using namespace arma;

class MvClus {
public:
	MvClus(const vector<mat> &);
	virtual ~MvClus();

	/** run clustering*/
	virtual void clustering() = 0;

	/* retrieve cluster assignment after clustering has been performed */
	virtual const uvec & getCluster() = 0;

protected:
	const vector<mat> & m_datasets;
	uint8_t m_nView; // number of views
	uint32_t m_nSample; // number of samples
	uvec m_nFeat; // number of features in each view
};

#endif /* SRC_MVCLUS_H_ */
