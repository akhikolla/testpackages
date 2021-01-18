/*
 * MvBiClus.h
 *
 *  Created on: Aug 20, 2015
 *      Author: Javon
 */

#ifndef SRC_MVBICLUS_H_
#define SRC_MVBICLUS_H_

#include "MvClus.h"

class MvBiClus: public MvClus {
public:
	MvBiClus(const vector<mat> &);
	virtual ~MvBiClus();

	/** obtain feature clusters after clustering is performed */
	virtual const vector<uvec> & getFeatClus() = 0;
};

#endif /* SRC_MVBICLUS_H_ */
