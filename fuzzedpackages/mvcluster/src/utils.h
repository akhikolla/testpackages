/*
 * utils.h
 *
 *  Created on: Aug 27, 2015
 *      Author: Javon
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <RcppArmadillo.h>

/** count number of none zero components in a vector */
uint32_t countNoneZero(const vec& v);

double selectionByMedianOfMedians(const vector<double>& values, uint32_t k);


#endif /* SRC_UTILS_H_ */
