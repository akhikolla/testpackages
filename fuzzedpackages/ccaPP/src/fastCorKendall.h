/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#ifndef _ccaPP_FASTCORKENDALL_H
#define _ccaPP_FASTCORKENDALL_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include "utils.h"

using namespace arma;

// functions to be used within C++
double fastCorKendall(const vec& x, const vec& y, const uword& n);

#endif
