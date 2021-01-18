// #ifndef RIEMBASE_RIEMFACTORY_H
// #define RIEMBASE_RIEMFACTORY_H
#pragma once

#include <RcppArmadillo.h>

/* For Each Manifold, We have individual headers.*/
#include "manifolds/riem_01euclidean.h"
#include "manifolds/riem_02sphere.h"
#include "manifolds/riem_03spd.h"
#include "manifolds/riem_04stiefel.h"
#include "manifolds/riem_05grassmann.h"


/* Methods/Procedures under the names */
#include "procedures/riemfunc_dist.h"
#include "procedures/riemfunc_nearest.h"
#include "procedures/riemfunc_log.h"
#include "procedures/riemfunc_exp.h"
#include "procedures/riemfunc_norm.h"
#include "procedures/riemfunc_equiv.h"
#include "procedures/riemfunc_invequiv.h"
#include "procedures/riemfunc_extdist.h"

// #endif

