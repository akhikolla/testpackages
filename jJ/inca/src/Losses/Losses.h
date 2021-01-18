#ifndef LOSSES_HPP
#define LOSSES_HPP

// Losses based on L1-norm
#define IS_L1BASED(X) ((X) < 10)
#include "L1.h"
#define LOSS_L1 1
#include "aL1.h"
#define LOSS_aL1 2
#include "rL1.h"
#define LOSS_rL1 3
#include "LB1.h"
#define LOSS_LB1 4
#include "rB1.h"
#define LOSS_rB1 5
//#include "Lasso1.h"
#define LOSS_RBLASSO1 6

// Losses based on L2-norm
#define IS_L2BASED(X) ((X) > 9 && (X) < 100)
#include "L2.h"
#define LOSS_L2 10
#include "aL2.h"
#define LOSS_aL2 20
#include "rL2.h"
#define LOSS_rL2 30
#include "LB2.h"
#define LOSS_LB2 40
#include "rB2.h"
#define LOSS_rB2 50
//#include "Lasso2.h"
#define LOSS_RBLASSO2 60

// Tools for losses based on boundaries
#define IS_LB(X) (((X) == LOSS_LB1) || ((X) == LOSS_LB2))
#define IS_RB(X) (((X) == LOSS_rB1) || ((X) == LOSS_rB2))
#define IS_RBLASSO(X) (((X) == LOSS_RBLASSO1) || ((X) == LOSS_RBLASSO2))
#define OFF_WHEN_BND(X, Y) ((X) || IS_LB((Y)) || IS_RB((Y)) || IS_RBLASSO((Y)))
#endif // LOSSES_HPP

#ifndef SIGN_TOOLS
#define SIGN_TOOLS
//Tools for sign of int, double, ...
#define i_sign(x) ((x) > 0 ? 1 : ((x) < 0 ? -1 : 0))
#define d_sign(x) ((x) > 0.0 ? 1.0 : ((x) < 0.0 ? -1.0 : 0.0))
#endif // SIGN_TOOLS
