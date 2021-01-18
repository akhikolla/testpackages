#ifndef ACA_H
#define ACA_H
#include <functional>
#include <RcppEigen.h>
#include "tlr.h"
int aca(std::function<double(int,int)> f, int sz, double epsl, TLRNode &tlrNode);
#endif
