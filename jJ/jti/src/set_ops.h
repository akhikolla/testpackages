#ifndef SETOPS_H
#define SETOPS_H

#include "jti_types.h"

// For outlier_utils
VS   set_intersect(VS &v1, VS &v2);
VI   int_set_intersect(VI &v1, VI &v2);
bool set_in(std::string & a, VS &b);
bool set_issubeq(VS &a, VS &b);
bool set_any(std::vector<bool> &v);

#endif
