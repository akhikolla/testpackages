
#ifndef ___BOUND_FUNC_H___
#define ___BOUND_FUNC_H___

#include "adfunc.h"
#include <boost/bind.hpp>

typedef Mat_1 ADmat;

extern boost::function<ADmat (const ADmat&)> bound_func;


#endif
