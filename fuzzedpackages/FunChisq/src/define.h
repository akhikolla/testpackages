//
//  define.h
//  eft-mhg
//
//  Created by Hua Zhong on 6/28/15.
//  Copyright (c) 2015 New Mexico State University. All rights reserved.
//

#if defined _WIN32
typedef double mydouble;
//typedef long double mydouble;
#else
typedef long double mydouble;
#endif

#include <vector>
#include <string>

mydouble funchisq(const std::vector<std::vector<int> > & O, mydouble & estimate,
                  const std::string index_kind);

mydouble funchisq(const std::vector<std::vector<int> > & O, const std::vector<int> & rowsums,
                  const std::vector<int> & colsums, int n);
