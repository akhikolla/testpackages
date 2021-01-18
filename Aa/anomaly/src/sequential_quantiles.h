#ifndef ___SEQUENTIAL_QUANTILES_H___
#define ___SEQUENTIAL_QUANTILES_H___

#include <tuple>
#include <vector>

std::tuple<std::vector<double>,std::vector<double>> sequential_ests(const std::vector<double>&, int, int, std::tuple<double,double>, std::tuple<double,double>, std::tuple<double,double>);


#endif
