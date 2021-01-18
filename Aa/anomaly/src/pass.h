
#ifndef ___PASS_H___
#define ___PASS_H___


#include <list>
#include <vector>
#include <tuple>
#include <thread>
#include <future>

std::tuple<std::list<std::tuple<int,int> >, std::list<double> >
pass(const std::list<std::vector<double> >&,
     const int&,
     const int&,
     const int&,
     const double&,
     std::shared_future<void>);

#endif
