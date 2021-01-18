

#include <Rcpp.h>

#include <list>
#include <vector>
#include <tuple>

#include <iostream>

#include "pelt.h"

#include <thread>
#include <iostream>
#include <assert.h>
#include <chrono>
#include <future>

#include "user.interrupt.h"

auto cost_func_factory = [](const std::string& str)
{
  if(str == "mean") return(norm_mean);
  if(str == "var") return(norm_var);
  if(str == "meanvar") return(norm_meanvar);
  return(norm_mean); // default
};


// [[Rcpp::export]]
Rcpp::List rcppeigen_peltuv(const std::vector<double>& data, const std::string& str_cost_func, const double& alpha,const int& m)
{

  // Create a std::promise object
  std::promise<void> exitSignal;
  // fetch std::future object associated with promise
  // std::future<void> futureObj = exitSignal.get_future();
  std::shared_future<void> futureObj = exitSignal.get_future();
  peltargs args(data,cost_func_factory(str_cost_func),alpha,m,futureObj);
  auto future = std::async(std::launch::async,pelt,args);
  /*
  auto future = std::async(std::launch::async,
			   [](auto arg1,auto arg2,auto arg3,auto arg4,auto arg5) {return pelt(arg1,arg2,arg3,arg4,std::move(arg5));},
			   data,cost_func_factory(str_cost_func),alpha,m,std::move(futureObj));
  */
  while(std::future_status::ready != future.wait_for(std::chrono::milliseconds(0)))
    {
      if(check_user_interrupt())
	{
	  exitSignal.set_value();
	  auto result = future.get(); // wait for it to tidy up
	  throw(std::exception());
	}
    }
  
  auto result = future.get();
  // auto result = pelt(data,cost_func_factory(str_cost_func),alpha,m);
  int n = result.size();
  // get the lastchangelik
  auto it_result = result.begin();
  std::list<double> lastchangelike;
  while(it_result != result.end())
    {
      lastchangelike.push_back(std::get<0>(*it_result));
      it_result++;
    }

  return Rcpp::List::create(Rcpp::Named("cpts")=std::get<1>(result[n-1]),
			    Rcpp::Named("lastchangelike")=lastchangelike);
  
}


// [[Rcpp::export]]
std::vector<double> rcppeigen_tail_costs(const std::vector<double>& data, const std::string& str_cost_func,const int& m)
{
    auto cost_func = cost_func_factory(str_cost_func);
    auto stats = sumstats(data);
    int n = data.size();
    std::vector<double> res(n);
    for(int i = 0; i < n; i++)
      {
	res[i] = cost_func(stats,i+1,n,m);
      }
    return res;
}







