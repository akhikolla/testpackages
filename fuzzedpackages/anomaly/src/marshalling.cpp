

#include <Rcpp.h>

#include <vector>
#include <list>
#include <tuple>
#include <algorithm>
#include <thread>
#include <chrono>
#include <future>


#include "MeanVarAnomaly.h"
#include "MeanAnomaly.h"
#include "RobustMeanAnomaly.h"
#include "RobustMeanAnomalyMV.h"
#include "sequential_quantiles.h"
#include "recursive_anomalies.h"
#include "recursive_mvanomalies.h"
#include "MultivariateCOD.h"
#include "MeanAnomalyMV.h"
#include "pass.h"





//[[Rcpp::export]]
std::vector<int> marshall_MeanVarAnomaly(SEXP a,
				      SEXP b,
				      SEXP c,
				      SEXP d,
				      SEXP e,
				      SEXP f,
				      SEXP g)
{
  return MeanVarAnomaly(a,b,c,d,e,f,g);
}



//[[Rcpp::export]]
std::vector<int> marshall_MeanAnomaly(SEXP a,
				      SEXP b,
				      SEXP c,
				      SEXP d,
				      SEXP e,
				      SEXP f,
				      SEXP g)
{
  return MeanAnomaly(a,b,c,d,e,f,g);
}

//[[Rcpp::export]]
std::vector<int> marshall_RobustMeanAnomaly(SEXP a,
				            SEXP b,
				            SEXP c,
				            SEXP d,
				            SEXP e,
				            SEXP f,
				            SEXP g)
{
  return RobustMeanAnomaly(a,b,c,d,e,f,g);
}


//[[Rcpp::export]]
std::vector<int> marshall_recursive_anomalies(SEXP a,
					      SEXP b,
					      SEXP c)
{
  return recursive_anomalies(a,b,c);
}

//[[Rcpp::export]]
std::vector<int> marshall_recursive_mvanomalies(SEXP a,
						SEXP b,
						SEXP c,
						SEXP d,
						SEXP e,
						SEXP f,
						SEXP g)
						
{
  return recursive_mvanomalies(a,b,c,d,e,f,g);
}


//[[Rcpp::export]]
std::vector<int> marshall_MeanVarAnomalyMV(SEXP a,
					   SEXP b,
					   SEXP c,
					   SEXP d,
					   SEXP e,
					   SEXP f,
					   SEXP g,
					   SEXP h,
					   SEXP i
					   )
{
  return MeanVarAnomalyMV(a,b,c,d,e,f,g,h,i);
}

//[[Rcpp::export]]
std::vector<int> marshall_RobustMeanAnomalyMV(SEXP a,
					SEXP b,
					SEXP c,
					SEXP d,
					SEXP e,
					SEXP f,
					SEXP g,
					SEXP h,
					SEXP i)
{
  return RobustMeanAnomalyMV(a,b,c,d,e,f,g,h,i);
}

//[[Rcpp::export]]
std::vector<int> marshall_MeanAnomalyMV(SEXP a,
					SEXP b,
					SEXP c,
					SEXP d,
					SEXP e,
					SEXP f,
					SEXP g,
					SEXP h,
					SEXP i)
{
  return MeanAnomalyMV(a,b,c,d,e,f,g,h,i);
}


// [[Rcpp::export]]
Rcpp::List marshall_sequential_ests(const std::vector<double>& data, int n, int burnin, double lqs, double lqf0, double meds, double medf0, double uqs, double uqf0)
{

  std::tuple<double,double> lqinit = std::make_tuple(lqs, lqf0);
  std::tuple<double,double> medinit = std::make_tuple(meds, medf0);
  std::tuple<double,double> uqinit = std::make_tuple(uqs, uqf0);
  
  std::tuple<std::vector<double>,std::vector<double>> result = sequential_ests(data, n, burnin, lqinit, medinit, uqinit);

  return Rcpp::List::create(Rcpp::Named("mu")=std::get<0>(result),
			    Rcpp::Named("sigma")=std::get<1>(result));
}



// [[Rcpp::export]]
std::list<std::vector<double> > marshall_pass(const std::list<std::vector<double> >& Xi,
					      const int& Lmax,
					      const int& Lmin,
					      const int& alpha,
					      const double& lambda)
{
  // Create a std::promise object
  std::promise<void> exitSignal;
  // fetch std::future object associated with promise
  std::shared_future<void> futureObj = exitSignal.get_future();

  // start worker thread
  auto future = std::async(std::launch::async,pass,Xi,Lmax,Lmin,alpha,lambda,std::move(futureObj));

  // check for user interrupt
  try
    {
      while(std::future_status::ready != future.wait_for(std::chrono::milliseconds(0)))
	{
	  Rcpp::checkUserInterrupt();
	}
    }  
  catch(std::bad_alloc &e)
    {
      exitSignal.set_value();
      auto result = future.get(); // wait for it to tidy up
      Rcpp::stop("insufficient memory");
    }
  catch(...)
    {
      exitSignal.set_value();
      auto result = future.get(); // wait for it to tidy up
      Rcpp::stop("user interrupt");
    }
  
  auto result = future.get();
  auto cpts = std::get<0>(result);
  auto xstar = std::get<1>(result);
  std::list<std::vector<double> > marshalled_result(cpts.size());
  if(marshalled_result.size() > 0)
    {
      transform(cpts.begin(),cpts.end(),xstar.begin(),marshalled_result.begin(),[](auto& cpt,auto& x)
		{
		  std::vector<double> entry(3);
		  entry[0] = (double)std::get<0>(cpt);
		  entry[1] = (double)std::get<1>(cpt);
		  entry[2] = x;
		  return(entry);
		});
    }
  
  return marshalled_result;  
}



