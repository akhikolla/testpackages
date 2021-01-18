


#include <tuple>
#include <list>
#include <vector>
#include <algorithm>
#include <limits>
#include <numeric>
#include <cmath>
#include "pelt.h"
#include <thread>
#include <future>


peltargs::peltargs(const std::vector<double>& _data,
		   double(_cost_func)(const std::vector<std::vector<double> >&, const int&, const int&,const int&),
		   const double& _alpha)
{
  data = _data;
  cost_func = _cost_func;
  alpha = _alpha;
  m = 2;
  futureobj = std::shared_future<void>();
}

peltargs::peltargs(const std::vector<double>& _data,
		   double(_cost_func)(const std::vector<std::vector<double> >&, const int&, const int&,const int&))
{
  data = _data;
  cost_func = _cost_func;
  alpha = 2*log((double)data.size());
  m = 2;
  futureobj = std::shared_future<void>();
}

peltargs::peltargs(const std::vector<double>& _data,
		   const double& _alpha)
{
  data = _data;
  cost_func = norm_mean;
  alpha = _alpha;
  m = 2;
  futureobj = std::shared_future<void>();
}

peltargs::peltargs(const std::vector<double>& _data)
{
  data = _data;
  cost_func = norm_mean;
  alpha = 2*log((double)data.size());
  m = 2;
  futureobj = std::shared_future<void>();
}

peltargs::peltargs(const std::vector<double>& _data,
		   double(_cost_func)(const std::vector<std::vector<double> >&, const int&, const int&,const int&),
		   const double& _alpha,
		   const int& _m)
{
  data = _data;
  cost_func = _cost_func;
  alpha = _alpha;
  m = _m;
  futureobj = std::shared_future<void>();
}

peltargs::peltargs(const std::vector<double>& _data,
		   double(_cost_func)(const std::vector<std::vector<double> >&,const int&, const int&,const int&),
		   const int& _m)
{
  data = _data;
  cost_func = _cost_func;
  alpha = 2*log((double)data.size());
  m = _m;
  futureobj = std::shared_future<void>();
}

peltargs::peltargs(const std::vector<double>& _data,
		   const double& _alpha,
		   const int& _m)
{
  data = _data;
  cost_func = norm_mean;
  alpha = _alpha;
  m = _m;
  futureobj = std::shared_future<void>();
}

peltargs::peltargs(const std::vector<double>& _data,
		   const int& _m)
{
  data = _data;
  cost_func = norm_mean;
  alpha = 2*log((double)data.size());
  m = _m;
  futureobj = std::shared_future<void>();
}



// asynchronous constructors
peltargs::peltargs(const std::vector<double>& _data,
		   double(_cost_func)(const std::vector<std::vector<double> >&, const int&, const int&,const int&),
		   const double& _alpha,
		   std::shared_future<void> _futureobj)
{
  data = _data;
  cost_func = _cost_func;
  alpha = _alpha;
  m = 2;
  futureobj = _futureobj;
}

peltargs::peltargs(const std::vector<double>& _data,
		   double(_cost_func)(const std::vector<std::vector<double> >&, const int&, const int&,const int&),
		   std::shared_future<void> _futureobj)
{
  data = _data;
  cost_func = _cost_func;
  alpha = 2*log((double)data.size());
  m = 2;
  futureobj = _futureobj;
}

peltargs::peltargs(const std::vector<double>& _data,
		   const double& _alpha,
		   std::shared_future<void> _futureobj)
{
  data = _data;
  cost_func = norm_mean;
  alpha = _alpha;
  m = 2;
  futureobj = _futureobj;
}

peltargs::peltargs(const std::vector<double>& _data,
		   std::shared_future<void> _futureobj)
{
  data = _data;
  cost_func = norm_mean;
  alpha = 2*log((double)data.size());
  m = 2;
  futureobj = _futureobj;
}

peltargs::peltargs(const std::vector<double>& _data,
		   double(_cost_func)(const std::vector<std::vector<double> >&, const int&, const int&,const int&),
		   const double& _alpha,
		   const int& _m,
		   std::shared_future<void> _futureobj)
{
  data = _data;
  cost_func = _cost_func;
  alpha = _alpha;
  m = _m;
  futureobj = _futureobj;
}

peltargs::peltargs(const std::vector<double>& _data,
		   double(_cost_func)(const std::vector<std::vector<double> >&,const int&, const int&,const int&),
		   const int& _m,
		   std::shared_future<void> _futureobj)
{
  data = _data;
  cost_func = _cost_func;
  alpha = 2*log((double)data.size());
  m = _m;
  futureobj = _futureobj;
}

peltargs::peltargs(const std::vector<double>& _data,
		   const double& _alpha,
		   const int& _m,
		   std::shared_future<void> _futureobj)
{
  data = _data;
  cost_func = norm_mean;
  alpha = _alpha;
  m = _m;
  futureobj = _futureobj;
}

peltargs::peltargs(const std::vector<double>& _data,
		   const int& _m,
		   std::shared_future<void> _futureobj)
{
  data = _data;
  cost_func = norm_mean;
  alpha = 2*log((double)data.size());
  m = _m;
  futureobj = _futureobj;
}






std::vector<std::vector<double> > sumstats(const std::vector<double>& data)
{
  int n = data.size();
  std::vector<std::vector<double> > result(3);
  std::vector<double> temp(n+1);
  temp[0] = 0.0;
  std::partial_sum(data.begin(),data.end(),temp.begin()+1,std::plus<double>());
  result[0] = temp;
  double mu = temp[n]/n;
  std::transform(data.begin(),data.end(),data.begin(),temp.begin()+1,std::multiplies<double>());
  std::partial_sum(temp.begin(),temp.end(),temp.begin(),std::plus<double>());
  result[1] = temp;
  std::copy(data.begin(),data.end(),temp.begin()+1);
  std::transform(temp.begin()+1, temp.end(), temp.begin()+1, [&](double x){return(x-mu);});
  std::transform(temp.begin(),temp.end(),temp.begin(),temp.begin(),std::multiplies<double>());
  std::partial_sum(temp.begin(),temp.end(),temp.begin(),std::plus<double>());
  result[2] = temp;
  return(result);
}


double norm_mean(const std::vector<std::vector<double> >& X, const int& i, const int& j, const int& m)
{
  if((j-i) < m) return(std::numeric_limits<double>::infinity());
  double sx =  X[0][j] - X[0][i-1];
  double sx2 = X[1][j] - X[1][i-1];
  return(sx2-sx*sx/(j-i+1));
}

double norm_var(const std::vector<std::vector<double> >& X, const int& i, const int& j, const int& m)
{
  if((j-i) < m) return(std::numeric_limits<double>::infinity());
  unsigned int n = j-i+1;
  double sxmu2 = X[2][j] - X[2][i-1];
  if(sxmu2 <= 0.0) {sxmu2 = 0.00000000001;} 
  return n*(log(2.0*M_PI) +  log((sxmu2)/n)  + 1);
}

double norm_meanvar(const std::vector<std::vector<double> >& X, const int& i, const int& j,const int& m)
{
  if((j-i) < m) return(std::numeric_limits<double>::infinity());
  unsigned int n = j-i+1;
  double sigsq = ((X[1][j] - X[1][i-1]) - pow(X[0][j] - X[0][i-1],2)/n)/n;
  if(sigsq <= 0) {sigsq = 0.00000000001;}
  return n*(log(2.0*M_PI) +  log(sigsq) + 1);
}



costfunc_type cost_function(double(*cost_func)(const std::vector<std::vector<double> >&, const int&, const int&,const int&), const std::vector<double>& data,const int& m)
{
  auto X = sumstats(data);
  return([=](const int& i, const int& j) mutable -> double { return cost_func(X,i,j,m); });
}


auto mv_cost_function(std::list<double(*)(const std::vector<std::vector<double> >&, const int&, const int&,const int&)> lcostfuncs,
		      const std::list<std::vector<double> >& ldata,const int& m)
{
  std::list<std::vector<std::vector<double> > > lsumstats;
  std::list<std::vector<double> >::const_iterator it_ldata = ldata.begin();
  while(it_ldata != ldata.end())
    {
      lsumstats.push_back(sumstats(*it_ldata));
      it_ldata++;
    }
  
  return([=](const int& i, const int& j){
      double total = 0.0;
      std::list<std::vector<std::vector<double> > >::const_iterator it_lsumstats = lsumstats.begin();
      std::list<double(*)(const std::vector<std::vector<double> >&, const int&, const int&,const int&)>::const_iterator it_lcostfuncs = lcostfuncs.begin();
      while(it_lsumstats != lsumstats.end())
	{
	  total += (*it_lcostfuncs)(*it_lsumstats,i,j,m);
	  it_lcostfuncs++;
	  it_lsumstats++;
	}
      return(total);
    });  
}


std::vector<std::tuple<double,std::list<int> > > pelt_algo(costfunc_type C, const double& alpha, const int& n,const int& m,std::shared_future<void> futureobj)
{
  std::list<int> R{0};
  std::vector<std::tuple<double,std::list<int> > > F(n+1);
  F[0] =  std::make_tuple(-alpha,std::list<int>());
  for(int i = 1; i < n+1; i++)
    {
      
      if(i % 100 == 0 && futureobj.valid() && !(futureobj.wait_for(std::chrono::milliseconds(0)) == std::future_status::timeout))
	{
	  return(F); // interrupted - send back what we have got so far
	}
      
      std::vector<double> costs(R.size());
      std::transform(R.begin(),R.end(),costs.begin(),[&](auto k) {return(alpha+C(k+1,i)+std::get<0>(F[k])); });
      auto pos = std::min_element(costs.begin(),costs.end());
      auto it = R.begin();
      std::advance(it,std::distance(costs.begin(),pos));
      auto cpts = std::get<1>(F[*it]);
      cpts.push_back(*it);
      F[i] = make_tuple(*pos,cpts);
      R.remove_if([&](auto j) {return(i-j > m  &&  std::get<0>(F[j]) + C(j+1,i)  > std::get<0>(F[i]));});
      R.push_back(i);
    }
  std::get<1>(F[n]).push_back(n);
  return(F);
}


std::tuple<double,std::list<int> > amoc_algo(costfunc_type C, const double& alpha, const int& n)
{
  std::vector<double> costs(n,0);
  costs[0] = C(1,n);
  for(int k = 1; k < n; k++)
    {
      costs[k] = C(1,k) + C(k+1,n) + alpha; 
    }
  auto pos = std::min_element(costs.begin(),costs.end());
  auto cpts = std::list<int>();
  cpts.push_back(0);
  auto dist = std::distance(costs.begin(),pos);
  if(dist != 0)
    {
      cpts.push_back(dist);
    }
  cpts.push_back(n);
  return std::make_tuple(costs[dist],cpts);
}



// NOTE - retrun type used to be std::list<int> ie just the changepoints 


std::vector<std::tuple<double,std::list<int> > > pelt(const peltargs& args)
{
  return pelt_algo(cost_function(args.cost_func,args.data,args.m),args.alpha,args.data.size(),args.m,std::move(args.futureobj));
}


std::vector<std::tuple<double,std::list<int> > > peltmv(const std::list<std::vector<double> >& ldata,
							std::list<double(*)(const std::vector<std::vector<double> >&, const int&, const int&,const int&)> lcost_func,
							const double& alpha,
							const int& m,
							std::shared_future<void> futureobj)
{
  int n = ldata.begin()->size();
  return pelt_algo(mv_cost_function(lcost_func,ldata,m),alpha,n,m,std::move(futureobj)); 
}

std::tuple<double,std::list<int> > amocmv(const std::list<std::vector<double> >& ldata,
							std::list<double(*)(const std::vector<std::vector<double> >&, const int&, const int&,const int&)> lcost_func,
							const double& alpha,
							const int& m)
{
  int n = ldata.begin()->size();
  return amoc_algo(mv_cost_function(lcost_func,ldata,m),alpha,n); 
}





