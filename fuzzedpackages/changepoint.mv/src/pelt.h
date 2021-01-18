

#ifndef ___PELT_H___
#define ___PELT_H___


#include <list>
#include <vector>
#include <tuple>
#include <functional>
#include <future>

// used to wrap arguments for the pelt algorithm into a single argument. Makes using default values easy (via the constructor) without needing overloaded function
// definitions which would make using asynch much more difficult
class peltargs
{
  typedef  double(*f_ptr_type)(const std::vector<std::vector<double> >&, const int&, const int&,const int&);
  // arguments
 public:
  std::vector<double> data;
  f_ptr_type cost_func;
  double alpha;
  int m;
  std::shared_future<void> futureobj;

  
  // standard constructors
 public:
  peltargs(const std::vector<double>&,
	   double(*)(const std::vector<std::vector<double> >&,const int&, const int&,const int&),
	   const int&);
    
  peltargs(const std::vector<double>&,
	   double(*)(const std::vector<std::vector<double> >&, const int&, const int&,const int&),
	   const double&);
  
  peltargs(const std::vector<double>&,
	   double(*)(const std::vector<std::vector<double> >&, const int&, const int&, const int&));
  
  peltargs(const std::vector<double>&,
	   const double&);
  
  peltargs(const std::vector<double>&);
  
  peltargs(const std::vector<double>&,
	   double(*)(const std::vector<std::vector<double> >&,const int&, const int&,const int&),
	   const double&,
	   const int&);
  
  peltargs(const std::vector<double>&,
	   const double&,
	   const int&);
  
  peltargs(const std::vector<double>&,
	   const int&);

  // constructors for asynchronous use   
  peltargs(const std::vector<double>&,
	   double(*)(const std::vector<std::vector<double> >&, const int&, const int&,const int&),
	   const double&,
	   std::shared_future<void>);
  
  peltargs(const std::vector<double>&,
	   double(*)(const std::vector<std::vector<double> >&, const int&, const int&, const int&),
	   std::shared_future<void> _futureobj);
  
  peltargs(const std::vector<double>&,
	   const double&,
	   std::shared_future<void> _futureobj);
  
  peltargs(const std::vector<double>&,
	   std::shared_future<void> _futureobj);
  
  peltargs(const std::vector<double>&,
	   double(*)(const std::vector<std::vector<double> >&,const int&, const int&,const int&),
	   const double&,
	   const int&,
	   std::shared_future<void>);
    
  peltargs(const std::vector<double>&,
	   const double&,const int&,
	   std::shared_future<void>);
  
  peltargs(const std::vector<double>&,
	   const int&,
	   std::shared_future<void>);

  peltargs(const std::vector<double>&,
	   double(*)(const std::vector<std::vector<double> >&,const int&, const int&,const int&),
	   const int&,
	   std::shared_future<void>);

  
};



typedef std::function<double(const int&,const int&)> costfunc_type;


std::vector<std::tuple<double,std::list<int> > > pelt(const peltargs&);


std::vector<std::tuple<double,std::list<int> > > peltmv(const std::list<std::vector<double> >&,
							std::list<double(*)(const std::vector<std::vector<double> >&, const int&, const int&,const int&)>,
							const double&,
							const int&,
							std::shared_future<void> futureobj = std::shared_future<void>());

std::tuple<double,std::list<int> > amocmv(const std::list<std::vector<double> >&,
					  std::list<double(*)(const std::vector<std::vector<double> >&, const int&, const int&,const int&)>,
					  const double&,
					  const int&);


// cost functions
double norm_mean(const std::vector<std::vector<double> >&, const int&, const int&, const int&);
double norm_var(const std::vector<std::vector<double> >&, const int&, const int&, const int&);
double norm_meanvar(const std::vector<std::vector<double> >&, const int&, const int&, const int&);


// utilities
std::vector<std::vector<double> > sumstats(const std::vector<double>&);
costfunc_type cost_function(double(*)(const std::vector<std::vector<double> >&, const int&, const int&,const int&), const std::vector<double>&,const int&);


#endif
