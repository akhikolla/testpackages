#include <Rcpp.h>
#include <thread>
#include <future>
#include <iostream>
#include "boost/iterator/counting_iterator.hpp"
#include "tbb/task_group.h"

#define GCC_VERSION (__GNUC__ * 10000 \
             + __GNUC_MINOR__ * 100   \
             + __GNUC_PATCHLEVEL__)                            \
             
#if GCC_VERSION == 40603
  #undef USE_C11_THREADS
#else         
  #define USE_C11_THREADS
#endif
 
#ifdef USE_C11_THREADS
  #include "ThreadPool.h"
#endif 


//#include "complexvec.h"
//#include "RcppParallel.h"

namespace mytype { // define data type
  typedef std::complex<double> ComplexNumber;
  //typedef Complex2d ComplexNumber;
  typedef std::vector<ComplexNumber> ComplexVector;  
} // namespace mytype

// define some constants
const mytype::ComplexNumber one(1.0,0.0), two(2.0,0.0), zero(0.0,0.0), tiny(1e-16,0.0), huge(1e16,0.0);

/////////////////////////////////////////
// Inline opertators  for complex number
/////////////////////////////////////////

//inline double real(const Complex2d& x) {
//  return(x.extract(0));
//}

//inline double imag(const Complex2d& x) {
//  return(x.extract(1));
//}

inline std::complex<double> operator*(const std::complex<double>& x, const std::complex<double>& y) {
                return {
                        x.real()*y.real() - x.imag()*y.imag(),
                        x.real()*y.imag() + x.imag()*y.real()
                };
              }
              
inline std::complex<double> operator/(const std::complex<double>& x, const std::complex<double>& y) {
                const double denominator = y.real() * y.real() + y.imag() * y.imag();
                return {
                        (x.real()*y.real() + x.imag()*y.imag())/denominator,
                        (-x.real()*y.imag() + x.imag()*y.real())/denominator
                };
              }
              
inline std::complex<double> operator/(const double& x, const std::complex<double>& y) {
                const double denominator = y.real() * y.real() + y.imag() * y.imag();
                return {
                        x*y.real()/denominator,
                        -x*y.imag()/denominator
                };
              }

////////////////////////////////////////////////
// Inline function for rearrange data structure
////////////////////////////////////////////////

inline int Trimat(int i, int j) {
  return(i + (j+1)*j/2); 
} // i <= j
  
inline int get_phi(int i, int j, int k, int Bp1) {
  return(i*Bp1*Bp1 + j*Bp1 + k);
} // phi(a,b,m)

////////////////////////////////////////
// Struct Levin for series acceleration
////////////////////////////////////////

struct Levin {
    
    std::vector<double> numer,denom;
    int n, ncv;
    bool cnvgd;
    double small;
    double eps, lastval, lasteps;
    
    Levin(double epss) {
        n = 0;
        ncv = 0;
        cnvgd = 0;
        eps = epss;
        lastval = 0.0;
        small = 1e-8;
    }
    
    double next(double sum, double omega, double beta) {
        
        double fact,ratio,term,val;
        if ((sum==0)&&(omega==0)) return(0);
        term = 1.0/(beta+n);
        denom.push_back(term/omega);
        numer.push_back(sum*denom[n]);
        if (n > 0) {
            ratio = (beta+n)*term;
            for (int j=1; j<=n; j++) {
                fact = (n-j+beta)*term;
                numer[n-j] = numer[n-j+1]-fact*numer[n-j];
                denom[n-j] = denom[n-j+1]-fact*denom[n-j];
                term = term*ratio;
            }
        }
        n++;
        val = std::abs(denom[0]) < small ? lastval : numer[0]/denom[0];
        if (std::isnan(val)) val = 0;
        //Rcpp::Rcout << "denom = " << denom[0] << ", numer = " << numer[0] << std::endl;
        lasteps = std::abs(val-lastval);
        if (lasteps <= eps) ++ncv;
        //if((ncv>0) && (lasteps > eps)) ncv = 0;
        if (ncv >= 5) cnvgd = 1; 
        lastval = val;
        
        return val;
    }
};

/////////////////////
// Declare functions
/////////////////////

void lentz_plus_invBk1dBk_Cpp(const int Bp1, const std::vector<double>& xvec, 
  const mytype::ComplexVector& yvec, const mytype::ComplexVector& inv_Bk1dBk, 
  mytype::ComplexVector& lentz_plus_invBk1dBk);
  
void inv_Bk1dBk_Cpp(const int Bp1, const std::vector<double>& xvec, const mytype::ComplexVector& yvec, 
    mytype::ComplexVector& inv_Bk1dBk);
    
void BidBj_Cpp(const int Bp1, const std::vector<double>& xvec, const mytype::ComplexVector& yvec, 
    const mytype::ComplexVector& inv_Bk1dBk, mytype::ComplexVector& BidBj);

std::vector<double> prod_mu2_Cpp(const int a, const int A, const int Bp1, const std::vector<double>& mat);
std::vector<double> prod_lambda2_Cpp(const int a, const int A, const int Bp1, const std::vector<double>& mat);

void phi_Cpp (const mytype::ComplexNumber s, const int a0, const int b0, const std::vector<double>& lambda2, 
    const std::vector<double>& mu2, const int A, const int Bp1, const int maxdepth, 
    mytype::ComplexVector& phi, const std::deque<std::vector<double>>& prod_mu2, 
    const std::deque<std::vector<double>>& prod_lambda2, const std::deque<std::vector<double>>& xvec, 
    const std::deque<std::vector<double>>& yvec_minus_s, mytype::ComplexVector& yvec, 
    mytype::ComplexVector& lentz_plus_invBk1dBk, mytype::ComplexVector& inv_Bk1dBk, 
    mytype::ComplexVector& BidBj);
    
void bbd_lt_Cpp(const mytype::ComplexNumber s, const int a0, const int b0, const std::vector<double>& lambda1, 
    const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, 
    const int A, const int Bp1, const int maxdepth, mytype::ComplexVector& phi, 
    const std::deque<std::vector<double>>& prod_mu2, const std::deque<std::vector<double>>& prod_lambda2, 
    const std::deque<std::vector<double>>& xvec, const std::deque<std::vector<double>>& yvec_minus_s, 
    mytype::ComplexVector& yvec, mytype::ComplexVector& lentz_plus_invBk1dBk, 
    mytype::ComplexVector& inv_Bk1dBk, mytype::ComplexVector& BidBj, 
    mytype::ComplexVector& f);
    
std::vector<double> bbd_lt_invert_Cpp(double t, const int a0, const int b0, 
    const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, 
    const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, 
    const int A, const int Bp1, const int nblocks, const double tol, const int computeMode, 
    const int nThreads, const int maxdepth);
    
std::vector<double> SIR_Cpp(const double t, const double alpha, const double beta,
    const long int S0, const long int I0, const int Ap1, const int Bp1, const int direction, 
    const int nblocks, const double tol, const int computeMode, const int nThreads);
    
std::vector<double> bb_lt_invert_Cpp(double t, const std::vector<double>& lambda1, 
    const std::vector<double>& lambda2,
    const int Ap1, const int Bp1, const int direction, const int nblocks, 
    const double tol, const int computeMode, const int nThreads);
    
void bb_lt_Cpp(const mytype::ComplexNumber s, const std::vector<double>& lambda1, 
    const std::vector<double>& lambda2, const int Ap1, const int Bp1, const int direction,
    const std::vector<double>& yvec, mytype::ComplexVector& f);

/////////////////
// Generic loops
/////////////////

//#define DEBUG

namespace loops {
  
    struct STL {  
    
    	template <class InputIt, class UnaryFunction>
    	UnaryFunction for_each(InputIt first, InputIt last, UnaryFunction f) const {
#ifdef DEBUG      
      		Rcpp::Rcout << "Thread: " << *first << " to " << *last << std::endl;
#endif      	
			return std::for_each(first, last, f);        		
    	}   
    	   	
    	inline size_t id(size_t t) const {
    		return 0;
    	} 
    	
    	inline size_t private_size() const {
    		return 1;
    	}
    }; // STL
    
    struct AbstractC11Thread {
    
    	AbstractC11Thread(int t, int w) : nThreads(t), chunkSize(w / t) { }
    	
    	inline size_t id(size_t w)  {
    		return std::min(w / chunkSize,nThreads - 1);
    	}
    	
    	inline size_t private_size() {
    		return nThreads;
    	}
    	
      size_t nThreads;
      int chunkSize;
    }; // AbstractC11Thread

    struct TbbThreads : public AbstractC11Thread {
      
      TbbThreads(int t, int w) : AbstractC11Thread(t, w) { }
      
      template <class InputIt, class UnaryFunction>
      inline UnaryFunction for_each(InputIt begin, InputIt end, UnaryFunction function) {
        const int minSize = 0;
        
        if (nThreads > 1 && std::distance(begin, end) >= minSize) {	
          
          tbb::task_group workers;
          
          size_t chunkSize = this->chunkSize;
          
           size_t start = 0;
          for (int i = 0; i < nThreads - 1; ++i, start += chunkSize) {
            
            workers.run([begin, end, start, chunkSize, function]() {
              std::for_each(begin + start, begin + start + chunkSize, function);
            });
          }
          
          auto rtn = std::for_each(begin + start, end, function);
          
          workers.wait();
          
          return rtn;
        } else {				
          return std::for_each(begin, end, function);
        }            
        
        
        
        
        
        return function;
      }
      
    };
    
#ifdef USE_C11_THREADS    
        
    struct C11Threads : public AbstractC11Thread {
    	
    	//using AbstractC11Thread::AbstractC11Thread; // inherit constructor
    	C11Threads(int t, int w) : AbstractC11Thread(t, w) { }
    	
		template <class InputIt, class UnaryFunction>
		inline UnaryFunction for_each(InputIt begin, InputIt end, UnaryFunction function) const {
	
			const int minSize = 0;

			if (nThreads > 1 && std::distance(begin, end) >= minSize) {				  
				std::vector<std::thread> workers(nThreads - 1);
				// size_t chunkSize = std::distance(begin, end) / nThreads;
				size_t start = 0;
				for (int i = 0; i < nThreads - 1; ++i, start += chunkSize) {
					workers[i] = std::thread(
						std::for_each<InputIt, UnaryFunction>,
						begin + start, 
						begin + start + chunkSize, 
						function);
#ifdef DEBUG
					Rcpp::Rcout << "Thread #" << i << ": " << *(begin + start) << " to " << *(begin + start + chunkSize) << std::endl;
#endif					
				}
#ifdef DEBUG
					Rcpp::Rcout << "Thread #" << (nThreads - 1) << ": " << *(begin + start) << " to " << *(end) << std::endl;
#endif				
    
				auto rtn = std::for_each(begin + start, end, function);

				for (int i = 0; i < nThreads - 1; ++i) {
					workers[i].join();
				}
				return rtn;
			} else {				
				return std::for_each(begin, end, function);
			}                  
		}         
    }; // C11Threads
    
    
    struct C11Async : public AbstractC11Thread {
      
    	//using AbstractC11Thread::AbstractC11Thread; // inherit constructor
    	C11Async(int t, int w) : AbstractC11Thread(t, w) { }    	
      
		template <class InputIt, class UnaryFunction>
		inline UnaryFunction for_each(InputIt begin, InputIt end, UnaryFunction function) const {
	
			const int minSize = 0;

			if (nThreads > 1 && std::distance(begin, end) >= minSize) {				  
				std::vector<std::future<UnaryFunction>> workers(nThreads - 1);
				// size_t chunkSize = std::distance(begin, end) / nThreads;
				size_t start = 0;
				for (int i = 0; i < nThreads - 1; ++i, start += chunkSize) {
					workers[i] = std::async(
            std::launch::async,
						std::for_each<InputIt, UnaryFunction>,
						begin + start, 
						begin + start + chunkSize, 
						function);
#ifdef DEBUG
					Rcpp::Rcout << "Thread #" << i << ": " << *(begin + start) << " to " << *(begin + start + chunkSize) << std::endl;
#endif					
				}
#ifdef DEBUG
					Rcpp::Rcout << "Thread #" << (nThreads - 1) << ": " << *(begin + start) << " to " << *(end) << std::endl;
#endif				
        
				auto rtn = std::for_each(begin + start, end, function);

				for (int i = 0; i < nThreads - 1; ++i) {
					workers[i].get();
				}
				return rtn;
			} else {				
				return std::for_each(begin, end, function);
			}                  
		}         
    }; // C11Async
    
    struct C11ThreadPool : public AbstractC11Thread {
    
    C11ThreadPool(int t, int w) : AbstractC11Thread(t,w), pool(t) {}
    virtual ~C11ThreadPool() { };

	  //using AbstractC11Thread::AbstractC11Thread; // inherit constructor   
    ThreadPool pool;
    
    template <class InputIt, class UnaryFunction>
  	inline UnaryFunction for_each(InputIt begin, InputIt end, UnaryFunction function) {
	
      std::vector< std::future<void> > results;
      
			//size_t chunkSize = std::distance(begin, end) / nThreads;
			size_t start = 0;
			for (int i = 0; i < nThreads - 1; ++i, start += chunkSize) {
        results.emplace_back(
  					pool.enqueue([=] {
							std::for_each(
								begin + start, 
								begin + start + chunkSize,
								function);					
						})
					);
			}
			results.emplace_back(
  				pool.enqueue([=] {
						std::for_each(begin + start, end, function);
					})
				);
			
      for (auto&& result: results) result.get();	
      return function;
		  }
    }; // C11ThreadPool
    
#endif // USE_C11_THREADS     
    
 //   mutex;
//    ostream stream;
    
//   template <typename InputIt, typename UnaryFunction>
//    struct WrapWorker : public ::RcppParallel::Worker {
      
//      WrapWorker(InputIt begin, InputIt end, UnaryFunction function) : 
//        begin(begin), end(end), function(function) { }
        
//        void operator()(std::size_t i, std::size_t j) {
          // grab lock
  //        stream << i << ":" << j;
          // release lock
//            std::for_each(begin + i, begin + j, function);
          //  std::cout << i << ":" << j << std::endl;
//        }

//        InputIt begin;
//        InputIt end;
//        UnaryFunction function;
//    }; // WrapWorker
    
//    struct RcppThreads : public AbstractC11Thread {
//      using AbstractC11Thread::AbstractC11Thread; // inherit constructor
      
//      template <class InputIt, class UnaryFunction>  
//      inline UnaryFunction for_each(InputIt begin, InputIt end, UnaryFunction function) {
//            auto worker = WrapWorker<InputIt, UnaryFunction>(begin, end, function);
//            ::RcppParallel::parallelFor(0, std::distance(begin, end), worker);
//            return function;
//          }
      
//    }; // RcppThreads
    
} // namespace loops

////////////////
// Loops unroll
////////////////

//#define ROUND_DOWN(x, s) ((x) & ~((s)-1))

//namespace unroll {
  
//  template <class InputIt, class UnaryFunction>
//  UnaryFunction for_each_1(InputIt begin, InputIt end, UnaryFunction f) {
//    for (; begin != end; ++begin) {
//        f(*begin);
//    }
//    return f;
//  }

//  template <class InputIt, class UnaryFunction>
//  UnaryFunction for_each_2(InputIt begin, InputIt end, UnaryFunction f) {      
//    InputIt endr = end - std::distance(begin, end) % 2;
//    for (; begin != endr; ++begin) {
//        f(*begin);
//        f(*(++begin));
//    }
    
//    for (; begin != end; ++begin) {
//      f(*begin);
//    }
    
//    return f;
//  }
  
//  template <class InputIt, class UnaryFunction>
//  UnaryFunction for_each_4(InputIt begin, InputIt end, UnaryFunction f) {      
//    InputIt endr = end - std::distance(begin, end) % 4;
//    for (; begin != endr; ++begin) {
//        f(*begin);
//        f(*(++begin));
//        f(*(++begin));
//        f(*(++begin)); 
//    }
    
//    for (; begin != end; ++begin) {
//      f(*begin);
//    }
    
//    return f;
//  }

//template <typename Iterator, typename Function>
//void vectorized_for_each(Iterator begin, Iterator end, Function vector_function) {
  // Assumes that std::distance(begin, end) % 2 == 0
//  while (begin != end) {
//    vector_function(*begin); // SSE, AVX, OpenCL
//    begin += 2;
//  }
//}

//  template <class InputIt, class UnaryFunction>
//  UnaryFunction vectorized_for_each(InputIt begin, InputIt end, UnaryFunction vec_f) {      
    // Assumes that std::distance(begin, end) % 2 == 0
//    for (; begin != end; begin += 2) {
//        vec_f(*begin);
//    }    
//    return vec_f;
//  }

//} // namespace unroll


