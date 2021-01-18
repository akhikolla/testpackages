#include "bbd.h"

template <class ParallelizationScheme>
std::vector<double> bb_lt_invert_Cpp_impl(double t,
    const std::vector<double>& lambda1, const std::vector<double>& lambda2,
    const int Ap1, const int Bp1, const int direction,
    const int nblocks, const double tol,
    ParallelizationScheme& scheme) {
  
//  auto start = std::chrono::steady_clock::now();
  
  const double double_PI  = 3.141592653589793238463, AA = 20.0;
  const int matsize = Ap1*Bp1;
  int kmax = nblocks;

  std::vector<mytype::ComplexVector> ig;
  std::vector<double> res(matsize), yvec(matsize);
  
  for (int i=0; i<matsize; ++i) {
    yvec[i] = lambda1[i] + lambda2[i];
  }
  
/////////////////////////////////////////////////////////////
// The following code computes the inverse Laplace transform
// Algorithm from Abate and Whitt using a Riemann sum
// Levin tranform is used to accelerate the convergence
/////////////////////////////////////////////////////////////

  ig.resize(kmax);
  
  scheme.for_each( boost::make_counting_iterator(0), boost::make_counting_iterator(kmax),
    [&](int w) {
      mytype::ComplexNumber s(AA/(2*t),double_PI*(w+1)/t);
      ig[w].resize(Ap1*Bp1);
      bb_lt_Cpp(s,lambda1,lambda2,Ap1,Bp1,direction,yvec,ig[w]);
    });
    
  mytype::ComplexNumber s(AA/(2*t),0.0);
  mytype::ComplexVector psum0(matsize);
  bb_lt_Cpp(s,lambda1,lambda2,Ap1,Bp1,direction,yvec,psum0);
  
  std::for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(matsize),
    [&](int i) {
      Levin levin(tol); // A struct for Levin transform
      double term = 1e16, sdiff = 1e16;
      int k = 1;
      double psum = real(psum0[i])/(2*t);
      double sk,sk1;
      while ((std::abs(sdiff) > 1e-16)||(std::abs(term)>1e-3)) {
        double sgn = (k%2 == 0) ? 1.0 : -1.0;
        term = sgn*real(ig[k-1][i])/t;
        psum += term;
        double omega = k*term;
        sk = levin.next(psum,omega,1.0);
        if (k>1) sdiff = sk - sk1;
        k++;
        sk1 = sk;
        if (k > kmax) {
          ig.resize(kmax+nblocks);
                      
          scheme.for_each( boost::make_counting_iterator(0), boost::make_counting_iterator(nblocks),
            [&](int w) {
              mytype::ComplexNumber s(AA/(2*t),double_PI*(w+kmax+1)/t);
              ig[w+kmax].resize(matsize);
              bb_lt_Cpp(s,lambda1,lambda2,Ap1,Bp1,direction,yvec,ig[w+kmax]);
            });
            
          kmax += nblocks;
        }
      }
      
      res[i] = sk1*exp(AA/2);
    });
    
//  auto end = std::chrono::steady_clock::now();  
//  
//  using TimingUnits = std::chrono::microseconds;
//  Rcpp::Rcout << "Time: " << std::chrono::duration_cast<TimingUnits>(end - start).count() << std::endl;  
  
  return(std::move(res));
}

// [[Rcpp::export]]
std::vector<double> bb_lt_invert_Cpp(double t, const std::vector<double>& lambda1, 
    const std::vector<double>& lambda2,
    const int Ap1, const int Bp1, const int direction, const int nblocks, 
    const double tol, const int computeMode, const int nThreads) {
      
    switch(computeMode) {  // Run-time selection on compute_mode    
    
      case 4: {
        // std::cerr << "TbbThreads" << std::endl;
        loops::TbbThreads loopTbbThreads(nThreads, nblocks);
        return bb_lt_invert_Cpp_impl(t, lambda1, lambda2, Ap1, Bp1, direction,
                    nblocks, tol, loopTbbThreads);
      }
    
#ifdef USE_C11_THREADS    
    
      case 1: {
        // std::cerr << "C11Threads" << std::endl;
        loops::C11Threads loopC11Threads(nThreads, nblocks);
        return bb_lt_invert_Cpp_impl(t, lambda1, lambda2, Ap1, Bp1, direction,
                    nblocks, tol, loopC11Threads);
      }
           
      case 2: {
        // std::cerr << "C11ThreadPool" << std::endl;        
        loops::C11ThreadPool loopC11ThreadPool(nThreads, nblocks);
        return bb_lt_invert_Cpp_impl(t, lambda1, lambda2, Ap1, Bp1, direction,
                    nblocks, tol, loopC11ThreadPool);      
      }
      
      case 3: {
        // std::cerr << "C11Async" << std::endl;        
        loops::C11Async loopC11Async(nThreads, nblocks);
        return bb_lt_invert_Cpp_impl(t, lambda1, lambda2, Ap1, Bp1, direction,
                    nblocks, tol, loopC11Async);      
      }
#endif // USE_C11_THREADS
      
//      case 4: {        
//        loops::RcppThreads loopRcppThreads(nThreads, nblocks);
//        return bbd_lt_invert_Cpp_impl(t, a0, b0, lambda1, lambda2, mu2, gamma, x, y, A, Bp1,
//                    maxdepth, nblocks, tol, loopRcppThreads);      
//      }
        
      default: {
        // std::cerr << "STL" << std::endl;
        loops::STL loopSTL; 
        return bb_lt_invert_Cpp_impl(t, lambda1, lambda2, Ap1, Bp1, direction,
                    nblocks, tol, loopSTL);                
      }                   
    }
                
}


