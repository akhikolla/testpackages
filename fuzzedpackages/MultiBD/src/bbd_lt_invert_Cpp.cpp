#include "bbd.h"

template <class ParallelizationScheme>
std::vector<double> bbd_lt_invert_Cpp_impl(double t, const int a0, const int b0, 
    const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, 
    const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, 
    const int A, const int Bp1, const int maxdepth, 
    const int nblocks, const double tol,
    ParallelizationScheme& scheme) {
  
//  auto start = std::chrono::steady_clock::now();  
  
  const double double_PI  = 3.141592653589793238463, AA = 20.0;
  const int dimsq = Bp1*Bp1;
  const int matsize = (A-a0+1)*Bp1;
  int kmax = nblocks;

  std::vector<mytype::ComplexVector> ig;
  std::deque<std::vector<double>> prod_mu2, prod_lambda2, xvec, yvec_minus_s;
  std::vector<double> res(matsize);
  
  const size_t size = scheme.private_size();  
  std::vector<mytype::ComplexVector> phi(size), yvec(size), lentz_plus_invBk1dBk(size), inv_Bk1dBk(size),BidBj(size);

  for (int i = 0; i < size; ++i) {
    phi[i].resize((A-a0+1)*dimsq);
    yvec[i].resize(Bp1 + maxdepth);
    lentz_plus_invBk1dBk[i].resize(Bp1);
    inv_Bk1dBk[i].resize(Bp1);
    BidBj[i].resize(Bp1*(Bp1+1)/2);
  }
    
  for (int a=0; a<(A-a0+1); ++a) {
    prod_mu2.push_back(prod_mu2_Cpp(a+1,A-a0,Bp1,mu2));
    prod_lambda2.push_back(prod_lambda2_Cpp(a+1,A-a0,Bp1,lambda2));
    std::vector<double> tmpx(Bp1 + maxdepth), tmpy(Bp1 + maxdepth);
    std::copy_n(&x[a*(Bp1 + maxdepth)],Bp1 + maxdepth,tmpx.begin());
    std::copy_n(&y[a*(Bp1 + maxdepth)],Bp1 + maxdepth,tmpy.begin());
    xvec.push_back(tmpx);
    yvec_minus_s.push_back(tmpy);
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
      ig[w].resize((A+1-a0)*Bp1);
      bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,A,Bp1,maxdepth,phi[scheme.id(w)],prod_mu2,prod_lambda2,xvec,yvec_minus_s,
          yvec[scheme.id(w)],lentz_plus_invBk1dBk[scheme.id(w)],inv_Bk1dBk[scheme.id(w)],BidBj[scheme.id(w)],ig[w]);
    });
  
  mytype::ComplexNumber s(AA/(2*t),0.0);
  mytype::ComplexVector psum0(matsize);
  bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,A,Bp1,maxdepth,phi[0],prod_mu2,prod_lambda2,
      xvec,yvec_minus_s,yvec[0],lentz_plus_invBk1dBk[0],inv_Bk1dBk[0],BidBj[0],psum0);
  
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
//        Rcpp::Rcout << "psum = " << psum << ", omega = " << omega << std::endl;
        sk = levin.next(psum,omega,1.0);
        if (k>1) sdiff = sk - sk1;
        k++;
        sk1 = sk;
//        Rcpp::Rcout << "sk1 = " << sk1 << std::endl;
        if (k > kmax) {
          ig.resize(kmax+nblocks);
                      
          scheme.for_each( boost::make_counting_iterator(0), boost::make_counting_iterator(nblocks),
            [&](int w) {
              mytype::ComplexNumber s(AA/(2*t),double_PI*(w+kmax+1)/t);
              ig[w+kmax].resize(matsize);
              bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,A,Bp1,maxdepth,phi[scheme.id(w)],
                prod_mu2,prod_lambda2,xvec,yvec_minus_s,yvec[scheme.id(w)],lentz_plus_invBk1dBk[scheme.id(w)],
                inv_Bk1dBk[scheme.id(w)],BidBj[scheme.id(w)],ig[w+kmax]);
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
//  Rcpp::Rcout << "kmax: " << kmax << std::endl;  
  
  return(std::move(res));
}

// [[Rcpp::export]]
std::vector<double> bbd_lt_invert_Cpp(double t, const int a0, const int b0, 
    const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, 
    const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, 
    const int A, const int Bp1, const int nblocks, const double tol, const int computeMode, 
    const int nThreads, const int maxdepth) {
            
    switch(computeMode) {  // Run-time selection on compute_mode    
    
      case 4: {
        // std::cerr << "TbbThreads" << std::endl;
        loops::TbbThreads loopTbbThreads(nThreads, nblocks);
        return bbd_lt_invert_Cpp_impl(t, a0, b0, lambda1, lambda2, mu2, gamma, x, y, A, Bp1,
                                      maxdepth, nblocks, tol, loopTbbThreads);
      }
    
#ifdef USE_C11_THREADS    
    
      case 1: {
        // std::cerr << "C11Threads" << std::endl;
        loops::C11Threads loopC11Threads(nThreads, nblocks);
        return bbd_lt_invert_Cpp_impl(t, a0, b0, lambda1, lambda2, mu2, gamma, x, y, A, Bp1,
                    maxdepth, nblocks, tol, loopC11Threads);         
      }
           
      case 2: {
        // std::cerr << "C11ThreadPool" << std::endl;        
        loops::C11ThreadPool loopC11ThreadPool(nThreads, nblocks);
        return bbd_lt_invert_Cpp_impl(t, a0, b0, lambda1, lambda2, mu2, gamma, x, y, A, Bp1,
                    maxdepth, nblocks, tol, loopC11ThreadPool);      
      }
      
      case 3: {
        // std::cerr << "C11Async" << std::endl;        
        loops::C11Async loopC11Async(nThreads, nblocks);
        return bbd_lt_invert_Cpp_impl(t, a0, b0, lambda1, lambda2, mu2, gamma, x, y, A, Bp1,
                    maxdepth, nblocks, tol, loopC11Async);      
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
        return bbd_lt_invert_Cpp_impl(t, a0, b0, lambda1, lambda2, mu2, gamma, x, y, A, Bp1,
                    maxdepth, nblocks, tol, loopSTL);        
      }                   
    }    
}


