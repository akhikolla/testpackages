// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef __cplusplus
extern "C" {
#endif
 
#include <sys/time.h>

#ifdef __cplusplus
}
#endif

#include <cstddef>

// [[Rcpp::export]]
double fasttimer(  ){
    static int lastcall = 0 ;
    
    // alternates difference in time between these two variables (=> first call will return -ve value)
    static struct timeval t1, t2;
    double elapsedTime;

    if (lastcall == 0) {
      // start timer
      gettimeofday(&t1, NULL);
      lastcall = 1 ;
    } else {
      // stop timer
      gettimeofday(&t2, NULL);
      lastcall = 0 ;
    }
  
    if (lastcall == 0) {
      // compute and print the elapsed time in millisec
      elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
      elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
      // cout << elapsedTime << " ms.\n";
    } else {
      elapsedTime = (t1.tv_sec - t2.tv_sec) * 1000.0;      // sec to ms
      elapsedTime += (t1.tv_usec - t2.tv_usec) / 1000.0;   // us to ms
    }

    if (elapsedTime < 0.0)
        elapsedTime = 0.0 ;
    return elapsedTime ;
}


