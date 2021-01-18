#include <Rcpp.h>
// Force compilation mode to C++11
// [[Rcpp::plugins(cpp11)]]
#include <chrono>
using namespace Rcpp;


//' @title Measure elapsed time with C++11 chrono library
//' @description Returns a time indicator that can be used to accurately measure elapsed time. The C++11 clock used is \code{chrono::high_resolution_clock}.
//' @return A double with the number of nanoseconds elapsed since a fixed epoch.
//' @examples
//' # Measure 1 second sleep
//' initT<-get_chronotime()
//' Sys.sleep(1)
//' measT<-(get_chronotime()-initT)*1e-9
//' cat("1 second passed in ",measT," seconds.\n")
//' @export
// [[Rcpp::export]]
double get_chronotime() {
  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = std::chrono::time_point<Clock>;
  using Ns = std::chrono::nanoseconds;
  const TimePoint time_point_now = Clock::now();
  const Clock::duration since_epoch = time_point_now.time_since_epoch();
//  double fs = std::chrono::duration_cast<Ns>(since_epoch).count();
  // std::cout << fs << " ms\n";
  return std::chrono::duration_cast<Ns>(since_epoch).count();
}
