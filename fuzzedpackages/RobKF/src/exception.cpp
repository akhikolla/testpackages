
#include <string>
#include "Rcpp.h"

using namespace Rcpp;

void throw_exception(const std::string& msg)
{
  stop(msg);
  return;
}
