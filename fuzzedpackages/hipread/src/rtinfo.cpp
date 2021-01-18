#include <Rcpp.h>
using namespace Rcpp;

#include "rtinfo.h"

RtInfo::RtInfo(List rt_info, std::vector<std::string> rectypes_) : rectypes(rectypes_) {
  start = as<int>(rt_info["start"]);
  width = as<int>(rt_info["width"]);
  verbose_warning = as<bool>(rt_info["verbose_warning"]);
  hierarchical = width > 0;
}

bool RtInfo::getRtIndex(const char* line_start, const char* line_end, size_t& out) {
  if (!hierarchical) {
    out = 0;
    return true;
  }
  if (line_start + start + width > line_end) {
    Rcpp::stop("rectype variable cannot be longer than line.");
  }
  std::string rt(line_start + start, line_start + start + width);

  long rt_pos = std::distance(
    rectypes.begin(),
    std::find(rectypes.begin(), rectypes.end(), rt)
  );
  if (rt_pos < 0) Rcpp::stop("Could not parse rectype");
  if (static_cast<size_t>(rt_pos) == rectypes.size()) {
    if (verbose_warning) Rcpp::warning("Data has unknown record type '" + rt + "'");
    return false;
  }
  out = static_cast<size_t>(rt_pos);
  return true;
}

size_t RtInfo::getNumRts() {
  return rectypes.size();
}
