#ifndef HIPREAD_RTINFO_H_
#define HIPREAD_RTINFO_H_

#include <Rcpp.h>

class RtInfo {
private:
  int start;
  int width;
  std::vector<std::string> rectypes;
  bool hierarchical;
  bool verbose_warning;

public:
  RtInfo(Rcpp::List rt_info, std::vector<std::string> rectypes_);
  bool getRtIndex(const char* line_start, const char* line_end, size_t& out);
  size_t getNumRts();
};

#endif
