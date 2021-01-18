#include <fstream>
#include "column.h"
#include "Progress.h"
#include "datasource.h"
#include "varinfo.h"
#include "rtinfo.h"
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
RObject read_freqs(
    CharacterVector filename,
    CharacterVector var_names,
    List rt_info_,
    List var_pos_info_,
    bool isGzipped,
    bool progress
) {
  const int PROGRESS_TICK = 131072;
  List rt_info = as<List>(rt_info_);
  List var_pos_info = as<List>(var_pos_info_);

  DataSourcePtr data = newDataSource(as<std::string>(filename[0]), isGzipped);

  Progress ProgressBar = Progress();

  RtInfo rts(rt_info, var_pos_info.names());
  VarInfo vars(var_pos_info, rts.getNumRts());

  std::vector<size_t> vnum_per_rt = vars.get_num_vars_rectype();
  std::vector<std::vector<size_t> > vpos_per_rt = vars.get_var_pos_rectype();
  std::vector<std::vector<int> > start_per_rt = vars.get_var_starts_rectype();
  std::vector<std::vector<int> > width_per_rt = vars.get_var_widths_rectype();
  std::vector<int> max_ends_per_rt = vars.get_max_ends_rectype();

  std::vector<std::map<std::string, int> > out;
  for (int i = 0; i < var_names.size(); ++i) {
    out.push_back(std::map<std::string, int>());
  }

  int i = 0;
  const char* line_start;
  const char* line_end;
  while (!data->isDone()) {
    data->getLine(line_start, line_end);

    if (line_end - line_start == 0 ||
        (line_end - line_start == 1 && std::string(line_start, line_end) == "\r")) {
      if (data->isDone()) {
        break;
      } else {
        continue;
      }
    }

    size_t rt_index;
    bool rt_found = rts.getRtIndex(line_start, line_end, rt_index);
    if (!rt_found) {
      continue;
    }

    // Check if raw line is long enough
    if (line_end - line_start < max_ends_per_rt[rt_index]) {
      Rcpp::stop("Line is too short for rectype.");
    }

    // Loop through vars in rectype and add to out
    for (size_t j = 0; j < vnum_per_rt[rt_index]; j++) {
      const char *x_start = line_start + start_per_rt[rt_index][j];
      const char *x_end = x_start + width_per_rt[rt_index][j];
      std::string x(x_start, static_cast<size_t>(x_end - x_start));
      size_t cur_var_pos = vpos_per_rt[rt_index][j];

      if (out[cur_var_pos].count(x) == 0) {
        out[cur_var_pos][x] = 1;
      } else {
        out[cur_var_pos][x] += 1;
      }
    }

    if (i % PROGRESS_TICK == 0) {
      Rcpp::checkUserInterrupt();
      if (progress) {
        ProgressBar.show(data->progress_info());
      }
    }
    ++i;
  }
  if (progress) {
    ProgressBar.show(data->progress_info());
  }
  ProgressBar.stop();
  List wrapped = wrap(out);
  wrapped.names() = var_names;
  return wrapped;
}
