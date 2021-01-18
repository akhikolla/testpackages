#include <fstream>
#include "column.h"
#include "datasource.h"
#include "varinfo.h"
#include "rtinfo.h"
#include "iconv.h"
#include "hipread_types.h"
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
XPtrDataSource start_yield(
  CharacterVector filename,
  bool isGzipped,
  int skip
) {
  XPtrDataSource data = newXptrDataSource(as<std::string>(filename[0]), isGzipped);
  data->skipLines(skip);
  return data;
}

// [[Rcpp::export]]
void reset_yield(XPtrDataSource data, int skip) {
  data->reset();
  data->skipLines(skip);
}

// [[Rcpp::export]]
bool yield_is_done(XPtrDataSource data) {
  return data->isDone();
}

// [[Rcpp::export]]
RObject next_yield_long(
    XPtrDataSource data,
    CharacterVector var_names,
    CharacterVector var_types,
    List rt_info_,
    List var_pos_info_,
    List var_opts_,
    int yield_size,
    CharacterVector encoding
) {
  if (data->isDone()) return R_NilValue;

  List rt_info = as<List>(rt_info_);
  List var_pos_info = as<List>(var_pos_info_);
  List var_opts = as<List>(var_opts_);
  Iconv pEncoder_(as<std::string>(encoding));

  RtInfo rts(rt_info, var_pos_info.names());
  VarInfo vars(var_pos_info, rts.getNumRts());

  std::vector<size_t> vnum_per_rt = vars.get_num_vars_rectype();
  std::vector<std::vector<size_t> > vpos_per_rt = vars.get_var_pos_rectype();
  std::vector<std::vector<int> > start_per_rt = vars.get_var_starts_rectype();
  std::vector<std::vector<int> > width_per_rt = vars.get_var_widths_rectype();
  std::vector<int> max_ends_per_rt = vars.get_max_ends_rectype();

  std::vector<ColumnPtr> yield = createAllColumns(var_types, var_opts, &pEncoder_);
  resizeAllColumns(yield, yield_size);

  int i;
  const char* line_start;
  const char* line_end;
  for (i = 0; i < yield_size; ++i) {
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
      break;
    }

    // Check if raw line is long enough
    if (line_end - line_start < max_ends_per_rt[rt_index]) {
      Rcpp::stop("Line is too short for rectype.");
    }

    // Loop through vars in rectype and add to out
    for (size_t j = 0; j < vnum_per_rt[rt_index]; j++) {
      const char *x_start = line_start + start_per_rt[rt_index][j];
      const char *x_end = x_start + width_per_rt[rt_index][j];

      size_t cur_var_pos = vpos_per_rt[rt_index][j];

      yield[cur_var_pos]->setValue(i, x_start, x_end);
    }
  }

  resizeAllColumns(yield, i);
  return columnsToDf(yield, var_names, i);
}

// [[Rcpp::export]]
RObject next_yield_list(
    XPtrDataSource data,
    List var_names_,
    List var_types_,
    List rt_info_,
    List var_pos_info_,
    List var_opts_,
    int yield_size,
    CharacterVector encoding
) {
  if (data->isDone()) return R_NilValue;

  List var_names = as<List>(var_names_);
  List var_types = as<List>(var_types_);
  List rt_info = as<List>(rt_info_);
  List var_pos_info = as<List>(var_pos_info_);
  List var_opts = as<List>(var_opts_);
  Iconv pEncoder_(as<std::string>(encoding));

  RtInfo rts(rt_info, var_pos_info.names());
  VarInfo vars(var_pos_info, rts.getNumRts());

  std::vector<size_t> vnum_per_rt = vars.get_num_vars_rectype();
  std::vector<std::vector<int> > start_per_rt = vars.get_var_starts_rectype();
  std::vector<std::vector<int> > width_per_rt = vars.get_var_widths_rectype();
  std::vector<int> max_ends_per_rt = vars.get_max_ends_rectype();


  std::vector<std::vector<ColumnPtr> > yield;
  std::vector<int> cur_pos_rt;
  for (size_t i = 0; i < rts.getNumRts(); ++i) {
    yield.push_back(createAllColumns(
        var_types[static_cast<R_xlen_t>(i)], var_opts[static_cast<R_xlen_t>(i)], &pEncoder_
    ));
    resizeAllColumns(yield[i], yield_size); // TODO: Could try to be smarter about allocation
    cur_pos_rt.push_back(-1);
  }

  int i;
  const char* line_start;
  const char* line_end;
  for (i = 0; i < yield_size; ++i) {
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
      break;
    }
    cur_pos_rt[rt_index]++;
    // Check if raw line is long enough
    if (line_end - line_start < max_ends_per_rt[rt_index]) {
      Rcpp::stop("Line is too short for rectype.");
    }

    // Loop through vars in rectype and add to out
    for (size_t j = 0; j < vnum_per_rt[rt_index]; j++) {
      const char *x_start = line_start + start_per_rt[rt_index][j];
      const char *x_end = x_start + width_per_rt[rt_index][j];

      yield[rt_index][j]->setValue(cur_pos_rt[rt_index], x_start, x_end);
    }
  }

  List list_yield;
  for (size_t j = 0; j < rts.getNumRts(); ++j) {
    resizeAllColumns(yield[j], cur_pos_rt[j] + 1);
    list_yield.push_back(columnsToDf(yield[j], var_names[static_cast<R_xlen_t>(j)], cur_pos_rt[j] + 1));
  }
  list_yield.names() = var_pos_info.names();

  return list_yield;
}
