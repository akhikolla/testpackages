#include <fstream>
#include "column.h"
#include "Progress.h"
#include "datasource.h"
#include "varinfo.h"
#include "rtinfo.h"
#include "iconv.h"
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;


static Function R6method(Environment env, const std::string& method) {
  return as<Function>(env[method]);
}
static bool isTrue(SEXP x) {
  if (!(TYPEOF(x) == LGLSXP && Rf_length(x) == 1)) {
    stop("`continue()` must return a length 1 logical vector");
  }
  return LOGICAL(x)[0] == TRUE;
}

// [[Rcpp::export]]
void read_chunked_long(
    CharacterVector filename,
    Environment callback,
    int chunksize,
    CharacterVector var_names,
    CharacterVector var_types,
    List rt_info_,
    List var_pos_info_,
    List var_opts_,
    int skip,
    bool isGzipped,
    CharacterVector encoding,
    bool progress
) {
  List rt_info = as<List>(rt_info_);
  List var_pos_info = as<List>(var_pos_info_);
  List var_opts = as<List>(var_opts_);
  Iconv pEncoder_(as<std::string>(encoding));

  DataSourcePtr data = newDataSource(as<std::string>(filename[0]), isGzipped);
  data->skipLines(skip);

  Progress ProgressBar = Progress();

  RtInfo rts(rt_info, var_pos_info.names());
  VarInfo vars(var_pos_info, rts.getNumRts());

  std::vector<size_t> vnum_per_rt = vars.get_num_vars_rectype();
  std::vector<std::vector<size_t> > vpos_per_rt = vars.get_var_pos_rectype();
  std::vector<std::vector<int> > start_per_rt = vars.get_var_starts_rectype();
  std::vector<std::vector<int> > width_per_rt = vars.get_var_widths_rectype();
  std::vector<int> max_ends_per_rt = vars.get_max_ends_rectype();

  int chunk_start = 1;
  while (isTrue(R6method(callback, "continue")()) && !data->isDone()) {
    std::vector<ColumnPtr> chunk = createAllColumns(var_types, var_opts, &pEncoder_);
    resizeAllColumns(chunk, chunksize);

    int i;
    const char* line_start;
    const char* line_end;
    for (i = 0; i < chunksize; ++i) {
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

        chunk[cur_var_pos]->setValue(i, x_start, x_end);
      }
    }

    resizeAllColumns(chunk, i);
    RObject chunk_df = columnsToDf(chunk, var_names, i);
    R6method(callback, "receive")(chunk_df, chunk_start);
    chunk_start += i;

    Rcpp::checkUserInterrupt();
    if (progress) {
      ProgressBar.show(data->progress_info());
    }
  }
  if (progress) {
    ProgressBar.show(data->progress_info());
  }
  ProgressBar.stop();
}

// [[Rcpp::export]]
void read_chunked_list(
    CharacterVector filename,
    Environment callback,
    int chunksize,
    List var_names_,
    List var_types_,
    List rt_info_,
    List var_pos_info_,
    List var_opts_,
    int skip,
    bool isGzipped,
    CharacterVector encoding,
    bool progress
) {
  List var_names = as<List>(var_names_);
  List var_types = as<List>(var_types_);
  List rt_info = as<List>(rt_info_);
  List var_pos_info = as<List>(var_pos_info_);
  List var_opts = as<List>(var_opts_);
  Iconv pEncoder_(as<std::string>(encoding));

  DataSourcePtr data = newDataSource(as<std::string>(filename[0]), isGzipped);
  data->skipLines(skip);

  Progress ProgressBar = Progress();

  RtInfo rts(rt_info, var_pos_info.names());
  VarInfo vars(var_pos_info, rts.getNumRts());

  std::vector<size_t> vnum_per_rt = vars.get_num_vars_rectype();
  std::vector<std::vector<int> > start_per_rt = vars.get_var_starts_rectype();
  std::vector<std::vector<int> > width_per_rt = vars.get_var_widths_rectype();
  std::vector<int> max_ends_per_rt = vars.get_max_ends_rectype();

  int chunk_start = 1;
  while (isTrue(R6method(callback, "continue")()) && !data->isDone()) {
    std::vector<std::vector<ColumnPtr> > chunk;
    std::vector<int> cur_pos_rt;
    for (size_t i = 0; i < rts.getNumRts(); ++i) {
      chunk.push_back(createAllColumns(
          var_types[static_cast<R_xlen_t>(i)], var_opts[static_cast<R_xlen_t>(i)], &pEncoder_
      ));
      resizeAllColumns(chunk[i], chunksize); // TODO: Could try to be smarter about allocation
      cur_pos_rt.push_back(-1);
    }

    int i;
    const char* line_start;
    const char* line_end;
    for (i = 0; i < chunksize; ++i) {
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

        chunk[rt_index][j]->setValue(cur_pos_rt[rt_index], x_start, x_end);
      }
    }

    List list_chunk;
    for (size_t j = 0; j < rts.getNumRts(); ++j) {
      resizeAllColumns(chunk[j], cur_pos_rt[j] + 1);
      list_chunk.push_back(columnsToDf(chunk[j], var_names[static_cast<R_xlen_t>(j)], cur_pos_rt[j] + 1));
    }
    list_chunk.names() = var_pos_info.names();
    R6method(callback, "receive")(list_chunk, chunk_start);
    chunk_start += i;

    Rcpp::checkUserInterrupt();
    if (progress) {
      ProgressBar.show(data->progress_info());
    }
  }
  if (progress) {
    ProgressBar.show(data->progress_info());
  }
  ProgressBar.stop();
}
