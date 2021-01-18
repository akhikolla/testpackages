#include "varinfo.h"
#include <Rcpp.h>
using namespace Rcpp;

VarInfo::VarInfo(List var_pos_info, size_t num_rt) {
  for (int i = 0; i < static_cast<int>(num_rt); i++) {
    starts.push_back(as<List>(var_pos_info[i])["start"]);
    widths.push_back(as<List>(var_pos_info[i])["width"]);
    var_pos.push_back(as<List>(var_pos_info[i])["var_pos"]);
    num_vars_rectype.push_back(starts[static_cast<size_t>(i)].size());
    max_ends.push_back(as<IntegerVector>(as<List>(var_pos_info[i])["max_end"])[0]);
  }
}

int VarInfo::get_start(size_t rt_index, size_t col_num) {
  return starts[rt_index][col_num];
}

int VarInfo::get_width(size_t rt_index, size_t col_num) {
  return widths[rt_index][col_num];
}

size_t VarInfo::get_var_pos(size_t rt_index, size_t col_num) {
  return var_pos[rt_index][col_num];
}

size_t VarInfo::get_num_vars(size_t rt_index) {
  return num_vars_rectype[rt_index];
}

int VarInfo::get_max_end(size_t rt_index) {
  return max_ends[rt_index];
}

std::vector<std::vector<int> > VarInfo::get_var_starts_rectype() {
  return starts;
}

std::vector<std::vector<int> > VarInfo::get_var_widths_rectype() {
  return widths;
}

std::vector <size_t> VarInfo::get_num_vars_rectype() {
  return num_vars_rectype;
}

std::vector<std::vector<size_t> > VarInfo::get_var_pos_rectype() {
  return var_pos;
}

std::vector <int> VarInfo::get_max_ends_rectype() {
  return max_ends;
}
