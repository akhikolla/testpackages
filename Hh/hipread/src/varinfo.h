#ifndef HIPREAD_VARINFO_H_
#define HIPREAD_VARINFO_H_

#include <Rcpp.h>

class VarInfo {
private:
  std::vector<std::vector<int> > starts;
  std::vector<std::vector<int> > widths;
  std::vector<std::vector<size_t> > var_pos;
  std::vector<size_t> num_vars_rectype;
  std::vector<int> max_ends;

public:
  VarInfo(Rcpp::List var_pos_info, size_t num_rt);
  int get_start(size_t rt_index, size_t col_num);
  int get_width(size_t rt_index, size_t col_num);
  size_t get_var_pos(size_t rt_index, size_t col_num);
  size_t get_num_vars(size_t rt_index);
  int get_max_end(size_t rt_index);
  std::vector<std::vector<int> > get_var_starts_rectype();
  std::vector<std::vector<int> > get_var_widths_rectype();
  std::vector <size_t> get_num_vars_rectype();
  std::vector<std::vector<size_t> > get_var_pos_rectype();
  std::vector <int> get_max_ends_rectype();
};

#endif
