#ifndef _haptools_PRINT_H
#define _haptools_PRINT_H

#include "common.h"

void print_vector(std::vector<int> other);
std::string sprint_vector(std::vector<int> other);
void print_alpha(const Rcpp::NumericVector alpha, int G);
void print_save_gs(const Rcpp::IntegerVector save_gs, int G);

#endif

