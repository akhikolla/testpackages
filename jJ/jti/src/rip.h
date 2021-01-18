#ifndef RIP_H
#define RIP_H

#include "jti_types.h"
#include "set_ops.h"     // For set_any

Rcpp::List mcs(Rcpp::List & adj, std::string start_node, bool check);
VVS        perfect_cliques(VVS & x);
Rcpp::List perfect_separators(VVS & x);
Rcpp::List rip(Rcpp::List & adj, std::string start_node);

#endif
