#include "print.h"

void print_vector(std::vector<int> other) {
  std::ostringstream oss;

  if (!other.empty())
  {
    std::copy(other.begin(), other.end() - 1, std::ostream_iterator<int>(oss, ","));
    oss << other.back();
  }
  
  Rcpp::Rcout << "(" << oss.str() << ")";
}

std::string sprint_vector(std::vector<int> other) {
  std::ostringstream oss;

  if (!other.empty())
  {
    std::copy(other.begin(), other.end() - 1, std::ostream_iterator<int>(oss, ","));
    oss << other.back();
  }
  
  return oss.str();
}

void print_alpha(const Rcpp::NumericVector alpha, int G) {
  int printed = 0;
  int skipped = 0;
  double last_alpha = -1;

  Rcpp::Rcout << "(";

  for (int i = 0; i < G; i++) {
    if (alpha[i] == last_alpha) {
      skipped++;
      continue;
    }

    if (skipped == 0 && i > 0) {
      Rcpp::Rcout << ", ";
    }
    
    last_alpha = alpha(i);

    printed++;
    
    if (skipped > 0) {
      Rcpp::Rcout << " x " << (skipped + 1);
      printed++;
      skipped = 0;
    }
    
    Rcpp::Rcout << alpha(i);
  }

  if (skipped > 0) {
    Rcpp::Rcout << " x " << (skipped + 1);
  }
    
  Rcpp::Rcout << ")" << std::endl;
}

void print_save_gs(const Rcpp::IntegerVector save_gs, int G) {
  int last_i = -1;
  int range_count = 0;
  
  for (int i = 0; i < G; i++) {
    if (save_gs(i) == 1) {
      last_i = i;
      range_count = 1;
      break;
    }
  }

  if (last_i == -1) {
    Rcpp::Rcout << "NONE" << std::endl;
    return;
  }  

  Rcpp::Rcout << "Generations: ";
  
  for (int i = last_i + 1; i < G; i++) {
    if (save_gs(i) == 0) {
      if (range_count > 0) {
        if (range_count == 1) {
          Rcpp::Rcout << (last_i + 1) << " ";
        } else {
          Rcpp::Rcout << (last_i + 1) << "-" << ((last_i + 1) + range_count - 1) << " ";
        }
      }
      
      range_count = 0;
      last_i = -1;
    }
    
    if (save_gs(i) == 1) {
      if (last_i == -1) {
        last_i = i;
        range_count = 1;
      } else {
        range_count++;
      }
    }
  }

  if (range_count > 0) {
    if (range_count == 1) {
      Rcpp::Rcout << (last_i + 1) << " ";
    } else {
      Rcpp::Rcout << (last_i + 1) << "-" << ((last_i + 1) + range_count - 1) << " ";
    }
  }
    
  Rcpp::Rcout << std::endl;
}

