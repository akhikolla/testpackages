#ifndef partial_sort_visitor_h
#define partial_sort_visitor_h

#include <lvec_interface.h>
#include <algorithm>
#include <vector>

class psort_visitor : public ldat::lvec_visitor {
  public:
    psort_visitor(std::vector<ldat::vec::vecsize> pivots) : pivots_(pivots) {
      // check pivots
      if (pivots_.size() == 0)
        throw Rcpp::exception("No pivots given");
      std::sort(pivots_.begin(), pivots_.end(), std::greater<ldat::vec::vecsize>());
    }

    template<class T>
    class compare {
      public:
        bool operator()(const T& lhs, const T& rhs) {
          if (ldat::is_nan(lhs)) return false;
          if (ldat::is_nan(rhs)) return true;
          return lhs < rhs;
        }
    };

    template<typename T>
    void visit_template(ldat::lvec<T>& vec) {
      ldat::lvec_iterator<T> p = vec.end();
      for (auto piv = pivots_.begin(); piv != pivots_.end(); ++piv) {
        if (((*piv) >= vec.size()) || ((*piv) < 0)) 
          throw Rcpp::exception("Pivots out of range.");
        ldat::lvec_iterator<T> q = vec.begin() + (*piv);
        std::nth_element(vec.begin(), q, p, compare<T>());
        p = vec.begin() + (*piv);
      }
    }

    void visit(ldat::lvec<double>& vec) {
      return visit_template(vec);
    }

    void visit(ldat::lvec<int>& vec) {
      return visit_template(vec);
    }

    void visit(ldat::lvec<ldat::boolean>& vec) {
      return visit_template(vec);
    }

    void visit(ldat::lvec<std::string>& vec) {
      return visit_template(vec);
    }

  private:
    std::vector<ldat::vec::vecsize> pivots_;
};


#endif
