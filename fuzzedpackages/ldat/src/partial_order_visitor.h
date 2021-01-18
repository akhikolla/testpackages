#ifndef partial_order_visitor_h
#define partial_order_visitor_h

#include <lvec_interface.h>
#include <algorithm>
#include <vector>

class porder_visitor : public ldat::lvec_visitor {
  public:
    porder_visitor(std::vector<ldat::vec::vecsize> pivots) : result_(0), pivots_(pivots) {
      // check pivots
      if (pivots_.size() == 0)
        throw Rcpp::exception("No pivots given");
      std::sort(pivots_.begin(), pivots_.end(), std::greater<ldat::vec::vecsize>());
    }

    template<class T>
    class compare {
      public:
        compare(const ldat::lvec<T>& vec) : vec_(vec) {};

        bool operator()(ldat::vec::vecsize lhs, ldat::vec::vecsize rhs) {
          T val_lhs = vec_.get(lhs-1);
          if (ldat::is_nan(val_lhs)) return false;
          T val_rhs = vec_.get(rhs-1);
          if (ldat::is_nan(val_rhs)) return true;
          return val_lhs < val_rhs;
        }

      private:
        const ldat::lvec<T>& vec_;
    };

    template<typename T>
    void visit_template(ldat::lvec<T>& vec) {
      // initialise result with current order
      ldat::vec::vecsize size = vec.size();
      std::unique_ptr<ldat::lvec<double> > result(new ldat::lvec<double>(size));
      for (ldat::vec::vecsize i = 0; i < size; ++i) result->set(i, i+1);
      // we are now going to sort the current order; but when comparing
      // we use the values pointed to by the order vector
      ldat::lvec_iterator<double> p = result->end();
      for (auto piv = pivots_.begin(); piv != pivots_.end(); ++piv) {
        if ((*piv) >= result->size()) 
          throw Rcpp::exception("Pivots out of range.");
        ldat::lvec_iterator<double> q = result->begin() + (*piv);
        std::nth_element(result->begin(), q, p, compare<T>(vec));
        p = result->begin() + (*piv);
      }
      // copy temporary result to result
      if (result_) delete result_;
      result_ = result.release();
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

    ldat::vec* result() {
      return result_;
    }

  private:
    ldat::lvec<double>* result_;
    std::vector<ldat::vec::vecsize> pivots_;
};


#endif
