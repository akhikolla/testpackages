#ifndef match_visitor_h
#define match_visitor_h

#include <lvec_interface.h>
#include <algorithm>
#include <vector>

class match_visitor : public ldat::lvec_visitor {
  public:
    match_visitor(ldat::vec* order, ldat::vec* tab, ldat::vec* order_tab, bool na_incomparable = false) : 
        order_(order), tab_(tab), order_tab_(order_tab), result_(0), 
        na_incomparable_(na_incomparable) {
    }

    template<class T>
    class compare_lt {
      public:
        bool operator()(const T& lhs, const T& rhs) {
          if (ldat::is_nan(lhs)) return false;
          if (ldat::is_nan(rhs)) return true;
          return lhs < rhs;
        }
    };

    template<class T>
    class compare_eq {
      public:
        bool operator()(const T& lhs, const T& rhs, bool na_incomparable) {
          // we need the manually handle missing values, because for double na's
          // are coded as nan and therefore regular comparison doesn't work.
          if (na_incomparable) {
            if (ldat::is_nan(lhs) || ldat::is_nan(rhs)) return false;
          } else {
            if (ldat::is_nan(lhs) && ldat::is_nan(rhs)) return true;
            if (ldat::is_nan(lhs) || ldat::is_nan(rhs)) return false;
          }
          return lhs == rhs;
        }
    };

    template<typename T>
    void visit_template(ldat::lvec<T>& vec) {
      ldat::vec::vecsize size = vec.size();
      std::unique_ptr<ldat::lvec<double> > result(new ldat::lvec<double>(size));

      compare_lt<T> less_than;
      compare_eq<T> equal;
      
      if (size > 0 && tab_->size() > 0) {
      
        ldat::vec::vecsize j = 0;
        ldat::vec::vecsize index_tab = order_tab_->get_of_type(j, double())-1.0;
        T el_tab = tab_->get_of_type(index_tab, ldat::base_type(T()));
        for (ldat::vec::vecsize i = 0; i != size; ++i) {
          ldat::vec::vecsize index_vec = order_->get_of_type(i, double())-1.0;
          T el = vec.get(index_vec);
          while (less_than(el_tab, el) && j < (tab_->size()-1)) {
            ++j;
            index_tab = order_tab_->get_of_type(j, double())-1.0;
            el_tab = tab_->get_of_type(index_tab, ldat::base_type(T()));
  
          }
          if (equal(el, el_tab, na_incomparable_)) {
            result->set(index_vec, index_tab + 1.0);
          } else {
            result->set(index_vec, ldat::na<double>());
          }
        }
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
    ldat::vec* order_;
    ldat::vec* tab_;
    ldat::vec* order_tab_;
    ldat::lvec<double>* result_;
    bool na_incomparable_;
};

#endif
