#ifndef compare_visitor_h
#define compare_visitor_h

#include <lvec_interface.h>
#include <algorithm>
#include <vector>

class compare_visitor : public ldat::lvec_visitor {
  public:

    compare_visitor() : i_(0), j_(0) {
    }

    void set_ij(ldat::vec::vecsize i, ldat::vec::vecsize j) {
      i_ = i;
      j_ = j;
    }

    int compare(ldat::vec* vec, ldat::vec::vecsize i, ldat::vec::vecsize j) {
      set_ij(i, j);
      vec->visit(this);
      return result_;
    }

    template<typename T>
    void visit_template(const ldat::lvec<T>& vec) {
        // -1 i_ > j_
        //  0 i_ == j_
        //  1 i_ < j_
        T val_lhs = vec.get(i_);
        if (ldat::is_nan(val_lhs)) {
          result_ = -1;
          return;
        }
        T val_rhs = vec.get(j_);
        if (ldat::is_nan(val_rhs)) {
          result_ = 1;
          return;
        }
        result_ = -1;
        if (val_lhs == val_rhs) {
          result_ = 0;
        } else if (val_lhs < val_rhs) {
          result_ = 1;
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

    int result() const {
      return result_;
    }
    
  private:
    ldat::vec::vecsize i_, j_;
    int result_;
};

class ldat_compare {
  public:
    ldat_compare(const std::vector<ldat::vec*> vecs) : vecs_(vecs) {};

    bool operator()(ldat::vec::vecsize lhs, ldat::vec::vecsize rhs) {
      compare_visitor visitor;
      for (std::vector<ldat::vec*>::const_iterator p = vecs_.begin(); p != vecs_.end(); ++p) {
        int c = visitor.compare(*p, lhs-1, rhs-1);
        if (c != 0) return c > 0;
      }
      return false;
    }

  private:
    const std::vector<ldat::vec*> vecs_;
};


#endif
