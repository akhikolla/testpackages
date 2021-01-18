/// \file f_leading_ones_ruggedness3.hpp
/// \brief cpp file for class f_leading_ones_ruggedness3.
///
/// A detailed file description.
///
/// \author Furong Ye
/// \date 2019-06-27
#ifndef _F_LEADING_ONES_RUGGEDNESSTHREE_H
#define _F_LEADING_ONES_RUGGEDNESSTHREE_H

#include "IOHprofiler_problem.h"
#include "wmodels.hpp"

class LeadingOnes_Ruggedness3 : public IOHprofiler_problem<int> {
public:
  LeadingOnes_Ruggedness3(int instance_id = DEFAULT_INSTANCE, int dimension = DEFAULT_DIMENSION) {

    IOHprofiler_set_instance_id(instance_id);
    IOHprofiler_set_problem_name("LeadingOnes_Ruggedness3");
    IOHprofiler_set_problem_type("pseudo_Boolean_problem");
    IOHprofiler_set_number_of_objectives(1);
    IOHprofiler_set_lowerbound(0);
    IOHprofiler_set_upperbound(1);
    IOHprofiler_set_best_variables(1);
    IOHprofiler_set_number_of_variables(dimension);
  }
  
  ~LeadingOnes_Ruggedness3() {}

  std::vector<double> info;
  void prepare_problem() {
    info = ruggedness3(IOHprofiler_get_number_of_variables());
  }

  double internal_evaluate(const std::vector<int> &x) {
    
    int n = x.size();
    int result = 0;
    for (int i = 0; i != n; ++i) {
      if(x[i] == 1) {
        result = i + 1;
      } else {
        break;
      }
    }
    result = this->info[(int)(result + 0.5)];
    return (double)result;
  }

  static LeadingOnes_Ruggedness3 * createInstance(int instance_id = DEFAULT_INSTANCE, int dimension = DEFAULT_DIMENSION) {
    return new LeadingOnes_Ruggedness3(instance_id, dimension);
  }
};

#endif
