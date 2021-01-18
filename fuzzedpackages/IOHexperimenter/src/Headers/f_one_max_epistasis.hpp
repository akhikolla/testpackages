/// \file f_one_max_epistasis.hpp
/// \brief cpp file for class f_one_max_epistasis.
///
/// This file implements a OneMax problem with epistasis transformation method from w-model.
/// The parameter v is chosen as 4.
///
/// \author Furong Ye
/// \date 2019-06-27
#ifndef _F_ONE_MAX_EPISTASIS_H
#define _F_ONE_MAX_EPISTASIS_H

#include "IOHprofiler_problem.h"
#include "wmodels.hpp"

class OneMax_Epistasis : public IOHprofiler_problem<int> {
public:
  OneMax_Epistasis(int instance_id = DEFAULT_INSTANCE, int dimension = DEFAULT_DIMENSION) {
    IOHprofiler_set_instance_id(instance_id);
    IOHprofiler_set_problem_name("OneMax_Epistasis");
    IOHprofiler_set_problem_type("pseudo_Boolean_problem");
    IOHprofiler_set_number_of_objectives(1);
    IOHprofiler_set_lowerbound(0);
    IOHprofiler_set_upperbound(1);
    IOHprofiler_set_number_of_variables(dimension);
  }

  ~OneMax_Epistasis() {}

  void customize_optimal() {
    IOHprofiler_set_optimal(IOHprofiler_get_number_of_variables());
  }

  double internal_evaluate(const std::vector<int> &x) {
    
    std::vector<int> new_variables = epistasis(x,4);
    int n = new_variables.size();
    int result = 0;
    for (int i = 0; i != n; ++i) {
      result += new_variables[i];
    }
    return (double)result;
  }

  static OneMax_Epistasis * createInstance(int instance_id = DEFAULT_INSTANCE, int dimension = DEFAULT_DIMENSION) {
    return new OneMax_Epistasis(instance_id, dimension);
  }
};

#endif
