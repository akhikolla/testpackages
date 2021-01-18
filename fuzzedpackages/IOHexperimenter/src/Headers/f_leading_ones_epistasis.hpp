/// \file f_leading_ones_epistasis.hpp
/// \brief cpp file for class f_leading_ones_epistasis.
///
/// This file implements a LeadingOnes problem with epistasis transformation method from w-model.
/// The parameter v is chosen as 4.
///
/// \author Furong Ye
/// \date 2019-06-27
#ifndef _F_LEADING_ONES_EPISTASIS_H
#define _F_LEADING_ONES_EPISTASIS_H

#include "IOHprofiler_problem.h"
#include "wmodels.hpp"

class LeadingOnes_Epistasis : public IOHprofiler_problem<int> {
public:
  LeadingOnes_Epistasis(int instance_id = DEFAULT_INSTANCE, int dimension = DEFAULT_DIMENSION) {    
    IOHprofiler_set_instance_id(instance_id);
    IOHprofiler_set_problem_name("LeadingOnes_Epistasis");
    IOHprofiler_set_problem_type("pseudo_Boolean_problem");
    IOHprofiler_set_number_of_objectives(1);
    IOHprofiler_set_lowerbound(0);
    IOHprofiler_set_upperbound(1);
    IOHprofiler_set_number_of_variables(dimension);
  }

  ~LeadingOnes_Epistasis() {}

  void customize_optimal() {
    IOHprofiler_set_optimal(IOHprofiler_get_number_of_variables());
  }
  
  double internal_evaluate(const std::vector<int> &x) {

    std::vector<int> new_variables = epistasis(x,4);
    int n = new_variables.size();
    int result = 0;
    for (int i = 0; i != n; ++i) {
      if(new_variables[i] == 1) {
        result = i + 1;
      } else {
        break;
      }
    }
    return (double)result;
  }
  
  static LeadingOnes_Epistasis * createInstance(int instance_id = DEFAULT_INSTANCE, int dimension = DEFAULT_DIMENSION) {
    return new LeadingOnes_Epistasis(instance_id, dimension);
  }
};

#endif
