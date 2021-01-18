//
//  LearnOutputHandling.h
//  Rmixmod
//
//  Created by Rémi Lebret on 06/02/12.
//  Copyright (c) 2012 CNRS. All rights reserved.
//

#ifndef Rmixmod_LearnOutputHandling_h
#define Rmixmod_LearnOutputHandling_h

#include "OutputHandling.h"
#include "mixmod/Utilities/mixmod.h"

/** base class for handling the mixmod outputs and filling them in
 *  a Rcpp list 
 */
class LearnOutputHandling : public OutputHandling
{
public:
  
  /** Default constructor
   *  @param cMOutput the model estimated by mixmod
   *  @param output the Rcpp list to fill in
   **/
  LearnOutputHandling( XEM::LearnModelOutput* lMOutput
                     , Rcpp::S4& xem
                     , const XEM::DataType dataType
                     , Rcpp::CharacterVector const & Rcriterion
                     , std::vector<int64_t> labels
                     );
  LearnOutputHandling( XEM::LearnModelOutput* lMOutput
                     , Rcpp::S4& xem
                     , const XEM::DataType dataType
                     , std::vector <XEM::CriterionName> const & iCriterion
                     , std::vector<int64_t> labels
                     );
  
  /** destructor */
  ~LearnOutputHandling();
  
};

#endif
