//
//  OutputHandling.h
//  Rmixmod
//
//  Created by Rémi Lebret on 06/02/12.
//  Copyright (c) 2012 CNRS. All rights reserved.
//

#ifndef Rmixmod_OutputHandling_h
#define Rmixmod_OutputHandling_h

#include <Rcpp.h>
#include <stdint.h>
#include "mixmod/Utilities/Util.h"
#include "mixmod/Kernel/Parameter/BinaryParameter.h"
#include "mixmod/Kernel/Parameter/GaussianEDDAParameter.h"
#include "mixmod/Kernel/Parameter/CompositeParameter.h"

/** base class for handling the mixmod outputs and filling them in
 *  a Rcpp list */
class OutputHandling
{
  
public:
  /** Default constructor
   *  @param cMOutput the model estimated by mixmod
   *  @param output the Rcpp list to fill in
   **/
  OutputHandling( XEM::ModelOutput* MOutput, Rcpp::S4& xem, const XEM::DataType dataType);
  OutputHandling( XEM::Parameter* par, Rcpp::S4& xem, const XEM::DataType dataType, int nbCluster);  
  
  /** destructor */
  ~OutputHandling();
  
private:
  /** set gaussian paramaters */
  void setGaussianParameter(XEM::GaussianEDDAParameter const * parArg=nullptr);
  
  /** set multinomial parameters */
  void setMultinomialParameter(XEM::BinaryParameter const * parArg=nullptr);
  
  /** set composite parameters */
  void setCompositeParameter(XEM::CompositeParameter *parArg=nullptr);

protected:
  /** A pointer on the MIXMOD output */
  XEM::ModelOutput* MOutput_;
  
  /** A reference on the R output */
  Rcpp::S4& xem_;
  
  /** Number of clusters*/
  int nbCluster_;
  
  /** Number of variables */
  int64_t nbVariable_;
};


#endif
