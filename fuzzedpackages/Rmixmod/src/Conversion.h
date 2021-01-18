/*
 * Conversion.h
 *
 *  Created on: 8 juin 2011
 *      Author: aude l
 */



#ifndef CONVERSIONRCPP_H_
#define CONVERSIONRCPP_H_

#include <Rcpp.h>
#include <vector>
#include <stdint.h>

#include "mixmod/Utilities/mixmod.h"

namespace Conversion
{
  /**
   * Convert a C matrix in a Rcpp Matrix
   * @param nbSample numbor of rows of the matrix
   * @param pdDimension number of column of the matrix
   * @param matrix the matrix to convert
   * @return the Rcpp Matrix
   */
  Rcpp::NumericMatrix CMatrixToRcppMatrix(int64_t nbSample, int64_t pdDimension, double** matrix);
  
  /**
   * Convert a matrix vector of vector in a Rcpp Matrix
   * @param source the matrix to convert
   * @return the Rcpp Matrix
   */
  Rcpp::NumericMatrix XEMMatrixToRcppMatrix(std::vector<std::vector<double> > const& source);
  /**
   * Convert a c vector in a Rcpp Vector
   * @param dim the dimension of the vector
   * @param vector the C vector
   * @return the Rcpp Vector
   */
  Rcpp::NumericMatrix CMatrixToRcppMatrixForInt(int64_t nbSample, int64_t pbDimension, int64_t** matrix);
  /**
   * Convert a c vector in a Rcpp Vector
   * @param dim the dimension of the vector
   * @param vector the C vector
   * @return the Rcpp Vector
   */
  Rcpp::NumericVector CVectorToRcppVector(int64_t dim, double* vector);
  /**
   * Create
   * @param nbCluster number of cluster
   * @param labels array with the label of the individuals
   * @return A Rcpp matrix with the partition of in a binary form
   */
  Rcpp::NumericVector CVectorToRcppVectorForInt(int64_t dim, int64_t* vector);

  /**
   * Create
   * @param nbCluster number of cluster
   * @param labels array with the label of the individuals
   * @return A Rcpp matrix with the partition of in a binary form
   */
  Rcpp::NumericVector VectorToRcppVectorForInt( std::vector<int64_t> const & vector);

  Rcpp::NumericMatrix LabelToPartition( int64_t nbCluster, std::vector<int64_t> const & labels);
  /**
   *  convert a numeric dataset to a XEM Gaussian data set
   * @param data the input data set
   * @return a pointer on a XEMGaussianData set
   */
  XEM::GaussianData * DataToXemGaussianData(Rcpp::NumericMatrix& data);
  /**
   *  convert a numeric dataset to a XEM Binary data set
   * @param data the input data set
   * @param factor the input factor
   * @return a pointer on a XEMBinaryData set
   */
  XEM::BinaryData * DataToXemBinaryData(Rcpp::NumericMatrix& data, Rcpp::NumericVector& factor);
  /**
   *  convert a numeric dataset to a XEM composite data set
   * @param data the input data set
   * @param factor the input factor
   * @return a pointer on a XEMCompositeData set
   */
  XEM::CompositeData * DataToXemCompositeData(Rcpp::NumericMatrix& data, Rcpp::NumericVector& factor);
  
  /*
   * convert a R numeric vector to a C array
   * @param in the R object
   * @return a pointer to a double* array
   */
  double * RcppVectorToCArray(Rcpp::NumericVector& in);
  
  /*
   * convert a R numeric matrix to a C 2D array
   * @param in the R object
   * @return a pointer to a double** array
   */
  double ** RcppMatrixToC2DArray(Rcpp::NumericMatrix& in);
  
  /*
   * convert a R numeric matrix to a C 2D array
   * @param in the R object
   * @return a pointer to a double** array
   */
  double *** RcppListOfMatrixToC3DArray(Rcpp::List& in);
}


#endif /* CONVERSIONRCPP_H_ */
