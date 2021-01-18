/*
 * Conversion.cpp
 *
 *  Created on: 8 juin 2011
 *      Author: aude
 */

#include "Conversion.h"

#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/CompositeData.h"

namespace Conversion
{


//To spin off double** in Rcpp::NumericMatrix
Rcpp::NumericMatrix CMatrixToRcppMatrix(int64_t nbSample, int64_t pbDimension, double** matrix)
{
  //NumericMatrix for matrix
  Rcpp::NumericMatrix matrixOutput(nbSample,pbDimension);
  for(int i=0;i<nbSample;i++)
  {
    for(int j=0;j<pbDimension;j++)
    {
      matrixOutput(i,j) = matrix[i][j];
    }
  }
  return matrixOutput;
}

//To spin off vector of vector in Rcpp::NumericMatrix
Rcpp::NumericMatrix XEMMatrixToRcppMatrix(std::vector<std::vector<double> > const& source)
{
  int64_t nbSample = source.size();
  int64_t pbDimension = source[0].size();

  //NumericMatrix for matrix
  Rcpp::NumericMatrix matrixOutput(nbSample,pbDimension);
  for (int64_t i=0; i<nbSample; i++)
  {
    for (int64_t j=0; j<pbDimension; j++)
    {
      matrixOutput(i,j) = source[i][j];
    }
  }
 return matrixOutput;
}

//To spin off int** in Rcpp::NumericMatrix
Rcpp::NumericMatrix CMatrixToRcppMatrixForInt(int64_t nbSample, int64_t pbDimension, int64_t** matrix)
{
    //NumericMatrix for matrix
    Rcpp::NumericMatrix matrixOutput(nbSample,pbDimension);
    for(int i=0;i<nbSample;i++)
    {
      for(int j=0;j<pbDimension;j++)
      {
        matrixOutput(i,j) = matrix[i][j];
      }
    }
    return matrixOutput;
}

//To spin off double* in Rcpp::NumericVector
Rcpp::NumericVector CVectorToRcppVector(int64_t dim, double* vector)
{
  // NumericVector for vector
  Rcpp::NumericVector vectorOutput(dim);
  for(int i=0;i<dim;i++)
  {
      vectorOutput(i) = vector[i];
  }
  return vectorOutput;
}

//To spin off int* in Rcpp::NumericVector
Rcpp::NumericVector CVectorToRcppVectorForInt(int64_t dim, int64_t* vector)
{
  // NumericVector for vector
  Rcpp::NumericVector vectorOutput(dim);
  for(int i=0;i<dim;i++)
  {
      vectorOutput(i) = vector[i];
  }
  return vectorOutput;
}

//To spin off int* in Rcpp::NumericVector
Rcpp::NumericVector VectorToRcppVectorForInt( std::vector<int64_t> const & data)
{
  // get dimension
  std::vector<int64_t>::size_type dim = data.size();
  //
  Rcpp::NumericVector vectorOutput(dim);
  for(unsigned int i=0;i<dim;i++)
  {
    vectorOutput[i] = (double)data[i];
  }
  return vectorOutput;
}

/*To spin off information for partition */
Rcpp::NumericMatrix LabelToPartition( int64_t nbCluster, std::vector<int64_t> const & labels)
{
  int nbSample = labels.size() ;

  // Matrix for partition
  Rcpp::NumericMatrix partitionOutput(nbSample,(int)nbCluster);
  for(int i=0;i<nbSample;i++)
  {
    int ind = labels[i]-1;
    for(int j=0;j<nbCluster;j++)
    {
      (j == ind) ? partitionOutput(i,j) = 1 : partitionOutput(i,j) = 0;
    }
  }
  return partitionOutput;
}

/*
 *  convert a gaussian dataset to a Xem gaussian data set
 * @param data the input data set
 * @return a pointer on a GaussianData set
 */
XEM::GaussianData * DataToXemGaussianData(Rcpp::NumericMatrix& data)
{
  // wrap int variables to int64_t variables
  int64_t nbSample64 = data.nrow(), nbVariable64 = data.ncol();

  double ** matrix = new double * [nbSample64];
  for (int i=0; i<nbSample64; i++)
  {
   matrix[i] = new double[nbVariable64];
   for (int j=0; j<nbVariable64; j++)
     matrix[i][j] = data(i,j);
  }
  // create XEMGaussianData
  XEM::GaussianData*  gData = new XEM::GaussianData(nbSample64, nbVariable64, matrix);

  // release memory
  for (int64_t i=0; i<nbSample64; i++)
  { delete [] matrix[i];}
  delete [] matrix;
  matrix = 0;

  return gData;
}

/*
 *  convert a binary dataset to a Xem binary data set
 * @param data the input data set
 * @return a pointer on a XEMBinaryData set
 */
XEM::BinaryData * DataToXemBinaryData(Rcpp::NumericMatrix& data, Rcpp::NumericVector& factor)
{
  // wrap int variables to int64_t variables
  int64_t nbSample64 = data.nrow(), nbVariable64 = data.ncol();

  int64_t ** matrix = new int64_t * [nbSample64];
  for (int i=0; i<nbSample64; i++)
  {
   matrix[i] = new int64_t[nbVariable64];
   for (int j=0; j<nbVariable64; j++)
   {
     matrix[i][j] = (int64_t)data(i,j);
   }
  }

  // compute the modalities of the data
  std::vector<int64_t> nbModality(nbVariable64);
  for(int j=0;j<nbVariable64;j++)
  { nbModality[j]=factor[j]; }
  
  // create XEMBinaryData
  XEM::BinaryData*  bData = new XEM::BinaryData(nbSample64, nbVariable64, nbModality, matrix);

  // release memory
  for (int64_t i=0; i<nbSample64; i++)
  { delete [] matrix[i];}
  delete [] matrix;
  matrix = 0;

  return bData;
}

XEM::CompositeData * DataToXemCompositeData(Rcpp::NumericMatrix& data, Rcpp::NumericVector& factor)
{
  int64_t nbSample = data.nrow();
  int64_t totalColumns = data.ncol();
  int64_t nbcolbinary =0,nbcolgaussian=0;
  for (int i = 0; i < totalColumns; ++i) {
    if(factor(i)>0)
      nbcolbinary++;
    else
      nbcolgaussian++;
  }

  int64_t ** bmatrix = new int64_t * [nbSample];
  double ** gmatrix = new double * [nbSample];
  std::vector<int64_t> nbModality(nbcolbinary);
  int64_t bcols=0,gcols=0;
  for (int i=0; i<nbSample; i++)
  {
   bmatrix[i] = new int64_t[nbcolbinary];
   gmatrix[i] = new double[nbcolgaussian];
   for (int j=0; j<totalColumns; j++)
   {
     if(factor(j)>0)
     {
       bmatrix[i][bcols] = (int64_t)data(i,j);
       nbModality[bcols++] = factor(j);
     }
     else
       gmatrix[i][gcols++] = (double)data(i,j);
   }
   bcols=0;gcols=0;
  }

  XEM::BinaryData*  bData = new XEM::BinaryData(nbSample, nbcolbinary, nbModality, bmatrix);
  XEM::GaussianData*  gData = new XEM::GaussianData(nbSample, nbcolgaussian, gmatrix);
  XEM::CompositeData * cData = new XEM::CompositeData(bData,gData);

  // release memory
  for (int64_t i=0; i<nbSample; i++)
  {
    delete [] gmatrix[i];
    delete [] bmatrix[i];

  }
  delete [] gmatrix;
  delete [] bmatrix;
  gmatrix = NULL;
  bmatrix = NULL;

  return cData;

}

  
/*
 * convert a R numeric vector to a C array
 * @param in the R object
 * @return a pointer to a double* array
 */
double * RcppVectorToCArray(Rcpp::NumericVector& in)
{
  // get the vector dimension
  int size = in.size();
  // allocate vector
  double * out = new double[size];
  for (int i=0; i<size; i++)
  {
    out[i] = in[i];
  }
  // return C vector
  return out;
}

  
/*
 * convert a R numeric matrix to a C 2D array
 * @param in the R object
 * @return a pointer to a double** array
 */
double ** RcppMatrixToC2DArray(Rcpp::NumericMatrix& in)
{
  // get the vector dimension
  int nrow = in.nrow();
  int ncol = in.ncol();
  
  // allocate columns
  double ** out = new double * [nrow];
  for (int i=0; i<nrow; i++)
  {
    // allocate rows
    out[i] = new double[ncol];
    // get values
    for (int j=0; j<ncol; j++)
    {
      out[i][j] = in(i,j);
    }
  }
  // return C matrix
  return out;
}
  
/*
 * convert a R numeric matrix to a C 2D array
 * @param in the R object
 * @return a pointer to a double** array
 */
double *** RcppListOfMatrixToC3DArray(Rcpp::List& in)
{
  // get the vector dimension
  int size = in.size();
  
  // allocate first dimension
  double *** out = new double ** [size];
  
  // loop over list
  for (int k=0; k<size; k++) 
  {
    Rcpp::NumericMatrix mat = SEXP(in[k]);
    
    // get dimensions
    int nrow = mat.nrow();
    int ncol = mat.ncol();
    
    // allocate rows
    out[k] = new double * [nrow];
    
    for (int i=0; i<nrow; i++)
    {
      // allocate rows
      out[k][i] = new double[ncol];
      // get values
      for (int j=0; j<ncol; j++)
      {
        out[k][i][j] = mat(i,j);
      }
    }
  }

  // return C matrix
  return out;
}


}


















