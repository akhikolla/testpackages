#ifndef MATRIX2D_H
#define MATRIX2D_H

#pragma once

#include "_Vector.h"
#include <iostream>

class Matrix2D
{
public:
  /** constructors **/
  Matrix2D();
  Matrix2D(int m, int n);
//  Matrix2D(double m, double n,...);
  Matrix2D(const Matrix2D & reducedMat,int* validRows, int reducedM, int reducedN); 	

  /** deconstructors **/
  ~Matrix2D();

  /** set values of matrix **/
  bool setRow(int i, Vector& v);
  bool setColumn(int j, Vector& v);
  bool setRow(int i, double* v); 
  bool setColumn(int j, double* v);
  bool setColumn(int j, double p, Vector& v);
  bool setValue(int i, int j, double v);

  /** get values of matrix **/
  Vector getRow(int i);
  Vector *getColumn(int j);
  double getValue(int i, int j);

  int getNumberColumns();
  int getNumberRows();

  void deleteRow(int i, Matrix2D* matrix);
  void deleteRowAndColumn(int i, int j, Matrix2D* matrix);
  void deleteColumn(int j, Matrix2D* matrix);

  /** print **/
  void print() const;


  /** calculate determinante **/
  double determinant(int* validRows, int reducedM, int reducedN);
  double determinant2(int *validRows, int reducedM, int reducedN);

  /** calculate variance of column vectors **/
  double getVariance();

  /** transpose matrix **/
  Matrix2D & transpose()const; // the result is the new A^T matrix
  
  /** mult matrix  **/ 
  Matrix2D & mult(const Matrix2D & other_mat) const;
  
  /** set the Matrix to Identity - need to be a sqaure Matrix **/
  void loadIdentity();

  /* The following procedure calculates an balanced matrix which has 
   * the same eigenvalues as the original one.
   * */
  Matrix2D & balance() const;
  
  /* The following procedure generates a copy of this matrix */
  Matrix2D & getACopy()const;
  
  /* The following procedure returns a new Matrix which is in 
   * Hessenberg Form and has the same eigenvalues as this one */
  Matrix2D & reduceToHessenbergForm(bool doBalancing = false) const;
  
  void printReducedToMapleFile (char * filename, char * matrix_variable_name , int * validRows, int reducedM, int reducedN);
  void printToMapleFile(char * filename,char * matrix_variable_name);
  
  /* the following procedure calculates the determinat of an reduced Matrix using gaussian elimination 
   * you may set the options between 0 ... 2 
   * piv_option 0:   normal pivoting (the first element in the pivot_col neq 0) 
   * piv_option 1:   column pivoting (the max.  element in the pivot_col)
   * piv_option 2:   total pivoting  (the max.  element of the whole submatrix) 
   * */
  double reduced_determinant_pivoting(int * validRows,int reducedM, int reducedN,int option = 0);

  /* the following procedure swaps two rows between >= start_col_index and < end_col_index */
  void swap_rows(int row_1 , int row_2, int start_col_index , int end_col_index);
  
  /* swaps two columns between >= start_row_index and < end_row_index */
  void swap_columns(int col_1, int col_2, int start_row_index, int end_row_index);
  
  /* pivots the submatrix (r+1 ... m) X (c+1 ... n) with the pivot element a_(r,c) */
  void do_pivoting( int r , int c);
  //double do_pivoting( int r , int c);
  
  /* the following procedure the determinat of an Matrix using gaussian elimination 
   * you may set the options between 0 ... 2 
   * piv_option 0:   normal pivoting (the first element in the pivot_col neq 0) 
   * piv_option 1:   column pivoting (the max.  element in the pivot_col)
   * piv_option 2:   total pivoting  (the max.  element of the whole submatrix) 
   * 
   * if copydata = false your matrix will be replaced by an pivoted one
   * */
  double determinant_pivoting(bool copydata = false, int piv_option = 0);
  
private:
  /** handle values of reduced matrix **/
  //inline double getValueOfReducedMat(int i, int j, int reducedM, int reducedN);
  
  
  /** size of the matrix **/
  int m; //rows
  int n; //columns

  /** values of the matrix **/
  double** values;
};

#endif

