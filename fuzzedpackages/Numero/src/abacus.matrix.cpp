/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
Matrix::Matrix() {
  this->buffer = new MatrixBuffer();
}

/*
 *
 */
Matrix::Matrix(const Matrix& t) {
  this->buffer = new MatrixBuffer(t.buffer);
}

/*
 *
 */
void
Matrix::operator=(const Matrix& t) {
  if(this == &t) return;
  MatrixBuffer* p = (MatrixBuffer*)buffer; delete p;
  this->buffer = new MatrixBuffer(t.buffer);
}

/*
 *
 */
Matrix::~Matrix() {
  MatrixBuffer* p = (MatrixBuffer*)(this->buffer);
  delete p;
}
