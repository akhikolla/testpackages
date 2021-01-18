//*--------------------------------------------------------------------------*//
//  File:               hypermatrix.h
//  Created by:         Pavlo Mozharovskyi
//  First released:     15.10.2016
//
// A binary hypermatrix that contains indexed cells only as bits.
//
//  Subsequent changes are listed below:
//  15.10.2016 (Pavlo Mozharovskyi): First version released.
//*--------------------------------------------------------------------------*//

#include "stdafx.h"

template<typename BlockType> class hypermatrix{
protected:
  int d, blockSize;
  int* ns;
  unsigned long long size;
  void clean(){
    d = 0;
    blockSize = 0;
    ns = 0;
    body = 0;
    size = 0;
  }
  hypermatrix(){
    clean();
  }
public:
  BlockType* body;
  void free(){
    if (body){
      delete[] ns;
      delete[] body;
      clean();
    }
  }
  int clear(){
    if (!body){
      return 1;
    }
    for (int i = 0; i < size; i++){
      body[i] = 0;
    }
    return 0;
  }
  ~hypermatrix(){
    free();
  }
  void getSize(int* d, int** ns){
    *d = this->d;
    *ns = new int[this->d];
    memcpy(*ns, this->ns, this->d * sizeof(int));
  }
};

template<typename BlockType> class binaryHypermatrix :
  public hypermatrix<BlockType>{
private:
  unsigned long long nNonzero;
 	binaryHypermatrix();
public:
	binaryHypermatrix(int, int*);
	bool setIfNotSet(int*);
	bool isSet(int*);
	unsigned long long sum();
	void clear(){
	  hypermatrix<BlockType>::clear();
	  nNonzero = 0;
	}
};

template<typename BlockType> binaryHypermatrix<BlockType>::
	binaryHypermatrix(int d, int* ns){
	this->d = d;
	this->ns = new int[d];
	memcpy(this->ns, ns, d * sizeof(int));
	this->blockSize = 8 * sizeof(BlockType);
	this->size = 1;
	for (int i = 0; i < d; i++){
	  this->size *= ns[i];
	}
	this->size = this->size / this->blockSize + 1;
	this->body = new BlockType[this->size]();
	this->clear();
	nNonzero = 0;
}

template<typename BlockType> bool binaryHypermatrix<BlockType>::
	setIfNotSet(int* index){
	unsigned long long bitNumber = 0;
	unsigned long long curOffsetStep = 1;
	for (int i = this->d - 1; i >= 0; i--){
		bitNumber += index[i] * curOffsetStep;
	  curOffsetStep *= this->ns[i];
	}
	unsigned long long blockNumber = bitNumber / this->blockSize;
	int bitInBlock = bitNumber % this->blockSize;
	BlockType mask = (BlockType)1 << bitInBlock;
	if (this->body[blockNumber] & mask){
		return false;
	}else{
		this->body[blockNumber] |= mask;
	  nNonzero++;
		return true;
	}
}

template<typename BlockType> bool binaryHypermatrix<BlockType>::
  isSet(int* index){
  unsigned long long bitNumber = 0;
  unsigned long long curOffsetStep = 1;
  for (int i = this->d - 1; i >= 0; i--){
    bitNumber += index[i] * curOffsetStep;
    curOffsetStep *= this->ns[i];
  }
  unsigned long long blockNumber = bitNumber / this->blockSize;
  int bitInBlock = bitNumber % this->blockSize;
  BlockType mask = (BlockType)1 << bitInBlock;
  if (this->body[blockNumber] & mask){
    return true;
  }else{
    return false;
  }
}

template<typename BlockType> unsigned long long binaryHypermatrix<BlockType>::
  sum(){
    return nNonzero;
}

template<typename UnitType> class typeHypermatrix :
  public hypermatrix<UnitType>{
public:
  typeHypermatrix(int d, int* ns);
  int add(typeHypermatrix<UnitType> &matrix);
  int addStandardized(typeHypermatrix<UnitType> &matrix);
  int standardize();
  UnitType getTotal();
  UnitType getNNonzero();
};

template<typename UnitType> typeHypermatrix<UnitType>::
  typeHypermatrix(int d, int* ns) : hypermatrix<UnitType>(){
  // std::cout << "Constructor typeHypermatrix(int d, int* ns) called." << std::endl;
  this->d = d;
  this->ns = new int[this->d];
  memcpy(this->ns, ns, d * sizeof(int));
  this->blockSize = 8 * sizeof(UnitType);
  this->size = 1;
  for (int i = 0; i < this->d; i++){
    this->size *= ns[i];
  }
  this->body = new UnitType[this->size]();
  this->clear();
}

template<typename UnitType> int typeHypermatrix<UnitType>::
  add(typeHypermatrix<UnitType> &matrix){
  if (!this->body || this->size != matrix.size){
    return 1;
  }
  for (int i = 0; i < this->size; i++){
    this->body[i] += matrix.body[i];
  }
  return 0;
}

template<typename UnitType> int typeHypermatrix<UnitType>::
  addStandardized(typeHypermatrix<UnitType> &matrix){
    //Rcout << this->size << " " << matrix.size << std::endl;
    if (!this->body || this->size != matrix.size){
      return 1;
    }
    UnitType total = matrix.getTotal();
    for (int i = 0; i < this->size; i++){
      this->body[i] += matrix.body[i] / total;
    }
    return 0;
  }

template<typename UnitType> int typeHypermatrix<UnitType>::
  standardize(){
    //Rcout << this->size << " " << matrix.size << std::endl;
    UnitType total = getTotal();
    for (int i = 0; i < this->size; i++){
      this->body[i] /= total;
    }
    return 0;
  }

template<typename UnitType> UnitType typeHypermatrix<UnitType>::
  getTotal(){
  UnitType sum = 0;
  for (int i = 0; i < this->size; i++){
    sum += this->body[i];
  }
  return sum;
}

template<typename UnitType> UnitType typeHypermatrix<UnitType>::
  getNNonzero(){
  UnitType num = 0;
  for (int i = 0; i < this->size; i++){
    if (this->body[i] != 0)
    num++;
  }
  return num;
}
