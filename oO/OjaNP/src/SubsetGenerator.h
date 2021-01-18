#ifndef SUBSETGENERATOR_H
#define SUBSETGENERATOR_H

#pragma once

#include "Matrix2D.h"
#include "MersenneTwister.h"

extern int FIXED_RANDOM_SUBSET;


class SubsetGenerator
{
public:
  virtual ~SubsetGenerator(){};

  virtual bool* getNextSubset(Matrix2D* m, Matrix2D& points)=0;
  virtual bool* getNextSubset(Matrix2D* m, Matrix2D* m2, Matrix2D& points)=0;

  virtual void reset()=0;
  virtual void reset(Matrix2D* m, Matrix2D& points, Vector& mu)=0;
  virtual void reset(Matrix2D* m, Matrix2D* m2, Matrix2D& points, Vector& mu)=0;

  virtual bool isFinished()=0;

};

class AllSubsets:public SubsetGenerator
{
public:
  
  AllSubsets(int numberOfElements, int sizeOfSubsample);
  ~AllSubsets();

  bool* getNextSubset(Matrix2D* m, Matrix2D& points);
  bool* getNextSubset(Matrix2D* m, Matrix2D* m2, Matrix2D& points);

  void reset();
  void reset(Matrix2D* m, Matrix2D& points, Vector& mu);
  void reset(Matrix2D* m, Matrix2D* m2, Matrix2D& points, Vector& mu);

  bool isFinished()
  {
	  return (finished);
  };

private:
  bool finished;

  enum States
  {
    _1a = 0,
    _1b = 1,
    _1c = 2,
    _2a = 3,
    _2b = 4,
    _2c = 5,
    _start = 6
  };
  States state;

  int sizeOfSubsample;
  int numberOfElements;

  int* p;
  int n;
  int h;


  bool* contains;
  int* selectedElements; // = -1, not selected - != -1 position in matrix
};

class RandomSubsets:public SubsetGenerator
{
public:
  RandomSubsets(int _numberOfElements, int _sizeOfSubsample,int _seed, int _numberOfSubsets)
  {
    sizeOfSubsample=_sizeOfSubsample;
    numberOfElements=_numberOfElements;
	seed=_seed;
	r=new MTRand(seed);
	contains=new bool[numberOfElements];
	numberOfSubsets=_numberOfSubsets;
  };
  
  ~RandomSubsets();
  
  void reset();
  void reset(Matrix2D* m, Matrix2D& points, Vector& mu);
  void reset(Matrix2D* m, Matrix2D* m2, Matrix2D& points, Vector& mu);

  bool* getNextSubset(Matrix2D* m, Matrix2D& points);
  bool* getNextSubset(Matrix2D* m, Matrix2D* m2, Matrix2D& points);
  
  bool isFinished()
  {
	  return (count>=numberOfSubsets);
  };

private:
	int seed;
	MTRand *r;
    int sizeOfSubsample;
    int numberOfElements;
	int numberOfSubsets;
	int count;

	void getRandomSubset(bool* subset);

  bool* contains;
};

#endif

