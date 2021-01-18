#ifndef VECTOR_H
#define VECTOR_H

#pragma once

#include <iostream>

class Vector
{
public:
  /** constructors **/
  Vector();
  Vector(int n);
  Vector(int n, double* v); 
  Vector(double n, ...);

  /** deconstructors **/
  ~Vector();

  /** set value **/
  bool setValue(int i, double v);
  bool setValues(Vector& v);

  bool addValue(int i, double v);

  /** get value **/
  double getValue(int i);
  int getSize();
  void normalize();
  double getLength();
  void setLength(double l);
  
  /** print **/
  void print();

  Vector& operator+=(Vector& other)
  {
    for(int i = 0; i < n; i++)
      values[i] += other.getValue(i);
    return *this;
  }

  //copy Konstruktor
  Vector(const Vector & t)
  {
	n=t.n;
	values=new double[n];
	for(int i=0;i<n;i++)
		values[i]=t.values[i];
  }

  //Zuweisungs konstruktor
  Vector & operator=(const Vector &t)
  {
	  if (this!=&t)
	  {
		if (n!=t.n)
		{
			delete []values;
			n=t.n;
			values=new double[n];
		}
		for(int i=0;i<n;i++)
			values[i]=t.values[i];
	  }
      return *this;
  };


private:
  /** size of vector **/
  int n;

  /** values of vector **/
  double* values;
};

#endif

