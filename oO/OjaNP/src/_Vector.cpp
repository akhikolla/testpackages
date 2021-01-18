#include "stdafx.h"
#include "_Vector.h"

#include <stdarg.h>
#include <math.h>

using namespace std;

Vector::Vector()
{
  this->n = 0;
  values = new double[1];
}

Vector::Vector(int n)
{
  this->n = n;
  values = new double[n];
  for(int i = 0; i < n; i++)
    values[i] = 0;
}

Vector::Vector(int n, double* v)
{
  this->n = n;
  values = new double[n];
  for(int i = 0; i < n; i++)
    values[i] = v[i];
}

Vector::Vector(double n, ...)
{
  this->n = (int)n;
  values = new double[this->n];

  va_list params;
  va_start(params, n);              // Aufruf an Initialisierungmakro

  for (int i = 0 ; i < n ; i++)
  {
    double v = va_arg(params, double); // Extrahiere einen Parameter
    values[i] = v;
  }
  va_end(params);                        // Schliessmakro;
}

Vector::~Vector()
{
  delete[] values;
}

bool Vector::setValue(int i, double v)
{
  if( i >= n)
  {
    //cout << "Vector::setValue - index too large!" << endl;
    return false;
  }

  values[i] = v;
  return true;
}

bool Vector::setValues(Vector& v)
{
  int dim = (int)v.getSize();
  if(dim != n)
  {
    //cout << "Vector::setValues - sizes do not match!" << endl;
    return false;
  }

  for(int i = 0; i < n; i++)
    values[i] = v.getValue(i);
  return true;
}

bool Vector::addValue(int i, double v)
{
  if( i >= n)
  {
    //cout << "Vector::addValue - index too large!" << endl;
    return false;
  }

  values[i] += v;
  return true;
}

double Vector::getValue(int i)
{
  if(i >= n)
  {
    //cout << "Vector::getValue - index too large!" << endl;
  }

  return values[i];
}

int Vector::getSize()
{
  return n;
}

double Vector::getLength()
{
  double l = 0;
  for(int i = 0; i < n; i++)
  {
    l += values[i] * values[i];
  }
  return sqrt(l);
}

void Vector::print()
{
  //cout << "(";
  //for(int i = 0; i < n - 1; i++)
    //cout << values[i] << ", ";
  
 // cout << values[n-1] << ")";// << endl;
}

void Vector::normalize()
{
  double l = getLength();
  for(int i = 0; i < n; i++)
    values[i] /= l;
}

void Vector::setLength(double l)
{
  normalize();
  for(int i = 0; i < n; i++)
    values[i] *= l;
}

