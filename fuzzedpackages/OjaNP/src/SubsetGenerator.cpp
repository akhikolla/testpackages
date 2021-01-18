#include "stdafx.h"
#include "SubsetGenerator.h"

#include <iostream>
using namespace std;

AllSubsets::AllSubsets(int numberOfElements, int sizeOfSubsample)
{
  this->sizeOfSubsample = sizeOfSubsample;
  this->numberOfElements = numberOfElements;
  contains = new bool[numberOfElements];
  p = new int[numberOfElements];
  selectedElements = new int[numberOfElements];

  reset();
}

AllSubsets::~AllSubsets()
{
  delete[] contains;
  delete[] p;
  delete[] selectedElements;
}


void AllSubsets::reset(Matrix2D* m, Matrix2D& points, Vector& mu)
{
  reset();

  //create first subset matrix
  for(int j = 0; j <= sizeOfSubsample; j++)
    m->setValue(0, j, 1.0);

  for(int i = 1; i <= sizeOfSubsample; i++)
  {
    m->setValue(i, 0, mu.getValue(i-1));
  }

  for(int j = 0; j < sizeOfSubsample - 1; j++)
  {
    selectedElements[j] = j;
    for(int i = 0; i < sizeOfSubsample; i++)
    {
      m->setValue(i+1, j+1, points.getValue(i,j));
    }
  }

  selectedElements[numberOfElements - 1] = sizeOfSubsample - 1;
  for(int i = 0; i < sizeOfSubsample; i++)
    m->setValue(i+1, sizeOfSubsample, points.getValue(i, numberOfElements - 1));
}

void AllSubsets::reset(Matrix2D* m, Matrix2D* m2, Matrix2D& points, Vector& mu)
{
  reset();
  reset(m, points, mu);

  //create first subset matrix
  for(int j = 0; j < sizeOfSubsample; j++)
    m2->setValue(0, j, 1.0);

  for(int j = 0; j < sizeOfSubsample - 1; j++)
  {
    for(int i = 0; i < sizeOfSubsample; i++)
    {
      m2->setValue(i+1, j, points.getValue(i,j));
    }
  }

  for(int i = 0; i < sizeOfSubsample; i++)
    m2->setValue(i+1, sizeOfSubsample - 1, points.getValue(i, numberOfElements - 1));
}

void AllSubsets::reset()
{
  state = AllSubsets::_1a;
  finished = false;
  for(int i = 0; i < numberOfElements; i++)
  {
    p[i] = 0;
  }

  n = numberOfElements;
  h = sizeOfSubsample;

  /*
   * The algorithm doesn't work properly, if n = h or h = 0.
   * Therefore these cases have to be handled first.
   */
  if (sizeOfSubsample >= numberOfElements)
  {
		for (int i = 0; i < numberOfElements; i++)
			contains[i] = true;
  }
  else if (sizeOfSubsample <= 0)
  {
    for (int i = 0; i < numberOfElements; i++)
			contains[i] = false;
  }
  else
  {
    contains[n-1] = true;
    for (int i = h - 1; i < n - 1; i++)
    {
      contains[i] = false;
    }
    for (int i = 0; i < h - 1; i++)
    {
      contains[i] = true;
    }
  }

  if((n == h) || (h == 0))
    state = AllSubsets::_start;
  else if ((h > n) || (h < 0))
    finished = true;

  for(int i = 0; i < numberOfElements; i++)
    selectedElements[i] = -1;
}


bool* AllSubsets::getNextSubset(Matrix2D* m, Matrix2D& points)
{
  while(true)
  {
    switch (state)
    {
    case AllSubsets::_start:
      finished = true;
      return contains;
      break;
    case AllSubsets::_1a:
      n--;
      h--;
      p[n-1] = 0;

      if (h > 0)
      {
        state = AllSubsets::_2a;
      }
      else if (h == 0)
      {
        state = AllSubsets::_1b;
        return contains;
      }
      break;
    case AllSubsets::_1b:
      h++;
      contains[(n + 1) - 1] = false;
      if (h > 1)
      {
        contains[(h - 1) - 1] = true;
        for(int i = 0; i < sizeOfSubsample; i++)
        {
          m->setValue(i+1, selectedElements[(n+1)-1]+1, points.getValue(i, (h-1)-1));
        }
        selectedElements[(h-1)-1] = selectedElements[(n+1)-1];
        selectedElements[(n+1)-1] = -1;
      }
      else
      {
        contains[n - 1] = true;
        for(int i = 0; i < sizeOfSubsample; i++)
        {
          m->setValue(i+1, selectedElements[(n+1)-1]+1, points.getValue(i, (h-1)-1));
        }
        selectedElements[n-1] = selectedElements[(n+1)-1];
        selectedElements[(n+1)-1] = -1;
      }
      p[n - 1] = p[(n + 1) - 1] + 1;
      if (n > h)
      {
        state = AllSubsets::_1a;
      }
      else if (n == h)
      {
        state = AllSubsets::_1c;
      }
      break;
    case AllSubsets::_1c:
      n += p[n - 1];
      if (n == numberOfElements)
      {
			  /* we are done and have produced all subsamples! */
        state = AllSubsets::_start;
        finished = true;
      }
      else
      {
        state = AllSubsets::_2c;
      }
      return contains;
      break;
    case AllSubsets::_2a:
      p[h - 1] = p[n - 1] + n - h;
      n = h;
      state = _2b;
      return contains;
      break;
    case AllSubsets::_2b:
      contains[(n + 1) - 1] = true;
      if (h > 1)
      {
        contains[(h - 1) - 1] = false;
        for(int i = 0; i < sizeOfSubsample; i++)
        {
          m->setValue(i+1, selectedElements[(h - 1) - 1]+1, points.getValue(i, (n + 1) - 1));
        }
        selectedElements[(n + 1) - 1] = selectedElements[(h - 1) - 1];
        selectedElements[(h - 1) - 1] = -1;
      }
      else
      {
        contains[n - 1] = false;
        for(int i = 0; i < sizeOfSubsample; i++)
        {
          m->setValue(i+1, selectedElements[n - 1]+1, points.getValue(i, (n + 1) - 1));  
        }
        selectedElements[(n + 1) - 1] = selectedElements[n - 1];
        selectedElements[n - 1] = -1;
      }
      p[(n + 1) - 1] = p[n - 1] - 1;
      p[n - 1] = 0;
      h--;

      if (h > 0)
      {
        state = AllSubsets::_1a;
      }
      else if (h == 0)
      {
        state = AllSubsets::_2c;
        return contains;
      }
      break;
    case AllSubsets::_2c:
      n++;
      h++;
      if (contains[(n + 1) - 1])
      {
        state = AllSubsets::_1b;
      }
      else
      {
        state = AllSubsets::_2b;
      }
      break;
    default:
      break;
    } // switch (state)
  }
}

bool* AllSubsets::getNextSubset(Matrix2D* m, Matrix2D* m2, Matrix2D& points)
{

  while(true)
  {
    switch (state)
    {
    case AllSubsets::_start:
      finished = true;
      return contains;
      break;
    case AllSubsets::_1a:
      n--;
      h--;
      p[n-1] = 0;

      if (h > 0)
      {
        state = AllSubsets::_2a;
      }
      else if (h == 0)
      {
        state = AllSubsets::_1b;
        return contains;
      }
      break;
    case AllSubsets::_1b:
      h++;
      contains[(n + 1) - 1] = false;
      if (h > 1)
      {
        contains[(h - 1) - 1] = true;
        for(int i = 0; i < sizeOfSubsample; i++)
        {
          m->setValue(i+1, selectedElements[(n+1)-1]+1, points.getValue(i, (h-1)-1));
          m2->setValue(i+1, selectedElements[(n+1)-1], points.getValue(i, (h-1)-1));
        }
        selectedElements[(h-1)-1] = selectedElements[(n+1)-1];
        selectedElements[(n+1)-1] = -1;
      }
      else
      {
        contains[n - 1] = true;
        for(int i = 0; i < sizeOfSubsample; i++)
        {
          m->setValue(i+1, selectedElements[(n+1)-1]+1, points.getValue(i, (h-1)-1));
          m2->setValue(i+1, selectedElements[(n+1)-1], points.getValue(i, (h-1)-1));
        }
        selectedElements[n-1] = selectedElements[(n+1)-1];
        selectedElements[(n+1)-1] = -1;
      }
      p[n - 1] = p[(n + 1) - 1] + 1;
      if (n > h)
      {
        state = AllSubsets::_1a;
      }
      else if (n == h)
      {
        state = AllSubsets::_1c;
      }
      break;
    case AllSubsets::_1c:
      n += p[n - 1];
      if (n == numberOfElements)
      {
			  /* we are done and have produced all subsamples! */
        state = AllSubsets::_start;
        finished = true;
      }
      else
      {
        state = AllSubsets::_2c;
      }
      return contains;
      break;
    case AllSubsets::_2a:
      p[h - 1] = p[n - 1] + n - h;
      n = h;
      state = _2b;
      return contains;
      break;
    case AllSubsets::_2b:
      contains[(n + 1) - 1] = true;
      if (h > 1)
      {
        contains[(h - 1) - 1] = false;
        for(int i = 0; i < sizeOfSubsample; i++)
        {
          m->setValue(i+1, selectedElements[(h - 1) - 1]+1, points.getValue(i, (n + 1) - 1));
          m2->setValue(i+1, selectedElements[(h - 1) - 1], points.getValue(i, (n + 1) - 1));
        }
        selectedElements[(n + 1) - 1] = selectedElements[(h - 1) - 1];
        selectedElements[(h - 1) - 1] = -1;
      }
      else
      {
        contains[n - 1] = false;
        for(int i = 0; i < sizeOfSubsample; i++)
        {
          m->setValue(i+1, selectedElements[n - 1]+1, points.getValue(i, (n + 1) - 1));
          m2->setValue(i+1, selectedElements[n - 1], points.getValue(i, (n + 1) - 1));  
        }
        selectedElements[(n + 1) - 1] = selectedElements[n - 1];
        selectedElements[n - 1] = -1;
      }
      p[(n + 1) - 1] = p[n - 1] - 1;
      p[n - 1] = 0;
      h--;

      if (h > 0)
      {
        state = AllSubsets::_1a;
      }
      else if (h == 0)
      {
        state = AllSubsets::_2c;
        return contains;
      }
      break;
    case AllSubsets::_2c:
      n++;
      h++;
      if (contains[(n + 1) - 1])
      {
        state = AllSubsets::_1b;
      }
      else
      {
        state = AllSubsets::_2b;
      }
      break;
    default:
      break;
    } // switch (state)
  }
}


void RandomSubsets::reset(Matrix2D* m, Matrix2D& points, Vector& mu)
{
  reset();

  //create first subset matrix
  for(int j = 0; j <= sizeOfSubsample; j++)
    m->setValue(0, j, 1.0);

  for(int i = 1; i <= sizeOfSubsample; i++)
  {
    m->setValue(i, 0, mu.getValue(i-1));
  }
 
  getNextSubset(m,points);
  count=0;
}

void RandomSubsets::reset(Matrix2D* m, Matrix2D* m2, Matrix2D& points, Vector& mu)
{
  reset(m, points, mu);
  reset();

  //create first subset matrix
  for(int j = 0; j < sizeOfSubsample; j++)
    m2->setValue(0, j, 1.0);

  for(int j = 0; j < sizeOfSubsample - 1; j++)
  {
    for(int i = 0; i < sizeOfSubsample; i++)
    {
      m2->setValue(i+1, j, points.getValue(i,j));
    }
  }

  getNextSubset(m,m2,points);
  count=0;
}

void RandomSubsets::reset()
{
  count=0;
  delete r;
  r=new MTRand(seed);
}

bool* RandomSubsets::getNextSubset(Matrix2D* m, Matrix2D& points)
{
	getRandomSubset(contains);

	int j = 0;
	for(int k = 0; k < numberOfElements; k++)
	{
		if(contains[k]) // besser auf selectedElements umsteigen!
		{
			for(int i = 0; i < sizeOfSubsample; i++)
				m->setValue(i+1, j+1, points.getValue(i,k)); // thorsten j->k
			j++;
			if(j == sizeOfSubsample)
				break;
		}

	}
	return contains; //???
}

bool* RandomSubsets::getNextSubset(Matrix2D* m, Matrix2D* m2, Matrix2D& points)
{
	getRandomSubset(contains);

	int j = 0;
	for(int k = 0; k < numberOfElements; k++)
	{
		if(contains[k]) // besser auf selectedElements umsteigen!
		{
			for(int i = 0; i < sizeOfSubsample   ; i++)
			{
				m->setValue(i+1, j+1, points.getValue(i,k));  //thorsten  j -> k
				m2->setValue(i+1, j, points.getValue(i,k));  //thorsten j -> k
			}
			j++;
			if(j == sizeOfSubsample)
				break;
		}
	}
	return contains; //???
}

void RandomSubsets::getRandomSubset(bool* subset)
{
	for(int i = 0; i < numberOfElements; i++)
		subset[i] = false;

	int nr;
	for(int i = 0; i < sizeOfSubsample; i++) 
	{
		do
		{
			nr = (int)r->rand(numberOfElements);
		}while(nr == numberOfElements || subset[nr]);
		subset[nr] = true;
	}
	count++;
}

RandomSubsets::~RandomSubsets()
{
delete r;
delete[] contains;
}
