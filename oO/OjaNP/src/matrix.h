/* $Id: matrix.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <iostream>
#include <stdexcept>

using namespace std;

// THIS FILE IS TAKEN FROM TCL-LITE WITHOUT PERMISSION !!!!
// ========================================================

namespace Matrix
{
	
class matrix_error : public logic_error
{
    public:
        matrix_error (const string& what_arg) : logic_error( what_arg) {}
};

template <class T> class matrix
{
public:
    
   matrix (const matrix<T> & m);
   matrix (size_t row = 6, size_t col = 6);

    
   ~matrix ();

    
   matrix<T> & operator = (const matrix<T> & m) throw () ;

    
   size_t RowNo () const { return _m->Row; }
   size_t ColNo () const { return _m->Col; }

    
   T& operator () (size_t row, size_t col) ;//throw (matrix_error) ;

    
   matrix<T>  operator + () throw ()  { return *this; }
   matrix<T>  operator - () throw () ;

    
   matrix<T> & operator += (const matrix<T> & m) ;//throw (matrix_error) ;
   matrix<T> & operator -= (const matrix<T> & m) ;//throw (matrix_error) ;
   matrix<T> & operator *= (const matrix<T> & m) ;//throw (matrix_error) ;
   matrix<T> & operator *= (const T& c) throw () ;
   matrix<T> & operator /= (const T& c) throw () ;
   matrix<T> & operator ^= (const size_t& pow) ;//throw (matrix_error) ;

    
   void Null (const size_t& row, const size_t& col) throw () ;
   void Null () throw () ;
   void Unit (const size_t& row) throw () ;
   void Unit () throw () ;
   void SetSize (size_t row, size_t col) throw () ;

    
   matrix<T>  Solve (const matrix<T> & v) const ;//throw (matrix_error) ;
   matrix<T>  Adj () const ;//throw (matrix_error) ;
   T Det () const ;//throw (matrix_error) ;
   T Det_destroyable() ;//throw (matrix_error);
   T Norm () throw () ;
   T Cofact (size_t row, size_t col) const ;//throw (matrix_error) ;
   T Cond () throw () ;

    
   bool IsSquare () throw ()  { return (_m->Row == _m->Col); } 
   bool IsSingular () throw () ;
   bool IsDiagonal () throw () ;
   bool IsScalar () throw () ;
   bool IsUnit () throw () ;
   bool IsNull () throw () ;
   bool IsSymmetric () throw () ;
   bool IsSkewSymmetric () throw () ;
   bool IsUpperTriangular () throw () ;
   bool IsLowerTriangular () throw () ;

    struct base_mat
    {
        T **Val;
        size_t Row, Col, RowSiz, ColSiz;
        int Refcnt;

        base_mat (size_t row, size_t col, T** v)
        {
            Row = row; RowSiz = row;
            Col = col; ColSiz = col;
            Refcnt = 1;

            Val = new T* [row];
            size_t rowlen = col * sizeof(T);

            for (size_t i=0; i < row; i++)
            { 
	        Val[i] = new T [col];
                if (v) memcpy( Val[i], v[i], rowlen);
            }
        }
        ~base_mat ()
        {
            for (size_t i=0; i < RowSiz; i++)
                delete [] Val[i];
            delete [] Val;
        }
    };
    base_mat *_m;

    void clone ();
    void realloc (size_t row, size_t col);
    int pivot (size_t row);
};


template <class T> 
matrix<T> ::matrix (size_t row, size_t col)
{
  _m = new base_mat( row, col, 0);
}

 
template <class T> 
matrix<T> ::matrix (const matrix<T> & m)
{
    _m = m._m;
    _m->Refcnt++;
}

 
template <class T>  inline void
matrix<T> ::clone ()
{
    _m->Refcnt--;
    _m = new base_mat( _m->Row, _m->Col, _m->Val);
}

 
template <class T>  inline
matrix<T> ::~matrix ()
{
   if (--_m->Refcnt == 0) delete _m;
}

 
template <class T>  inline matrix<T> &
matrix<T> ::operator = (const matrix<T> & m) throw () 
{
    m._m->Refcnt++;
    if (--_m->Refcnt == 0) delete _m;
    _m = m._m;
    return *this;
}

 
template <class T>  inline void 
matrix<T> ::realloc (size_t row, size_t col)
{
   if (row == _m->RowSiz && col == _m->ColSiz)
   {
      _m->Row = _m->RowSiz;
      _m->Col = _m->ColSiz;
      return;
   }

   base_mat *m1 = new base_mat( row, col, ((void *)0) );
   size_t colSize = ((( _m->Col ) < ( col )) ? ( _m->Col ) : ( col ))  * sizeof(T);
   size_t minRow = ((( _m->Row ) < ( row )) ? ( _m->Row ) : ( row )) ;

   for (size_t i=0; i < minRow; i++)
      memcpy( m1->Val[i], _m->Val[i], colSize);

   if (--_m->Refcnt == 0) 
       delete _m;
   _m = m1;

   return;
}

 
template <class T>  inline void
matrix<T> ::SetSize (size_t row, size_t col) throw () 
{
   size_t i,j;
   size_t oldRow = _m->Row;
   size_t oldCol = _m->Col;

   if (row != _m->RowSiz || col != _m->ColSiz)
      realloc( row, col);

   for (i=oldRow; i < row; i++)
      for (j=0; j < col; j++)
         _m->Val[i][j] = T(0);

   for (i=0; i < row; i++)			
      for (j=oldCol; j < col; j++)
         _m->Val[i][j] = T(0);

   return;
}

 
template <class T>  inline T&
matrix<T> ::operator () (size_t row, size_t col) //throw (matrix_error) 
{
   if (row >= _m->Row || col >= _m->Col)
      throw matrix_error(   "matrixT::operator(): Index out of range!" ); ;
   if (_m->Refcnt > 1) clone();
   return _m->Val[row][col];
}

 
template <class T>  inline istream&
operator >> (istream& istrm, matrix<T> & m)
{
    if (m._m->Refcnt > 1)
    {
        m._m->Refcnt--;
        m._m = new typename matrix<T>::base_mat( m._m->Row, m._m->Col, ((void *)0) );
    }
    for (size_t i=0; i < m._m->Row; i++)
       for (size_t j=0; j < m._m->Col; j++)
          istrm >> m._m->Val[i][j];
    return istrm;
}

 
template <class T>  inline ostream&
operator << (ostream &ostrm, const matrix<T> & m)
{
   for (size_t i=0; i < m._m->Row; i++)
   {
       for (size_t j=0; j < m._m->Col; j++)
          ostrm << m._m->Val[i][j] << '\t';
       ostrm << endl;
   }
   return ostrm;
}


 
template <class T>  inline bool
operator == (const matrix<T> & m1, const matrix<T> & m2) throw () 
{
   if (m1.RowNo() != m2.RowNo() || m1.ColNo() != m2.ColNo())
      return false;

   for (size_t i=0; i < m1.RowNo(); i++)
      for (size_t j=0; j < m1.ColNo(); j++)
         if (m1._m->Val[i][j] != m2._m->Val[i][j])
            return false;

   return true;
}

 
template <class T>  inline bool
operator != (const matrix<T> & m1, const matrix<T> & m2) throw () 
{
    return (m1 == m2) ? false : true;
}

 
template <class T>  inline matrix<T> &
matrix<T> ::operator += (const matrix<T> & m) //throw (matrix_error) 
{
   if (_m->Row != m._m->Row || _m->Col != m._m->Col)
      throw matrix_error(   "matrixT::operator+= : Inconsistent matrix sizes in addition!" ); ;
   if (_m->Refcnt > 1) clone();
   for (size_t i=0; i < m._m->Row; i++)
      for (size_t j=0; j < m._m->Col; j++)
         _m->Val[i][j] += m._m->Val[i][j];
   return *this;
}

 
template <class T>  inline matrix<T> &
matrix<T> ::operator -= (const matrix<T> & m) //throw (matrix_error) 
{
   if (_m->Row != m._m->Row || _m->Col != m._m->Col)
      throw matrix_error(   "matrixT::operator-= : Inconsistent matrix sizes in subtraction!" ); ;
   if (_m->Refcnt > 1) clone();
   for (size_t i=0; i < m._m->Row; i++)
      for (size_t j=0; j < m._m->Col; j++)
         _m->Val[i][j] -= m._m->Val[i][j];
   return *this;
}

 
template <class T>  inline matrix<T> &
matrix<T> ::operator *= (const T& c) throw () 
{
    if (_m->Refcnt > 1) clone();
    for (size_t i=0; i < _m->Row; i++)
        for (size_t j=0; j < _m->Col; j++)
            _m->Val[i][j] *= c;
    return *this;
}

 
template <class T>  inline matrix<T> &
matrix<T> ::operator *= (const matrix<T> & m) //throw (matrix_error) 
{
   if (_m->Col != m._m->Row)
      throw matrix_error(   "matrixT::operator*= : Inconsistent matrix sizes in multiplication!" ); ;

   *this = *this * m;

   return *this;
}

 
template <class T>  inline matrix<T> &
matrix<T> ::operator /= (const T& c) throw () 
{
    if (_m->Refcnt > 1) clone();
    for (size_t i=0; i < _m->Row; i++)
        for (size_t j=0; j < _m->Col; j++)
            _m->Val[i][j] /= c;

    return *this;
}

 
template <class T>  inline matrix<T> &
matrix<T> ::operator ^= (const size_t& pow) //throw (matrix_error) 
{
	matrix<T>  temp(*this);

	for (size_t i=2; i <= pow; i++)
      *this = *this * temp;

	return *this;
}

 
template <class T>  inline matrix<T> 
matrix<T> ::operator - () throw () 
{
   matrix<T>  temp(_m->Row,_m->Col);

   for (size_t i=0; i < _m->Row; i++)
      for (size_t j=0; j < _m->Col; j++)
         temp._m->Val[i][j] = - _m->Val[i][j];

   return temp;
}

 
template <class T>  inline matrix<T> 
operator + (const matrix<T> & m1, const matrix<T> & m2) //throw (matrix_error) 
{
   if (m1._m->Row != m2._m->Row || m1._m->Col != m2._m->Col)
      throw matrix_error(   "matrixT::operator+: Inconsistent matrix size in addition!" ); ;

   matrix<T>  temp(m1._m->Row,m1._m->Col);

   for (size_t i=0; i < m1._m->Row; i++)
      for (size_t j=0; j < m1._m->Col; j++)
         temp._m->Val[i][j] = m1._m->Val[i][j] + m2._m->Val[i][j];

   return temp;
}

 
template <class T>  inline matrix<T> 
operator - (const matrix<T> & m1, const matrix<T> & m2) //throw (matrix_error) 
{
   if (m1._m->Row != m2._m->Row || m1._m->Col != m2._m->Col)
     throw matrix_error(   "matrixT::operator-: Inconsistent matrix size in subtraction!" ); ;

   matrix<T>  temp(m1._m->Row,m1._m->Col);

   for (size_t i=0; i < m1._m->Row; i++)
      for (size_t j=0; j < m1._m->Col; j++)
         temp._m->Val[i][j] = m1._m->Val[i][j] - m2._m->Val[i][j];

   return temp;
}

 
template <class T>  inline matrix<T> 
operator * (const matrix<T> & m, const T& no) throw () 
{
   matrix<T>  temp(m._m->Row,m._m->Col);

   for (size_t i=0; i < m._m->Row; i++)
      for (size_t j=0; j < m._m->Col; j++)
         temp._m->Val[i][j] = no * m._m->Val[i][j];

   return temp;
}


 
template <class T>  inline matrix<T> 
operator * (const T& no, const matrix<T> & m) throw () 
{
   return (m * no);
}

 
template <class T>  inline matrix<T> 
operator * (const matrix<T> & m1, const matrix<T> & m2) //throw (matrix_error) 
{
   if (m1._m->Col != m2._m->Row)
      throw matrix_error(   "matrixT::operator*: Inconsistent matrix size in multiplication!" ); ;

   matrix<T>  temp(m1._m->Row,m2._m->Col);

   for (size_t i=0; i < m1._m->Row; i++)
      for (size_t j=0; j < m2._m->Col; j++)
      {
         temp._m->Val[i][j] = T(0);
         for (size_t k=0; k < m1._m->Col; k++)
            temp._m->Val[i][j] += m1._m->Val[i][k] * m2._m->Val[k][j];
      }
   return temp;
}

 
template <class T>  inline matrix<T> 
operator / (const matrix<T> & m, const T& no) throw () 
{
    return (m * (1.0 / no));
}


 
template <class T>  inline matrix<T> 
operator / (const T& no, const matrix<T> & m) //throw (matrix_error) 
{
    return (!m * no);
}

 
template <class T>  inline matrix<T> 
operator / (const matrix<T> & m1, const matrix<T> & m2) //throw (matrix_error) 
{
    return (m1 * !m2);
}

 
template <class T>  inline matrix<T> 
operator ^ (const matrix<T> & m, const size_t& pow) //throw (matrix_error) 
{
   matrix<T>  temp(m);

   for (size_t i=2; i <= pow; i++)
      temp = temp * m;

   return temp;
}

 
template <class T>  inline matrix<T> 
operator  ~ (const matrix<T> & m) throw () 
{
   matrix<T>  temp(m._m->Col,m._m->Row);

   for (size_t i=0; i < m._m->Row; i++)
      for (size_t j=0; j < m._m->Col; j++)
	    temp._m->Val[j][i] = m._m->Val[i][j];

   return temp;
}

 
template <class T>  inline matrix<T> 
operator ! (matrix<T>  m) //throw (matrix_error) 
{
   size_t i,j,k;
   T a1,a2,*rowptr;

   if (m._m->Row != m._m->Col)
      throw matrix_error(   "matrixT::operator!: Inversion of a non-square matrix" ); ;

   matrix<T>  temp(m._m->Row,m._m->Col);
   if (m._m->Refcnt > 1) m.clone();


   temp.Unit();
   for (k=0; k < m._m->Row; k++)
   {
      int indx = m.pivot(k);
      if (indx == -1)
         throw matrix_error(   "matrixT::operator!: Inversion of a singular matrix" ); ;

      if (indx != 0)
      {
         rowptr = temp._m->Val[k];
         temp._m->Val[k] = temp._m->Val[indx];
         temp._m->Val[indx] = rowptr;
      }
      a1 = m._m->Val[k][k];
      for (j=0; j < m._m->Row; j++)
      {
         m._m->Val[k][j] /= a1;
         temp._m->Val[k][j] /= a1;
      }
      for (i=0; i < m._m->Row; i++)
         if (i != k)
         {
            a2 = m._m->Val[i][k];
            for (j=0; j < m._m->Row; j++)
            {
               m._m->Val[i][j] -= a2 * m._m->Val[k][j];
               temp._m->Val[i][j] -= a2 * temp._m->Val[k][j];
            }
         }
   }
   return temp;
}

 
template <class T>  inline matrix<T> 
matrix<T> ::Solve (const matrix<T> & v) const //throw (matrix_error) 
{
   size_t i,j,k;
   T a1;

   if (!(_m->Row == _m->Col && _m->Col == v._m->Row))
      throw matrix_error(   "matrixT::Solve():Inconsistent matrices!" ); ;

   matrix<T>  temp(_m->Row,_m->Col+v._m->Col);
   for (i=0; i < _m->Row; i++)
   {
      for (j=0; j < _m->Col; j++)
         temp._m->Val[i][j] = _m->Val[i][j];
      for (k=0; k < v._m->Col; k++)
         temp._m->Val[i][_m->Col+k] = v._m->Val[i][k];
   }
   for (k=0; k < _m->Row; k++)
   {
      int indx = temp.pivot(k);
      if (indx == -1)
         throw matrix_error(   "matrixT::Solve(): Singular matrix!" ); ;

      a1 = temp._m->Val[k][k];
      for (j=k; j < temp._m->Col; j++)
         temp._m->Val[k][j] /= a1;

      for (i=k+1; i < _m->Row; i++)
      {
         a1 = temp._m->Val[i][k];
         for (j=k; j < temp._m->Col; j++)
	   temp._m->Val[i][j] -= a1 * temp._m->Val[k][j];
      }
   }
   matrix<T>  s(v._m->Row,v._m->Col);
   for (k=0; k < v._m->Col; k++)
      for (int m=int(_m->Row)-1; m >= 0; m--)
      {
         s._m->Val[m][k] = temp._m->Val[m][_m->Col+k];
         for (j=m+1; j < _m->Col; j++)
            s._m->Val[m][k] -= temp._m->Val[m][j] * s._m->Val[j][k];
      }
   return s;
}

 
template <class T>  inline void
matrix<T> ::Null (const size_t& row, const size_t& col) throw () 
{
    if (row != _m->Row || col != _m->Col)
        realloc( row,col);

    if (_m->Refcnt > 1) 
        clone();

    for (size_t i=0; i < _m->Row; i++)
        for (size_t j=0; j < _m->Col; j++)
            _m->Val[i][j] = T(0);
    return;
}

 
template <class T>  inline void
matrix<T> ::Null() throw () 
{
    if (_m->Refcnt > 1) clone();   
    for (size_t i=0; i < _m->Row; i++)
        for (size_t j=0; j < _m->Col; j++)
 	        _m->Val[i][j] = T(0);
    return;
}

 
template <class T>  inline void
matrix<T> ::Unit (const size_t& row) throw () 
{
    if (row != _m->Row || row != _m->Col)
        realloc( row, row);
        
    if (_m->Refcnt > 1) 
        clone();

    for (size_t i=0; i < _m->Row; i++)
        for (size_t j=0; j < _m->Col; j++)
            _m->Val[i][j] = i == j ? T(1) : T(0);
    return;
}

 
template <class T>  inline void
matrix<T> ::Unit () throw () 
{
    if (_m->Refcnt > 1) clone();   
    size_t row = ((( _m->Row ) < ( _m->Col )) ? ( _m->Row ) : ( _m->Col )) ;
    _m->Row = _m->Col = row;

    for (size_t i=0; i < _m->Row; i++)
        for (size_t j=0; j < _m->Col; j++)
            _m->Val[i][j] = i == j ? T(1) : T(0);
    return;
}

 
template <class T>  inline int
matrix<T> ::pivot (size_t row)
{
  int k = int(row);
  double amax,temp;

  amax = -1;
  for (size_t i=row; i < _m->Row; i++)
    if ( (temp = abs( _m->Val[i][row])) > amax)
     {
       if(temp!=0.0){
         amax = temp;
         k = i;
       }
     }
  if (_m->Val[k][row] == T(0))
     return -1;
  if (k != int(row))
  {
     T* rowptr = _m->Val[k];
     _m->Val[k] = _m->Val[row];
     _m->Val[row] = rowptr;
     return k;
  }
  return 0;
}

 
template <class T>  inline T
matrix<T> ::Det () const //throw (matrix_error) 
{
   size_t i,j,k;
   T piv,detVal = T(1);

   if (_m->Row != _m->Col)
      throw matrix_error(   "matrixT::Det(): Determinant a non-square matrix!" ); ;
   
   matrix<T>  temp(*this);
   if (temp._m->Refcnt > 1) temp.clone();

   for (k=0; k < _m->Row; k++)
   {
      int indx = temp.pivot(k);
      if (indx == -1)
         return 0;
      if (indx != 0)
         detVal = - detVal;
      detVal = detVal * temp._m->Val[k][k];
      for (i=k+1; i < _m->Row; i++)
      {
         piv = temp._m->Val[i][k] / temp._m->Val[k][k];
         for (j=k+1; j < _m->Row; j++)
            temp._m->Val[i][j] -= piv * temp._m->Val[k][j];
      }
   }
   return detVal;
}

/**
* Calculates det without making a temp matrix, destroying the source one. To be used if the given matrix is already a copy
*/
template <class T>  inline T
matrix<T> ::Det_destroyable() //throw (matrix_error)
{
	size_t i, j, k;
	T piv, detVal = T(1);

	if (_m->Row != _m->Col)
		throw matrix_error("matrixT::Det(): Determinant a non-square matrix!");

	for (k = 0; k < _m->Row; k++)
	{
		int indx = pivot(k);
		if (indx == -1)
			return 0;
		if (indx != 0)
			detVal = -detVal;
		detVal = detVal * _m->Val[k][k];
		for (i = k + 1; i < _m->Row; i++)
		{
			piv = _m->Val[i][k] / _m->Val[k][k];
			for (j = k + 1; j < _m->Row; j++)
				_m->Val[i][j] -= piv * _m->Val[k][j];
		}
	}
	return detVal;
}
 
template <class T>  inline T
matrix<T> ::Norm () throw () 
{
   T retVal = T(0);

   for (size_t i=0; i < _m->Row; i++)
      for (size_t j=0; j < _m->Col; j++)
         retVal += _m->Val[i][j] * _m->Val[i][j];
   retVal = sqrt( retVal);

   return retVal;
}

 
template <class T>  inline T
matrix<T> ::Cond () throw () 
{
   matrix<T>  inv = ! (*this);
   return (Norm() * inv.Norm());
}

 
template <class T>  inline T
matrix<T> ::Cofact (size_t row, size_t col) const //throw (matrix_error) 
{
   size_t i,i1,j,j1;

   if (_m->Row != _m->Col)
      throw matrix_error(   "matrixT::Cofact(): Cofactor of a non-square matrix!" ); ;

   if (row > _m->Row || col > _m->Col)
      throw matrix_error(   "matrixT::Cofact(): Index out of range!" ); ;

   matrix<T>  temp (_m->Row-1,_m->Col-1);

   for (i=i1=0; i < _m->Row; i++)
   {
      if (i == row)
        continue;
      for (j=j1=0; j < _m->Col; j++)
      {
      	 if (j == col)
            continue;
      	 temp._m->Val[i1][j1] = _m->Val[i][j];
         j1++;
      }
      i1++;
   }
   T  cof = temp.Det_destroyable();
   if ((row+col)%2 == 1)
      cof = -cof;

   return cof;
}
 
template <class T>  inline matrix<T> 
matrix<T> ::Adj () const //throw (matrix_error) 
{
   if (_m->Row != _m->Col)
      throw matrix_error(   "matrixT::Adj(): Adjoin of a non-square matrix." ); ;

   matrix<T>  temp(_m->Row,_m->Col);

   for (size_t i=0; i < _m->Row; i++)
      for (size_t j=0; j < _m->Col; j++)
          temp._m->Val[j][i] = Cofact(i,j);
   return temp;
}

 
template <class T>  inline bool
matrix<T> ::IsSingular () throw () 
{
   if (_m->Row != _m->Col)
      return false;
   return (Det() == T(0));
}

 
template <class T>  inline bool
matrix<T> ::IsDiagonal () throw () 
{
   if (_m->Row != _m->Col)
      return false;
   for (size_t i=0; i < _m->Row; i++)
     for (size_t j=0; j < _m->Col; j++)
        if (i != j && _m->Val[i][j] != T(0))
          return false;
   return true;
}

 
template <class T>  inline bool
matrix<T> ::IsScalar () throw () 
{
   if (!IsDiagonal())
     return false;
   T v = _m->Val[0][0];
   for (size_t i=1; i < _m->Row; i++)
     if (_m->Val[i][i] != v)
        return false;
   return true;
}

 
template <class T>  inline bool
matrix<T> ::IsUnit () throw () 
{
   if (IsScalar() && _m->Val[0][0] == T(1))
     return true;
   return false;
}

 
template <class T>  inline bool
matrix<T> ::IsNull () throw () 
{
   for (size_t i=0; i < _m->Row; i++)
      for (size_t j=0; j < _m->Col; j++)
         if (_m->Val[i][j] != T(0))
            return false;
   return true;
}

 
template <class T>  inline bool
matrix<T> ::IsSymmetric () throw () 
{
   if (_m->Row != _m->Col)
      return false;
   for (size_t i=0; i < _m->Row; i++)
      for (size_t j=0; j < _m->Col; j++)
         if (_m->Val[i][j] != _m->Val[j][i])
            return false;
   return true;
}
           
 
template <class T>  inline bool
matrix<T> ::IsSkewSymmetric () throw () 
{
   if (_m->Row != _m->Col)
      return false;
   for (size_t i=0; i < _m->Row; i++)
      for (size_t j=0; j < _m->Col; j++)
         if (_m->Val[i][j] != -_m->Val[j][i])
            return false;
   return true;
}
   
 
template <class T>  inline bool
matrix<T> ::IsUpperTriangular () throw () 
{
   if (_m->Row != _m->Col)
      return false;
   for (size_t i=1; i < _m->Row; i++)
      for (size_t j=0; j < i-1; j++)
         if (_m->Val[i][j] != T(0))
            return false;
   return true;
}

 
template <class T>  inline bool
matrix<T> ::IsLowerTriangular () throw () 
{
   if (_m->Row != _m->Col)
      return false;

   for (size_t j=1; j < _m->Col; j++)
      for (size_t i=0; i < j-1; i++)
         if (_m->Val[i][j] != T(0))
            return false;

   return true;
}
}
