/**
 * @authors Hmamouche Youssef
 **/

#include "../inst/include/struct.h"
using namespace std;

namespace Struct
{
/***********************************************************************************/
double CVDouble::Mean () const //throw (Exception)
  {
    if (this->size () == 0)
         throw Exception ("Vector of size null");
    double Sum = 0;
    for (auto it = this->begin() ; it!= this->end() ; ++it)
        Sum += (*it);
    return Sum / this->size();
  }
/***********************************************************************************/
bool CVDouble::Contains (unsigned & x)
{
    if (this->size () == 0) return false;
    bool test = false;
    for (auto Val = this->begin() ; Val != this->end() ; ++Val)
        if (*Val == x)
        {
            test = true;
            break;
        }
    return test;
}
/***********************************************************************************/
double CVDouble::CMean () const //throw (Exception)
{
    if (this->size () == 0)
         throw Exception ("Vector of size null");

    double m = 0.0;
    int c = 0;
    for (auto it = this->begin() ; it!= this->end() ; ++it)
    {
        if ( ! std::isnan (*it) )
        {
            m += (*it);
            ++c;
        }
    }
    return m / c;
}

/***********************************************************************************/

double CVDouble::StdDev () const //throw (Exception)
{
    if (0 == this->size ())
               throw Exception ("Vector of size null");
    double Sum = 0.0, Mean = this->Mean();
    for( auto it = this->begin() ; it != this->end() ; ++it)
    {
        Sum += pow (*it - Mean, 2);
    }
    return sqrt(Sum / (this->size ()));
}

/***********************************************************************************/

double CVDouble::Min () const //throw (Exception)
{
    if (0 == this->size ())
               throw Exception ("Vector of size null");
    double min = *(this->begin());
     for (auto Val = this->begin() ; Val != this->end() ; ++Val)
         if ( *Val < min)
             min = *Val;
     return min;
}

/***********************************************************************************/

double CVDouble::Max () const //throw (Exception)
{
    if (0 == this->size ())
               throw Exception ("Vector of size null");
    double max = *(this->begin());
     for (auto Val = this->begin() ; Val != this->end() ; ++Val)
         if (*Val > max)
             max = *Val;
     return max;
}

/***********************************************************************************/

void   CVDouble::Standardise () //throw (Exception)
{
    if (0 == this->size ())
               throw Exception ("Vector of size null");
    double Mean = this->Mean (), sigma = this->StdDev();
    if (Mean != 0 && sigma > 1E-9)
      for (auto it = this->begin () ; it != this->end () ; ++it)
         (*it) = ((*it) - Mean);

}

/***********************************************************************************/
void   CVDouble::Normalise () //throw (Exception)
{
    if (0 == this->size ())
               throw Exception ("Vector of size null");
    double Min = this->Min ();
    double Max = this->Max ();
    if (abs (Max) > abs (Min))
      for (auto it = this->begin () ; it != this->end () ; ++it)
            (*it) = (*it) / Max;
    if (abs (Max) < abs (Min))
      for (auto it = this->begin () ; it != this->end () ; ++it)
            (*it) = (*it) / Min;


}

/***********************************************************************************/

void CVDouble::Add  (double & m)
{
        for (auto it = this->begin () ; it != this->end () ; ++it)
            (*it) += m;
}

/***********************************************************************************/

void CMatDouble::Init_Mat( const vector< vector <double> > & M)
  {
      unsigned int i,j;
      this->clear ();
      unsigned long Nblign = M.size ();
      unsigned long NbColumn = M[0].size ();

      this->resize (NbColumn);
      for( j = 0 ; j < NbColumn ; ++j)
        {
          (*this) [j] = CVDouble (Nblign);
          for( i = 0 ; i < Nblign  ; ++i)
              (*this) [j][i] = M[i][j];
      }
}
/****************************************************/
vector<vector<double> > CMatDouble::to_Mat()
{
    unsigned long m = (*this).size(), n = (*this)[0].size ();
    vector<vector<double>>  T (n);
    for (unsigned long i = 0 ; i < n ; ++i)
    {
        T[i]. resize (m);
         for (unsigned long j = 0 ; j < m ; ++j)
             T[i][j] = (*this)[j][i];
    }
    return T;
}
/***********************************************************************************/

// Check if a vector contains a nan value
bool CVDouble::NBR_NAN() const
{
    bool test = 1;
    for (const auto & val : (*this) )
        if (std::isnan (val) )
        {
            test = 0;
            break;
        }
     return test;
}

/***********************************************************************************/

// Interpolation x(t) = x(t-1)
void CMatDouble::Interpol () //throw (Exception)
{
   if (this->size () == 0)
       throw Exception ("Vector of size null");
   double moy ;
   CVDouble T ;
   for (auto it = this->begin () ; it!= this->end () ; ++it)
   {
       if(it->NBR_NAN() == 0)
       {
           moy = it->CMean ();

           if (std::isnan(*(it->begin ())))
               *(it->begin ()) = moy;

           for (auto val = it->begin () + 1 ; val!= it->end () ; ++val)
               if (std::isnan(*val))
                   *val = *(val - 1);
           }
       }
   }
/***********************************************************************************/

void CMatDouble::Standardise ()
{
    for( auto it = this->begin() ; it!= this->end() ; ++it)
        it->Standardise ();
}

/***********************************************************************************/
CMatDouble CMatDouble::Normalise ()
{
    if (this->size () == 0)
        throw Exception ("Matrix of size null");

    double min, max;
    CMatDouble minMax (2);

    for (auto Vect = this->begin () ; Vect != this->end () ; ++Vect)
    {
        min = Vect->Min ();
        max = Vect->Max ();
        minMax[0].push_back  (min);
        minMax[1].push_back  (max);

        for (auto it = Vect->begin () ; it != Vect->end () ; ++it)
            (*it) = ((*it) - min) / (max - min);
    }

    return (minMax);
}

/***********************************************************************************/
void CMatDouble::Denormalising (const CMatDouble & minMax)
{
    int i = 0;
    for (auto Vect = this->begin () ; Vect != this->end () ; ++Vect)
    {
        for (auto it = Vect->begin () ; it != Vect->end () ; ++it)
        {
            (*it) = (*it) * (minMax[1][i] - minMax[0][i]) + minMax[0][i];
        }
        ++i;
    }
}
/***********************************************************************************/

CMatDouble Trans (const CMatDouble & M)
{
    unsigned int n = M.size(), m = M[0].size ();
    CMatDouble T (m);
    for (unsigned int i = 0 ; i < m ; ++i)
    {
        T[i] = CVDouble (n);
         for (unsigned int j = 0 ; j < n ; ++j)
             T[i][j] = M[j][i];
    }
        return T;
}

/***********************************************************************************/

void permute( CVDouble &X , CVDouble &Y)
{
    double a;
    for ( unsigned int k = 0 ; k < X.size() ; ++k)
    {
        a = X[k];
        X[k] = Y[k];
        Y[k] = a;
    }

}

/***********************************************************************************/

bool Trig (CMatDouble & A ,  CMatDouble & B)
{
    double pivot = 1, coef=0, max ;
    unsigned int i,j,k,l, indice_max;
    unsigned int n = A.size();
    vector<double> T, P;

// Pivot de Gauss
    // Reorganisation de la matrice
    for ( k = 0 ; k < n-1 ; ++k)
    {
        T.resize(0);
        P.resize(0);
        indice_max = k, max = A[k][k];
        for( l = k+1 ; l < n ; l++)
        {
            if (abs(A[l][k]) > abs(max) )
            {
                max = A[l][k];
                indice_max = l;
            }
        }
        if(indice_max > k)
            {
            permute(B[k] , B[indice_max]);
            permute(A[k] , A[indice_max]);
            }
     // Affectation du pivot
      pivot = A[k][k];

        if( pivot == 0 )
          {
             return 0;
          }

// Triangularisation de la matrice A
        for ( i = k+1 ; i < n ; ++i)
           {
             coef = A[i][k];
              for ( j = 0 ; j < n ; ++j)
                 {
                  A[i][j] -= (coef * A[k][j] / pivot);
                  B[i][j] -= (B[k][j] * coef / pivot);
                  }
           }

        } // Matrice Triangulaire

        if( A[n-1][n-1] == 0 )
          {
            return 0;
          }

    B = Trans(B);
    return 1;
}

/***********************************************************************************/

double Det ( const CMatDouble & M)
{
    CMatDouble A = M;
    double pivot = 1, coef = 0, max, det = 1 ;
    unsigned int i, j, k, l, indice_max;
    unsigned int n = A.size ();

// Pivot de Gauss
    for (k = 0 ; k < n-1 ; ++k)
    {
        indice_max = k, max = A[k][k];
        for (l = k+1 ; l < n ; l++)
        {
            if (abs (A[l][k]) > abs (max) )
            {
                max = A[l][k];
                indice_max = l;
            }
        }
        if(indice_max > k)
            {
            permute(A[k] , A[indice_max]);
            }
        pivot = A[k][k];
        if (pivot == 0)
            return 0.0;
        det *= pivot;

// Triangularisation de la matrice A
        for ( i = k+1 ; i < n ; ++i)
           {
             coef = A[i][k];
              for ( j = 0 ; j < n ; ++j)
                  A[i][j] -= (coef * A[k][j] / pivot);
           }
        } // Matrice Triangulaire
    det *= A[n-1][n-1];
    return det;
}

/***********************************************************************************/

// Resolution du système linéaire A*X = B
void Resolve(const CMatDouble  & A , const CVDouble  & B, CVDouble & X)
{
    double somme = 0;
    int ind = 0, j;
    unsigned int n = B.size();
    X.clear ();
    X.resize (n);
// Resolution du système linéaire
    X[n - 1] = B[n-1]/A[n - 1][n - 1];
    for (ind = n-2 ; ind >= 0 ; ind--)
    {
        somme = B[ind];
        for (j = ind + 1 ; j < (int) n ; ++j)
            somme -= A[ind][j] * X[j];
        X[ind]  =  somme / A[ind][ind];
    }
}

/***********************************************************************************/

// Matrice Inverse
bool Inverse ( const CMatDouble & B, CMatDouble & M) //throw (Exception)
{
    bool a;
    unsigned int i , j , n = B.size();
    CMatDouble I(n), A = B;

    M.clear ();
    M.resize (n);

    for (i = 0 ; i < n ; ++i)
        {
            M[i] = CVDouble (n);
            I[i] = CVDouble (n);
            for (j = 0 ; j < n ; ++j)
            {
                if ( i == j) I[i][j] = 1;
                else I[i][j] = 0;
            }
        }
        a = Trig(A,I);
        if( a == 1)
         {
            for (i = 0 ; i < n ; ++i)
            {
                Resolve(A , I[i], M[i]);
             }
         }
    else
        {
            //throw Exception ("Singular matrix");
            return 0;
        }

        I.clear ();

    return 1;
}

/************************************************************************/

// Fonction de repartition empirique
double F (const CVDouble &T,const double &x, const int &n)
{
    double F = 0;

    for (int i = 0; i < n; ++i)
        if (T[i] <= x)
            F = F + 1.0;
     return F / n;
}
/************************************************************************/

// 1er quartile
double Quartil_1 (const CVDouble &T)
{
    CVDouble Tab;
    double Q1 = 0;

    for (double a:T)
        if (!std::isnan (a))
        {
            Tab.reserve (1);
            Tab.push_back(a);
        }
    for (double val:Tab)
        if (F(Tab, val, T.size ()) >= 0.25)
        {
            Q1 = val;
            break;
        }
    Tab.clear ();
    return Q1;
}
/************************************************************************/

// 3eme quartile
double Quartil_3 (const CVDouble &T)
{
    CVDouble Tab;
    double Q3 = 0;

    for (double a:T)
        if (!std::isnan (a))
        {
            //Tab.reserve (1);
            Tab.push_back (a);
        }
    for (double val:Tab)
        if (F (Tab, val, T.size ()) >= 0.75)
        {
            Q3 = val;
            break;
        }
    Tab.clear ();
    return Q3;
}
/************************************************************************/

/* detects outliers in a numeric column
 * using box-plot method
 * and replace them with NULL
 */
void boxPlotOutliersDetection (CMatDouble & M,
                                 unsigned fstd)
{
    double Q1, Q3, LU, LF;
    for (auto cVec = M.begin () ; cVec != M.end () ; ++cVec)
     {
        Q1 =  Quartil_1 (*cVec);
        Q3 =  Quartil_3 (*cVec);

        LU = Q1 - (fstd * (Q3 - Q1));
        LF = Q3 + (fstd * (Q3 - Q1));

        for (auto cVal = cVec->begin () ; cVal != cVec->end () ; ++cVal)
           {
            if (!std::isnan (*cVal))
                if ( (*cVal) < LU || (*cVal) > LF)
                {
                    (*cVal) = NAN;
                }
            }
       }
}
/************************************************************************/

/* detects outliers in a numeric column
* using algebraic method
* and replace them with NULL
*/
void algebraicOutlier (CMatDouble & M, unsigned fstd)
{

    double avg, stdev;
for (auto cVec = M.begin () ; cVec != M.end () ; ++cVec)
    {
    avg = cVec->Mean ();
    stdev = cVec->StdDev ();
    for (auto cVal = cVec->begin () ; cVal != cVec->end () ; ++cVal)
       {
        if (!std::isnan (*cVal))
            if ( abs ((*cVal - avg) / stdev) >= (double) fstd)
                    (*cVal) =  NAN;
            }
     }
}
/************************************************************************/
} // namespace Struct
