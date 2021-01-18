/**
 * @authors Hmamouche Youssef
 **/

#include "../inst/include/matrixOperations.h"
#include "../inst/include/operateurs.h"
#include "../inst/include/tests.h"

using namespace std;
using namespace Struct;


namespace MatrixOperations {
bool regression (const CMatDouble & Mn, const CVDouble & Y,  CVDouble & pBeta ) //throw (Exception)
{
    CMatDouble A, B, M = Mn;
    CVDouble  Beta, Teta;
    pBeta.clear ();
    pBeta.resize (M.size ());
    bool test = true;

    MultCVDouble (Trans (M), M, B);

    test = Inverse (B , A);

    if (!test)
        return 0;

    MultCVDouble (Trans (M), Y, Teta);
    MultCVDouble (A, Teta, Beta);

    pBeta = Beta;

    Beta.clear();
    Teta.clear ();
    A.clear();
    M.clear();
    B.clear();

    return 1;
}
/***********************************************************************************/
// Building the  matrix of lagged variables of the Var(p) Model

void P_Part (CVDouble & V , CMatDouble & Present, CMatDouble & M, unsigned int  p)
{
    unsigned int  N = V.size ();
    unsigned int i;
    unsigned int j;
    CVDouble P;
    P.reserve(N - p);
    for (i = p ; i < N ; i++)
    {
        P.push_back (V[i]) ;
    }
    Present.push_back(P);
    P.clear ();

    for (j = 1 ; j <= p ; j++)
    {
        for (i = p ; i < N ; i++)
        {
            P.push_back (V[i - j]) ;
        }
        M.push_back (P);
        P.clear ();
    }
}

/*********************************************************************************/
// Building the  matrix of lagged variables of the Var(p) Model

void Pr_Part (CVDouble & V, CMatDouble & M, unsigned int  p)
{
    unsigned int  N = V.size ();
    CVDouble P ;

    for (unsigned j = 1 ; j <= p ; j++)
    {
        for (unsigned i = p + 1 ; i <= N ; i++)
        {
            P.push_back (V[i - j]) ;
        }
        M.push_back (P);
        P.clear ();
    }
}

/***********************************************************************************/

void Diff (CVDouble & V) {
    for (auto it = V.end () -1 ; it != V.begin () ; --it)
        *it -= *(it-1) ;
    V.erase (V.begin ());
}
/***********************************************************/
CVDouble VECbivar (CMatDouble  M, unsigned lag, bool d /* = false */) //throw (Exception)
{
    CMatDouble initials;
    CMatDouble residus;
    CMatDouble Alpha;
    vector <unsigned> orders, testCoint;
    CMatDouble Test;
    unsigned nlins, i, maxOrder = 0;
    CMatDouble B, A, Present, Beta, pCible;
    CVDouble  pBeta, one, trend, Aecorrection, Becorrection, Conf;
    bool test;

    unsigned nbreCol = M.size(), p = lag;
    unsigned nbreLn = M[0].size();
    double SSE;

    Conf.resize(nbreCol);
    Alpha.resize(nbreCol);
    orders.reserve (nbreCol);

    for (auto Vect:M)
        orders.push_back (order (Vect, p));

    if (d)
    {
      for (auto & in:orders){
        if (in > maxOrder)
            maxOrder = in;
        }

      initials.resize (nbreCol);

      for (unsigned i(0) ; i < nbreCol ; ++i)
      {
        initials[i].resize (maxOrder);

        for (unsigned m (0) ; m < orders[i] ; ++m)
        {
          for (auto it = M[i].end () -1;it != M[i].begin ();--it){
              *it -= *(it-1) ;}

          initials[i][m] = M[i][0];
          M[i].erase (M[i].begin ());
        }

        if (orders[i] < maxOrder)
        for (unsigned l(orders[i]) ; l < maxOrder ; ++l)
          {
              initials[i][l] = M[i][0];
              M[i].erase (M[i].begin ());
          }
      }

      // Compute the number of lines after differentiation
      nlins = nbreLn - maxOrder;
    }
    else {
        nlins = nbreLn;
      }

    for (i = p ; i < nlins ; i++) {
        one.push_back (1);
    }
    A.reserve (nbreCol * lag + 2);
    B.reserve (nbreCol * lag + 2);
    B.push_back (one);
    A.push_back (one);

    pCible.resize (nbreCol);

    for (auto vect:M)
    {
        P_Part (vect, Present, B, lag);
        Pr_Part (vect, A, lag);
    }

    for (i = 0 ; i < nbreCol ; ++i)
    {
       test = regression (B , Present[i], pBeta);
       if (!test)
            throw Exception ("Singular Matrix");

        MultCVDouble (B, pBeta, pCible[i]);
        pBeta.clear();
    } // END i

    for (i = 0 ; i < nbreCol ; ++i)
    {
      SSE = 0;
      for (unsigned m = 0 ; m < Present[i].size() ; m++){
          SSE += pow (pCible[i][m] - Present[i][m], 2);
        }
      Conf[i] = SSE;
    }

    one.clear ();
    trend.clear ();
    pCible.clear ();
    Beta.clear ();
    pBeta. clear ();
    A.clear ();
    B.clear ();
    Present.clear ();

    orders.clear();
    initials.clear();

    return Conf;
}
}
