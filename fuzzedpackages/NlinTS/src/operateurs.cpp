/**
 * @authors Hmamouche Youssef
 **/

#include "../inst/include/operateurs.h"

using namespace std;
using namespace Struct;

/******** Produit vecteur--vecteur *********/
void MultCVDouble (const CVDouble & A, const CVDouble & B, CVDouble & Res) //throw (Exception)
{
    unsigned int i, n = A.size();
    Res.clear();
    Res.resize (n);
    
    for (i = 0 ; i < n ; ++i)
        Res[i]  = (A[i] * B[i]) ;
        }

/********* Produit Matric--vecteur *********/
void MultCVDouble (const CMatDouble & A, const CVDouble & B, CVDouble & Res) //throw (Exception)
{
    unsigned int i , j, n = A[0].size (), m = B.size ();
    Res.clear();
    Res.resize(n);
    
    for (i = 0 ; i < n ; ++i)
    {
        for (j = 0 ; j < m ; ++j)
        {
            Res[i] += (A[j][i] * B[j]) ;
        }
    }
}

/*************** Produit Matrice--Matrice  ***************/
void MultCVDouble (const CMatDouble & A, const CMatDouble & B, CMatDouble & Res) //throw (Exception)
{
    unsigned int  n = B.size();
    Res.clear ();
    Res.resize (n);
    for (unsigned int i = 0 ; i < n ; ++i)
    {
        MultCVDouble (A, B[i], Res[i]);
    }
}
