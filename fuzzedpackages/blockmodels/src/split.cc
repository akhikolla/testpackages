#include "split.h"

SEXP split(SEXP sexp_vars)
{
    Rcpp::List vars = sexp_vars;
    mat coordinates = vars["coordinates"];
    mat mass = vars["mass"];

    mat centroids(2,coordinates.n_cols);

    mat classif = zeros<mat>(coordinates.n_rows,1);

    centroids.row(0) = coordinates.row(0);
    
    {
        unsigned int k=1;
        while(all(coordinates.row(0)-coordinates.row(k)==0))
        {
            k++;
            if(k>=coordinates.n_rows)
                return Rcpp::wrap(classif);
        }
        centroids.row(1) = coordinates.row(k);
    }

    mat old_classif;

    do
    {
        old_classif = classif;

        for(unsigned int i=0;i<classif.n_rows;i++)
        {
            double d0 = norm(rowvec(coordinates.row(i)-centroids.row(0)),2);
            double d1 = norm(rowvec(coordinates.row(i)-centroids.row(1)),2);

            classif(i) = (d0<d1)?0:1;
        }

        centroids.fill(0);

        colvec S(2);
        S.fill(0);
        for(unsigned int i=0;i<classif.n_rows;i++)
        {
            centroids.row(classif(i)) += mass(i)*coordinates.row(i);
            S(classif(i))+=mass(i);

        }

        centroids.row(0) /= 1.0*S(0);
        centroids.row(1) /= 1.0*S(1);

    } while(accu(classif-old_classif!=0)!=0);

    return Rcpp::wrap(classif);
}



