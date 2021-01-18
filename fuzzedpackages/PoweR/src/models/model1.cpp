#include <R.h>
#include "Rmath.h"

extern "C" {


  void model1 (int *xlen, double *x) {

    // Les x_i's contiennent en fait les epsilon_i's

    int i, n=*(xlen+0);


    double *Y;
    Y = new double [n];


    // On fabrique les Y_i's en utilisant la formule Y_i=f(xvec,thetavec)+epsilon_i:


    for (i=0;i<n;i++) {

      *(Y+i) = 0.0 + *(x+i); // Attention!! Ici j'ai pris un modèle très simple où f==0

    }


    // On estime les paramètres thetavec par thetavecchap:

       // Là c'est vite fait puisque dans ce modèle volontairement simplifié, il n'y a pas de paramètres !

    // On fabrique les résidus en utilisant la formule epschap_i=Y_i-f(xvec,thetavecchap)
    // Ces résidus sont renvoyés dans le pointeur x


    for (i=0;i<n;i++) {

      *(x+i) = *(Y+i) - 0.0;

    }

// We free the unused array of pointers. Then we return.
    delete[] Y;
    return;
    
  }
  
}
