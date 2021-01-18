#define TMB_LIB_INIT R_init_mumm
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data Section*/
  DATA_MATRIX(X);
  DATA_MATRIX(Xnu);
  DATA_SPARSE_MATRIX(Z);
  DATA_MATRIX(rfac);
  DATA_MATRIX(ffac);
  DATA_VECTOR(y);
  DATA_IVECTOR(npar); /* The number of levels in all of the random effects belonging to a */
  DATA_IVECTOR(nlevelsf); /* The number of levels in all of the fixed effects in the multiplicative terms (nu) */
  DATA_IVECTOR(nlevelsr); /* The number of levels in all of the random effects in the multiplicative terms* (b) */
  DATA_IVECTOR(sizenu);
  DATA_IVECTOR(indexIna); /* The index for the part of a that is correlated with b */
  DATA_SCALAR(indexInSiga);  /* The index in sigma_a for the part of a that is correlated with b */


  using CppAD::Integer;

  /* Parameter Section */
  PARAMETER_VECTOR(beta); /*related to X */
  PARAMETER_VECTOR(a);      /*related to Z */
  PARAMETER_VECTOR(b);      /*random scalling coef*/
  PARAMETER_VECTOR(nu);

  PARAMETER_VECTOR(log_sigma_a);
  PARAMETER(log_sigma_b);
  PARAMETER(log_sigma);
  PARAMETER(transf_rho);

  Type sigma = exp(log_sigma);
  vector<Type> sigma_a = exp(log_sigma_a);
  Type sigma_b = exp(log_sigma_b);
  Type rho = transf_rho / sqrt(1. + transf_rho*transf_rho);


  /* Including the estimates and standard errors on natural scale */
  ADREPORT(sigma);
  ADREPORT(sigma_a);
  ADREPORT(sigma_b);
  ADREPORT(rho);

  Type nll = 0;
  int nobs = y.size();

  /* The linear part of the model */
  vector<Type> linmodel = X*beta + Z*a + Xnu*nu;


  /* The multiplicative part of the model */

  vector<Type> mult(nobs);      /*"nu"*/
  int jumpf = 0;
  int jumpr = 0;
  int indx = 0;

  /*Going through all of the multiplicative terms*/

  /*The first multiplicative term j = 0 */

  vector<Type> nuj(nlevelsf[0]);
  vector<Type> bj(nlevelsr[0]);
  int start = nlevelsf[0]-sizenu[0];

  for (int i= 0; i<nuj.size(); i++) {   /*to overcome identifiability issues*/
    nuj[i] = 0;
  }


  /*Build nu for the j'th multiplicative term*/
  for (int l = start; l<nlevelsf[0]; l++) {
    indx = l-start+jumpf;
    nuj[l] = nu[indx];
  }

  indx = 0;
  /*Build b for the j'th multiplicative term*/
  for (int l=0; l<nlevelsr[0]; l++) {
    indx = l + jumpr;
    bj[l] = b[indx];
  }

  vector<Type> nu2 = nuj - nuj.sum()/nuj.size();

  /*Going through all of the observations */
  for (int i = 0; i<nobs; i++){

    mult[i] += bj[Integer(rfac(i,0))-1]*nu2[Integer(ffac(i,0))-1];

  }
  jumpf += sizenu[0];
  jumpr += nlevelsr[0];


  /* The rest of the multiplicative terms - Not in use yet! */

  for (int j = 1; j<nlevelsf.size(); j++){
    vector<Type> nuj(nlevelsf[j]);
    vector<Type> bj(nlevelsr[j]);
    indx = 0;
    int start = nlevelsf[j]-sizenu[j];

    for (int i= 0; i<nuj.size(); i++) {   /*to overcome identifiability issues*/
      nuj[i] = 0;
    }

    /*Build nu for the j'th multiplicative term*/
    for (int l = start; l<nlevelsf[j]; l++) {
      indx = l-start+jumpf;
      nuj[l] = nu[indx];
    }


    indx = 0;
    /*Build b for the j'th multiplicative term*/
    for (int l=0; l<nlevelsr[j]; l++) {
      indx = l + jumpr;
      bj[l] = b[indx];
    }

    vector<Type> nu2 = nuj - nuj.sum()/nuj.size();

    /*Going through all of the observations */
    for (int i = 0; i<nobs; i++){

      mult[i] += bj[Integer(rfac(i,j))-1]*nu2[Integer(ffac(i,j))-1];

    }
    jumpf += sizenu[j];
    jumpr += nlevelsr[j];
  }



  /* The negative joint log-likelihood function */

  for(int i=0; i<nobs; i++){

    nll -= dnorm(y[i],
                 linmodel[i]+
                   mult[i],sigma,true);

  }

  int index_count = 0;
  int index = 0;



  vector<Type> acor (indexIna[1]-indexIna[0]+1);
  vector<Type> anotcor (a.size()-acor.size());
  matrix<Type> Sigma(2,2);


  //if no part of a is related to b
  if (indexInSiga == 0) {

    anotcor = a;
    index_count = 0;

    /* Going through all of the random effects in b */
     //for(int i=0; i<sigma_b.size(); i++){

     for(int j=0; j<nlevelsr[0]; j++){
     index = j+index_count;
     nll -= dnorm(b[index], Type(0), sigma_b, true);
     }


     //index_count += nlevelsr[i];
     //}


  }
  else {

    //Making a subvector of a, with the elements that are correlated with b

    for (unsigned i=0; i<acor.size(); i++){
      acor[i]=a[indexIna[0]-1+i];
    }
    //Making a subvector of a, with the elements that are NOT correlated with b

    for (unsigned i=0; i<(indexIna[0]-1); i++){
      anotcor[i]=a[i];
    }
    for (unsigned i=indexIna[1]; i<a.size(); i++){
      anotcor[i-indexIna[1]+indexIna[0]-1]=a[i];
    }


    Sigma(0,0)=(sigma_a[Integer(indexInSiga)-1])*(sigma_a[Integer(indexInSiga)-1]);
    Sigma(1,1)=sigma_b*sigma_b;
    Sigma(0,1)=rho*sigma_a[Integer(indexInSiga)-1]*sigma_b;
    Sigma(1,0)=Sigma(0,1);


  }

  index_count = 0;

  /* Going through all of the random effects in a */
  for(int i=0; i<sigma_a.size(); i++){

    if (i == (indexInSiga-1)) {
      for(int k=0; k<acor.size(); k++) {
        vector<Type> ab(2); ab(0)=acor(k); ab(1)=b(k);
        nll += density::MVNORM(Sigma)(ab);
      }
    }
    else {
      for(int j=0; j<npar[i]; j++){
        index = j+index_count;
        nll -= dnorm(anotcor[index], Type(0), sigma_a[i], true);
      }
      index_count += npar[i];
    }
  }




  return nll;
}

