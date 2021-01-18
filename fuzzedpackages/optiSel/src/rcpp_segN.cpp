
#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericMatrix rcpp_segN(std::string pathNative, int NFileN, int NC, const arma::ivec& ArmaIndexN, const arma::vec& ArmaNkb) {
  int m, i, j;
  char str2[100];
  FILE *fN;
  Rcpp::NumericMatrix ArmasegN(NC, NC);
  
  int M = ArmaNkb.n_elem;
  size_t bufsize = 2*NFileN;  
  char* Line = (char*)malloc(bufsize*sizeof(char));
  if(Line == NULL){error_return("Memory allocation failed.");};
  
  double** fsegN = (double**)calloc(NC,sizeof(double*));            /*  NCxNC - matrix */
  int* indexN    = (int*)calloc(NC,sizeof(int));                    /*     NC - vector */
  int* Nat       = (int*)calloc(NC,sizeof(int));                    /*     NC - vector */
  double* Nkb    = (double*)calloc(ArmaNkb.n_elem, sizeof(double)); /*    M+1 - vector */
  if(fsegN   == NULL){error_return("Memory allocation failed.");};
  if(indexN  == NULL){error_return("Memory allocation failed.");};
  if(Nat     == NULL){error_return("Memory allocation failed.");};
  if(Nkb     == NULL){error_return("Memory allocation failed.");};
  
  for(m=0;m<M;m++){
    Nkb[m]  = ArmaNkb.at(m);
  }
  
  for(i=0; i<NC;i++){
    indexN[i] = ArmaIndexN.at(i);
    fsegN[i]  = (double*)calloc(i+1, sizeof(double));
    if(fsegN[i] == NULL){error_return("Memory allocation failed.");};
  }
  
  /* ******* Main part ******** */
  fN = fopen(pathNative.c_str(),"r");
  if(fN == NULL){error_return("File opening failed.");}; 
  while(fgetc(fN)!='\n'){}
  
  m=0;
  while(fscanf(fN, "%s ", str2)>0){
    if(fgets(Line, 2*NFileN, fN)!=NULL){
      for(i=0; i<NC;i++){
        Nat[i] = ((Line[2*indexN[i]]=='1')?1:0);
      }
      for(i=0; i<NC;i++){
        if(Nat[i]>0){
          for(j=0; j<i+1; j++){
            if(Nat[j]>0){
              fsegN[i][j] += Nkb[m];
            }
          }
        }
      }
      m=m+1;
    }
  }
  fclose(fN);	
  
  Rprintf("M=%d\n",m);
  
  for(i=0; i<NC;i++){
    for(j=0; j<i+1; j++){
      ArmasegN.at(j,i) = fsegN[i][j];
      ArmasegN.at(i,j) = fsegN[i][j];
    }
  }

  for(i=0; i<NC;i++){free(fsegN[i]);}
  free(Nat); 
  free(fsegN);
  free(Nkb);
  free(indexN);
  free(Line);
  
  return ArmasegN;
}
