
#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

Rcpp::NumericMatrix rcpp_segBreedComp(std::vector<std::string> pathNative, int Nfile, int N, const arma::ivec& ArmaIndexN, const arma::ivec& MatChr, const arma::vec& Armakb) {
  int m, i, iB, nB, cB, iSNP;
  unsigned int chr;
  int M = Armakb.n_elem;
  char str[100], merge[2];
  FILE *fN;

  size_t bufsize = 2*Nfile;  
  char* Line = (char*)malloc(bufsize*sizeof(char));
  if(Line == NULL){error_return("Memory allocation failed.");};
   
  int* indexN        = (int*)calloc(N,sizeof(int));                    /*     N - vector */
  double* kb         = (double*)calloc(Armakb.n_elem, sizeof(double)); /*    MatChr+1 - vector */
  double** BreedCont = (double**)calloc(255, sizeof(double*));  
  int* hasCont       = (int*)calloc(255,sizeof(int));                    /*     vector */
  
  if(indexN   == NULL){error_return("Memory allocation failed.");};
  if(hasCont  == NULL){error_return("Memory allocation failed.");};
  if(kb       == NULL){error_return("Memory allocation failed.");};
  if(BreedCont== NULL){error_return("Memory allocation failed.");};
  for(i=0;i<255;i++){
    hasCont[i]   = 0;
    BreedCont[i] = (double*)calloc(N, sizeof(double));
    if(BreedCont[i]== NULL){error_return("Memory allocation failed.");};
  }
  
  for(m=0;m<M;m++){kb[m]     = Armakb.at(m);}
  for(i=0;i<N;i++){indexN[i] = ArmaIndexN.at(i);}
  merge[1] = '\0';
  
  /* ******* Main part ******** */
  iSNP = 0;
  for(chr=0;chr<pathNative.size();chr++){
    fN = fopen(pathNative[chr].c_str(),"r");
    if(fN== NULL){error_return("File opening failed.");}; 
    while(fgetc(fN)!='\n'){}
    m=0;
    while(fscanf(fN, "%s ", str)>0){
      if(fgets(Line, 2*Nfile, fN)!=NULL){
        for(i=0; i<N; i++){
          iB = (int) Line[2*indexN[i]];
          hasCont[iB] = 1;
          BreedCont[iB][i] += kb[iSNP];
        }
      m = m + 1;
      iSNP = iSNP + 1;
      }
    }
    if(m!=MatChr.at(chr)){Rprintf("Numbers of marker in map and file are not equal at chromosome %d\n",chr);}
    fclose(fN);
  }
  
 /* Rprintf("MatChr=%d\n",m); */
  nB = 0;
  for(i=0;i<255;i++){
    nB = nB + hasCont[i];
  }
  Rcpp::NumericMatrix RcppBreedCont(N, nB);
  Rcpp::CharacterVector cnamesBreedCont(nB);
  
  cB = 0;
  for(iB=0; iB<255; iB++){
    if(hasCont[iB]){
      merge[0]= (char) iB;
      cnamesBreedCont(cB) = merge;
      if(merge[0]=='1'){
        cnamesBreedCont(cB) = "native";
      }
      for(i=0; i<N;i++){
        RcppBreedCont.at(i, cB) = BreedCont[iB][i];
      }
      cB = cB + 1;
    }
    free(BreedCont[iB]);
  }
  colnames(RcppBreedCont) = cnamesBreedCont;
  free(BreedCont);
  free(kb);
  free(indexN);
  free(hasCont);
  
  return RcppBreedCont;
}
