
#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericMatrix rcpp_segIBDandN(std::string pathThisBreed, std::string pathNative, int NFileC, int NFileN, const arma::ivec& ArmaIndexC, const arma::ivec& ArmaIndexN, int NC, int minSNP, double minL, const arma::vec& ArmaPos, const arma::vec& Armakb, double a, std::string stdsymB, int skip, int cskip) {
  /* ***** initialize variables ****** */
  int m, i, j, r, r2, rK, endoffile, gleich;
  double L;
  char str1[100];
  FILE *fC, *fN;
  char symB = stdsymB.at(0);
  Rcpp::NumericMatrix confROH(NC, NC);
  int K  = (minSNP<=60)?(minSNP/2):(30);
  int M  = Armakb.n_elem - 1;
  
  size_t bufsize = 2*(NFileC+NFileN);  
  char* Line = (char*)malloc(bufsize*sizeof(char));
  if(Line == NULL){error_return("Memory allocation failed.");};
  
  int** Nat         = (int**)calloc(M,sizeof(int*));                   /*  M xNC - matrix */
  double** fROH     = (double**)calloc(NC,sizeof(double*));            /*  NCxNC - matrix */
  int** thisROH     = (int**)calloc(NC,sizeof(int*));                  /*  NCxNC - matrix */
  double** lSEG     = (double**)calloc(NC,sizeof(double*));            /*  NCxNC - matrix */
  int* currAllelesC = (int*)calloc(NC,sizeof(int));                    /*     NC - vector */
  int* prevAllelesC = (int*)calloc(NC,sizeof(int));                    /*     NC - vector */
  int* indexC       = (int*)calloc(NC,sizeof(int));                    /*     NC - vector */
  int* indexN       = (int*)calloc(NC,sizeof(int));                    /*     NC - vector */
  double* Pos       = (double*)calloc(ArmaPos.n_elem, sizeof(double)); /*    M+1 - vector */
  double* kb        = (double*)calloc(Armakb.n_elem, sizeof(double));  /*    M+1 - vector */
  if(fROH        == NULL){error_return("Memory allocation failed.");};
  if(Nat         == NULL){error_return("Memory allocation failed.");};
  if(thisROH     == NULL){error_return("Memory allocation failed.");};
  if(lSEG        == NULL){error_return("Memory allocation failed.");};
  if(currAllelesC== NULL){error_return("Memory allocation failed.");};
  if(prevAllelesC== NULL){error_return("Memory allocation failed.");};
  if(indexC      == NULL){error_return("Memory allocation failed.");};
  if(indexN      == NULL){error_return("Memory allocation failed.");};
  if(Pos         == NULL){error_return("Memory allocation failed.");};
  if(kb          == NULL){error_return("Memory allocation failed.");};
  
  for(m=0;m<M+1;m++){
    Pos[m] = ArmaPos.at(m);
    kb[m]  = Armakb.at(m);
  }
  
  for(i=0; i<NC;i++){
    indexC[i] = ArmaIndexC.at(i);
    indexN[i] = ArmaIndexN.at(i);
    fROH[i]   = (double*)calloc(i+1,sizeof(double));
    thisROH[i]= (int*)calloc(i+1,sizeof(int));
    lSEG[i]   = (double*)calloc(i+1,sizeof(double));
    if(fROH[i]   == NULL){error_return("Memory allocation failed.");};
    if(thisROH[i]== NULL){error_return("Memory allocation failed.");};
    if(lSEG[i]   == NULL){error_return("Memory allocation failed.");};
  }
  
  /* ******* Main part ******** */
  fC = fopen(pathThisBreed.c_str(),"r");
  fN = fopen(pathNative.c_str(),"r");
  if(fC == NULL){error_return("File opening failed.");};
  if(fN == NULL){error_return("File opening failed.");};
  while(fgetc(fN)!='\n'){}
  for(i=0;i<skip+1;i++){
    while(fgetc(fC)!='\n'){}
  }
  endoffile=0;
  m=0;
  while(!endoffile){
    /* *** Determine previous alleles and current alleles (K at a time) *** */
    /* ***           and native alleles for candidates                  *** */
    for(i=0; i<NC;i++){
      prevAllelesC[i] = currAllelesC[i];
      currAllelesC[i] = 0;
    }
    rK=0;
    while(rK<K){
      for(i=0; i<cskip; i++){
        endoffile = fscanf(fC, "%s ", str1)<1;
        if(endoffile){break;}
      }
      if(endoffile){break;}
      endoffile = fscanf(fN, "%s ", str1)<1;
      if(endoffile){break;}
      endoffile = fgets(Line,2*NFileC,fC)==NULL;
      if(endoffile){break;}
      for(i=0; i<NC;i++){
        if(Line[2*indexC[i]]==symB){currAllelesC[i]= currAllelesC[i] | (1u<<rK);}
      }
      endoffile = fgets(Line, 2*NFileN, fN)==NULL;
      if(endoffile){break;}
      Nat[m+rK] = (int*)calloc(NC,sizeof(int));
      for(i=0; i<NC;i++){
        Nat[m+rK][i] = ((Line[2*indexN[i]]=='1')?1:0);
      }
      rK++;
    }
    if(endoffile){Rprintf("M=%d\n",m+rK);}
    if(rK==0){break;}
    
    for(i=0; i<NC;i++){
      for(j=0; j<i+1; j++){
        if(currAllelesC[i]==currAllelesC[j]){
          if(prevAllelesC[i]==prevAllelesC[j] && m>0){ /* ROH verlängern */
            thisROH[i][j] += rK;
            for(r2=0;r2<rK;r2++){if(Nat[m+r2][i]*Nat[m+r2][j]>0){lSEG[i][j] += kb[m+r2+1]-kb[m+r2];}} /* !!!!! */
          }else{  /* neuer ROH */
            thisROH[i][j] = rK;
            for(r2=0;r2<rK;r2++){if(Nat[m+r2][i]*Nat[m+r2][j]>0){lSEG[i][j] += kb[m+r2+1]-kb[m+r2];}} /* !!!!! */
            if(m>0){
              gleich = ~(prevAllelesC[i] ^ prevAllelesC[j]);
              r = K-1;
              while(r>=0 && ((gleich>>r)&1u)){
                thisROH[i][j] += 1;
                if(Nat[m-K+r][i]*Nat[m-K+r][j]>0){lSEG[i][j] += kb[m-K+r+1]-kb[m-K+r];} /* !!!!! */
                r--;
              }
            }
          }
        }else{
          if(prevAllelesC[i]==prevAllelesC[j] && m>0){ /* ROH beenden */
            gleich = ~(currAllelesC[i] ^ currAllelesC[j]);
            r = 0;
            while(r<K && ((gleich>>r)&1u)){
              thisROH[i][j] += 1;
              if(Nat[m+r][i]*Nat[m+r][j]>0){lSEG[i][j] += kb[m+r+1]-kb[m+r];} /* !!!!! */
              r++;
            }
            
            if(thisROH[i][j]>=minSNP){
              L = Pos[m+r]-Pos[m+r-thisROH[i][j]];
              if(L>=minL){
                fROH[i][j] += (L*L/(a+L*L))*lSEG[i][j];
                }
              }
            thisROH[i][j] = 0;
            lSEG[i][j] = 0.0;
          }
        }
      }
    }
    m=m+rK;
  }
  fclose(fC);	
  fclose(fN);	
  
  
  for(i=0; i<NC;i++){
    for(j=0; j<i+1; j++){
      if(thisROH[i][j]>=minSNP){
        L = Pos[M]-Pos[M-thisROH[i][j]];
        if(L>=minL){
          fROH[i][j] += (L*L/(a+L*L))*lSEG[i][j];
          }
        }
      confROH.at(j,i) = fROH[i][j];
      confROH.at(i,j) = fROH[i][j];
    }
  }

  for(m=0; m<M;m++){
    free(Nat[m]);
  }
  free(Nat); 
  
  for(i=0; i<NC;i++){
    free(fROH[i]);
    free(thisROH[i]);
    free(lSEG[i]);
  }
  free(fROH);
  free(thisROH);
  free(lSEG);
  free(kb);
  free(Pos);
  free(indexC);
  free(indexN);
  free(prevAllelesC);
  free(currAllelesC);
  free(Line);
  
  return confROH;
}
