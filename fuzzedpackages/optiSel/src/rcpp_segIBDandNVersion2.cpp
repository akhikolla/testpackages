
#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericMatrix rcpp_segIBDandNVersion2(std::string pathThisBreed, int NFileC, int NC, const arma::ivec& ArmaIndexC, const arma::mat& ArmaNat, int minSNP, double minL, const arma::vec& ArmaPos, const arma::vec& Armakb, double a, std::string stdsymB, int skip, int cskip) {
  int m, m2, i, j, r, rK, endoffile, gleich;
  double L, w, lSEG ;
  char str1[100];
  char symB = stdsymB.at(0);
  FILE *fC;
  Rcpp::NumericMatrix confROH(NC, NC);
  int K  = (minSNP<=60)?(minSNP/2):(30);
  int M  = Armakb.n_elem - 1;
  
  size_t bufsize = 2*NFileC;  
  char* Line = (char*)malloc(bufsize*sizeof(char));
  if(Line == NULL){error_return("Memory allocation failed.");};
  
  int** Nat         = (int**)calloc(NC,sizeof(int*));
  double** fROH     = (double**)calloc(NC,sizeof(double*));
  int** thisROH     = (int**)calloc(NC,sizeof(int*));
  int* currAllelesC = (int*)calloc(NC,sizeof(int));
  int* prevAllelesC = (int*)calloc(NC,sizeof(int));
  int* indexC       = (int*)calloc(NC,sizeof(int));
  double* Pos       = (double*)calloc(ArmaPos.n_elem, sizeof(double));
  double* kb        = (double*)calloc(Armakb.n_elem, sizeof(double));
  
  if(Nat          == NULL){error_return("Memory allocation failed.");};
  if(fROH         == NULL){error_return("Memory allocation failed.");};
  if(thisROH      == NULL){error_return("Memory allocation failed.");};
  if(currAllelesC == NULL){error_return("Memory allocation failed.");};
  if(prevAllelesC == NULL){error_return("Memory allocation failed.");};
  if(indexC       == NULL){error_return("Memory allocation failed.");};
  if(Pos          == NULL){error_return("Memory allocation failed.");};
  if(kb           == NULL){error_return("Memory allocation failed.");};
  
  for(m=0;m<M+1;m++){
    Pos[m] = ArmaPos.at(m);
    kb[m]  = Armakb.at(m);
  }
  
  for(i=0; i<NC;i++){
    indexC[i] = ArmaIndexC.at(i);
    fROH[i]   = (double*)calloc(i+1, sizeof(double));
    thisROH[i]=    (int*)calloc(i+1, sizeof(int));
    Nat[i]    =    (int*)calloc(M,   sizeof(int));
    if(fROH[i]    == NULL){error_return("Memory allocation failed.");};
    if(thisROH[i] == NULL){error_return("Memory allocation failed.");};
    if(Nat[i]     == NULL){error_return("Memory allocation failed.");};
    for(m=0; m<M;m++){
      Nat[i][m] = ArmaNat.at(m,i);
    }
  }
  
  
  fC = fopen(pathThisBreed.c_str(),"r");
  if(fC == NULL){error_return("File opening failed.");}; 
  for(i=0;i<skip+1;i++){
    while(fgetc(fC)!='\n'){}
  }
  
  endoffile=0;
  m=0;
  while(!endoffile){
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
      endoffile = fgets(Line,2*NFileC,fC)==NULL;
      if(endoffile){break;}
      for(i=0; i<NC;i++){
        if(Line[2*indexC[i]]==symB){currAllelesC[i]= currAllelesC[i] | (1u<<rK);}
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
          }else{  /* neuer ROH */
            thisROH[i][j] = rK;
            if(m>0){
              gleich = ~(prevAllelesC[i] ^ prevAllelesC[j]);
              r = K-1;
              while(r>=0 && ((gleich>>r)&1u)){
                thisROH[i][j] += 1;
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
              r++;
            }
            
            if(thisROH[i][j]>=minSNP){
              L = Pos[m+r]-Pos[m+r-thisROH[i][j]];
              if(L>=minL){
                w = L*L/(a+L*L);
                lSEG = 0.0;
                for(m2=m+r-thisROH[i][j];m2<m+r;m2++){
                  if(Nat[i][m2]*Nat[j][m2]>0){lSEG += kb[m2+1]-kb[m2];}
                }
                fROH[i][j] += w*lSEG;
                }
              }
            thisROH[i][j] = 0;
          }
        }
      }
    }
    m=m+rK;
  }

  fclose(fC);	
  
  
  for(i=0; i<NC;i++){
    for(j=0; j<i+1; j++){
      if(thisROH[i][j]>=minSNP){
        L = Pos[M]-Pos[M-thisROH[i][j]];
        if(L>=minL){
          w = L*L/(a+L*L);
          lSEG = 0.0;
          for(m2=M-thisROH[i][j];m2<M;m2++){
            if(Nat[i][m2]*Nat[j][m2]>0){lSEG += kb[m2+1]-kb[m2];}
          }
          fROH[i][j] += w*lSEG;
          }
        }
      confROH.at(j,i) = fROH[i][j];
      confROH.at(i,j) = fROH[i][j];
    }
  }

  for(i=0; i<NC;i++){
    free(Nat[i]);
    free(fROH[i]);
    free(thisROH[i]);
  }
  free(Nat);
  free(fROH);
  free(thisROH);
  free(kb);
  free(Pos);
  free(indexC);
  free(currAllelesC);
  free(prevAllelesC);
  free(Line);
  
  return confROH;
}
