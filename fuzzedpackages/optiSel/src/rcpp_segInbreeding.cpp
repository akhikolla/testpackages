
#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericVector rcpp_segInbreeding(std::string path1, std::string path2, int NFile1, int NFile2, const arma::ivec& ArmaIndex1, const arma::ivec& ArmaIndex2, int N1, int N2, int M, int minSNP, double minL, const arma::vec& ArmacM, const arma::vec& Armakb, double a, std::string stdsymB, int skip, int cskip) {
  int m, i, i0, j, rK, r, endoffile, gleich;
  double L;
  char str1[100];
  char symB = stdsymB.at(0);
  FILE *f1, *f2;
  int N  = N1 + N2;
  int K  = (minSNP<=60)?(minSNP/2):(30);
  Rcpp::NumericVector ArmasegInbr(N/2);
  
  size_t bufsize = 2*(NFile1+NFile2);  
  char* Line = (char*)malloc(bufsize*sizeof(char));
  if(Line == NULL){error_return("Memory allocation failed.");};
  
  
  double* fROH   = (double*)calloc(N/2,sizeof(double));
  int* thisROH   = (int*)calloc(N/2,sizeof(int));
  int* thisAllel = (int*)calloc(N,sizeof(int));
  int* prevAllel = (int*)calloc(N,sizeof(int));
  double* cM     = (double*)calloc(M+1,sizeof(double));
  double* kb     = (double*)calloc(M+1,sizeof(double));
  int* index1    = (int*)calloc(N1,sizeof(int));          /*     N1 - vector */
  int* index2    = (int*)calloc(N2,sizeof(int));          /*     N2 - vector */
  if(fROH      == NULL){error_return("Memory allocation failed.");};
  if(thisROH   == NULL){error_return("Memory allocation failed.");};
  if(thisAllel == NULL){error_return("Memory allocation failed.");};
  if(prevAllel == NULL){error_return("Memory allocation failed.");};
  if(cM        == NULL){error_return("Memory allocation failed.");};
  if(kb        == NULL){error_return("Memory allocation failed.");};
  if(index1    == NULL){error_return("Memory allocation failed.");};
  if(index2    == NULL){error_return("Memory allocation failed.");};
  
  for(i=0;i<N1;i++){index1[i]=ArmaIndex1.at(i);}
  for(i=0;i<N2;i++){index2[i]=ArmaIndex2.at(i);}
  
  for(m=0;m<M+1;m++){
    cM[m]=ArmacM.at(m);
    kb[m]=Armakb.at(m);
  }
  
  f1 = fopen(path1.c_str(),"r");
  if(f1 == NULL){error_return("File opening failed.");};
  for(i=0;i<skip+1;i++){
    while(fgetc(f1)!='\n'){}
  }
  
  if(N2>0){
    f2 = fopen(path2.c_str(),"r");
    if(f2 == NULL){error_return("File opening failed.");};
    for(i=0;i<skip+1;i++){
      while(fgetc(f2)!='\n'){}
    }
  }else{f2 = f1; /* avoid warnings */}
  
  endoffile=0;
  m=0;
  while(!endoffile){
    for(i=0; i<N;i++){
      prevAllel[i] = thisAllel[i];
      thisAllel[i] = 0;
      }
    rK=0;
    while(rK<K){
      for(i=0; i<cskip; i++){
        endoffile = fscanf(f1, "%s ", str1)<1;
        if(endoffile){break;}
      }
      if(endoffile){break;}
      endoffile = fgets(Line, 2*NFile1, f1)==NULL;
      if(endoffile){break;}
      for(i=0; i<N1;i++){
        if(Line[2*index1[i]]==symB){thisAllel[i]= thisAllel[i] | (1u<<rK);}
      }
      rK++;
    }
    if(N2>0){
      rK=0;
      while(rK<K){
        for(i=0; i<cskip; i++){
          endoffile = fscanf(f2, "%s ", str1)<1;
          if(endoffile){break;}
        }
        if(endoffile){break;}
        endoffile = fgets(Line,2*NFile2,f2)==NULL;
        if(endoffile){break;}
        for(i=0; i<N2;i++){
          if(Line[2*index2[i]]==symB){thisAllel[i+N1]= thisAllel[i+N1] | (1u<<rK);}
        }
        rK++;
      }
    }
    if(endoffile){Rprintf("M=%d\n",m+rK);}
    if(rK==0){break;}
    
    for(i0=0; i0<N/2;i0++){
      i = 2*i0;
      j = 2*i0+1;
      if(thisAllel[i]==thisAllel[j]){
          if(prevAllel[i]==prevAllel[j] && m>0){ /* ROH verlängern */
            thisROH[i0] += rK;
          }else{  /* neuer ROH */
            thisROH[i0] = rK;
            if(m>0){
              gleich = ~(prevAllel[i] ^ prevAllel[j]);
              r = K-1;
              while(r>=0 && ((gleich>>r)&1u)){
                thisROH[i0] += 1;
                r--;
              }
            }
          }
        }else{
          if(prevAllel[i]==prevAllel[j] && m>0){ /* ROH beenden */
            gleich = ~(thisAllel[i] ^ thisAllel[j]);
            r = 0;
            while(r<rK && ((gleich>>r)&1u)){
              thisROH[i0] += 1;
              r++;
            }
            
            if(thisROH[i0]>=minSNP){
              L = cM[m+r]-cM[m+r-thisROH[i0]];
              if(L>=minL){fROH[i0] += (L*L/(a+L*L))*(kb[m+r]-kb[m+r-thisROH[i0]]);}
              }
            thisROH[i0] = 0;
          }
        }
    }
    m=m+rK;
  }

  fclose(f1);	
  if(N2>0){fclose(f2);}	
  
  for(i0=0; i0<N/2;i0++){
    if(thisROH[i0]>=minSNP){
        L = cM[M]-cM[M-thisROH[i0]];
        if(L>=minL){fROH[i0] += (L*L/(a+L*L))*(kb[M]-kb[M-thisROH[i0]]);}
        }
    ArmasegInbr.at(i0) = fROH[i0];
  }
  
  free(fROH);
  free(thisROH);
  free(cM);
  free(kb);
  free(thisAllel);
  free(prevAllel);
  free(index1);
  free(index2);
  free(Line);
  
  return ArmasegInbr;
}
