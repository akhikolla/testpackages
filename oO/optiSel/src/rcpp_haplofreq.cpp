
#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::List rcpp_haplofreq(std::string pathThisBreed,std::string pathRefBreeds, std::string pathFreq, std::string pathOrig, std::vector< std::string > MarkerName, std::string stdBreedSymbol, const arma::ivec& ArmaIndexC, const arma::imat& ArmaIndexR, int NFileC,int NFileR, int NC, const arma::ivec& ArmaNR, int minSNP, double minL, double ubFreq, const arma::vec& ArmaPos, std::string  stdsymB, int skip, int cskip, int getFreq, int getOrig) {
  /* ***** initialize variables ****** */
  int m, mx, ms, i, j, r, rK, b, gleich, BreedIndex, endoffile;
  double L, maxFreq, thisFreq;
  char str1[100];
  FILE *fC, *fR, *fFreq, *fOrig;
  char merge[2];
  char symB = stdsymB.at(0);
  int M  = ArmaPos.n_elem - 1;
  
  int K            = (minSNP<=60)?(minSNP/2):(30);
  int B            = ArmaIndexR.n_cols;
  int saveOrig     = (pathOrig.length()>0);
  int saveFreq     = (pathFreq.length()>0);
  int returnResult = ((saveOrig || saveFreq)?0:1);
  int returnOrig   = getOrig && returnResult;
  int returnFreq   = getFreq && returnResult;
  arma::mat ArmaHapFreq(M*returnFreq, NC*returnFreq);
  Rcpp::CharacterMatrix RcppHapOrig(M*returnOrig, NC*returnOrig);
    
  char* BreedSymbol= (char*)malloc((B+1)*sizeof(char));   /*      B+1 - vector */
  char* wLine      = (char*)malloc((2*NC)*sizeof(char));  /*      2NC - vector */
  char* fLine      = (char*)malloc((6*NC)*sizeof(char));  /*      6NC - vector */
  char* smaxFreq   = (char*)malloc(6*sizeof(char));       /*        6 - vector */
  if(BreedSymbol == NULL){error_return("Memory allocation failed.");};
  if(wLine       == NULL){error_return("Memory allocation failed.");};
  if(fLine       == NULL){error_return("Memory allocation failed.");};
  if(smaxFreq    == NULL){error_return("Memory allocation failed.");};

  size_t bufsize = 2*(NFileC+NFileR);  
  char* rLine = (char*)malloc(bufsize*sizeof(char));
  if(rLine == NULL){error_return("Memory allocation failed.");};
  
  int*** HapCount  = (int***)calloc(M, sizeof(int**));    /* M xNCxB  - matrix */
  int*** thisROH   = (int***)calloc(NC,sizeof(int**));    /* NCxB xNR - matrix */
  double* Pos      = (double*)calloc(M+1, sizeof(double));/*    (M+1) - vector */
  int* saveLimit   = (int*)calloc(M, sizeof(int));        /*        M - vector */
  int* NR          = (int*)calloc(B, sizeof(int));        /*        B - vector */
  if(HapCount  == NULL){error_return("Memory allocation failed.");};
  if(thisROH   == NULL){error_return("Memory allocation failed.");};
  if(Pos       == NULL){error_return("Memory allocation failed.");};
  if(saveLimit == NULL){error_return("Memory allocation failed.");};
  if(NR        == NULL){error_return("Memory allocation failed.");};
  
  int*  currAllelesC = (int*) calloc(NC, sizeof(int));      /*       NC - vector */
  int** currAllelesR = (int**)calloc(B, sizeof(int*));      /*        B - matrix */
  int*  prevAllelesC = (int*) calloc(NC, sizeof(int));      /*       NC - vector */
  int** prevAllelesR = (int**)calloc(B, sizeof(int*));      /*        B - matrix */
  int*  indexC       = (int*) calloc(NC, sizeof(int));      /*       NC - vector */
  int** indexR       = (int**)calloc(B, sizeof(int*));      /*     BxNR - matrix */
  if(currAllelesC== NULL){error_return("Memory allocation failed.");};	
  if(currAllelesR== NULL){error_return("Memory allocation failed.");};
  if(prevAllelesC== NULL){error_return("Memory allocation failed.");};
  if(prevAllelesR== NULL){error_return("Memory allocation failed.");};
  if(indexC      == NULL){error_return("Memory allocation failed.");};
  if(indexR      == NULL){error_return("Memory allocation failed.");};

  strcpy(BreedSymbol, stdBreedSymbol.c_str());
  
  for(b=0;b<B;b++){
    NR[b]          = ArmaNR.at(b);
    currAllelesR[b]= (int*)calloc(NR[b], sizeof(int));
    prevAllelesR[b]= (int*)calloc(NR[b], sizeof(int));
    indexR[b]      = (int*)calloc(NR[b], sizeof(int));
    if(currAllelesR[b]== NULL){error_return("Memory allocation failed.");};
    if(prevAllelesR[b]== NULL){error_return("Memory allocation failed.");};
    if(indexR[b]      == NULL){error_return("Memory allocation failed.");};
    for(j=0;j<NR[b];j++){indexR[b][j] = ArmaIndexR.at(j,b);}
  }
  
  for(m=0; m<M+1;   m++){Pos[m]=ArmaPos.at(m);}
  for(i=0; i<NC;    i++){indexC[i]=ArmaIndexC.at(i);}
  for(i=0; i<6*NC-1;i++){fLine[i]=' ';} 
  for(i=0; i<2*NC-1;i++){wLine[i]=' ';} 
  wLine[2*NC-1]= '\0';
  fLine[6*NC-1]= '\0';
  merge[1]     = '\0';
  smaxFreq[5]  = '\0';
  
  for(i=0; i<NC;i++){
    thisROH[i]= (int**)calloc(B, sizeof(int*));
    if(thisROH[i]== NULL){error_return("Memory allocation failed.");};
    for(b=0;b<B;b++){
      thisROH[i][b]= (int*)calloc(NR[b], sizeof(int));
      if(thisROH[i][b]== NULL){error_return("Memory allocation failed.");};
    }
  }
  
  /* ***** For every m get saveLimit  ****** */
  
  ms = -1;
  m  = 0;
  while(m<M){
    while(Pos[ms+1]+minL<Pos[m]){ms=ms+1;}
    saveLimit[m] = ((ms<m-minSNP)?(ms):(m-minSNP));
    m=m+1;
  }
  
  /* ******* Main part ******** */
  
  fC = fopen(pathThisBreed.c_str(),"r"); if(fC== NULL){error_return("File opening failed.");};	 
  fR = fopen(pathRefBreeds.c_str(),"r"); if(fR== NULL){error_return("File opening failed.");};	 
  if(saveFreq){fFreq = fopen(pathFreq.c_str(), "a"); if(fFreq== NULL){error_return("File opening failed.");};}else{fFreq=fC;/* avoid warning */}
  if(saveOrig){fOrig = fopen(pathOrig.c_str(), "a"); if(fOrig== NULL){error_return("File opening failed.");};}else{fOrig=fC;/* avoid warning */}
  
  for(i=0;i<skip+1;i++){
    while(fgetc(fC)!='\n'){}
    while(fgetc(fR)!='\n'){}
   }
  endoffile=0;
  m=0;
  ms=0;
  while(!endoffile){
    /* *** Determine previous alleles and current alleles (K at a time) *** */
    /* ***         for candidates and reference individuals             *** */
    for(i=0; i<NC;i++){
      prevAllelesC[i] = currAllelesC[i];
      currAllelesC[i] = 0;
    }
    for(b=0;b<B;b++){
      for(j=0; j<NR[b];j++){
        prevAllelesR[b][j] = currAllelesR[b][j];
        currAllelesR[b][j] = 0;
      }
    }
    rK=0;
    while(rK<K){
      for(i=0; i<cskip; i++){
        endoffile = fscanf(fC, "%s ", str1)<1;
        if(endoffile){break;}
        endoffile = fscanf(fR, "%s ", str1)<1;
        if(endoffile){break;}
      }
      if(endoffile){break;}
      
      endoffile = fgets(rLine,2*NFileC,fC)==NULL;
      if(endoffile){break;}
      for(i=0; i<NC;i++){
        if(rLine[2*indexC[i]]==symB){currAllelesC[i]= currAllelesC[i] | (1u<<rK);}
      }
      endoffile = fgets(rLine,2*NFileR,fR)==NULL;
      if(endoffile){break;}
      for(b=0;b<B;b++){
        for(j=0; j<NR[b];j++){
          if(rLine[2*indexR[b][j]]==symB){currAllelesR[b][j]= currAllelesR[b][j] | (1u<<rK);}
        }
      }
      HapCount[m+rK] = (int**)calloc(NC, sizeof(int*));
      if(HapCount[m+rK]== NULL){error_return("Memory allocation failed.");};
      for(i=0;i<NC;i++){
        HapCount[m+rK][i] = (int*)calloc(B, sizeof(int));
        if(HapCount[m+rK][i]== NULL){error_return("Memory allocation failed.");};
      }      
      rK++;
    }
    /* if(endoffile){break;}*/
    if(rK==0){break;}
    /* *** Variable thisROH[i][b][j] contains the number of markers included in the current    *** */
    /* *** ROH that have not yet been used to increase haplotype counts                        *** */
    /* *** - Increase length of the current ROH for each pair of haplotypes if possible        *** */
    /* *** - Increase Haplotype Count when length of a ROH exceeds the minimum length          *** */
    /* *** - Decrease length of the current ROH when the corresponding marker has been counted *** */
    
    for(i=0; i<NC;i++){
      for(b=0; b<B;b++){
        for(j=0; j<NR[b]; j++){
          if(currAllelesC[i]==currAllelesR[b][j]){
            if(prevAllelesC[i]==prevAllelesR[b][j] && m>0){ /* ROH verlängern */
              thisROH[i][b][j] += rK;
              /* ******************************************* */
              while(thisROH[i][b][j]>minSNP && (Pos[m+rK]-Pos[m+rK-thisROH[i][b][j]+1])>minL){
                HapCount[m+rK-thisROH[i][b][j]][i][b] += 1;
                thisROH[i][b][j] -= 1;
              }
              /* ******************************************* */
              
            }else{  /* neuer ROH */
              thisROH[i][b][j] = rK;
              if(m>0){
                gleich = ~(prevAllelesC[i] ^ prevAllelesR[b][j]);
                r = K-1;
                while(r>=0 && ((gleich>>r)&1u)){
                  thisROH[i][b][j] += 1;
                  r--;
                }
              }
            }
          }else{
            if(prevAllelesC[i]==prevAllelesR[b][j] && m>0){ /* ROH beenden */
              gleich = ~(currAllelesC[i] ^ currAllelesR[b][j]);
              r = 0;
              while(r<K && ((gleich>>r)&1u)){
                thisROH[i][b][j] += 1;
                r++;
              } 
              if(thisROH[i][b][j]>=minSNP){
                L = Pos[m+r]-Pos[m+r-thisROH[i][b][j]];
                if(L>=minL){
                  /* ******************************** */
                  for(mx=m+r-thisROH[i][b][j];mx<m+r;mx++){
                    HapCount[mx][i][b] += 1;
                    }
                  /* ******************************** */
                  }
              }
              thisROH[i][b][j] = 0;
            }
          }
        }
      }
    }
    
    
    /* ***  Save Haplotype Frequences until   *** */
    /* ***  the current saveLimit is reached  *** */
    
    /* Rprintf("m=%d, rK=%d, ms=%d, saveL=%d, M=%d\n",m,rK, saveLimit[m],ms,M);*/
    while(ms<saveLimit[m+rK-1]){
      if(saveOrig){fputs(MarkerName.at(ms).c_str(), fOrig);}
      if(saveFreq){fputs(MarkerName.at(ms).c_str(), fFreq);}
      for(i=0; i<NC; i++){
        maxFreq = 0.0;
        BreedIndex = -1;
        for(b=0;b<B;b++){
          thisFreq = ((double)HapCount[ms][i][b])/((double)NR[b]);
          if(thisFreq>maxFreq){maxFreq=thisFreq;BreedIndex=b;}
        }
        if(maxFreq<ubFreq){BreedIndex=-1;}
        free(HapCount[ms][i]);
        if(returnFreq){
          ArmaHapFreq(ms, i) = maxFreq;
        }        
        if(returnOrig){
          if(BreedIndex<0){
            RcppHapOrig(ms, i) = "1";
          }else{
            merge[0] = BreedSymbol[BreedIndex];
            RcppHapOrig(ms, i) = merge;
          }
        }
        if(saveFreq){
          sprintf(smaxFreq, "%0.3f", maxFreq);
          fLine[6*i+0]=smaxFreq[0];
          fLine[6*i+1]=smaxFreq[1];
          fLine[6*i+2]=smaxFreq[2];
          fLine[6*i+3]=smaxFreq[3];
          fLine[6*i+4]=smaxFreq[4];
        }
        if(saveOrig){
          if(BreedIndex<0){wLine[2*i]='1';}else{wLine[2*i]=BreedSymbol[BreedIndex];}
        }
      }
      free(HapCount[ms]);
      if(saveOrig){
        fputs(wLine, fOrig);
        fputs("\n",fOrig);
      }
      if(saveFreq){
        fputs(fLine, fFreq);
        fputs("\n",fFreq);
      }      
      ms++;
    }
    m=m+rK;
  }
  fclose(fC);
  fclose(fR);
  
  
  /* *** Update Haplotype counts for the remaining part of the Chromosome ** */
  
  for(i=0; i<NC;i++){
    for(b=0;b<B;b++){
      for(j=0; j<NR[b]; j++){
        if(thisROH[i][b][j]>=minSNP){
          L = Pos[M]-Pos[M-thisROH[i][b][j]];
          if(L>=minL){
            /* ************************** */
            for(mx=M-thisROH[i][b][j];mx<M;mx++){
              HapCount[mx][i][b] += 1;
              }
            /* ************************** */
            }
          }
        }
      }
    }
  
  /* *** Save Haplotype frequencies for the remaining part of the chromosome ** */
  
  while(ms<M){
    if(saveOrig){fputs(MarkerName.at(ms).c_str(), fOrig);}
    if(saveFreq){fputs(MarkerName.at(ms).c_str(), fFreq);}
    for(i=0; i<NC; i++){
      maxFreq = 0.0;
      BreedIndex = -1;
      for(b=0;b<B;b++){
        thisFreq = ((double)HapCount[ms][i][b])/((double)NR[b]);
        if(thisFreq>maxFreq){maxFreq=thisFreq;BreedIndex=b;}
      }
      if(maxFreq<ubFreq){BreedIndex=-1;}
      free(HapCount[ms][i]);
      if(returnFreq){
        ArmaHapFreq(ms, i) = maxFreq;
      }
      if(returnOrig){
        if(BreedIndex<0){
          RcppHapOrig(ms, i) = "1";
        }else{
          merge[0] = BreedSymbol[BreedIndex];
          RcppHapOrig(ms, i) = merge;
        }
      }
      if(saveFreq){
        sprintf(smaxFreq, "%0.3f", maxFreq);
        fLine[6*i+0]=smaxFreq[0];
        fLine[6*i+1]=smaxFreq[1];
        fLine[6*i+2]=smaxFreq[2];
        fLine[6*i+3]=smaxFreq[3];
        fLine[6*i+4]=smaxFreq[4];
      }
      if(saveOrig){
        if(BreedIndex<0){wLine[2*i]='1';}else{wLine[2*i]=BreedSymbol[BreedIndex];}
      }
    }
    free(HapCount[ms]);
    if(saveOrig){
      fputs(wLine, fOrig); 
      fputs("\n",fOrig);
    }
    if(saveFreq){
      fputs(fLine, fFreq);
      fputs("\n",fFreq);
    }
    ms=ms+1;
  }
  if(saveOrig){    fclose(fOrig);}
  if(saveFreq){fclose(fFreq);}
  
  /* *** Release memory  *** */
  
  for(i=0; i<NC;i++){
    for(b=0;b<B;b++){
      free(thisROH[i][b]);
    }
    free(thisROH[i]);
  }
  free(thisROH);

  
  for(b=0;b<B;b++){
    free(indexR[b]);
    free(currAllelesR[b]);
    free(prevAllelesR[b]);
  }
  free(indexR);
  free(currAllelesR);
  free(prevAllelesR);
  
  free(HapCount);
  free(Pos);
  free(indexC);
  free(currAllelesC);
  free(prevAllelesC);
  free(BreedSymbol);
  free(NR);
  free(smaxFreq);
  free(rLine);
  free(wLine);
  free(fLine);
  free(saveLimit);
  
  return Rcpp::List::create(Rcpp::Named("freq")=ArmaHapFreq, Rcpp::Named("match")=RcppHapOrig);
}
