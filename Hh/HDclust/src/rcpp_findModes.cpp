/*==========================================================================================*/
  /*                                                                                          */
  /* Copyright (C) [Dec 2015]-[April 2017] Jia Li, Department of Statistics,                  */
  /* The Pennsylvania State University, USA - All Rights Reserved                             */
  /*                                                                                          */
  /* Unauthorized copying of this file, via any medium is strictly prohibited                 */
  /*                                                                                          */
  /* Proprietary and CONFIDENTIAL                                                             */
  /*                                                                                          */
  /* NOTICE: All information contained herein is, and remains the property of The             */
  /* Pennsylvania State University. The intellectual and technical concepts                   */
  /* contained herein are proprietary to The Pennsylvania State University and may            */
  /* be covered by U.S. and Foreign Patents, patents in process, and are protected            */
  /* by trade secret or copyright law. Dissemination of this information or                   */
  /* reproduction of this material is strictly forbidden unless prior written                 */
  /* permission is obtained from Jia Li at The Pennsylvania State University. If              */
  /* you obtained this code from other sources, please write to Jia Li.                       */
  /*                                                                                          */
  /*                                                                                          */
  /* The software is a part of the package for                                                */
  /* Clustering with Hidden Markov Models on Variable Blocks                                  */
  /*                                                                                          */
  /* Written by Jia Li <jiali@stat.psu.edu>, April 7, 2017                                    */ 
  /*                                                                                          */
  /*==========================================================================================*/
  #include <Rcpp.h>
  #include "utils_rcpp.h"
  #include "hmm.h"
  
  #ifdef _OPENMP
  #include <omp.h>
  #endif
  
  #include <string.h>
  
  
  // [[Rcpp::plugins(openmp)]]

using namespace Rcpp;


// [[Rcpp::export]]
S4 rcpp_findModes(NumericMatrix dataTranspose, S4 HmmVb, IntegerVector nthread)
{
  //double *wt=NULL;
  
  /*----------------------------------------------------------------*/
    /*---------------- Parse input data ------------------------------*/
    /*----------------------------------------------------------------*/
  int nseq = dataTranspose.ncol(); // num datapoints
  int dim = dataTranspose.nrow(); // data dimensionality
  
  // print parsed data
  //Rcout << "nseq = " << nseq << "\n";
  //Rcout << "dim = " << dim << "\n";
  
  double **u; // 2d representation of the data. Each row - new datapoint
  
  u = (double **)calloc(nseq,sizeof(double *));
  
  for (int i = 0; i < nseq; i++){
    u[i] = dataTranspose.begin() + i*dim;
  }
  
  /*----------------------------------------------------------------*/
    /*---------------- Parse model -----------------------------------*/
    /*----------------------------------------------------------------*/
  S4 VbStructure = HmmVb.slot("VbStructure");
  DIAGCOV = as<int>(HmmVb.slot("diagCov"));
  
  /*----------------------------------------------------------------*/
    /*---------------- Parse variable block structure ----------------*/
    /*----------------------------------------------------------------*/
    CondChain *md=NULL;
  
  if (!Rf_isNull(VbStructure)){
    md=(CondChain *)calloc(1,sizeof(CondChain));
    
    parseVbStructure(VbStructure, md);
    
    //~ Rcout << "nb = " << md->nb << "\n";
    //~ 
      //~ Rcout << "bdim:\n";
    //~ for (int i = 0; i < md->nb; i++)
      //~ Rcout << md->bdim[i] << " ";
    //~ Rcout << "\n";
    //~ 
      //~ Rcout << "numst:\n";
    //~ for (int i = 0; i < md->nb; i++)
      //~ Rcout << md->numst[i] << " ";
    //~ Rcout << "\n";
    //~ 
      //~ Rcout << "var:\n";
    //~ for (int i = 0; i < md->nb; i++){
      //~ for (int j = 0; j < md->bdim[i]; j++)
        //~ Rcout << md->var[i][j] << " ";
      //~ Rcout << "\n";
      //~}
  } else {
    Rcout << "VbStructure is NULL!\n";
    return S4();  
  }
  
  /*----------------------------------------------------------------*/
    /*---------------- Parse HMM chain -------------------------------*/
    /*----------------------------------------------------------------*/
    List HmmChain = HmmVb.slot("HmmChain");
  
  if (!Rf_isNull(HmmChain)){
    if (static_cast<int>(HmmChain.size()) != md->nb){
      Rcout << "number of Hmm models in HmmChain doesn't match with nb!\n";
      return S4();  
    }
    
    parseHmmChain(HmmChain, md);
    
  } else {
    Rcout << "HmmChain is NULL!\n";
    return S4();  
  }
  
  /*----------------------------------------------------------------*/
    /*----------------- Estimate Viterbi paths  ---------------------------------*/
    /*----------------------------------------------------------------*/
  ordervar(u,nseq,md->nb,md->bdim,md->var);
  
  int **optst;
  optst=(int **)calloc(nseq,sizeof(int *));
  for (int i=0;i<nseq;i++) optst[i]=(int *)calloc(md->nb,sizeof(int));
  
  
  #ifdef _OPENMP
  int num_threads = as<int>(nthread);
  omp_set_num_threads(num_threads);
  
  if (num_threads > 1)
    Rcout << "Number of threads used: " << num_threads << std::endl;
  #else
  int num_threads = as<int>(nthread);
  
  if (num_threads > 1)
    Rcout << "More than one thread is used with OpenMP not configured. Make sure your compiler supports OpenMP and reinstall the package" << std::endl;
  #endif
  
  #pragma omp parallel
  {
    double *meritbuf=(double *)calloc(md->maxnumst,sizeof(double));
    
    #pragma omp for
    for (int i=0;i<nseq;i++)
      viterbi(md,u[i],optst[i],NULL,meritbuf);
    
    free(meritbuf);
  }
  
  
  
  //---------------------------------------------------------//
    //----------- Find unique Viterbi sequences ----------------//
    //---------------------------------------------------------//
    
  int *clsid, *ct;
  
  clsid=(int *)calloc(nseq,sizeof(int));
  
  int **newoptst, newnseq, *vseqid;
  
  FindDifSeq(optst, nseq, md->nb, &newoptst, &newnseq, clsid);
  
  vseqid=(int *)calloc(newnseq,sizeof(int));
  
  for (int i=0; i<newnseq; i++) vseqid[i] = i;
  
  
  //=============================//
  // Compute modes               //
  //=============================//
  CompMode *cpm;
  cpm=(CompMode *)calloc(newnseq,sizeof(CompMode));
  
  double **mode = (double **)calloc(newnseq,sizeof(double *));
  
  for (int i=0;i<newnseq;i++) {
    SetCompMode(cpm+i,newoptst[i], md);
    cpm[i].logpdf=bwmem(md, cpm[i].mu, cpm[i].mode);
    mode[i]=cpm[i].mode;
  }

  double *sigmadat;
  sigmadat=(double *)calloc(dim,sizeof(double));


  //The overall deviation
  DataSigma(u,sigmadat,dim,nseq);

  //--------------------------------------------------------------------//
  //Second round alignment with reference based on newly computed modes //
  //--------------------------------------------------------------------//
    
  List newClust = wrapClust(mode, sigmadat, md->nb, dim, 
                            newnseq, newnseq, newoptst, vseqid);
  
    
    
  ct=(int *)calloc(newnseq,sizeof(int));
  for (int i=0;i<newnseq;i++) ct[i]=0;
  for (int i=0;i<nseq;i++) ct[clsid[i]]++;
  
  S4 HmmVbClust("HMMVBclust");
  
  HmmVbClust.slot("clustParam") = newClust;
  HmmVbClust.slot("clsid") = IntegerVector(clsid, clsid + nseq);
  HmmVbClust.slot("size") = IntegerVector(ct, ct + newnseq);
  
  // free memory
  for (int i=0; i<newnseq; i++) free(newoptst[i]);
  free(newoptst);
  
  free(vseqid);
  free(clsid);
  free(ct);
  free(sigmadat);
  
  free(mode);
  for (int i=0; i<newnseq; i++) freeCompMode(cpm+i);
  free(cpm);
  
  freeccm(&md);
  free(md);  
  free(u);
    

  
  return HmmVbClust;
  }
  
  
  
  
