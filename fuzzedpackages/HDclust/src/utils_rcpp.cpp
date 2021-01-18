#include "utils_rcpp.h"

List wrapClust(double **mode, double *sigma, int nb, int dim, 
               int ncls, int ndseq, int **path, int *seqcls)
{
  NumericMatrix vmode(ncls, dim);
  
  for (int i = 0; i < ncls; i++)
    std::copy(mode[i], mode[i] + dim, vmode(i,_).begin());
  
  List vseq(ndseq);
  
  for (int i = 0; i < ndseq; i++)
    vseq[i] = IntegerVector(path[i], path[i] + nb);
  
  
  return List::create(Named("ncls") = ncls,
                      Named("ndseq") = ndseq,
                      Named("vseqid") = IntegerVector(seqcls, seqcls + ndseq),
                      Named("sigma") = NumericVector(sigma, sigma + dim),
                      Named("mode") = vmode,
                      Named("vseq") = vseq); 	
}

void freeClust(double **mode, double *sigma,
               int ncls, int ndseq, int **path, int *cls)
{
  free(cls);
  free(sigma);
  
  for (int i = 0; i < ncls; i++) free(mode[i]);
  free(mode);
  
  for (int i = 0; i < ndseq; i++) free(path[i]);
  free(path);
}


void parseVbStructure(const S4 &VbStructure, CondChain *md)
{	  
  md->nb = VbStructure.slot("nb");
  
  md->bdim = (int *)calloc(md->nb,sizeof(int));  
  md->numst = (int *)calloc(md->nb,sizeof(int));
  
  IntegerVector bstrBDim = VbStructure.slot("bdim");
  IntegerVector bstrNumStates = VbStructure.slot("numst");
  
  std::copy(bstrBDim.begin(), bstrBDim.end(), md->bdim);
  std::copy(bstrNumStates.begin(), bstrNumStates.end(), md->numst);
  
  int dim = 0;
  for (int i = 0; i < md->nb; i++) dim += md->bdim[i];
  md->dim = dim;
  
  int maxnumst = 0;
  for (int i = 0; i < md->nb; i++)
    if (md->numst[i] > maxnumst) maxnumst = md->numst[i];
    
    md->maxnumst = maxnumst;
    
    std::vector<IntegerVector> bstrVarOrder = VbStructure.slot("varorder");
    
    md->var = (int **)calloc(md->nb,sizeof(int *));
    for (int i=0; i < md->nb; i++){
      md->var[i] = (int *)calloc(md->bdim[i],sizeof(int));
      std::copy(bstrVarOrder[i].begin(), bstrVarOrder[i].end(), md->var[i]);
      
      for (int j = 0; j < md->bdim[i]; j++)
        md->var[i][j]--;
    }
    
    // Set up the redundant fields for convenience of computation
    md->cbdim=(int *)calloc(md->nb,sizeof(int));
    md->cbdim[0]=0;
    md->cnumst=(int *)calloc(md->nb,sizeof(int));
    md->cnumst[0]=0;
    for (int i = 0; i < md->nb-1; i++) {
      md->cbdim[i+1] = md->cbdim[i] + md->bdim[i];
      md->cnumst[i+1] = md->cnumst[i] + md->numst[i];
    }
}

void parseHmmChain(const List &HmmChain, CondChain *md)
{
  md->mds=(HmmModel **)calloc(md->nb,sizeof(HmmModel *));
  for (int i = 0; i < md->nb; i++) {
    md->mds[i]=(HmmModel *)calloc(1,sizeof(HmmModel));
    
    S4 Hmm = HmmChain[i];
    
    md->mds[i]->dim = Hmm.slot("dim");
    md->mds[i]->numst = Hmm.slot("numst");
    md->mds[i]->prenumst = Hmm.slot("prenumst");
    
    NumericVector a00 = Hmm.slot("a00");
    
    md->mds[i]->a00 = (double *)calloc(md->mds[i]->numst,sizeof(double));
    std::copy(a00.begin(), a00.end(), md->mds[i]->a00);
    
    NumericMatrix a = Hmm.slot("a");
    
    md->mds[i]->a = (double **)calloc(md->mds[i]->prenumst,sizeof(double *));
    
    for (int j = 0; j < md->mds[i]->prenumst; j++){
      md->mds[i]->a[j] = (double *)calloc(md->mds[i]->numst,sizeof(double));
      std::copy(a(j,_).begin(), a(j,_).end(), md->mds[i]->a[j]);
    }
    
    NumericMatrix mean = Hmm.slot("mean");
    std::vector<NumericMatrix> sigmaList = Hmm.slot("sigma");
    std::vector<NumericMatrix> sigmaInvList = Hmm.slot("sigmaInv");
    NumericVector sigmaDetLog = Hmm.slot("sigmaDetLog");
    
    
    md->mds[i]->stpdf=(GaussModel **)calloc(md->mds[i]->numst, sizeof(GaussModel *));
    
    for (int j = 0; j < md->mds[i]->numst; j++){
      md->mds[i]->stpdf[j]=(GaussModel *)calloc(1, sizeof(GaussModel));
      
      md->mds[i]->stpdf[j]->exist = 1;
      md->mds[i]->stpdf[j]->dim = Hmm.slot("dim");
      
      md->mds[i]->stpdf[j]->mean = (double *)calloc(md->mds[i]->stpdf[j]->dim,sizeof(double));
      std::copy(mean(j,_).begin(), mean(j,_).end(), md->mds[i]->stpdf[j]->mean);
      
      md->mds[i]->stpdf[j]->sigma_det_log = sigmaDetLog[j];
      
      md->mds[i]->stpdf[j]->sigma = (double **)calloc(md->mds[i]->stpdf[j]->dim,sizeof(double *));
      
      for (int k =0; k<md->mds[i]->stpdf[j]->dim; k++){
        md->mds[i]->stpdf[j]->sigma[k] = (double *)calloc(md->mds[i]->stpdf[j]->dim,sizeof(double));
        std::copy((sigmaList[j])(k,_).begin(), (sigmaList[j])(k,_).end(), md->mds[i]->stpdf[j]->sigma[k]);
      }
      
      md->mds[i]->stpdf[j]->sigma_inv = (double **)calloc(md->mds[i]->stpdf[j]->dim,sizeof(double *));
      
      for (int k = 0; k < md->mds[i]->stpdf[j]->dim; k++){
        md->mds[i]->stpdf[j]->sigma_inv[k] = (double *)calloc(md->mds[i]->stpdf[j]->dim,sizeof(double));
        std::copy((sigmaInvList[j])(k,_).begin(), (sigmaInvList[j])(k,_).end(), md->mds[i]->stpdf[j]->sigma_inv[k]);
      }
    }
  }    
}
