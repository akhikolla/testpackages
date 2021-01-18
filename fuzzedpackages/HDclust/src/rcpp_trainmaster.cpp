#include <Rcpp.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "hmm.h"


// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
S4 rcpp_trainHmmVb(NumericMatrix dataTranspose, const RObject &VbStructure,
                       const List &searchControl, const List &trainControl, IntegerVector nthread,
                       Function VB, Function HMM, Function HMMVB, bool bprint = true) {
  
  ///////////////////
  // parse data
  ///////////////////
  int nseq = dataTranspose.ncol(); // num datapoints
  int dim = dataTranspose.nrow(); // data dimensionality

  // print parsed data
  //~ Rcout << "nseq = " << nseq << "\n";
  //~ Rcout << "dim = " << dim << "\n";

  double **u; // 2d representation of the data. Each row - new datapoint

  u = (double **)calloc(nseq,sizeof(double *));

  for (int i = 0; i < nseq; i++){
    u[i] = dataTranspose.begin() + i*dim;
  }
  
  ///////////////////////////////
  // parse training parameters
  ///////////////////////////////
  int ninit0 = trainControl["ninit0"];
  int ninit1 = trainControl["ninit1"];
  int ninit2 = trainControl["ninit2"];

  float epsilon = trainControl["epsilon"];

  DIAGCOV = as<int>(trainControl["diagCov"]);

  //~ Rcout << "ninit0 = " << ninit0 << "\n";
  //~ Rcout << "ninit1 = " << ninit1 << "\n";
  //~ Rcout << "ninit2 = " << ninit2 << "\n";
  //~ Rcout << "epsilon = " << epsilon << "\n";
  //~ Rcout << "DIAGCOV = " << DIAGCOV << "\n";


  // parse search parameters
  int *Numst0=NULL; // number of states per dimension
  int maxDim;
  int minDim;
  int relaxsearch;

  int nperm=1; // number of permutations for searching variable block structure inputted with -n parameter or 1 if non were supplied
  int nperm0=0; // number of permutations for searching variable block structure in file with permutations
  int **vlist0=NULL; // array that contains permutations used for searching variable block structure

  int nb; // num variable blocks
  int *bdim; // dimensions of variable blocks
  int **var; // variable id in variable blocks
  int *numst; // number of states in blocks

  /////////////////////////////////////////////////
  // Parse variable block structure if provided
  // or search parameters for variable block
  // structure otherwise
  /////////////////////////////////////////////////
  
  IntegerVector bstrBDim;
  IntegerVector bstrNumStates;
  IntegerVector numstatesPerDim;
  
  std::vector<std::vector<int> > bstrVarOrder;
  std::vector<std::vector<int> > vecPerm;
  
  if (!Rf_isNull(VbStructure)){
    nb = VbStructure.slot("nb");
    bstrBDim = VbStructure.slot("bdim");
    bstrNumStates = VbStructure.slot("numst");
    bdim = bstrBDim.begin();
    numst = bstrNumStates.begin();

	List lvar = VbStructure.slot("varorder");

    bstrVarOrder = as<std::vector<std::vector<int> > > (lvar);


    var = (int **)calloc(nb,sizeof(int *));

    for (int i=0; i<nb; i++){
	  for (std::vector<int>::iterator it = bstrVarOrder[i].begin(); it != bstrVarOrder[i].end(); it++) *it -= 1;
      var[i] = &(*(bstrVarOrder[i].begin()));
	}
	
    //~ Rcout << "nb = " << nb << "\n";
//~ 
    //~ Rcout << "bdim:\n";
    //~ for (int i = 0; i < nb; i++)
      //~ Rcout << bdim[i] << " ";
    //~ Rcout << "\n";
//~ 
    //~ Rcout << "numst:\n";
    //~ for (int i = 0; i < nb; i++)
      //~ Rcout << numst[i] << " ";
    //~ Rcout << "\n";
//~ 
    //~ Rcout << "var:\n";
    //~ for (int i = 0; i < nb; i++){
      //~ for (int j = 0; j < bdim[i]; j++)
        //~ Rcout << var[i][j] << " ";
      //~ Rcout << "\n";
    //~ }


  } else{
    maxDim = searchControl["maxDim"];
    minDim = searchControl["minDim"];
    relaxsearch = as<int>(searchControl["relax"]);

    //~ Rcout << "maxDim = " << maxDim << "\n";
    //~ Rcout << "minDim = " << minDim << "\n";
    //~ Rcout << "relaxsearch = " << relaxsearch << "\n";

	// read states per dimension array if provided
    if (!Rf_isNull(searchControl["numstPerDim"])){
      numstatesPerDim = searchControl["numstPerDim"];
      Numst0 = numstatesPerDim.begin();

      //~ Rcout << "Numst0:\n";
      //~ for (int i = 0; i < maxDim; i++)
        //~ Rcout << Numst0[i] << " ";
      //~ Rcout << "\n";
    }


	// read permutations if provided
    if (!Rf_isNull(searchControl["perm"])){
	  List lperm = searchControl["perm"];
      
      vecPerm = as<std::vector<std::vector<int> > >(lperm);

      nperm0 = static_cast<int>(vecPerm.size());
      vlist0=(int **)calloc(nperm0,sizeof(int *));

      for (int i=0; i<nperm0; i++){
		for (std::vector<int>::iterator it = vecPerm[i].begin(); it != vecPerm[i].end(); it++) *it -= 1;
        vlist0[i] = &(*(vecPerm[i].begin()));
	  }
	  
      //~ Rcout << "nperm0 = " << nperm0 << "\n";
//~ 
      //~ Rcout << "vlist0:\n";
      //~ for (int i = 0; i < nperm0; i++){
        //~ for (int j = 0; j < dim; j++)
          //~ Rcout << vlist0[i][j] << " ";
        //~ Rcout << "\n";
      //~ }
    }
    else
      nperm = searchControl["nperm"];
    
    //~ Rcout << "nperm = " << nperm << "\n";
    
  }
  
   //~ Rcout << "var:\n";
    //~ for (int i = 0; i < nb; i++){
      //~ for (int j = 0; j < bdim[i]; j++)
        //~ Rcout << var[i][j] << " ";
      //~ Rcout << "\n";
    //~ }
  
  ///////////////////////////////////////////////////////////////
  // if variable block structure was provided, train HMM
  // if not, search for variable block structure first and then
  // train HMM
  ///////////////////////////////////////////////////////////////
  double *wt=NULL; // array with weights for sample points
  CondChain *md=NULL;
  double *loglikehd, lhsum;
  int randomseed=0;

  loglikehd=(double *)calloc(nseq,sizeof(double));
  
  
  #ifdef _OPENMP
  int num_threads = as<int>(nthread);
  omp_set_num_threads(num_threads);
  
  if ((bprint) && (num_threads > 1))
	Rcout << "\nNumber of threads used: " << num_threads << std::endl;
  #else
  int num_threads = as<int>(nthread);
  
  if ((bprint)&&(num_threads > 1))
	Rcout << "\nMore than one thread is used with OpenMP not configured. Make sure your compiler supports OpenMP and reinstall the package" << std::endl;
  #endif

  
  if (!Rf_isNull(VbStructure)) {//variable block structure specified
    hmmfit_minit(u, nseq, nb, bdim, var, numst, &md, loglikehd, &lhsum, (double)epsilon, wt,
                 ninit0, ninit1, ninit2, randomseed); //lhsum is loglikelihood
    lhsum-=(double)(computenp(nb, bdim,numst))*log((double)nseq)*0.5; //BIC
  }
  else {//variable block structure not given and will be searched
    hmmfit_vb(u, nseq, dim, &nb, &bdim, &var, nperm, nperm0, vlist0,
              &md, loglikehd, &lhsum, (double)epsilon, wt, ninit0, ninit1,ninit2, randomseed,
              Numst0, maxDim, minDim, relaxsearch); //output lhsum is BIC not loglikelihood
  }
  
  
  // now wrap results to HMMVB object
  // first wrap HmmChain - list of HMM models
  List HmmChain(nb);
  
  for (int i = 0; i < nb; i++){
    
    NumericVector va00 = NumericVector(md->mds[i]->a00, md->mds[i]->a00 + md->mds[i]->numst);
    
    NumericMatrix a(md->mds[i]->prenumst, md->mds[i]->numst);
    
    for (int j = 0; j < md->mds[i]->prenumst; j++)
      std::copy(md->mds[i]->a[j], md->mds[i]->a[j] + md->mds[i]->numst, a(j,_).begin());
    
    NumericMatrix mean(md->mds[i]->numst, md->mds[i]->dim);
    
    List sigmaList(md->mds[i]->numst);
    List sigmaInvList(md->mds[i]->numst);
    
    NumericVector sigmaDetLog(md->mds[i]->numst);
    
    for (int j = 0; j < md->mds[i]->numst; j++){
	  sigmaDetLog[j] = md->mds[i]->stpdf[j]->sigma_det_log;
		
      std::copy(md->mds[i]->stpdf[j]->mean, md->mds[i]->stpdf[j]->mean + md->mds[i]->dim, mean(j,_).begin());
      
      NumericMatrix sigma(md->mds[i]->dim, md->mds[i]->dim);
      NumericMatrix sigmaInv(md->mds[i]->dim, md->mds[i]->dim);
      
      for (int k = 0; k < md->mds[i]->dim; k++){
        std::copy(md->mds[i]->stpdf[j]->sigma[k], md->mds[i]->stpdf[j]->sigma[k] + md->mds[i]->dim, sigma(k,_).begin());
         std::copy(md->mds[i]->stpdf[j]->sigma_inv[k], md->mds[i]->stpdf[j]->sigma_inv[k] + md->mds[i]->dim, sigmaInv(k,_).begin());
      }
      sigmaList[j] = sigma;    
      sigmaInvList[j] = sigmaInv;    
    }
    
   
    
    S4 Hmm = HMM(Named("dim", md->mds[i]->dim), Named("numst", md->mds[i]->numst), Named("prenumst", md->mds[i]->prenumst),
				 Named("a00", va00), Named("a", a), Named("mean", mean), Named("sigma", sigmaList), Named("sigmaInv", sigmaInvList),
				 Named("sigmaDetLog",sigmaDetLog));
    
    HmmChain[i] = Hmm;
  }
  
  
  // now wrap variable block structure.
  // if it was provided in the input, simply return it 
  S4 NewVbStructure;
  
  if (Rf_isNull(VbStructure)){ 
    NumericVector vbdim = NumericVector(bdim, bdim + nb);
    
    NumericVector vnumst(nb);
    
    for (int i = 0; i < nb; i++)
		vnumst[i] = md->mds[i]->numst;
    
    List varorderList(nb);
    
    for (int i = 0; i < nb; i++){
	  for (int *ptr = var[i]; ptr != var[i]+bdim[i]; ptr++) *ptr += 1;
      NumericVector varorder(var[i], var[i]+bdim[i]);
      varorderList[i] = varorder;
    }
    
	
	NewVbStructure = VB(Named("nb", nb), Named("dim", dim), Named("bdim", vbdim),
						Named("numst", vnumst), Named("varorder", varorderList));
 
  }
  else
    NewVbStructure = VbStructure;

  // output loglikelihood
  NumericVector loglh(nseq);
  for (int i=0;i<nseq;i++){
      loglh[i]=loglikehd[i];
  }
    
  // free memory  
  free(u);
  free(loglikehd);
  freeccm(&md);
  free(md);
  
  if (!Rf_isNull(searchControl["perm"])) free(vlist0);
  
  if (Rf_isNull(VbStructure)){
  	free(bdim);
  	
  	for (int i = 0; i < nb; i++)
  		free(var[i]);
  }
  
  free(var);
  
  
  S4 HmmVb = HMMVB(Named("HmmChain", HmmChain), Named("VbStructure", NewVbStructure), Named("BIC", -2*lhsum), Named("diagCov", LogicalVector(wrap(DIAGCOV))), Named("Loglikehd",loglh));
  
  
  return(HmmVb);
}


