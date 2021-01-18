#include <Rcpp.h>
#include "hmm.h"

using namespace Rcpp;

List wrapClust(double **mode, double *sigma, int nb, int dim, 
               int ncls, int ndseq, int **path, int *seqcls);

void freeClust(double **mode, double *sigma,
               int ncls, int ndseq, int **path, int *cls);

void parseVbStructure(const S4 &VbStructure, CondChain *md);
  
void parseHmmChain(const List &HmmChain, CondChain *md);
