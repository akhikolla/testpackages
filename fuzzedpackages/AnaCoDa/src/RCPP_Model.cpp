#ifndef STANDALONE
#include "include/ROC/ROCModel.h"
#include "include/PA/PAModel.h"
#include "include/FONSE/FONSEModel.h"
#include "include/PANSE/PANSEModel.h"
#include <Rcpp.h>
using namespace Rcpp;

/*
void roc_finalizer(ROCModel* m)
{
	delete m;
}
*/

//' @name simulateGenome
//' @title simulateGenome
//' @description Method of Model class (access via model$<function name>, where model is an object initialized by initializeModelObject). Will simulate a version of the given genome using the current set of parameters stored in the Parameter object. This can be written to a FASTA file using genome$writeFasta(<filename>,simulated = TRUE).
//' @param genome a Genome object initialized by initializeGenomeObject

RCPP_EXPOSED_CLASS(ROCParameter)
RCPP_EXPOSED_CLASS(PAParameter)
RCPP_EXPOSED_CLASS(PANSEParameter)
RCPP_EXPOSED_CLASS(FONSEParameter)
RCPP_EXPOSED_CLASS(Parameter)
RCPP_EXPOSED_CLASS(Genome)
RCPP_MODULE(Model_mod)
{
	class_<Model>("Model")
		;

	class_<ROCModel>( "ROCModel" )
		.derives<Model>("Model")
		.constructor<bool, bool>()
		//.finalizer(&roc_finalizer)
		.method("CalculateProbabilitiesForCodons", &ROCModel::CalculateProbabilitiesForCodons,
		        "Calculated codon probabilities. Input is one element shorter than output")
  		.method("getParameter", &ROCModel::getParameter)
  		.method("setParameter", &ROCModel::setParameter)
  		.method("simulateGenome", &ROCModel::simulateGenome)
		;
	
	class_<PAModel>("PAModel")
		.derives<Model>("Model")
		.constructor<unsigned,bool, bool>()
  		.method("getParameter", &PAModel::getParameter)
		.method("setParameter", &PAModel::setParameter)
		.method("simulateGenome", &PAModel::simulateGenome)
		;
	
    class_<PANSEModel>("PANSEModel")
		.derives<Model>("Model")
		.constructor<unsigned,bool, bool>()
  		.method("getParameter", &PANSEModel::getParameter)
		.method("setParameter", &PANSEModel::setParameter)
		.method("simulateGenome", &PANSEModel::simulateGenome)
		;

	class_<FONSEModel>("FONSEModel")
		.derives<Model>("Model")
		.constructor<bool, bool>()
		.method("CalculateProbabilitiesForCodons", &FONSEModel::CalculateProbabilitiesForCodons,
		        "Calculated codon probabilities. Input is one element shorter than output")
  		.method("getParameter", &FONSEModel::getParameter)
		.method("setParameter", &FONSEModel::setParameter)
		.method("simulateGenome", &FONSEModel::simulateGenome)
		;
}
#endif
