#include "rcpp_flan_module.h"

RCPP_MODULE(flan_module) {
  class_<FLAN_MutationModel>("FlanMutMod")
	.constructor<List>()
	.method("pflan",&FLAN_MutationModel::computeCumulativeFunction,"compute cumulative function")
	.method("dflan",&FLAN_MutationModel::computeProbability,"compute probability")
	.method("dflanda",&FLAN_MutationModel::computeProbability1DerivativeAlpha,"compute probability derivative wtt alpha")
	.method("dflandr",&FLAN_MutationModel::computeProbability1DerivativeRho,"compute probability derivative wtt rho")
	.method("dflangrad",&FLAN_MutationModel::computeProbability1DerivativesAlphaRho,"compute probability derivative wtt alpha and rho")
	.method("deduce.dflan",&FLAN_MutationModel::deduceProbability,"compute probability")
	.method("deduce.dflanda",&FLAN_MutationModel::deduceProbability1DerivativeAlpha,"compute probability derivative wtt alpha")
	.method("MutationGFEstimation",&FLAN_MutationModel::MutationGFEstimation,"estimate alpha with GF method")
	.method("CovGFEstimation",&FLAN_MutationModel::covariance,"standard deviation of GF method")
	.method("unbias.mutprob",&FLAN_MutationModel::unbiasPiEstimation,"unbias mutprob estimation")
  ;

  class_<FLAN_Sim>("FlanSim")
	.constructor<List>()
	.method("rflan",&FLAN_Sim::computeSamplesMutantsFinalsNumber,"compute sample mutants")
  ;


  class_<FLAN_SimInhomogeneous>("FlanIhSim")
	.constructor<List>()
	.method("rflan",&FLAN_SimInhomogeneous::computeSamplesMutantsFinalsNumber,"compute sample mutants")
  ;



  // class_<FLAN_SimClone>("FlanSimClone")
	// .constructor<double,double,List>()
	// .method("rclone",&FLAN_SimClone::computeSample,"compute Clone")
  // ;

  // class_<FLAN_SimInhomogeneousClone>("FlanSimInhomogeneousClone")
	// .constructor<double,double,List>()
	// .method("rcloneih",&FLAN_SimInhomogeneousClone::computeSample,"compute Inhomogeneous Clone")
  // ;


  class_<FLAN_ExponentialClone>("FlanExpClone")
  .constructor<List>()
	.method("dclone",&FLAN_ExponentialClone::computeProbability,"compute probability")
  // .method("integral",&FLAN_ExponentialClone::testIntegral,"compute integrale")
	.method("dclonedr",&FLAN_ExponentialClone::computeProbability1DerivativeRho,"compute probability")
	.method("pgf2",&FLAN_ExponentialClone::computeGeneratingFunction2,"compute generating function for several z")
  ;

  class_<FLAN_DiracClone>("FlanDirClone")
  .constructor<List>()
	.method("dclone",&FLAN_DiracClone::computeProbability,"compute probability")
	.method("dclonedr",&FLAN_DiracClone::computeProbability1DerivativeRho,"compute probability")
	.method("pgf2",&FLAN_DiracClone::computeGeneratingFunction2,"compute generating function for several z")
  ;


  class_<FLAN_InhomogeneousClone>("FlanInhClone")
  .constructor<List>()
	.method("dclone",&FLAN_InhomogeneousClone::computeProbability,"compute probability")
	.method("dclonedr",&FLAN_InhomogeneousClone::computeProbability1DerivativeRho,"compute probability")
// 	.method("pgf",&FLAN_InhomogeneousClone::computeGeneratingFunction,"compute generating function for several z")
	.method("pgf2",&FLAN_InhomogeneousClone::computeGeneratingFunction2,"compute generating function according to fitness")
	.method("pgfdr",&FLAN_InhomogeneousClone::computeGeneratingFunction1DerivativeRho,"compute generating function according to fitness")
  ;


  //
  // class_<MATH_Integration>("MathInt")
  // .constructor<List>
  //
}
