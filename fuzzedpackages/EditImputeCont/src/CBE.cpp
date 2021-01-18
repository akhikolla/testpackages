/* Version updates

EditImputCont_1.1.2 
- remove sparselib.c and sparselib.c

EditImputCont_1.1.3
- remove unncessary libraries from lp_solve, Newmat and Newran
	 newmat9.cpp, sparselib.c, sparselib.h, ufortify.h
*/


#include "CHeader.h"  
#include <Rcpp.h>
#include <R_ext/Utils.h>
#include "CData.h" 
#include "CFeasibilityMap.h"
#include "CParam.h"
#include "CBE.h"

////.constructor<Rcpp::NumericMatrix, Rcpp::NumericMatrix, Rcpp::NumericMatrix,Rcpp::NumericMatrix, Rcpp::NumericMatrix>()
RCPP_MODULE(cbei){
  using namespace Rcpp;
  using namespace R;
	
	class_<CBE>( "bei" )
    // expose the default constructor
    .constructor<Rcpp::NumericMatrix, Rcpp::NumericMatrix, Rcpp::NumericMatrix, int, int>()
    .method("InitializeSandD", &CBE::InitializeSandD, "Generate initial S and D matrix" )
    // .method("SetInitialSandD", &CBE::SetInitialSandD, "Set initial S and D matrix" )
    // .method("SetTrueS", &CBE::SetTrueS, "Set true S matrix" )
    .method("Iterate", &CBE::Iterate, "Run one iteration of MCMC algorithm" )
    .method("Run", &CBE::Run, "Run MCMC algorithm for given times" )
    // .method("Test", &CBE::test, "bug catcher" )
    .method("BuildFeasibilityMap", &CBE::BuildMap, "Build Feasibility Map" )
    //.method("SetOptionalData", &CBE::SetOptionalData, "Set Optional Data" )
    .method("Simulate_logUnif_case2", &CBE::Simulate_logUnif_case2, "Simulate logUnif values for case 2" )
    .property("RandomSeed", &CBE::GetRandomSeed, &CBE::SetRandomSeed, "Random Seed")
    .property("msg.level", &CBE::GetMsgLevel, &CBE::SetMsgLevel, "Message Level")
    // .property("useMap", &CBE::GetUseMap, &CBE::SetUseMap, "Use Map")
		.property("S.Compact.Initial", &CBE::Get_initialS_Compact, "Initial S matrix")
		// .property("Log.Uniform.Compact", &CBE::Get_logUnif_Compact, &CBE::SetS_logUnif_Compact, "log uniform for case 2 matrix")
		.property("D.Initial", &CBE::GetD,"D matrix")
		.method("K", &CBE::GetK, "Number of components")
    // .property("K", &CBE::GetK, &CBE::SetK, "Number of components")
    // .property("FeasibilityMap", &CBE::Get_FeasibilityMap, &CBE::Set_FeasibilityMap, "Feasibility map matrix")
		// .property( "ImputedX.Compact", &CBE::GetX_compact, "Retrieve compact version of the imputed data matrix" )
    .property( "Y.edited", &CBE::GetX, "Retrieve complete version of the imputed data matrix" )
		// .property( "S.Compact", &CBE::GetS_compact, "S compact matrix" )
    // .property( "alpha", &CBE::GetAlpha, "alpha" )
		// .property( "Phi", &CBE::GetPhi, "Phi" )
    // .property( "Sigma", &CBE::GetSigma, "Sigma" )
    // .property( "Mu", &CBE::GetMu, "Mu" )
    // .property( "Pi", &CBE::GetPi, "Pi" )
    // .property( "Z", &CBE::GetZ, "Z" )
		.property( "n.occ", &CBE::GetNZ, "n.occ" )
    // .property( "Accept", &CBE::GetAccept, "is_accept" )
    .property( "Prob.A", &CBE::GetProb_A, "Prob_A" )
		.property( "FaultyRecordID", &CBE::GetFaultyIndex, "Faulty record ID's" )
    // .property( "DrawFindS", &CBE::GetDrawFindS, "Draw Find S" )
    .property( "Y.input", &CBE::GetDobs, "Y.input" )
  ;
}       

CBE::CBE(Rcpp::NumericMatrix X_, Rcpp::NumericMatrix Edit_, Rcpp::NumericMatrix logB_, int nbalance_edit, int K) {
  Matrix D_obs = Rcpp2Mat2(X_);
  Matrix Edit = Rcpp2Mat2(Edit_);
  Matrix logB = Rcpp2Mat2(logB_);
  
  Data.SetData(D_obs,Edit,logB,nbalance_edit);
  Data.init();
  
  //srand (time(NULL));
  GetRNGstate();
  randseed = (float)unif_rand(); //to use double only later once verified
  //randseed = (float)((double)rand() / RAND_MAX); //to use double only later once verified
  PutRNGstate();
  urng = new MotherOfAll(randseed);
  Random::Set(*urng);
  randUnif = new Uniform;
  // hyper.init(5.0,0.25,0.25,0.25,0.25,(Data.n_var-Data.n_balance_edit)+1);
  hyper.init(1.0,0.25,0.25,0.25,0.25,(Data.n_var-Data.n_balance_edit)+1); // CHANGED 2/2/2015
  
  Param.K = -K; //do delayed inilization of componnts.
  
  is_S_and_D_initilized = false;
  is_feasiblity_Map_initilized = false;;
  is_log_unif_initilized = false;
  IterCount = 0;
	useMap = true;
  // useMap = false; // memory error if this is false when using large dataset
}
/*
CBE::CBE(Rcpp::NumericMatrix X_, Rcpp::NumericMatrix Edit_, Rcpp::NumericMatrix logB_,
  Rcpp::NumericMatrix S_, Rcpp::NumericMatrix D_) {

  Matrix D_obs = Rcpp2Mat(X_);
  Matrix Edit = Rcpp2Mat(Edit_);
  Matrix logB = Rcpp2Mat(logB_);
  Data.SetData(D_obs,Edit,logB,2);
  Data.init();
  
  
  Matrix temp = Rcpp2Mat(D_); 
  Data.D_initial = Data.D_Observed; 
  for (int i = 1; i <= Data.n_faulty; i++ ) {
    Data.D_initial.row(Data.Faulty2Original[i-1]) = temp.row(i);
  }
  temp = Rcpp2Mat(S_); 
  Data.initial_S_Mat = temp;
  
  srand (time(NULL));
  //randseed = (float)(rand() / RAND_MAX);
  randseed = 0.101;
  
  urng = new MotherOfAll(randseed);
  Random::Set(*urng);
  randUnif = new Uniform;
  
  hyper.init(5.0,0.25,0.25,0.25,0.25,(Data.n_var-Data.n_balance_edit)+1);
  
  FM.Build(Data);
  
  int n_simul = 100;  // To calculate logUnif_case2
  FM.Simulate_logUnif_case2(n_simul, *randUnif, Data);
  Param.initialize(Data,50);
  cout << "Param.Prob_A = " << Param.Prob_A << endl;
  for (int iter=1; iter<=55; iter++) {
     cout<< "iter: " <<  iter << endl;
    Param.iterate(iter, Data, FM, hyper, *randUnif);
  }
}
*/
//Destructor
CBE::~CBE(){
  delete urng;
}

void CBE::SetOptionalData(Rcpp::List OData) { //testing code
  //std::string method   = Rcpp::as<std::string>(OData["method"]);
  //double tolerance     = Rcpp::as<double>(OData["tolerance"]);
}

// Commented on 05/21/2015
// void CBE::test() {
//   FM.test(Data);
// }

void CBE::SetUseMap(bool use) {
  useMap = use;
}
bool CBE::GetUseMap() {
  return useMap;
}
void CBE::BuildMap() {
  if (is_log_unif_initilized) {
	Rprintf( "BuildFeasibilityMap must be called before Simulate_logUnif_case2.\n");
    return;
  }
	
	if ( Param.msg_level > 0 )	{
		Rprintf( "useMap\n");
	}
  FM.useMap = useMap;
	if ( Param.msg_level > 0 )	{
		Rprintf( "FM.Build\n");
	}
  FM.Build(Data);
	if ( Param.msg_level > 0)	{
		Rprintf( "end.FM.Build\n");
	}
  is_feasiblity_Map_initilized = true;
}

void CBE::Iterate() {
  if (!is_log_unif_initilized) {
    Simulate_logUnif_case2(100);
  }
  if (is_log_unif_initilized) {
    IterCount++;
    Param.iterate(IterCount, Data, FM, hyper, *randUnif,100);
  } else {
	Rprintf( "Model was not initilized properly\n");
  }
}

void CBE::Run(int iter) {
	
  if (!is_log_unif_initilized) {
    Simulate_logUnif_case2(100);
  }

  if (is_log_unif_initilized) {
    for (int i = 1; i <= iter; i++) {
      IterCount++;
      Param.iterate(IterCount, Data, FM, hyper, *randUnif,100);
    }
  } else {
	Rprintf( "Model was not initilized properly\n");
  }

}

void CBE::SetK(int ncomp) {
  if (Param.K == -1) { 
    Param.K = -ncomp; //record the number but delay the inilization 
  } else {
    if (Param.K != ncomp) {
	  Rprintf( "re-initilizing number of components to %d\n", ncomp);
      Param.initialize(Data,ncomp,FM, *randUnif,100);
    }
  }
}

int CBE::GetK() {
  return Param.K < 0?-Param.K:Param.K;
}

void CBE::SetMsgLevel(int level) {
  Param.msg_level = level;
}
int CBE::GetMsgLevel() {
  return Param.msg_level;
}
    

void CBE::SetS_logUnif_Compact(Rcpp::NumericMatrix X_) { //merge code with Simulate_logUnif_case2 later
  if (!is_S_and_D_initilized) {
    InitializeSandD();
  }
  if (is_S_and_D_initilized) {
    if (!is_feasiblity_Map_initilized) {
      BuildMap();
    }
  } else {
    return;
  }
  if (is_feasiblity_Map_initilized) {
    Data.logUnif_case2 = Rcpp2Mat(X_);
  }
  is_log_unif_initilized = true;
}

void CBE::Simulate_logUnif_case2(int nsim) {
  if (!is_S_and_D_initilized) {
    InitializeSandD();
  }
  if (is_S_and_D_initilized) {
    if (!is_feasiblity_Map_initilized) {
      BuildMap();
    }
  } else {
    return;
  }
  if (is_feasiblity_Map_initilized) {
    FM.Simulate_logUnif_case2(nsim, *randUnif, Data);
  }
  is_log_unif_initilized = true;
}
Rcpp::NumericMatrix CBE::Get_FeasibilityMap() {
   return Mat2Rcpp(FM.feasibleMap);
}

void CBE::Set_FeasibilityMap(Rcpp::NumericMatrix X_) { //check dim later
  if (is_log_unif_initilized) {
	Rprintf( "FeasibilityMap must be set before Simulate_logUnif_case2\n");
    return;
  }
  FM.feasibleMap = Rcpp2Mat(X_);
  is_feasiblity_Map_initilized = true;
}

Rcpp::NumericMatrix CBE::Get_initialS_Compact() {
   return Mat2Rcpp(Data.initial_S_Mat);
}

Rcpp::NumericMatrix CBE::Get_logUnif_Compact() {
  return Mat2Rcpp(Data.logUnif_case2);
}

Rcpp::NumericMatrix CBE::GetD() {
   return Mat2Rcpp(Data.D_initial);
}
Rcpp::NumericVector CBE::GetZ() {
  return Mat2Rcpp(Param.z_in);
}
Rcpp::NumericVector CBE::GetNZ() {
  Matrix n_z =  Param.n_z.t();
  return Mat2Rcpp(n_z);
}
Rcpp::NumericVector CBE::GetAccept() {
  return Mat2Rcpp(Param.is_accept);
}
    
Rcpp::NumericMatrix  CBE::GetX() {
  Matrix X_MI(Data.n_sample,Data.n_var);
  for (int i_sample=1; i_sample<=Data.n_sample; i_sample++){
    for (int i_var=1; i_var<=Data.n_var; i_var++){
      X_MI(i_sample,i_var) = exp(Param.Y_in(i_sample,i_var));	
		}
	}
  return Mat2Rcpp(X_MI);
}
Rcpp::NumericMatrix CBE::GetS_compact() {
  Matrix mat = Param.S_Mat.t();
  return Mat2Rcpp(mat);
}
double CBE::GetAlpha() {
  return Param.alpha;
}
Rcpp::NumericMatrix CBE::GetX_compact() {
  Matrix X_MI(Data.n_faulty,Data.n_var);
  for (int i=1; i<=Data.n_faulty; i++){
    int i_sample = Data.Faulty2Original[i-1];
		for (int i_var=1; i_var<=Data.n_var; i_var++){
      X_MI(i,i_var) = exp(Param.Y_in(i_sample,i_var));	
		}
	}
  return Mat2Rcpp(X_MI);
}

Rcpp::NumericMatrix CBE::Mat2Rcpp(Matrix &X_) { //optimize it later
  int n_col = X_.ncols();
  int n_row = X_.nrows();
  Rcpp::NumericMatrix mat(n_row, n_col);
  for (int i = 1; i<=n_row; i++) {
    for (int j = 1; j <= n_col; j++) {
      mat(i-1,j-1) = X_(i,j); 
    }
  }
  return mat;
}

Matrix CBE::Rcpp2Mat(Rcpp::NumericMatrix X_) { //optimize it later
  int n_col = X_.ncol();
  int n_row = X_.nrow();
  Matrix mat = Matrix(n_row,n_col);
  for (int i = 1; i<=n_row; i++) {
    for (int j = 1; j <= n_col; j++) {
      mat(i,j) = (float)X_(i-1,j-1); 
    }
  }
  return mat;
}

Matrix CBE::Rcpp2Mat2(Rcpp::NumericMatrix X_) { //optimize it later
  int n_col = X_.ncol();
  int n_row = X_.nrow();
  Matrix mat = Matrix(n_row,n_col);
  for (int i = 1; i<=n_row; i++) {
    for (int j = 1; j <= n_col; j++) {
      mat(i,j) = (float)X_(i-1,j-1); 
    }
  }
  return mat;
}


double CBE::GetRandomSeed() {
  return randseed;  
}

void CBE::SetRandomSeed(double seed) {
  delete urng;
  delete randUnif;
  randseed = (float)seed;
  urng = new MotherOfAll(randseed);
  Random::Set(*urng);
  randUnif = new Uniform;
}

Rcpp::NumericVector CBE::GetMu() {
  Matrix mat = Param.Mu;
  return Mat2Rcpp(mat);
}
Rcpp::NumericVector CBE::GetPi() {
  Matrix mat = Param.pi.t();
  return Mat2Rcpp(mat);
}

Rcpp::NumericVector CBE::GetSigma() {
  Rcpp::Dimension d(Param.n_var_independent,Param.n_var_independent,Param.K);
  Rcpp::NumericVector Sigma3D(d); 
  int s2 = Param.n_var_independent * Param.n_var_independent;
  for (int k = 1; k <= Param.K; k++) {
    Matrix sigma = Param.SIGMA[k-1];
    for (int i = 1; i <= Param.n_var_independent; i++) {
      for (int j = 1; j <= Param.n_var_independent; j++) {
        Sigma3D[(k-1) * s2 + (i-1) * Param.n_var_independent + (j-1)] = sigma(i,j);
      }
    }
  }
  return Sigma3D;
}
      
Rcpp::NumericMatrix CBE::GetPhi() {
  Matrix Phi = Param.Phi;
  return Mat2Rcpp(Phi);
}
double CBE::GetProb_A() {
  return Param.Prob_A;
}

void CBE::SetInitialSandD(Rcpp::NumericMatrix S_, Rcpp::NumericMatrix D_) {
  if (is_log_unif_initilized) {
	Rprintf( "InitializeSandD must be called before Simulate_logUnif_case2\n");
    return;
  }
  Set_initialS(S_);
  SetD(D_);
  is_S_and_D_initilized = true;
}

void CBE::Set_initialS(Rcpp::NumericMatrix X_) { 
  Matrix temp = Rcpp2Mat(X_); 
  if (temp.nrows() == Data.D_Observed.nrows()) {
    Data.initial_S_Mat = Matrix(Data.n_faulty,temp.ncols());
    for (int i = 1; i <= Data.n_faulty; i++ ) {
      Data.initial_S_Mat.row(i) = temp.row(Data.Faulty2Original[i-1]);
    }  
  } else {    //compact matrix, assuming dim matching for now
    Data.initial_S_Mat = temp;
  }
}

void CBE::SetD(Rcpp::NumericMatrix X_) {
  //initilize intial imputed data for records
  Matrix temp = Rcpp2Mat(X_); //full matrix
  if (temp.nrows() == Data.D_Observed.nrows()) {
     Data.D_initial = temp;
  } else {    //compact matrix
    Data.D_initial = Data.D_Observed; 
    for (int i = 1; i <= Data.n_faulty; i++ ) {
      Data.D_initial.row(Data.Faulty2Original[i-1]) = temp.row(i);
    }  
  }
}

void CBE::InitializeSandD() {
  if (is_log_unif_initilized) {
	Rprintf( "InitializeSandD must be called before Simulate_logUnif_case2\n");
    return;
  }
  FM.initilize_D_and_S(Data);
  if (!Data.InitialRecordValid()) {
    Data.initial_S_Mat = 0; //reset to zeros
	Rprintf( "Initial values can not be found!\n");
    return;
  }
  is_S_and_D_initilized = true;
}
Rcpp::IntegerVector CBE::GetFaultyIndex() {
  return Rcpp::wrap(Data.Faulty2Original);
}

void CBE::SetTrueS(Rcpp::NumericMatrix S_) {
  Matrix temp = Rcpp2Mat(S_); 
  if (temp.nrows() == Data.D_Observed.nrows()) {
    Data.True_S_Mat = Matrix(Data.n_faulty,temp.ncols());
    for (int i = 1; i <= Data.n_faulty; i++ ) {
      Data.True_S_Mat.row(i) = temp.row(Data.Faulty2Original[i-1]);
    }  
  } else {    //compact matrix, assuming dim matching for now
    Data.True_S_Mat = temp;
  }
  Data.True_S_Mat = Data.True_S_Mat.t(); //transposed internally 
  Data.has_trueS = true;
}
int CBE::GetDrawFindS() {
  if (Data.has_trueS) {
    return Data.count_Draw_find_s(Param.S_Mat);
  } else {
    return -1;
  }
}

Rcpp::NumericMatrix CBE::GetDobs() {
  return Mat2Rcpp(Data.D_Observed);
}
