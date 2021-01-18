
#if !defined(_CBE_H)
#define _CBE_H
class CBE {
  public:
    CBE(Rcpp::NumericMatrix X_, Rcpp::NumericMatrix Edit_, Rcpp::NumericMatrix logB_, int nbalance_edit, int K);
    /*CBE(Rcpp::NumericMatrix X_, Rcpp::NumericMatrix Edit_, Rcpp::NumericMatrix logB_,
  Rcpp::NumericMatrix S_, Rcpp::NumericMatrix D_);*/
	  ~CBE(); //destructor
    
    void InitializeSandD();
    void SetInitialSandD(Rcpp::NumericMatrix S_, Rcpp::NumericMatrix D_);
    
    void BuildMap();
    void Simulate_logUnif_case2(int nsim);
    
    void SetTrueS(Rcpp::NumericMatrix S_);
    
    double GetRandomSeed();
    void SetRandomSeed(double seed);
    
    Rcpp::NumericMatrix Get_initialS_Compact();
    Rcpp::NumericMatrix GetD();
    
    Rcpp::NumericMatrix Get_FeasibilityMap();
    void Set_FeasibilityMap(Rcpp::NumericMatrix X_);
    
    Rcpp::NumericMatrix Get_logUnif_Compact();
    void SetS_logUnif_Compact(Rcpp::NumericMatrix X_);
    
    void SetK(int ncomp);
    int GetK();
    
    void SetMsgLevel(int level);
    int GetMsgLevel();
    
    void SetUseMap(bool use);
    bool GetUseMap();
    
    void Iterate();
    // void test();
    void Run(int iter);
    
    Rcpp::NumericMatrix GetX();
    Rcpp::NumericMatrix GetX_compact();
    Rcpp::NumericMatrix GetS_compact();
    
    double GetAlpha();
    Rcpp::NumericMatrix GetPhi();
    Rcpp::NumericVector GetSigma();
    Rcpp::NumericVector GetMu();
    Rcpp::NumericVector GetPi();
    Rcpp::NumericVector GetZ();
    Rcpp::NumericVector GetNZ();
    Rcpp::NumericVector GetAccept();
    Rcpp::IntegerVector GetFaultyIndex();
    double GetProb_A();
      
    int GetDrawFindS();
    
    void SetOptionalData(Rcpp::List OData);
    
    Rcpp::NumericMatrix GetDobs();
    
  private:
    Matrix Rcpp2Mat(Rcpp::NumericMatrix X_);
    Matrix Rcpp2Mat2(Rcpp::NumericMatrix X_);
    Rcpp::NumericMatrix Mat2Rcpp(Matrix &X_);
    void Set_initialS(Rcpp::NumericMatrix X_);
    void SetD(Rcpp::NumericMatrix X_);
    
    CData Data;
    CFeasibilityMap FM;
    CParam Param;
    CHyperParam hyper;
    
    MotherOfAll *urng;
    float randseed; 
    
    Uniform *randUnif;
    
    bool is_S_and_D_initilized;
    bool is_feasiblity_Map_initilized;
    bool is_log_unif_initilized;
    int IterCount;
    bool useMap;
};

#endif  //_CBE_H
