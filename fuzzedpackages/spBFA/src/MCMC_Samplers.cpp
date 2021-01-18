#define ARMA_DONT_PRINT_ERRORS //So the cholesky warning is suppressed
#include <RcppArmadillo.h>
#include "MCMC_bfa_sp.h"

//Function to sample latent polya-gamma process using Gibbs sampling step------------------------------------------------
arma::mat SampleOmega(int f, datobj DatObj, para Para) {
  
  //Set data objects
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  int N = DatObj.N;
  int O = DatObj.O;
  arma::mat TrialsMat = DatObj.Trials(arma::span::all, arma::span(f, f), arma::span::all);
  
  //Set parameters
  arma::colvec Mean = Para.Mean;
  
  //Moments
  arma::cube MeanOut(N, 1, 1);
  MeanOut(arma::span::all, arma::span(0, 0), arma::span(0, 0)) = Mean;
  MeanOut = arma::reshape(MeanOut, M, O, Nu);
  arma::mat MeanMat = MeanOut(arma::span::all, arma::span(f, f), arma::span::all);

  //Sample latent Variable from full conditional
  arma::mat omega = arma::reshape(pgRcpp(arma::vectorise(TrialsMat), arma::vectorise(MeanMat)), M, Nu);
  return omega;
  
}



//Function to sample latent probit process using Gibbs sampling step------------------------------------------------
datobj SampleUpper(datobj DatObj, para Para, dataug DatAug) {
  
  //Set data objects
  arma::colvec YStar = DatObj.YStar;
  arma::colvec OneNu = DatObj.OneNu;
  arma::colvec OneO = DatObj.OneO;
  
  //Set parameters
  arma::colvec Mean = Para.Mean;
  arma::cube Cov = Para.Cov;
  
  //Set data augmentation objects
  arma::uvec WhichAbove = DatAug.WhichAbove;
  int NAbove = DatAug.NAbove;
  
  //Moments
  arma::colvec Mu = Mean(WhichAbove);
  arma::colvec SDFull = arma::sqrt(arma::vectorise(Cov));
  arma::colvec SD = SDFull(WhichAbove);
  
  //Sample latent Variable from full conditional
  for (arma::uword i = 0; i < NAbove; i++) {
    double Temp = rtnormRcpp(Mu(i), SD(i), false);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(Mu(i), SD(i), 0, arma::datum::inf);
    if (!arma::is_finite(Temp)) Rcpp::stop("infinte value sampled in Probit sampling step. Most likey cause for this error is that the data being used is inappropriate (i.e., to far from zero) for a Probit model. Consider scaling towards zero and re-running.");
    YStar(WhichAbove(i)) = Temp;
  }
  DatObj.YStar = YStar;
  return DatObj;
  
}



//Function to sample latent tobit process using Gibbs sampling step------------------------------------------------
datobj SampleLower(datobj DatObj, para Para, dataug DatAug) {
  
  //Set data objects
  arma::colvec YStar = DatObj.YStar;
  arma::colvec OneNu = DatObj.OneNu;
  arma::colvec OneO = DatObj.OneO;
  
  //Set parameters
  arma::colvec Mean = Para.Mean;
  arma::cube Cov = Para.Cov;
  
  //Set data augmentation objects
  int NBelow = DatAug.NBelow;
  arma::uvec WhichBelow = DatAug.WhichBelow;
  
  //Moments
  arma::colvec Mu = Mean(WhichBelow);
  arma::colvec SDFull = arma::sqrt(arma::vectorise(Cov));
  arma::colvec SD = SDFull(WhichBelow);
  
  //Sample latent Variable from full conditional
  for (arma::uword i = 0; i < NBelow; i++) {
    double Temp = rtnormRcpp(Mu(i), SD(i), true);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(Mu(i), SD(i), -arma::datum::inf, 0);
    if (!arma::is_finite(Temp)) Rcpp::stop("infinte value sampled in Tobit/Probit sampling step. Most likely cause for this error is that the data being used is inappropriate (i.e., to far from zero) for a Tobit/Probit model. Consider scaling towards zero and re-running.");
    YStar(WhichBelow(i)) = Temp;
  }
  DatObj.YStar = YStar;
  return DatObj;
  
}



//Function to sample latent process from its full conditional------------------------------------------------------
std::pair<datobj, para> SampleY(datobj DatObj, para Para, dataug DatAug) {
  
  //Set data objects
  arma::Col<int> FamilyInd = DatObj.FamilyInd;
  int N = DatObj.N;
  int M = DatObj.M;
  int O = DatObj.O;
  int Nu = DatObj.Nu;
  arma::cube Chi = DatObj.Chi;
  
  //Set parameter objects
  arma::cube Cov = Para.Cov;
  
  //Begin by updating all Tobit/Probit latent variables
  if (DatAug.NAbove > 0) DatObj = SampleUpper(DatObj, Para, DatAug);
  if (DatAug.NBelow > 0) DatObj = SampleLower(DatObj, Para, DatAug);
  arma::colvec YStar = DatObj.YStar;
  
  if (any(FamilyInd == 3)) {
    
    //Update Polya-Gamma latent variables
    arma::cube YStarOut(N, 1, 1);
    YStarOut(arma::span::all, arma::span(0, 0), arma::span(0, 0)) = YStar;
    YStarOut = arma::reshape(YStarOut, M, O, Nu);
  
    //Loop over latent dimensions
    for (arma::uword f = 0; f < O; f++) {
  
      //Observation specifics
      int FamInd = FamilyInd(f);
  
      //Polya-Gamma latent process updates
      if (FamInd == 3) {
        
        //Sample omega
        arma::mat chi = Chi(arma::span::all, arma::span(f, f), arma::span::all);
        arma::mat omega = SampleOmega(f, DatObj, Para); //M x Nu
        arma::mat omegaInv = arma::pow(omega, -1);
        
        //Update Polya-Gamma objects
        YStarOut(arma::span::all, arma::span(f, f), arma::span::all) = (chi % omegaInv);
        Cov(arma::span::all, arma::span(f, f), arma::span::all) = omegaInv;
  
      }
      
    //End loop over observation dimension
    }
    YStar = arma::vectorise(YStarOut);
  }
  
  //Save output
  DatObj.YStar = YStar;
  DatObj.YStarWide = arma::reshape(YStar, M * O, Nu);
  Para.Cov = Cov;
  return std::pair<datobj, para>(DatObj, Para);

}



//Function to sample sigma2(s_i) using a Gibbs sampler step---------------------------------------------------------------
para SampleSigma2(datobj DatObj, para Para, hypara HyPara) {
  
  //Set data objects
  arma::Col<int> FamilyInd = DatObj.FamilyInd;
  int Nu = DatObj.Nu;
  int O = DatObj.O;
  int M = DatObj.M;
  int C = DatObj.C;
  arma::mat YStarWide = DatObj.YStarWide;
  arma::mat EyeO = DatObj.EyeO;
  
  //Set parameter objects
  arma::colvec Mean = Para.Mean;
  arma::mat MeanMat = arma::reshape(Mean, M * O, Nu);
  arma::mat Sigma2 = Para.Sigma2;
  arma::cube Cov = Para.Cov;

  //Set hyperparameter objects
  double A = HyPara.A;
  double B = HyPara.B;
  
  //Shape is constant over locations
  double Shape = A + 0.5 * Nu;
  
  //Declarations
  arma::uword Index;
  
  //Loop over observations
  arma::uvec NotCount = find(FamilyInd != 3);
  for (arma::uword c = 0; c < (O - C); c++) {
    
    //Set the observation 
    int o = NotCount(c);
    
    //Loop over locations
    for (arma::uword i = 0; i < M; i++) {
    
      //Location specific objects
      Index = i + M * o;
      arma::rowvec Diff = YStarWide.row(Index) - MeanMat.row(Index);
  
      //Calculate rate
      double Resids = arma::as_scalar(Diff * arma::trans(Diff));
    
      //Sample sigma2Io
      double Rate = B + 0.5 * Resids;
      Sigma2(i, c) = rigammaRcpp(Shape, Rate);
      
    //End loop over locations  
    }
    
    //Update covariance
    Cov(arma::span::all, arma::span(o, o), arma::span::all) = arma::repmat(Sigma2.col(c), 1, Nu);

  //End loop over observations 
  }
  
  //Update parameters object
  Para.Sigma2 = Sigma2;
  Para.Cov = Cov;
  return Para;
  
}



//Function to sample new value of psi using a Metropolis sampler step-----------------------------------------------
std::pair<para, metrobj> SamplePsi(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj) {
  
  //Set data objects
  arma::mat TimeDist = DatObj.TimeDist;
  int Nu = DatObj.Nu;
  int TempCorInd = DatObj.TempCorInd;
  arma::mat EyeKbyNu = DatObj.EyeKbyNu;
  arma::colvec ZeroKbyNu = DatObj.ZeroKbyNu;
  arma::mat EyeNu = DatObj.EyeNu;
  
  //Set parameter objects
  double Psi = Para.Psi;
  arma::mat HPsi = Para.HPsi;
  arma::mat CholHPsi = Para.CholHPsi;
  arma::mat Upsilon = Para.Upsilon;
  arma::colvec Eta = Para.Eta;

  //Set hyperparameter objects
  double APsi = HyPara.APsi;
  double BPsi = HyPara.BPsi;
  double Gamma = HyPara.Gamma;
  double Beta = HyPara.Beta;
  
  //Set metropolis objects
  double MetropPsi = sqrt(MetrObj.MetropPsi);
  double AcceptancePsi = MetrObj.AcceptancePsi;
  
  //Transform current state to real line
  double BigDelta = log((Psi - APsi) / (BPsi - Psi));
  
  //Numerical fix for when the propopsal cholesky doesn't exist
  double PsiProposal, BigDeltaProposal;
  arma::mat CholHPsiProposal(Nu, Nu), HPsiProposal(Nu, Nu);
  bool Cholesky = false;
  while (!Cholesky) {
    
    //Sample a new Proposal
    BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropPsi));
    
    //Compute Phi Proposal
    PsiProposal = (BPsi * exp(BigDeltaProposal) + APsi) / (1 + exp(BigDeltaProposal));

    //Fix numerical issue where PsiProposal can equal APsi or BPsi
    // arma::vec PsiProposalVec(1), APsiVec(1), BPsiVec(1);
    // PsiProposalVec(0) = PsiProposal;
    // APsiVec(0) = APsi;
    // BPsiVec(0) = BPsi;
    // double TOL = 0.000001;
    // if ((rows_equal(PsiProposalVec, APsiVec, TOL)) || (rows_equal(PsiProposalVec, BPsiVec, TOL))) {
    //   if (rows_equal(PsiProposalVec, APsiVec, TOL)) PsiProposal *= 1.1; //doesn't work when APsi is negative
    //   if (rows_equal(PsiProposalVec, BPsiVec, TOL)) PsiProposal *= 0.99;
    //   BigDeltaProposal = log((PsiProposal - APsi) / (BPsi - PsiProposal));
    // }
    
    //Proposal temporal correlation
    HPsiProposal = H(PsiProposal, TempCorInd, TimeDist, Nu);
    Cholesky = arma::chol(CholHPsiProposal, HPsiProposal);
    
  }
  
  //Eta structure components
  arma::mat CholUpsilon = arma::chol(Upsilon);
  arma::mat RootiEta = arma::solve(arma::trimatu(arma::kron(CholHPsi, CholUpsilon)), EyeKbyNu);
  arma::mat RootiEtaProposal = arma::solve(arma::trimatu(arma::kron(CholHPsiProposal, CholUpsilon)), EyeKbyNu);
  double Component1A = lndMvn(Eta, ZeroKbyNu, RootiEtaProposal);
  double Component1B = lndMvn(Eta, ZeroKbyNu, RootiEta);
  double Component1 = Component1A - Component1B;
  
  //Prior component
  double Component2 = 0; //exponential
  if (TempCorInd == 1) { //ar1
    double Component2A = (Gamma - 1) * log(1 + PsiProposal) + (Beta - 1) * log(1 - PsiProposal);
    double Component2B = (Gamma - 1) * log(1 + Psi) + (Beta - 1) * log(1 - Psi);
    Component2 = Component2A - Component2B;
  }
  
  //Jacobian component 1
  double Component3A = BigDeltaProposal;
  double Component3B = BigDelta;
  double Component3 = Component3A - Component3B;
  
  //Jacobian component 2
  double Component4 = 2 * log((1 + exp(BigDelta)) / (1 + exp(BigDeltaProposal)));
  
  //Compute log acceptance ratio
  double LogR = Component1 + Component2 + Component3 + Component4;
  
  //Metropolis update
  double RandU = randuRcpp();
  if (log(RandU) < LogR) {
    
    //Keep Count of Acceptances
    AcceptancePsi++;
    MetrObj.AcceptancePsi = AcceptancePsi;
    
    //Update dependent parameters
    arma::mat P = arma::solve(arma::trimatu(CholHPsiProposal), EyeNu);
    
    //Update parameters object
    Para.Psi = PsiProposal;
    Para.HPsi = HPsiProposal;
    Para.CholHPsi = CholHPsiProposal;
    Para.HPsiInv = P * arma::trans(P);

  }
  
  //Return output object
  return std::pair<para, metrobj>(Para, MetrObj);
  
}



//Function to sample Upsilon using a Gibbs sampler step-------------------------------------------------------------------
para SampleUpsilon(datobj DatObj, para Para, hypara HyPara) {
  
  //Set data objects
  int Nu = DatObj.Nu;
  arma::mat EyeO = DatObj.EyeO;
  int K = DatObj.K;

  //Set parameters
  arma::mat BigPhi = Para.BigPhi;
  arma::mat HPsiInv = Para.HPsiInv;
  arma::mat Lambda = Para.Lambda;

  //Set hyperparameter objects
  double Zeta = HyPara.Zeta;
  arma::mat Omega = HyPara.Omega;
  
  //Compute SPhiPsi
  arma::mat SPhiPsi = BigPhi * HPsiInv * arma::trans(BigPhi);
  
  //Sample Upsilon
  double n = Zeta + Nu;
  arma::mat V = SPhiPsi + Omega;
  arma::mat Upsilon(K, K), UpsilonInv(K, K);
  if (K > 1) {
    UpsilonInv = rwishRcpp(n, CholInv(V));
    Upsilon = CholInv(UpsilonInv);
  } else {
    Upsilon(0, 0) = rigammaRcpp(0.5 * n, 0.5 * arma::as_scalar(V));
    UpsilonInv = 1 / Upsilon;
  }
  
  //Update parameters object
  Para.Upsilon = Upsilon;
  Para.UpsilonInv = UpsilonInv;
  return Para;
}



//Function to sample beta using a Gibbs sampler step---------------------------------------------------------------
para SampleBeta(datobj DatObj, para Para, hypara HyPara) {
  
  //Set data objects
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat YStarWide = DatObj.YStarWide;
  int Nu = DatObj.Nu;
  int P = DatObj.P;
  arma::mat X = DatObj.X;
  arma::Col<int> Indeces = DatObj.Indeces;
  
  //Set parameters
  arma::mat Lambda = Para.Lambda;
  arma::mat Sigma2 = Para.Sigma2;
  arma::cube Cov = Para.Cov;
  arma::colvec Eta = Para.Eta;
  arma::mat BigPhi = Para.BigPhi;
  
  //Set hyperparameters
  arma::mat SigmaBetaInv = HyPara.SigmaBetaInv;
  arma::colvec SigmaBetaInvMuBeta = HyPara.SigmaBetaInvMuBeta;

  //Compute moments
  arma::mat Sum1(P, P, arma::fill::zeros);
  arma::colvec Sum2(P, arma::fill::zeros);
  for (arma::uword t = 0; t < Nu; t++) {
    arma::mat SigmaTInv = arma::diagmat(arma::vectorise(1 / Cov.slice(t)));
    arma::mat XT = X.rows(find(Indeces == t));
    arma::mat tXTSigmaTInv = arma::trans(XT) * SigmaTInv;
    Sum1 += tXTSigmaTInv * XT;
    Sum2 += tXTSigmaTInv * (YStarWide.col(t) - Lambda * BigPhi.col(t));
  }
  
  //Sample Beta
  arma::mat CovBeta = CholInv(Sum1 + Nu * SigmaBetaInv);
  arma::colvec MeanBeta = CovBeta * (Sum2 + Nu * SigmaBetaInvMuBeta);
  arma::colvec Beta = rmvnormRcpp(1, MeanBeta, CovBeta);
  
  //Update parameters dependent on delta
  arma::colvec XBeta = X * Beta;
  arma::colvec Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;

  //Update parameters object
  Para.Beta = Beta;
  Para.Mean = Mean;
  Para.XBeta = XBeta;
  return Para;
  
}



//Function to sample eta using a Gibbs sampler step---------------------------------------------------------------
para SampleEta(datobj DatObj, para Para, hypara HyPara) {
  
  //Set data objects
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat EyeK = DatObj.EyeK;
  arma::mat YStarWide = DatObj.YStarWide;
  int K = DatObj.K;
  int Nu = DatObj.Nu;
  int M = DatObj.M;
  int O = DatObj.O;

  //Set parameters
  arma::mat BigPhi = Para.BigPhi;
  arma::mat Lambda = Para.Lambda;
  arma::cube Cov = Para.Cov;
  arma::mat HPsi = Para.HPsi;
  arma::mat UpsilonInv = Para.UpsilonInv;
  arma::colvec XBeta = Para.XBeta;
  arma::mat XBetaMat = arma::reshape(XBeta, M * O, Nu);

  //Declarations
  arma::vec SeqNu = arma::linspace<arma::vec>(0, Nu - 1, Nu);
  arma::uvec IndecesT(1), IndecesMinusT(Nu - 1);
  arma::colvec EtaT(K);
  
  //Loop over t
  for (arma::uword t = 0; t < Nu; t++) {
    
    //Conditional moments
    IndecesMinusT = find(SeqNu != t);
    IndecesT(0) = t;
    arma::mat BigPhiMinusT = BigPhi;
    BigPhiMinusT.shed_col(t);
    arma::rowvec HPlus = HPsi(IndecesT, IndecesMinusT) * CholInv(HPsi(IndecesMinusT, IndecesMinusT));
    double HStarInv = 1 / (arma::as_scalar(HPsi(IndecesT, IndecesT) - HPlus * HPsi(IndecesMinusT, IndecesT)));
    arma::colvec CondMuEta = arma::kron(HPlus, EyeK) * arma::vectorise(BigPhiMinusT);
    arma::mat CondPrecEta = HStarInv * UpsilonInv;
    
    //Sample EtaT
    arma::mat SigmaTInv = arma::diagmat(arma::vectorise(1 / Cov.slice(t)));
    arma::mat tLambdaSigmaInv = arma::trans(Lambda) * SigmaTInv;
    arma::mat CovEtaT = CholInv(tLambdaSigmaInv * Lambda + CondPrecEta);
    arma::colvec MeanEtaT = CovEtaT * (tLambdaSigmaInv * (YStarWide.col(t) - XBetaMat.col(t)) + CondPrecEta * CondMuEta);
    
    EtaT = rmvnormRcpp(1, MeanEtaT, CovEtaT);
    BigPhi.col(t) = EtaT;
    
  //End loop over visits
  }
  
  //Update parameters dependent on eta
  arma::colvec Eta = arma::vectorise(BigPhi);
  
  //Update parameters object
  Para.Eta = Eta;
  Para.BigPhi = BigPhi;
  Para.Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;
  return Para;
  
}



//Function to sample delta using a Gibbs sampler step---------------------------------------------------------------
para SampleDelta(datobj DatObj, para Para, hypara HyPara) {
  
  //Set data objects
  int K = DatObj.K;
  int L = DatObj.L;
  int LInf = DatObj.LInf;
  int GS = DatObj.GS;
  
  //Set parameter objects
  arma::mat Theta = Para.Theta;
  arma::colvec Delta = Para.Delta;
  arma::colvec LStarJ = Para.LStarJ;
  
  //Set hyperparameter objects
  double A1 = HyPara.A1;
  double A2 = HyPara.A2;
  
  //Upper bound for L
  int UpperL = L;
  
  //Gamma shrinkage prior
  if (GS == 1) {
    
    //Loop over K delta precision parameters
    double AH;
    for (arma::uword h = 0; h < K; h++) {
      
      //Remove hth delta
      arma::colvec DeltaMinusH = Delta;
      DeltaMinusH(h) = 1;
      
      //Asign hyperparameter
      if (h == 0) AH = A1;
      if (h > 0) AH = A2;
      
      //Obtain moments
      double Resids = 0;
      for (arma::uword j = h; j < K; j++) {
        arma::colvec ThetaJ = Theta.col(j);
        double tThetaTheta;
        if (LInf == 0) tThetaTheta = arma::as_scalar(arma::trans(ThetaJ) * ThetaJ);
        if (LInf == 1) {
          UpperL = LStarJ(j);
          arma::colvec ThetaJUpperL = ThetaJ(arma::span(0, UpperL));
          tThetaTheta = arma::as_scalar(arma::trans(ThetaJUpperL) * ThetaJUpperL);
        }
        Resids += arma::as_scalar(tThetaTheta * arma::prod(DeltaMinusH(arma::span(0, j))));
      }
      double Shape = AH;
      if (LInf == 0) Shape += 0.5 * (K - h) * UpperL;
      if (LInf == 1) Shape += 0.5 * arma::sum(LStarJ(arma::span(h, K - 1)) + 1);
      double Rate = 1 + 0.5 * Resids;
      
      //Sample Deltah
      Delta(h) = rgammaRcpp(Shape, Rate);
  
    //End loop over deltas  
    }
  
  //End shrinkage prior
  }
  
  //NO Gamma shrinkage prior
  if (GS == 0) {
    
    //Loop over K delta precision parameters
    for (arma::uword j = 0; j < K; j++) {
      
      //Obtain moments
      arma::colvec ThetaJ = Theta.col(j);
      double tThetaTheta;
      if (LInf == 0) tThetaTheta = arma::as_scalar(arma::trans(ThetaJ) * ThetaJ);
      if (LInf == 1) {
        UpperL = LStarJ(j);
        arma::colvec ThetaJUpperL = ThetaJ(arma::span(0, UpperL));
        tThetaTheta = arma::as_scalar(arma::trans(ThetaJUpperL) * ThetaJUpperL);
      }
      double Shape = A1;
      if (LInf == 0) Shape += 0.5 * UpperL;
      if (LInf == 1) Shape += 0.5 * LStarJ(j) + 1;
      double Rate = A2 + 0.5 * tThetaTheta;
      
      //Sample Deltah
      Delta(j) = rgammaRcpp(Shape, Rate);
      
    //End loop over deltas  
    }
    
  //End NO shrinkage prior
  }
  
  //Update parameters object
  Para.Delta = Delta;
  if (GS == 1) Para.Tau = arma::cumprod(Delta);
  if (GS == 0) Para.Tau = Delta;
  return Para;
  
}



//Function to sample new value of psi using a Metropolis sampler step-----------------------------------------------
std::pair<para, metrobj> SampleRho(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj) {
  
  //Set data objects
  int M = DatObj.M;
  int O = DatObj.O;
  arma::mat SpDist = DatObj.SpDist;
  int L = DatObj.L;
  int LInf = DatObj.LInf;
  int K = DatObj.K;
  arma::colvec ZeroOM = DatObj.ZeroOM;
  arma::mat EyeOM = DatObj.EyeOM;
  arma::mat EyeM = DatObj.EyeM;
  
  //Set parameter objects
  double Rho = Para.Rho;
  arma::mat CholKappa = Para.CholKappa;
  arma::mat SpCov = Para.SpCov;
  arma::mat CholSpCov = Para.CholSpCov;
  arma::colvec LStarJ = Para.LStarJ;
  arma::cube Alpha = Para.Alpha;
  
  //Set hyperparameter objects
  double ARho = HyPara.ARho;
  double BRho = HyPara.BRho;

  //Set metropolis objects
  double MetropRho = sqrt(MetrObj.MetropRho);
  double AcceptanceRho = MetrObj.AcceptanceRho;
  
  //Transform current state to real line
  double BigDelta = log((Rho - ARho) / (BRho - Rho));
  
  //Numerical fix for when the propopsal cholesky doesn't exist
  double RhoProposal, BigDeltaProposal;
  arma::mat CholSpCovProposal(M, M), SpCovProposal(M, M);
  bool Cholesky = false;
  while (!Cholesky) {
    
    //Sample a new Proposal
    BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropRho));
    
    //Compute Phi Proposal
    RhoProposal = (BRho * exp(BigDeltaProposal) + ARho) / (1 + exp(BigDeltaProposal));
    
    //Fix numerical issue where RhoProposal can equal ARho or BRho
    // arma::vec RhoProposalVec(1), ARhoVec(1), BRhoVec(1);
    // RhoProposalVec(0) = RhoProposal;
    // ARhoVec(0) = ARho;
    // BRhoVec(0) = BRho;
    // double TOL = 0.000001;
    // if ((rows_equal(RhoProposalVec, ARhoVec, TOL)) || (rows_equal(RhoProposalVec, BRhoVec, TOL))) {
    //   if (rows_equal(RhoProposalVec, ARhoVec, TOL)) RhoProposal *= 1.1; //doesn't work when ARho is negative
    //   if (rows_equal(RhoProposalVec, BRhoVec, TOL)) RhoProposal *= 0.99;
    //   BigDeltaProposal = log((RhoProposal - ARho) / (BRho - RhoProposal));
    // }
    
    //Proposal temporal correlation
    SpCovProposal = SpEXP(RhoProposal, SpDist, M);
    Cholesky = arma::chol(CholSpCovProposal, SpCovProposal);
    
  }
  
  //Upper bound for l
  int UpperL = L;
  
  //Alpha structure components
  double Component1A = 0, Component1B = 0;
  arma::mat RootiAlpha(O * M, O * M), RootiAlphaProposal(O * M, O * M);
  
  //Loop over columns of j
  for (arma::uword j = 0; j < K; j++) {
    
    //Objects that depend on j
    if (LInf == 1) UpperL = LStarJ(j);
    arma::mat AlphaJ = Alpha.slice(j);
    
    //Loop over mixture components
    for (arma::uword l = 0; l < UpperL; l++) {
      arma::colvec AlphaJL = AlphaJ.row(l).t();
      RootiAlpha = arma::solve(arma::trimatu(arma::kron(CholKappa, CholSpCov)), EyeOM);
      RootiAlphaProposal = arma::solve(arma::trimatu(arma::kron(CholKappa, CholSpCovProposal)), EyeOM);
      Component1A += lndMvn(AlphaJL, ZeroOM, RootiAlphaProposal);
      Component1B += lndMvn(AlphaJL, ZeroOM, RootiAlpha);
    }
  }
  double Component1 = Component1A - Component1B;

  //Jacobian component 1
  double Component2A = BigDeltaProposal;
  double Component2B = BigDelta;
  double Component2 = Component2A - Component2B;
  
  //Jacobian component 2
  double Component3 = 2 * log((1 + exp(BigDelta)) / (1 + exp(BigDeltaProposal)));
  
  //Compute log acceptance ratio
  double LogR = Component1 + Component2 + Component3;
  
  //Metropolis update
  double RandU = randuRcpp();
  if (log(RandU) < LogR) {
    
    //Keep Count of Acceptances
    AcceptanceRho++;
    MetrObj.AcceptanceRho = AcceptanceRho;
    
    //Update dependent parameters
    arma::mat P = arma::solve(arma::trimatu(CholSpCovProposal), EyeM);
    
    //Update parameters object
    Para.Rho = RhoProposal;
    Para.SpCov = SpCovProposal;
    Para.CholSpCov = CholSpCovProposal;
    Para.SpCovInv = P * arma::trans(P);
    
  }
  
  //Return output object
  return std::pair<para, metrobj>(Para, MetrObj);
  
}



//Function to sample kappa using a Gibbs sampler step---------------------------------------------------------------
para SampleKappa(datobj DatObj, para Para, hypara HyPara) {
  
  //Set data objects
  int K = DatObj.K;
  int M = DatObj.M;
  int O = DatObj.O;
  int L = DatObj.L;
  int LInf = DatObj.LInf;
  arma::mat EyeO = DatObj.EyeO;

  //Set parameter objects
  arma::cube Alpha = Para.Alpha;
  arma::colvec LStarJ = Para.LStarJ;
  arma::mat SpCovInv = Para.SpCovInv;
  
  //Set hyperparameter objects
  double SmallUpsilon = HyPara.SmallUpsilon;
  arma::mat BigTheta = HyPara.BigTheta;
  
  //Upper bound for L
  int UpperL = L;
  
  //Calculate moments
  arma::mat Resids(O, O, arma::fill::zeros);
  for (arma::uword j = 0; j < K; j++) {
    if (LInf == 1) UpperL = LStarJ(j);
    for (arma::uword l = 0; l < UpperL; l++) {
      arma::rowvec AlphaJL = Alpha.slice(j).row(l);
      arma::mat AlphaJLMat = arma::reshape(AlphaJL, M, O);
      Resids += arma::trans(AlphaJLMat) * SpCovInv * AlphaJLMat;
    }
  }
  double Shape = SmallUpsilon;
  if (LInf == 0) Shape += M * L * K;
  if (LInf == 1) Shape += M * arma::as_scalar(arma::sum(LStarJ + 1));
  arma::mat Rate = Resids + BigTheta;
  
  //Sample kappa
  arma::mat Kappa(O, O), KappaInv(O, O), CholKappa(O, O);
  if (O == 1) {
    Kappa(0, 0) = rigammaRcpp(0.5 * arma::as_scalar(Shape), 0.5 * arma::as_scalar(Rate));
    CholKappa = arma::sqrt(Kappa);
    KappaInv = 1 / Kappa;
  }
  if (O > 1) {
    KappaInv = rwishRcpp(Shape, CholInv(Rate));
    CholKappa = arma::solve(arma::trimatu(arma::chol(KappaInv)), EyeO);
    Kappa = arma::trans(CholKappa) * CholKappa;
  }
  
  //Update parameters object
  Para.Kappa = Kappa;
  Para.KappaInv = KappaInv;
  Para.CholKappa = CholKappa;
  return Para;
  
}



//Function to sample alphajl(s_i)'s using a Gibbs sampler step---------------------------------------------------------------
para SampleAlpha(datobj DatObj, para Para) {
  
  //Different samplers depending on whether a BNP model is used
  if (DatObj.CL == 1) {
  
    //Set data objects
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int L = DatObj.L;
    arma::mat EyeOM = DatObj.EyeOM;
    int LInf = DatObj.LInf;
    
    //Set parameter objects
    arma::mat KappaInv = Para.KappaInv;
    arma::cube Z = Para.Z;
    arma::cube Alpha = Para.Alpha;
    arma::colvec LStarJ = Para.LStarJ;
    arma::mat U = Para.U;
    arma::mat SpCovInv = Para.SpCovInv;
    
    //Upper bound for L
    int UpperL = L;
    
    //Covariance object
    arma::mat CovAlpha = CholInv(EyeOM + arma::kron(KappaInv, SpCovInv));
    
    //Loop over columns K
    for (arma::uword j = 0; j < K; j++) {
      
      //Get jth process objects
      arma::mat AlphaJ = Alpha.slice(j); 
      arma::mat ZJ = Z.slice(j);
      if (LInf == 1) UpperL = LStarJ(j);
      
      //Loop over clusters L
      for (arma::uword l = 0; l < UpperL; l++) {
        
        //Sample AlphaJL
        arma::colvec MeanAlpha = CovAlpha * arma::trans(ZJ.row(l));
        arma::colvec AlphaJL = rmvnormRcpp(1, MeanAlpha, CovAlpha);
        AlphaJ.row(l) = arma::trans(AlphaJL);  
          
      //End loop over clusters  
      }
        
      //Update Z
      Alpha.slice(j) = AlphaJ;
      
    //End loop over columns 
    }
    
    //Update Weights
    arma::cube Weights = GetWeights(Alpha, K, M, L, O);
    arma::cube logWeights = GetlogWeights(Alpha, K, M, L, O);
    if (LInf == 1) LStarJ = GetLStarJ(U, Weights, K, M, O);
    
    //Update parameters object
    Para.Alpha = Alpha;
    Para.logWeights = logWeights;
    Para.Weights = Weights;
    Para.LStarJ = LStarJ;
   
  //End sampler for BNP version 
  }

  //Sampler for non-BNP version
  if (DatObj.CL == 0) {
    
    //Set data objects
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int Nu = DatObj.Nu;
    arma::mat EyeNu = DatObj.EyeNu;
    arma::mat X = DatObj.X;
    arma::Col<int> Indeces = DatObj.Indeces;
    arma::mat YStarWide = DatObj.YStarWide;

    //Set parameter objects
    arma::mat KappaInv = Para.KappaInv;
    arma::mat SpCovInv = Para.SpCovInv;
    arma::mat BigPhi = Para.BigPhi;
    arma::cube Cov = Para.Cov;
    arma::colvec Beta = Para.Beta;
    arma::mat Lambda = Para.Lambda;
    arma::colvec Eta = Para.Eta;
    arma::colvec XBeta = Para.XBeta;

    //Covariance object
    arma::mat CovAlphaFixed = arma::kron(KappaInv, SpCovInv);
    
    //Loop over columns K
    arma::cube Alpha(1, M * O, K);
    for (arma::uword j = 0; j < K; j++) {
      
      //Factor specific objects
      arma::mat LambdaMinusJ = Lambda, BigPhiMinusJ = BigPhi;
      LambdaMinusJ.shed_col(j);
      BigPhiMinusJ.shed_row(j);
      
      //Loop over time
      arma::mat Sum1(M * O, M * O, arma::fill::zeros);
      arma::colvec Sum2(M * O, arma::fill::zeros);
      for (arma::uword t = 0; t < Nu; t++) {
        arma::mat SigmaTInv = arma::diagmat(arma::vectorise(1 / Cov.slice(t)));
        arma::mat XT = X.rows(find(Indeces == t));
        double EtaTJ = BigPhi(j, t);
        arma::colvec LambdaSum(M * O, arma::fill::zeros);
        for (arma::uword h = 0; h < (K - 1); h++) LambdaSum += LambdaMinusJ.col(h) * BigPhiMinusJ(h, t);
        arma::colvec MuTJ = (YStarWide.col(t) - XT * Beta - LambdaSum);
        Sum1 += SigmaTInv * (EtaTJ * EtaTJ);
        Sum2 += EtaTJ * SigmaTInv * MuTJ;
      }
      
      //Sample AlphaJ
      arma::mat CovAlpha = CholInv(Sum1 + CovAlphaFixed);
      arma::colvec MeanAlpha = CovAlpha * Sum2;
      arma::colvec AlphaJ = rmvnormRcpp(1, MeanAlpha, CovAlpha);
      Alpha(arma::span(0, 0), arma::span::all, arma::span(j, j)) = AlphaJ; 
      Lambda.col(j) = AlphaJ;
      
    //End loop over columns 
    }
    
    //Update parameters object
    Para.Alpha = Alpha;
    Para.Lambda = Lambda;
    Para.Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;

  //End non-BNP sampler  
  }
  return Para;
  
}



//Function to sample zjl(s_i)'s using a Gibbs sampler step---------------------------------------------------------------
para SampleZ(datobj DatObj, para Para) {
  
  //Set data objects
  int K = DatObj.K;
  int M = DatObj.M;
  int O = DatObj.O;
  int L = DatObj.L;
  int LInf = DatObj.LInf;
  
  //Set parameter objects
  arma::cube Alpha = Para.Alpha;
  arma::umat Xi = Para.Xi;
  arma::cube Z = Para.Z;
  arma::colvec LStarJ = Para.LStarJ;
  
  //Upper bound for L
  int UpperL = L;
  
  //Declarations
  arma::uword Index;
  
  //Loop over columns K
  for (arma::uword j = 0; j < K; j++) {
    
    //Get jth process objects
    arma::mat AlphaJ = Alpha.slice(j);
    if (LInf == 1) UpperL = LStarJ(j);
    
    //Loop over spatial observations
    arma::mat ZJ(L, M * O);
    for (arma::uword o = 0; o < O; o++) {
    
      //Loop over spatial locations
      for (arma::uword i = 0; i < M; i++) {
      
      //Location specific
      Index = i + M * o;
      arma::uword XiOIJ = Xi(Index, j);

      //Loop over clusters L
      for (arma::uword l = 0; l < UpperL; l++) {
        
        //Sample zjl(s_i)
        double ZOIJL;
        if (XiOIJ > l) ZOIJL = rtnormRcppMSM(AlphaJ(l, Index), 1, -arma::datum::inf, 0);
        if (XiOIJ == l) ZOIJL = rtnormRcppMSM(AlphaJ(l, Index), 1, 0, arma::datum::inf);
        ZJ(l, Index) = ZOIJL;
      
      //End loop over clusters  
      }
    
    //End loop over locations
    }
    
    //End loop over spatial observations
    }
    
    //Update Z
    Z.slice(j) = ZJ;
  
  //End loop over columns 
  }
  
  //Update parameters object
  Para.Z = Z;
  return Para;
  
}



//Function to sample xi's using a Gibbs sampler step---------------------------------------------------------------
para SampleXi(datobj DatObj, para Para) {
  
  //Set data objects
  int K = DatObj.K;
  int M = DatObj.M;
  int O = DatObj.O;
  int L = DatObj.L;
  int Nu = DatObj.Nu;
  int LInf = DatObj.LInf;
  arma::Col<int> SeqL = DatObj.SeqL;
  arma::mat YStarWide = DatObj.YStarWide;
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat EyeO = DatObj.EyeO;
  
  //Set parameter objects
  arma::umat Xi = Para.Xi;
  arma::cube logWeights = Para.logWeights;
  arma::mat Lambda = Para.Lambda;
  arma::mat BigPhi = Para.BigPhi;
  arma::cube Cov = Para.Cov;
  arma::mat Theta = Para.Theta;
  arma::colvec Eta = Para.Eta;
  arma::mat Upsilon = Para.Upsilon;
  arma::vec LStarJ = Para.LStarJ;
  arma::mat U = Para.U;
  arma::colvec XBeta = Para.XBeta;
  arma::mat XBetaMat = arma::reshape(XBetaMat, M * O, Nu);
    
  //Upper bound for L
  int UpperL = L;
  
  //Declarations
  arma::uword Index;

  //Loop over columns K
  for (arma::uword j = 0; j < K; j++) {
    
    //Get jth process objects
    arma::mat logWeightsJ = logWeights.slice(j);
    arma::colvec ThetaJ = Theta.col(j);
    if (LInf == 1) UpperL = LStarJ(j);
    arma::colvec UJ = U.col(j);
    
    //Loop over spatial observations
    for (arma::uword o = 0; o < O; o++) {
    
      arma::mat CovO = Cov(arma::span::all, arma::span(o, o), arma::span::all);
      
      //Loop over locations M
      for (arma::uword i = 0; i < M; i++) {
        
        //Get location specific objects
        Index = i + M * o;
        arma::colvec logWeightsOIJ = logWeightsJ.col(Index);
        arma::rowvec CovOI = CovO.row(i);
        arma::rowvec LambdaOI = Lambda.row(Index);
        double UOIJ = UJ(Index);
        
        //Loop over clusters L
        arma::vec logProbsRaw(UpperL);
        logProbsRaw.fill(-arma::datum::inf);
        for (arma::uword l = 0; l < UpperL; l++) {
          
          //For finite mixture model
          if (LInf == 0) {
            
            //Obtain the un-normalized (raw) probabilities on the log scale
            LambdaOI(j) = ThetaJ(l);
            arma::rowvec Resid = YStarWide.row(Index) - XBetaMat.row(Index) - LambdaOI * BigPhi;
            Resid = Resid % (1 / arma::sqrt(CovOI));
            double ResidQ = arma::as_scalar(Resid * arma::trans(Resid));
            double Likelihood = -0.5 * ResidQ;
            double logWeightsOIJL = logWeightsOIJ(l);
            logProbsRaw(l) = logWeightsOIJL + Likelihood;
            
          }
          
          //For infinite mixture model
          if (LInf == 1) {
            
            //Determine if it is non-zero
            bool Include = logWeightsOIJ(l) > log(UOIJ);
            
            //Obtain the un-normalized (raw) probabilities on the log scale
            if (Include) {
              LambdaOI(j) = ThetaJ(l);
              arma::rowvec Resid = YStarWide.row(Index) - XBetaMat.row(Index) - LambdaOI * BigPhi;
              Resid = Resid % (1 / arma::sqrt(CovOI));
              double ResidQ = arma::as_scalar(Resid * arma::trans(Resid));
              double Likelihood = -0.5 * ResidQ;
              logProbsRaw(l) = Likelihood;
            }
            
          }
        
        //End loop over clusters  
        }
      
        //Use log sum exponential trick to get normalized probabilities
        double Max = arma::max(logProbsRaw);
        double Delta = Max + log(arma::sum(arma::exp(logProbsRaw - Max)));
        arma::vec ProbsOIJ = arma::exp(logProbsRaw - Delta);
        
        //Sample a new label
        arma::uword XiOIJ;
        if (LInf == 0) XiOIJ = arma::as_scalar(sampleRcpp(SeqL, 1, true, ProbsOIJ));
        if (LInf == 1) {
          arma::vec Finite = logProbsRaw.elem(arma::find_finite(logProbsRaw));
          if (Finite.size() > 0) {
            arma::Col<int> SeqLStarJ(UpperL);
            for (arma::uword seq = 0; seq < UpperL; seq++) SeqLStarJ(seq) = seq;
            XiOIJ = arma::as_scalar(sampleRcpp(SeqLStarJ, 1, true, ProbsOIJ)); 
          }
          if (Finite.size() == 0) XiOIJ = Xi(Index, j);
        }
        Xi(Index, j) = XiOIJ;
        
        //Update Lambda
        Lambda(Index, j) = Theta(XiOIJ, j);
      
      //End loop over locations  
      }
     
    //End loop over observations 
    }
    
  //End loop over columns
  }
  
  //Update parameters object
  Para.Xi = Xi;
  Para.Lambda = Lambda;
  Para.Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;
  return Para;
  
}


//Function to sample theta_jl's using a Gibbs sampler step---------------------------------------------------------------
para SampleTheta(datobj DatObj, para Para) {
  
  //Set data objects
  int K = DatObj.K;
  int L = DatObj.L;
  int LInf = DatObj.LInf;
  int M = DatObj.M;
  int O = DatObj.O;
  int Nu = DatObj.Nu;
  arma::mat YStarWide = DatObj.YStarWide;
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat EyeO = DatObj.EyeO;
  arma::colvec OneO = DatObj.OneO;
  
  //Set parameter objects
  arma::mat BigPhi = Para.BigPhi;
  arma::umat Xi = Para.Xi;
  arma::mat Lambda = Para.Lambda;
  arma::cube Cov = Para.Cov;
  arma::mat CovMat = arma::reshape(arma::vectorise(Cov), M * O, Nu);
  arma::vec Tau = Para.Tau;
  arma::mat Theta = Para.Theta;
  arma::colvec Eta = Para.Eta;
  arma::colvec LStarJ = Para.LStarJ;
  arma::mat Upsilon = Para.Upsilon;
  arma::colvec XBeta = Para.XBeta;
  arma::mat XBetaMat = arma::reshape(XBetaMat, M * O, Nu);
  
  //L to loop over
  int UpperL = L;
  
  //Loop over columns K
  for (arma::uword j = 0; j < K; j++) {
    
    //Column specific objects
    arma::colvec EtaJ = arma::trans(BigPhi.row(j));
    arma::uvec XiJ = Xi.col(j);
    double TauJ = arma::as_scalar(Tau.row(j));
    arma::mat LambdaMinusJ = Lambda, BigPhiMinusJ = BigPhi;
    LambdaMinusJ.shed_col(j);
    BigPhiMinusJ.shed_row(j);
    if (LInf == 1) UpperL = LStarJ(j); 
    
    //Loop over clusters L
    for (arma::uword l = 0; l < UpperL; l++) {
        
      //Number of objects in cluster l
      arma::uvec WhichJL = find(XiJ == l);
      int NJL = WhichJL.size();
      
      //Declarations
      double ThetaJL;
      
      //Sample from prior if none belong to cluster l
      if (NJL == 0) ThetaJL = arma::as_scalar(rnormRcpp(1, 0, sqrt(1 / TauJ)));
      
      //Sample from posterior if at least 1 belongs to cluster l
      if (NJL > 0) {
        
        //Cluster l specific objects
        arma::mat CovInvL = 1 / CovMat.rows(WhichJL);
        arma::mat YJL = YStarWide.rows(WhichJL);
        arma::mat LambdaJL = LambdaMinusJ.rows(WhichJL);
        arma::mat XBetaComp = XBetaMat.rows(WhichJL);
        
        //Sample theta
        double ResidsJL = arma::as_scalar(arma::sum((YJL - XBetaComp - LambdaJL * BigPhiMinusJ) % CovInvL, 0) * EtaJ);
        double VarThetaJL = arma::as_scalar(1 / (arma::sum(CovInvL, 0) * (EtaJ % EtaJ) + TauJ));
        double MeanThetaJL = VarThetaJL * ResidsJL;
        ThetaJL = arma::as_scalar(rnormRcpp(1, MeanThetaJL, sqrt(VarThetaJL)));
        
      }

      //Update parameters
      Theta(l, j) = ThetaJL;
      arma::colvec LambdaJ = Lambda.col(j);
      arma::vec ThetaJLVec(1);
      ThetaJLVec(0) = ThetaJL;
      LambdaJ(WhichJL) = arma::repmat(ThetaJLVec, NJL, 1);
      Lambda.col(j) = LambdaJ;

    }
  }
  
  //Final updates
  arma::colvec Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;

  //Update parameters object
  Para.Theta = Theta;
  Para.Lambda = Lambda;
  Para.Mean = Mean;
  return Para;
  
}



//Function to sample Latent U using a Gibbs sampler step-------------------------------------------------------------------
para SampleU(datobj DatObj, para Para) {
  
  //Set data objects
  int M = DatObj.M;
  int K = DatObj.K;
  int O = DatObj.O;
  
  //Set parameters
  arma::cube Weights = Para.Weights;
  arma::umat Xi = Para.Xi;
  
  //Sample latent U
  arma::mat U(M * O, K);
  arma::uword Index;
  for (arma::uword j = 0; j < K; j++) {
    arma::mat WeightsJ = Weights.slice(j);
    for (arma::uword o = 0; o < O; o++) {
      for (arma::uword i = 0; i < M; i++) {
        Index = i + M * o;
        double WeightsOIJ = WeightsJ(Xi(Index, j), Index);
        U(Index, j) = randuRcpp() * WeightsOIJ;
      }
    }
  }
  
  //Update upper bounds LStarJ
  arma::colvec LStarJ = GetLStarJ(U, Weights, K, M, O);
  
  //Update parameters object
  Para.U = U;
  Para.LStarJ = LStarJ;
  return Para;
}
