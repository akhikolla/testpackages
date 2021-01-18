#include <RcppArmadillo.h>
#include "MCMC_STBDwDM.h"

//Function to sample delta using a Gibbs sampler step---------------------------------------------------------------
para SampleDelta(datobj DatObj, para Para, hypara HyPara) {

  //Set data objects
  arma::mat ZDelta = DatObj.ZDelta;
  arma::vec OneNu = DatObj.OneNu;

  //Set parameters
  arma::vec VecTheta = Para.VecTheta;
  arma::mat CovThetaInv = Para.CovThetaInv;

  //Set hyperparameter objects
  arma::mat OmegaDeltaInv = HyPara.OmegaDeltaInv;
  arma::vec OmegaDeltaInvMuDelta = HyPara.OmegaDeltaInvMuDelta;

  //Sample delta
  arma::mat CovDelta = CholInv(arma::trans(ZDelta) * CovThetaInv * ZDelta + OmegaDeltaInv);
  arma::vec MeanDelta = CovDelta * (arma::trans(ZDelta) * CovThetaInv * VecTheta + OmegaDeltaInvMuDelta);
  arma::vec Delta = rmvnormRcpp(1, MeanDelta, CovDelta);

  //Update parameters dependent on delta
  arma::mat MMat = Delta * arma::trans(OneNu);
  arma::vec MeanTheta = ZDelta * Delta;

  //Update parameters object
  Para.Delta = Delta;
  Para.MeanTheta = MeanTheta;
  Para.MMat = MMat;
  return Para;
}



//Function to sample new value of phi using a Metropolis sampler step-----------------------------------------------
std::pair<para, metrobj> SamplePhi(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj) {

  //Set data objects
  arma::mat TimeDist = DatObj.TimeDist;
  int TempCorInd = DatObj.TempCorInd;
  arma::mat EyeNTheta = DatObj.EyeNTheta;
  int Nu = DatObj.Nu;

  //Set hyperparameter objects
  double APhi = HyPara.APhi;
  double BPhi = HyPara.BPhi;

  //Set parameter objects
  arma::vec VecTheta = Para.VecTheta;
  arma::vec MeanTheta = Para.MeanTheta;
  arma::mat RootiTheta = Para.RootiTheta;
  double Phi = Para.Phi;
  arma::mat T = Para.T;
  arma::mat TInv = Para.TInv;

  //Set metropolis objects
  double MetropPhi = sqrt(MetrObj.MetropPhi);
  double AcceptancePhi = MetrObj.AcceptancePhi;

  //Transform current state to real line
  double BigDelta = log((Phi - APhi) / (BPhi - Phi ));

  //Sample a new Proposal
  double BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropPhi));

  //Compute Phi Proposal
  double PhiProposal = (BPhi * exp(BigDeltaProposal) + APhi) / (1 + exp(BigDeltaProposal));

  // //Fix numerical issue where PhiProposal can equal APhi or BPhi
  // arma::vec PhiProposalVec(1), APhiVec(1), BPhiVec(1);
  // PhiProposalVec = PhiProposal;
  // APhiVec = APhi;
  // BPhiVec = BPhi;
  // double TOL = 0.000001;
  // if ((rows_equal(PhiProposalVec, APhiVec, TOL)) || (rows_equal(PhiProposalVec, BPhiVec, TOL))) {
  //   if (rows_equal(PhiProposalVec, APhiVec, TOL)) PhiProposal *= 1.1;
  //   if (rows_equal(PhiProposalVec, BPhiVec, TOL)) PhiProposal *= 0.99;
  //   BigDeltaProposal = log( (PhiProposal - APhi) / (BPhi - PhiProposal) );
  // }

  //Proposal temporal correlation
  arma::mat SIGMAPhiProposal = SIGMA(PhiProposal, TempCorInd, TimeDist, Nu);

  //Theta structure components
  arma::mat CovThetaProposal = arma::kron(SIGMAPhiProposal, T);
  arma::mat RootiThetaProposal = GetRooti(CovThetaProposal, EyeNTheta);
  double Component1A = lndMvn(VecTheta, MeanTheta, RootiThetaProposal);
  double Component1B = lndMvn(VecTheta, MeanTheta, RootiTheta);
  double Component1 = Component1A - Component1B;

  //Prior components
  double Component2A = BigDeltaProposal;
  double Component2B = BigDelta;
  double Component2 = Component2A - Component2B;

  //Jacobian
  double Component3 = 2 * log( (1 + exp(BigDelta)) / (1 + exp(BigDeltaProposal)) );

  //Compute log acceptance ratio
  double LogR = Component1 + Component2 + Component3;

  //Metropolis update
  double RandU = randuRcpp();
  if (log(RandU) < LogR) {

    //Keep Count of Acceptances
    AcceptancePhi++;
    MetrObj.AcceptancePhi = AcceptancePhi;

    //Update parameters dependent on phi
    arma::mat SIGMAPhiInv = CholInv(SIGMAPhiProposal);
    arma::mat CovThetaInv = arma::kron(SIGMAPhiInv, TInv);

    //Update parameters object
    Para.Phi = PhiProposal;
    Para.SIGMAPhi = SIGMAPhiProposal;
    Para.SIGMAPhiInv = SIGMAPhiInv;
    Para.CovThetaInv = CovThetaInv;
    Para.RootiTheta = RootiThetaProposal;

  }

  //Return output object
  return std::pair<para, metrobj>(Para, MetrObj);

}



//Function to sample T using a Gibbs sampler step-------------------------------------------------------------------
para SampleT(datobj DatObj, para Para, hypara HyPara) {

  //Set data objects
  int Nu = DatObj.Nu;
  arma::vec OneNu = DatObj.OneNu;
  arma::mat EyeNTheta = DatObj.EyeNTheta;

  //Set parameters
  arma::mat Theta = Para.Theta;
  arma::mat SIGMAPhiInv = Para.SIGMAPhiInv;
  arma::mat SIGMAPhi = Para.SIGMAPhi;
  arma::mat MMat = Para.MMat;

  //Set hyperparameter objects
  double Xi = HyPara.Xi;
  arma::mat Psi = HyPara.Psi;

  //Compute STheta
  arma::mat ThetaDiff = Theta - MMat;
  arma::mat STheta = ThetaDiff * SIGMAPhiInv * arma::trans(ThetaDiff);

  //Sample T
  double n = Xi + Nu;
  arma::mat V = STheta + Psi;
  arma::mat TInv = rwishRcpp(n, Inv3(V));
  arma::mat T = Inv3(TInv);

  //Update parameters dependent on T
  arma::mat CovThetaInv = arma::kron(SIGMAPhiInv, TInv);
  arma::mat CovTheta = arma::kron(SIGMAPhi, T);
  arma::mat RootiTheta = GetRooti(CovTheta, EyeNTheta);

  //Update parameters object
  Para.T = T;
  Para.TInv = TInv;
  Para.CovThetaInv = CovThetaInv;
  Para.RootiTheta = RootiTheta;
  return Para;
}



//Function to sample new value of theta1 (i.e. Mu) using a Gibbs sampler step--------------------------------------
para SampleTheta1(datobj DatObj, para Para) {

  //Set data objects
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat EyeM = DatObj.EyeM;
  double Rho = DatObj.Rho;
  arma::vec OneM = DatObj.OneM;
  int Nu = DatObj.Nu;
  arma::mat YStarWide = DatObj.YStarWide;
  int M = DatObj.M;

  //Set parameter objects
  arma::mat T = Para.T;
  arma::mat SIGMAPhiInv = Para.SIGMAPhiInv;
  arma::mat MMat = Para.MMat;
  arma::mat Theta = Para.Theta;
  arma::cube WAlphas = Para.WAlphas;
  arma::vec Tau2 = Para.Tau2;

  //Get conditional moments
  arma::uvec j(1); j(0) = 0;
  arma::uvec k(2); k(0) = 1, k(1) = 2;
  arma::rowvec Tjk = T(j, k);
  arma::rowvec TPlus = Tjk * Inv2(T(k, k));
  double TStar = arma::as_scalar(T(j, j) - TPlus * arma::trans(Tjk));
  arma::mat CondCovInv = SIGMAPhiInv / TStar;
  arma::mat ThetaDiff = Theta.rows(k) - MMat.rows(k);
  arma::vec CondMean = arma::trans(MMat.rows(j)) + arma::kron(TPlus, EyeNu) * arma::vectorise(arma::trans(ThetaDiff));

  //Initialized likelihood objects
  arma::cube JointPrecisions = JointPrecisionCube(WAlphas, Tau2, EyeM, Rho, M, Nu); //Q(alpha_t)/tau_t^2
  arma::mat ZThetaT_Q_ZTheta = arma::eye(Nu, Nu);
  arma::vec ZThetaT_Q_YStar(Nu);
  arma::cube CubeTemp1A(M, 1, Nu), CubeTemp1B(1, 1, Nu), CubeTemp2A(1, M, Nu);

  //Compute covariance object
  CubeTemp1A = JointPrecisions.each_slice() * OneM;
  CubeTemp1B = arma::trans(OneM) * CubeTemp1A.each_slice();
  for (int i = 0; i < Nu; i++) ZThetaT_Q_ZTheta(i, i) = arma::as_scalar(CubeTemp1B.slice(i));

  //Compute mean object
  CubeTemp2A = arma::trans(OneM) * JointPrecisions.each_slice();
  for (int i = 0; i < Nu; i++)
    ZThetaT_Q_YStar(i) = arma::as_scalar(CubeTemp2A.slice(i) * YStarWide.col(i));

  //Get full conditional moments
  arma::mat CovTheta1 = CholInv(ZThetaT_Q_ZTheta + CondCovInv);
  arma::vec MeanTheta1 = CovTheta1 * (ZThetaT_Q_YStar + CondCovInv * CondMean);

  //Sample from full conditional
  arma::vec Mu = rmvnormRcpp(1, MeanTheta1, CovTheta1);

  //Update parameters dependent on theta1
  Theta.rows(j) = arma::trans(Mu);
  arma::vec VecTheta = arma::vectorise(Theta);

  //Return parameters object
  Para.Mu = Mu;
  Para.Theta = Theta;
  Para.VecTheta = VecTheta;
  return Para;

}



//Function to sample new value of theta2 (i.e. tau2) using a Metropolis sampler step--------------------------------
std::pair<para, metrobj> SampleTheta2(datobj DatObj, para Para, metrobj MetrObj) {

  //Set data objects
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat YStarWide = DatObj.YStarWide;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  double Rho = DatObj.Rho;
  int Nu = DatObj.Nu;
  int M = DatObj.M;

  //Set Metropolis Tuning Objects
  arma::vec MetropTheta2 = MetrObj.MetropTheta2;
  arma::vec AcceptanceTheta2 = MetrObj.AcceptanceTheta2;

  //Set parameter objects
  arma::mat T = Para.T;
  arma::mat SIGMAPhi = Para.SIGMAPhi;
  arma::mat MMat = Para.MMat;
  arma::mat Theta = Para.Theta;
  arma::vec Mu = Para.Mu;
  arma::vec Tau2 = Para.Tau2;
  arma::cube WAlphas = Para.WAlphas;
  arma::cube RootiLikelihoods = Para.RootiLikelihoods;
  arma::cube JointCovariances = Para.JointCovariances;

  //Get conditional moments
  arma::uvec j(1); j(0) = 1;
  arma::uvec k(2); k(0) = 0, k(1) = 2;
  arma::rowvec Tjk = T(j, k);
  arma::rowvec TPlus = Tjk * Inv2(T(k, k));
  double TStar = arma::as_scalar(T(j, j) - TPlus * arma::trans(Tjk));
  arma::mat CondRooti = GetRooti(TStar * SIGMAPhi, EyeNu);
  arma::mat ThetaDiff = Theta.rows(k) - MMat.rows(k);
  arma::vec CondMean = arma::trans(MMat.rows(j)) + arma::kron(TPlus, EyeNu) * arma::vectorise(arma::trans(ThetaDiff));

  //Initialize objects
  arma::uvec VisitInd(1);
  arma::mat WAlphaVisit(M, M), ThetaProposal(3, Nu), JointCovarianceProposal(M, M);
  arma::mat RootiLikelihoodVisit(M, M), RootiLikelihoodProposal(M, M);
  arma::vec YStarWideVisit(M), ProposalTheta(1), MeanLikelihood(M), RandU(1);
  double MuVisit, CurrentTheta, Tau2VisitProposal, LogR;
  double Component1A, Component1B, Component1;
  double Component2A, Component2B, Component2;
  double TuningSD;

  //Loop over visits
  for (int Visit = 0; Visit < Nu; Visit++) {

    //Visit specific objects
    VisitInd(0) = Visit;
    MuVisit = Mu(Visit);
    WAlphaVisit = WAlphas.slice(Visit);
    RootiLikelihoodVisit = RootiLikelihoods.slice(Visit);
    YStarWideVisit = YStarWide.col(Visit);
    CurrentTheta = arma::as_scalar(Theta(j, VisitInd));
    TuningSD = sqrt(MetropTheta2(Visit));

    //Numerical fix for when a propopsal is rounded to 0 or Inf
    Tau2VisitProposal = arma::datum::inf;
    ThetaProposal = Theta;
    while ( (!arma::is_finite(Tau2VisitProposal)) || (Tau2VisitProposal == 0) ) {

      //Sample proposal
      ProposalTheta = rnormRcpp(1, CurrentTheta, TuningSD);
      ThetaProposal(j, VisitInd) = ProposalTheta;
      Tau2VisitProposal = arma::as_scalar(arma::square(arma::exp(ProposalTheta))); //verified this step on 2/15/17

    }

    //Get proposal JointCovariance
    JointCovarianceProposal = JointCovarianceMatrix(WAlphaVisit, Tau2VisitProposal, EyeM, Rho, M);

    //Likelihood Component
    RootiLikelihoodProposal = GetRooti(JointCovarianceProposal, EyeM);
    MeanLikelihood = MuVisit * OneM;
    Component1A = lndMvn(YStarWideVisit, MeanLikelihood, RootiLikelihoodProposal);
    Component1B = lndMvn(YStarWideVisit, MeanLikelihood, RootiLikelihoodVisit);
    Component1 = Component1A - Component1B;

    //Prior components
    Component2A = lndMvn(arma::trans(ThetaProposal.rows(j)), CondMean, CondRooti);
    Component2B = lndMvn(arma::trans(Theta.rows(j)), CondMean, CondRooti);
    Component2 = Component2A - Component2B;

    //Log acceptance ratio
    LogR = Component1 + Component2;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceTheta2(Visit)++;

      //Update parameters output
      Tau2(Visit) = Tau2VisitProposal;
      Theta = ThetaProposal;
      RootiLikelihoods.slice(Visit) = RootiLikelihoodProposal;
      JointCovariances.slice(Visit) = JointCovarianceProposal;

    }

    //End loop over visits
  }

  //Update Metropolis object
  MetrObj.AcceptanceTheta2 = AcceptanceTheta2;

  //Update Para objects
  Para.Tau2 = Tau2;
  Para.Theta = Theta;
  arma::vec VecTheta = arma::vectorise(Theta);
  Para.VecTheta = VecTheta;
  Para.RootiLikelihoods = RootiLikelihoods;
  Para.JointCovariances = JointCovariances;

  //Return final object
  return std::pair<para, metrobj>(Para, MetrObj);

}



//Function to sample new value of theta3 (i.e. alpha) using a Metropolis sampler step------------------------------
std::pair<para, metrobj> SampleTheta3(datobj DatObj, para Para, metrobj MetrObj) {

  //Set data objects
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat YStarWide = DatObj.YStarWide;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  double Rho = DatObj.Rho;
  int Nu = DatObj.Nu;
  int M = DatObj.M;
  int WeightsInd = DatObj.WeightsInd;
  arma::vec Z = DatObj.Z;
  arma::umat AdjacentEdgesBoolean = DatObj.AdjacentEdgesBoolean;
  arma::mat W = DatObj.W;

  //Set Metropolis Tuning Objects
  arma::vec MetropTheta3 = MetrObj.MetropTheta3;
  arma::vec AcceptanceTheta3 = MetrObj.AcceptanceTheta3;

  //Set parameter objects
  arma::mat T = Para.T;
  arma::mat SIGMAPhi = Para.SIGMAPhi;
  arma::mat MMat = Para.MMat;
  arma::mat Theta = Para.Theta;
  arma::cube WAlphas = Para.WAlphas;
  arma::cube RootiLikelihoods = Para.RootiLikelihoods;
  arma::cube JointCovariances = Para.JointCovariances;
  arma::vec Mu = Para.Mu;
  arma::vec Tau2 = Para.Tau2;
  arma::vec Alpha = Para.Alpha;

  //Get conditional moments
  arma::uvec j(1); j(0) = 2;
  arma::uvec k(2); k(0) = 0, k(1) = 1;
  arma::rowvec Tjk = T(j, k);
  arma::rowvec TPlus = Tjk * Inv2(T(k, k));
  double TStar = arma::as_scalar(T(j, j) - TPlus * arma::trans(Tjk));
  arma::mat CondRooti = GetRooti(TStar * SIGMAPhi, EyeNu);
  arma::mat ThetaDiff = Theta.rows(k) - MMat.rows(k);
  arma::vec CondMean = arma::trans(MMat.rows(j)) + arma::kron(TPlus, EyeNu) * arma::vectorise(arma::trans(ThetaDiff));

  //Initialize objects
  arma::uvec VisitInd(1);
  arma::mat WAlphaProposal(M, M), ThetaProposal(3, Nu), JointCovarianceProposal(M, M);
  arma::mat RootiLikelihoodVisit(M, M), RootiLikelihoodProposal(M, M);
  arma::vec YStarWideVisit(M), ProposalTheta(1), MeanLikelihood(M), RandU(1);
  double MuVisit, Tau2Visit, CurrentTheta, AlphaVisitProposal, LogR;
  double Component1A, Component1B, Component1;
  double Component2A, Component2B, Component2;
  double TuningSD;

  //Loop over visits
  for (int Visit = 0; Visit < Nu; Visit++) {

    //Visit specific objects
    VisitInd(0) = Visit;
    MuVisit = Mu(Visit);
    Tau2Visit = Tau2(Visit);
    YStarWideVisit = YStarWide.col(Visit);
    RootiLikelihoodVisit = RootiLikelihoods.slice(Visit);
    CurrentTheta = arma::as_scalar(Theta(j, VisitInd));
    TuningSD = sqrt(MetropTheta3(Visit));

    //Numerical fix for when a propopsal is rounded to 0 or Inf
    AlphaVisitProposal = arma::datum::inf;
    ThetaProposal = Theta;
    while ( (!arma::is_finite(AlphaVisitProposal)) || (AlphaVisitProposal == 0) ) {

      //Sample proposal
      ProposalTheta = rnormRcpp(1, CurrentTheta, TuningSD);
      ThetaProposal(j, VisitInd) = ProposalTheta;
      AlphaVisitProposal = arma::as_scalar(arma::exp(ProposalTheta));

    }

    //Get proposal JointCovariance
    WAlphaProposal = WAlphaMatrix(AlphaVisitProposal, Z, AdjacentEdgesBoolean, W, M, WeightsInd);
    JointCovarianceProposal = JointCovarianceMatrix(WAlphaProposal, Tau2Visit, EyeM, Rho, M);

    //Likelihood Component
    RootiLikelihoodProposal = GetRooti(JointCovarianceProposal, EyeM);
    MeanLikelihood = MuVisit * OneM;
    Component1A = lndMvn(YStarWideVisit, MeanLikelihood, RootiLikelihoodProposal);
    Component1B = lndMvn(YStarWideVisit, MeanLikelihood, RootiLikelihoodVisit);
    Component1 = Component1A - Component1B;

    //Prior components
    Component2A = lndMvn(arma::trans(ThetaProposal.rows(j)), CondMean, CondRooti);
    Component2B = lndMvn(arma::trans(Theta.rows(j)), CondMean, CondRooti);
    Component2 = Component2A - Component2B;

    //Log acceptance ratio
    LogR = Component1 + Component2;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceTheta3(Visit)++;

      //Update parameters output
      Alpha(Visit) = AlphaVisitProposal;
      WAlphas.slice(Visit) = WAlphaProposal;
      RootiLikelihoods.slice(Visit) = RootiLikelihoodProposal;
      JointCovariances.slice(Visit) = JointCovarianceProposal;
      Theta = ThetaProposal;

    }

    //End loop over visits
  }

  //Update Metropolis object
  MetrObj.AcceptanceTheta3 = AcceptanceTheta3;

  //Update Para objects
  Para.Alpha = Alpha;
  Para.WAlphas = WAlphas;
  Para.Theta = Theta;
  arma::vec VecTheta = arma::vectorise(Theta);
  Para.VecTheta = VecTheta;
  Para.RootiLikelihoods = RootiLikelihoods;
  Para.JointCovariances = JointCovariances;

  //Return final object
  return std::pair<para, metrobj>(Para, MetrObj);

}



//Function to sample latent probit process using Gibbs sampling step------------------------------------------------
arma::mat SampleProbit(datobj DatObj, para Para, dataug DatAug) {

  //Set data objects
  double Rho = DatObj.Rho;
  arma::mat YStarWide = DatObj.YStarWide;
  arma::vec OneM = DatObj.OneM;

  //Set parameters
  arma::vec Mu = Para.Mu;
  arma::vec Tau2 = Para.Tau2;
  arma::cube WAlphas = Para.WAlphas;

  //Set data augmentation objects
  int NBelow = DatAug.NBelow;
  int NAbove = DatAug.NAbove;
  arma::mat TobitIndeces = DatAug.TobitIndeces;
  arma::mat ProbitIndeces = DatAug.ProbitIndeces;

  //Loop over latent observations with Y = 0
  arma::rowvec Row(2);
  int Location, Visit;
  for (int i = 0; i < NBelow; i++) {

    //Determine location and visit
    Row = TobitIndeces.row(i);
    Location = Row(0);
    Visit = Row(1);

    //Compute visit specific parameters
    double mu = Mu(Visit);
    double tau2 = Tau2(Visit);
    arma::mat WAlpha = WAlphas.slice(Visit);
    arma::colvec YStarVisit = YStarWide.col(Visit);

    //Compute conditional mean and variance
    arma::rowvec Wij = WAlpha.row(Location);
    double BigDelta = arma::as_scalar(Rho * (Wij * OneM) + (1 - Rho) );
    double ConditionalSD = sqrt(tau2 / BigDelta);
    double ConditionalMean = arma::as_scalar(Rho * (Wij * YStarVisit) + (1 - Rho) * mu) / BigDelta;

    //Sample latent Variable from full conditional
    double Temp = rtnormRcpp(ConditionalMean, ConditionalSD, true);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(ConditionalMean, ConditionalSD, -arma::datum::inf, 0);
    YStarWide(Location, Visit) = Temp;

    //End loop over truncated observations
  }

  //Loop over latent observations with Y = 1
  for (int i = 0; i < NAbove; i++) {

    //Determine location and visit
    Row = ProbitIndeces.row(i);
    Location = Row(0);
    Visit = Row(1);

    //Compute visit specific parameters
    double mu = Mu(Visit);
    double tau2 = Tau2(Visit);
    arma::mat WAlpha = WAlphas.slice(Visit);
    arma::colvec YStarVisit = YStarWide.col(Visit);

    //Compute conditional mean and variance
    arma::rowvec Wij = WAlpha.row(Location);
    double BigDelta = arma::as_scalar(Rho * (Wij * OneM) + (1 - Rho) );
    double ConditionalSD = sqrt(tau2 / BigDelta);
    double ConditionalMean = arma::as_scalar(Rho * (Wij * YStarVisit) + (1 - Rho) * mu) / BigDelta;

    //Sample latent Variable from full conditional
    double Temp = rtnormRcpp(ConditionalMean, ConditionalSD, true);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(ConditionalMean, ConditionalSD, 0, arma::datum::inf);
    YStarWide(Location, Visit) = Temp;

    //End loop over truncated observations
  }
  return YStarWide;

}



//Function to sample latent tobit process using Gibbs sampling step------------------------------------------------
arma::mat SampleTobit(datobj DatObj, para Para, dataug DatAug) {

  //Set data objects
  double Rho = DatObj.Rho;
  arma::mat YStarWide = DatObj.YStarWide;
  arma::vec OneM = DatObj.OneM;

  //Set parameters
  arma::vec Mu = Para.Mu;
  arma::vec Tau2 = Para.Tau2;
  arma::cube WAlphas = Para.WAlphas;

  //Set data augmentation objects
  int NBelow = DatAug.NBelow;
  arma::mat TobitIndeces = DatAug.TobitIndeces;

  //Loop over latent observations
  arma::rowvec Row(2);
  int Location, Visit;
  for (int i = 0; i < NBelow; i++) {

    //Determine location and visit
    Row = TobitIndeces.row(i);
    Location = Row(0);
    Visit = Row(1);

    //Compute visit specific parameters
    double mu = Mu(Visit);
    double tau2 = Tau2(Visit);
    arma::mat WAlpha = WAlphas.slice(Visit);
    arma::colvec YStarVisit = YStarWide.col(Visit);

    //Compute conditional mean and variance
    arma::rowvec Wij = WAlpha.row(Location);
    double BigDelta = arma::as_scalar(Rho * (Wij * OneM) + (1 - Rho) );
    double ConditionalSD = sqrt(tau2 / BigDelta);
    double ConditionalMean = arma::as_scalar(Rho * (Wij * YStarVisit) + (1 - Rho) * mu) / BigDelta;

    //Sample latent Variable from full conditional
    double Temp = rtnormRcpp(ConditionalMean, ConditionalSD, true);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(ConditionalMean, ConditionalSD, -arma::datum::inf, 0);
    // if (!arma::is_finite(Temp)) Rcpp::stop("infinte value sampled in Tobit sampling step. Most likey cause for this error is that the data being used is inappropriate (i.e., to far from zero) for a Tobit model. Consider scaling towards zero and re-running.");
    YStarWide(Location, Visit) = Temp;

    //End loop over truncated observations
  }
  return YStarWide;

}



//Function to sample latent process from its full conditional------------------------------------------------------
datobj SampleY(datobj DatObj, para Para, dataug DatAug) {

  //Set data objects
  int FamilyInd = DatObj.FamilyInd;
  int M = DatObj.M;
  int Nu = DatObj.Nu;

  //Sample latent process
  arma::mat YStarWide(M, Nu);
  if (FamilyInd == 1) YStarWide = SampleProbit(DatObj, Para, DatAug);
  if (FamilyInd == 2) YStarWide = SampleTobit(DatObj, Para, DatAug);

  //Save output
  DatObj.YStarWide = YStarWide;
  return DatObj;

}

