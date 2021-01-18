#include <RcppArmadillo.h>
#include "MCMC_spCP.h"

//Function to sample delta using a Gibbs sampler step---------------------------------------------------------------
para SampleDelta(datobj DatObj, para Para, hypara HyPara) {

  //Set data objects
  arma::mat Eye5 = DatObj.Eye5;
  arma::mat ZDelta = DatObj.ZDelta;
  arma::mat OneM = DatObj.OneM;

  //Set parameters
  arma::vec Phi = Para.Phi;
  arma::mat PhiPrec = Para.PhiPrec;
  arma::mat XTheta = Para.XTheta;

  //Set hyperparameter objects
  double Kappa2 = HyPara.Kappa2;

  //Sample delta
  arma::mat tZDelta = arma::trans(ZDelta);
  arma::mat CovDelta = CholInv(tZDelta * PhiPrec * ZDelta + Eye5 / Kappa2);
  arma::vec MeanDelta = CovDelta * (tZDelta * PhiPrec * Phi);
  arma::vec Delta = rmvnormRcpp(1, MeanDelta, CovDelta);

  //Update parameters dependent on delta
  arma::vec PhiMean = arma::kron(OneM, Delta);

  //Update parameters object
  Para.Delta = Delta;
  Para.PhiMean = PhiMean;
  return Para;
}



//Function to sample new value of alpha using a Metropolis sampler step-----------------------------------------------
std::pair<para, metrobj> SampleAlpha(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj) {

  //Set data objects
  int M = DatObj.M;
  int WeightsInd = DatObj.WeightsInd;
  arma::vec DMLong = DatObj.DMLong;
  arma::umat AdjacentEdgesBoolean = DatObj.AdjacentEdgesBoolean;
  arma::Mat<int> W = DatObj.W;
  arma::mat EyeM = DatObj.EyeM;
  arma::mat Eye5M = DatObj.Eye5M;
  double Rho = DatObj.Rho;

  //Set hyperparameter objects
  double AAlpha = HyPara.AAlpha;
  double BAlpha = HyPara.BAlpha;

  //Set parameter objects
  double Alpha = Para.Alpha;
  arma::mat Phi = Para.Phi;
  arma::vec PhiMean = Para.PhiMean;
  arma::mat Sigma = Para.Sigma;
  arma::mat QInv = Para.QInv;

  //Set metropolis objects
  double MetropAlpha = sqrt(MetrObj.MetropAlpha);
  double AcceptanceAlpha = MetrObj.AcceptanceAlpha;

  //Transform current state to real line
  double BigDelta = log((Alpha - AAlpha) / (BAlpha - Alpha ));

  //Sample a new Proposal
  double BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropAlpha));

  //Compute Phi Proposal
  double AlphaProposal = (BAlpha * exp(BigDeltaProposal) + AAlpha) / (1 + exp(BigDeltaProposal));

  //Fix numerical issue where AlphaProposal can equal AAlpha or BAlpha
  arma::vec AlphaProposalVec(1), AAlphaVec(1), BAlphaVec(1);
  AlphaProposalVec = AlphaProposal;
  AAlphaVec = AAlpha;
  BAlphaVec = BAlpha;
  double TOL = 0.000001;
  if ((rows_equal(AlphaProposalVec, AAlphaVec, TOL)) || (rows_equal(AlphaProposalVec, BAlphaVec, TOL))) {
    if (rows_equal(AlphaProposalVec, AAlphaVec, TOL)) AlphaProposal *= 1.1;
    if (rows_equal(AlphaProposalVec, BAlphaVec, TOL)) AlphaProposal *= 0.99;
    BigDeltaProposal = log( (AlphaProposal - AAlpha) / (BAlpha - AlphaProposal) );
  }

  //Proposal alpha objects
  arma::mat WAlphaProposal =  WAlphaFnc(AlphaProposal, DMLong, AdjacentEdgesBoolean, W, M, WeightsInd);
  arma::mat QProposal = QFnc(WAlphaProposal, EyeM, Rho, M);
  arma::mat QInvProposal = CholInv(QProposal);

  //Spatial components
  arma::mat CholSigma = arma::chol(Sigma);
  arma::mat RootiPhi = arma::solve(arma::trimatu(arma::kron(arma::chol(QInv), CholSigma)), Eye5M);
  arma::mat RootiPhiProposal = arma::solve(arma::trimatu(arma::kron(arma::chol(QInvProposal), CholSigma)), Eye5M);
  double Component1A = lndMvn(Phi, PhiMean, RootiPhiProposal);
  double Component1B = lndMvn(Phi, PhiMean, RootiPhi);
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
    AcceptanceAlpha++;
    MetrObj.AcceptanceAlpha = AcceptanceAlpha;

    //Update parameters object
    Para.Alpha = AlphaProposal;
    Para.WAlpha = WAlphaProposal;
    Para.QInv = QInvProposal;
    Para.Q = QProposal;
    Para.PhiPrec = arma::kron(QProposal, Para.SigmaInv);
    Para.PhiCov = arma::kron(QInvProposal, Sigma);

  }

  //Return output object
  return std::pair<para, metrobj>(Para, MetrObj);

}



//Function to sample T using a Gibbs sampler step-------------------------------------------------------------------
para SampleSigma(datobj DatObj, para Para, hypara HyPara) {

  //Set data objects
  int M = DatObj.M;
  arma::vec OneM = DatObj.OneM;

  //Set parameters
  arma::vec Delta = Para.Delta;
  arma::vec Phi = Para.Phi;
  arma::mat PhiMatrix = arma::reshape(Phi, 5, M);
  arma::mat Q = Para.Q;
  arma::mat QInv = Para.QInv;

  //Set hyperparameter objects
  double Xi = HyPara.Xi;
  arma::mat Psi = HyPara.Psi;

  //Compute SThetaKappa
  arma::mat MMatrix = Delta * arma::trans(OneM);
  arma::mat PhiDiff = PhiMatrix - MMatrix;
  arma::mat STheta = PhiDiff * Q * arma::trans(PhiDiff);

  //Sample T
  double n = Xi + M;
  arma::mat V = STheta + Psi;
  arma::mat SigmaInv = rwishRcpp(n, CholInv(V));
  arma::mat Sigma = CholInv(SigmaInv);

  //Update parameters object
  Para.Sigma = Sigma;
  Para.SigmaInv = SigmaInv;
  Para.PhiPrec = arma::kron(Q, SigmaInv);
  Para.PhiCov = arma::kron(QInv, Sigma);
  return Para;
}



//Function to sample new value of beta (i.e. Beta0 and Beta1) using a Gibbs sampler step--------------------------------------
para SampleBeta(datobj DatObj, para Para) {

  //Set data objects
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::vec YStar = DatObj.YStar;
  int M = DatObj.M;
  arma::umat PhiIndeces = DatObj.PhiIndeces;

  //Set parameter objects
  arma::vec Delta = Para.Delta;
  arma::mat Sigma = Para.Sigma;
  arma::mat XTheta = Para.XTheta;
  arma::mat OmegaInv = Para.OmegaInv;
  arma::vec Lambda = Para.Lambda;
  arma::vec Eta = Para.Eta;
  arma::mat Q = Para.Q;
  arma::vec Phi = Para.Phi;

  //Get conditional moments
  arma::uvec j(2); j(0) = 0, j(1) = 1;
  arma::uvec k(3); k(0) = 2, k(1) = 3, k(2) = 4;
  arma::mat Sigmajk = Sigma(j, k);
  arma::mat SigmaPlus = Sigmajk * Inv3(Sigma(k, k));
  arma::mat SigmaStar = Sigma(j, j) - SigmaPlus * arma::trans(Sigmajk);
  arma::mat CondPrec = arma::kron(Q, Inv2(SigmaStar));
  arma::vec Phi_k = Phi(arma::vectorise(PhiIndeces.rows(k)));
  arma::vec CondMean = arma::kron(OneM, Delta(j)) + arma::kron(EyeM, SigmaPlus) * (Phi_k - arma::kron(OneM, Delta(k)));

  //Get full conditional moments
  arma::mat tXThetaOmegaInv = arma::trans(XTheta) * OmegaInv;
  arma::mat CovBeta = CholInv(tXThetaOmegaInv * XTheta + CondPrec);
  arma::vec MeanBeta = CovBeta * (tXThetaOmegaInv * YStar + CondPrec * CondMean);

  //Sample from full conditional
  arma::vec Beta = rmvnormRcpp(1, MeanBeta, CovBeta);

  //Return parameters object
  Para.Beta = Beta;
  Para.Mu = XTheta * Beta;
  Para.Phi = CreatePhi(Beta, Lambda, Eta, M);
  return Para;

}



//Function to sample new value of lambda0 using a Metropolis sampler step--------------------------------
std::pair<para, metrobj> SampleLambda0(datobj DatObj, para Para, metrobj MetrObj) {

  //Set data objects
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  arma::umat PhiIndeces = DatObj.PhiIndeces;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::mat YStarWide = arma::trans(DatObj.YStarWide);
  arma::mat EyeNu = DatObj.EyeNu;

  //Set Metropolis Tuning Objects
  arma::vec MetropLambda0Vec = MetrObj.MetropLambda0Vec;
  arma::vec AcceptanceLambda0Vec = MetrObj.AcceptanceLambda0Vec;

  //Set parameter objects
  arma::vec Phi = Para.Phi;
  arma::vec Beta = Para.Beta;
  arma::vec Eta = Para.Eta;
  arma::mat Sigma = Para.Sigma;
  arma::vec Delta = Para.Delta;
  arma::mat QInv = Para.QInv;
  arma::vec Sigma2 = Para.Sigma2;
  arma::mat MuMatrix = arma::trans(arma::reshape(Para.Mu, M, Nu));
  arma::vec Lambda = Para.Lambda;
  arma::mat XTheta = Para.XTheta;

  //Get conditional moments
  arma::uvec j(1); j(0) = 2;
  arma::uvec k(4); k(0) = 0, k(1) = 1, k(2) = 3, k(3) = 4;
  arma::rowvec Sigmajk = Sigma(j, k);
  arma::rowvec SigmaPlus = Sigmajk * CholInv(Sigma(k, k));
  double SigmaStar = arma::as_scalar(Sigma(j, j) - SigmaPlus * arma::trans(Sigmajk));
  arma::vec Phi_k = Phi(arma::vectorise(PhiIndeces.rows(k)));
  arma::vec CondMean = arma::kron(OneM, Delta(j)) + arma::kron(EyeM, SigmaPlus) * (Phi_k - arma::kron(OneM, Delta(k)));
  arma::mat CondRooti = GetRooti(SigmaStar * QInv, EyeM);

  //Initialize objects
  arma::uvec LocInd(1), ColIndeces(2), RowIndeces(Nu);
  arma::colvec YStarWideLoc(Nu), MeanLikelihoodLoc(Nu);
  arma::colvec LambdaLoc(2), LambdaLocProposal(2), Lambda0(M), Lambda0Proposal(M), LambdaProposal(2 * M);
  arma::mat tLambda0(M, 2), tLambda0Proposal(M, 2), OmegaLocChol(Nu, Nu), OmegaLocProposalChol(Nu, Nu);
  arma::vec Sigma2LocProposal(Nu), RandU(1);
  arma::mat RootiLikelihood(Nu, Nu), RootiLikelihoodProposal(Nu, Nu);
  double Component1A, Component1B, Component1;
  double Component2A, Component2B, Component2;
  double TuningSD, LogR, Lambda0LocProposal;

  //Indeces
  arma::vec SeqNu(Nu), Seq2(2);
  for (int i = 0; i < Nu; i++) SeqNu(i) = i;
  for (int i = 0; i < 2; i++) Seq2(i) = i;

  //Loop over visits
  for (int Loc = 0; Loc < M; Loc++) {

    //Indeces
    RowIndeces = arma::conv_to<arma::uvec>::from(SeqNu * M + Loc);
    ColIndeces = arma::conv_to<arma::uvec>::from(Seq2 + Loc * 2);

    //Visit specific objects
    LocInd(0) = Loc;
    YStarWideLoc = YStarWide.col(Loc);
    LambdaLoc = Lambda(ColIndeces);
    tLambda0 = arma::trans(arma::reshape(Lambda, 2, M));
    Lambda0 = tLambda0.col(0);
    MeanLikelihoodLoc = MuMatrix.col(Loc);
    OmegaLocChol = arma::diagmat(sqrt(Sigma2(RowIndeces)));
    TuningSD = sqrt(MetropLambda0Vec(Loc));

    //Numerical fix for when a propopsal is rounded to -Inf or Inf
    Lambda0LocProposal = arma::datum::inf;
    LambdaProposal = Lambda;
    LambdaLocProposal = LambdaLoc;
    while ((!arma::is_finite(Lambda0LocProposal))) {

      //Sample proposal
      Lambda0LocProposal = arma::as_scalar(rnormRcpp(1, LambdaLoc(0), TuningSD));
      LambdaLocProposal(0) = Lambda0LocProposal;
      LambdaProposal.rows(ColIndeces) = LambdaLocProposal;
      tLambda0Proposal = arma::trans(arma::reshape(LambdaProposal, 2, M));
      Lambda0Proposal = tLambda0Proposal.col(0);
      Sigma2LocProposal = arma::exp(2 * (XTheta.rows(RowIndeces) * LambdaProposal));
      OmegaLocProposalChol = arma::diagmat(sqrt(Sigma2LocProposal));

    }

    //Likelihood Component
    RootiLikelihoodProposal = arma::solve(arma::trimatu(OmegaLocProposalChol), EyeNu);
    RootiLikelihood = arma::solve(arma::trimatu(OmegaLocChol), EyeNu);
    Component1A = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihoodProposal);
    Component1B = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihood);
    Component1 = Component1A - Component1B;

    //Prior components
    Component2A = lndMvn(Lambda0Proposal, CondMean, CondRooti);
    Component2B = lndMvn(Lambda0, CondMean, CondRooti);
    Component2 = Component2A - Component2B;

    //Log acceptance ratio
    LogR = Component1 + Component2;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceLambda0Vec(Loc)++;

      //Update parameters output
      Lambda = LambdaProposal;
      Sigma2.rows(RowIndeces) = Sigma2LocProposal;

    }

  //End loop over locations
  }

  //Update Metropolis object
  MetrObj.AcceptanceLambda0Vec = AcceptanceLambda0Vec;

  //Update Para objects
  Para.Lambda = Lambda;
  Para.Sigma2 = Sigma2;
  Para.Phi = CreatePhi(Beta, Lambda, Eta, M);
  Para.Omega = arma::diagmat(Sigma2);
  Para.OmegaInv = arma::diagmat(1 / Sigma2);

  //Return final object
  return std::pair<para, metrobj>(Para, MetrObj);

}



//Function to sample new value of lambda1 using a Metropolis sampler step------------------------------
std::pair<para, metrobj> SampleLambda1(datobj DatObj, para Para, metrobj MetrObj) {

  //Set data objects
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  arma::umat PhiIndeces = DatObj.PhiIndeces;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::mat YStarWide = arma::trans(DatObj.YStarWide);
  arma::mat EyeNu = DatObj.EyeNu;

  //Set Metropolis Tuning Objects
  arma::vec MetropLambda1Vec = MetrObj.MetropLambda1Vec;
  arma::vec AcceptanceLambda1Vec = MetrObj.AcceptanceLambda1Vec;

  //Set parameter objects
  arma::vec Phi = Para.Phi;
  arma::vec Beta = Para.Beta;
  arma::vec Eta = Para.Eta;
  arma::mat Sigma = Para.Sigma;
  arma::vec Delta = Para.Delta;
  arma::mat QInv = Para.QInv;
  arma::vec Sigma2 = Para.Sigma2;
  arma::mat MuMatrix = arma::trans(arma::reshape(Para.Mu, M, Nu));
  arma::vec Lambda = Para.Lambda;
  arma::mat XTheta = Para.XTheta;

  //Get conditional moments
  arma::uvec j(1); j(0) = 3;
  arma::uvec k(4); k(0) = 0, k(1) = 1, k(2) = 2, k(3) = 4;
  arma::rowvec Sigmajk = Sigma(j, k);
  arma::rowvec SigmaPlus = Sigmajk * CholInv(Sigma(k, k));
  double SigmaStar = arma::as_scalar(Sigma(j, j) - SigmaPlus * arma::trans(Sigmajk));
  arma::vec Phi_k = Phi(arma::vectorise(PhiIndeces.rows(k)));
  arma::vec CondMean = arma::kron(OneM, Delta(j)) + arma::kron(EyeM, SigmaPlus) * (Phi_k - arma::kron(OneM, Delta(k)));
  arma::mat CondRooti = GetRooti(SigmaStar * QInv, EyeM);

  //Initialize objects
  arma::uvec LocInd(1), ColIndeces(2), RowIndeces(Nu);
  arma::colvec YStarWideLoc(Nu), MeanLikelihoodLoc(Nu);
  arma::colvec LambdaLoc(2), LambdaLocProposal(2), Lambda1(M), Lambda1Proposal(M), LambdaProposal(2 * M);
  arma::mat tLambda1(M, 2), tLambda1Proposal(M, 2), OmegaLocChol(Nu, Nu), OmegaLocProposalChol(Nu, Nu);
  arma::vec Sigma2LocProposal(Nu), RandU(1);
  arma::mat RootiLikelihood(Nu, Nu), RootiLikelihoodProposal(Nu, Nu);
  double Component1A, Component1B, Component1;
  double Component2A, Component2B, Component2;
  double TuningSD, LogR, Lambda1LocProposal;

  //Indeces
  arma::vec SeqNu(Nu), Seq2(2);
  for (int i = 0; i < Nu; i++) SeqNu(i) = i;
  for (int i = 0; i < 2; i++) Seq2(i) = i;

  //Loop over visits
  for (int Loc = 0; Loc < M; Loc++) {

    //Indeces
    RowIndeces = arma::conv_to<arma::uvec>::from(SeqNu * M + Loc);
    ColIndeces = arma::conv_to<arma::uvec>::from(Seq2 + Loc * 2);

    //Visit specific objects
    LocInd(0) = Loc;
    YStarWideLoc = YStarWide.col(Loc);
    LambdaLoc = Lambda(ColIndeces);
    tLambda1 = arma::trans(arma::reshape(Lambda, 2, M));
    Lambda1 = tLambda1.col(1);
    MeanLikelihoodLoc = MuMatrix.col(Loc);
    OmegaLocChol = arma::diagmat(sqrt(Sigma2(RowIndeces)));
    TuningSD = sqrt(MetropLambda1Vec(Loc));

    //Numerical fix for when a propopsal is rounded to -Inf or Inf
    Lambda1LocProposal = arma::datum::inf;
    LambdaProposal = Lambda;
    LambdaLocProposal = LambdaLoc;
    while ((!arma::is_finite(Lambda1LocProposal))) {

      //Sample proposal
      Lambda1LocProposal = arma::as_scalar(rnormRcpp(1, LambdaLoc(1), TuningSD));
      LambdaLocProposal(1) = Lambda1LocProposal;
      LambdaProposal.rows(ColIndeces) = LambdaLocProposal;
      tLambda1Proposal = arma::trans(arma::reshape(LambdaProposal, 2, M));
      Lambda1Proposal = tLambda1Proposal.col(1);
      Sigma2LocProposal = arma::exp(2 * (XTheta.rows(RowIndeces) * LambdaProposal));
      OmegaLocProposalChol = arma::diagmat(sqrt(Sigma2LocProposal));

    }

    //Likelihood Component
    RootiLikelihoodProposal = arma::solve(arma::trimatu(OmegaLocProposalChol), EyeNu);
    RootiLikelihood = arma::solve(arma::trimatu(OmegaLocChol), EyeNu);
    Component1A = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihoodProposal);
    Component1B = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihood);
    Component1 = Component1A - Component1B;

    //Prior components
    Component2A = lndMvn(Lambda1Proposal, CondMean, CondRooti);
    Component2B = lndMvn(Lambda1, CondMean, CondRooti);
    Component2 = Component2A - Component2B;

    //Log acceptance ratio
    LogR = Component1 + Component2;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceLambda1Vec(Loc)++;

      //Update parameters output
      Lambda = LambdaProposal;
      Sigma2.rows(RowIndeces) = Sigma2LocProposal;

    }

    //End loop over locations
  }

  //Update Metropolis object
  MetrObj.AcceptanceLambda1Vec = AcceptanceLambda1Vec;

  //Update Para objects
  Para.Lambda = Lambda;
  Para.Sigma2 = Sigma2;
  Para.Phi = CreatePhi(Beta, Lambda, Eta, M);
  Para.Omega = arma::diagmat(Sigma2);
  Para.OmegaInv = arma::diagmat(1 / Sigma2);

  //Return final object
  return std::pair<para, metrobj>(Para, MetrObj);

}



//Function to sample new value of eta using a Metropolis sampler step------------------------------
std::pair<para, metrobj> SampleEta(datobj DatObj, para Para, metrobj MetrObj) {

  //Set data objects
  arma::umat PhiIndeces = DatObj.PhiIndeces;
  arma::mat EyeM = DatObj.EyeM;
  arma::colvec OneM = DatObj.OneM;
  arma::colvec OneNu = DatObj.OneNu;
  arma::mat YStarWide = arma::trans(DatObj.YStarWide);
  arma::mat EyeNu = DatObj.EyeNu;
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  int N = DatObj.N;
  arma::colvec Time = DatObj.Time;
  double tNu = DatObj.tNu;
  double t1 = DatObj.t1;

  //Set parameter objects
  arma::colvec Phi = Para.Phi;
  arma::colvec Beta = Para.Beta;
  arma::colvec Eta = Para.Eta;
  arma::mat Sigma = Para.Sigma;
  arma::colvec Delta = Para.Delta;
  arma::mat QInv = Para.QInv;
  arma::colvec Sigma2 = Para.Sigma2;
  arma::mat MuMatrix = arma::trans(reshape(Para.Mu, M, Nu));
  arma::colvec Lambda = Para.Lambda;
  arma::mat XTheta = Para.XTheta;
  arma::colvec Theta = Para.Theta;

  //Set Metropolis Tuning Objects
  arma::vec MetropEtaVec = MetrObj.MetropEtaVec;
  arma::vec AcceptanceEtaVec = MetrObj.AcceptanceEtaVec;

  //Get conditional moments
  arma::uvec j(1); j(0) = 4;
  arma::uvec k(4); k(0) = 0, k(1) = 1, k(2) = 2, k(3) = 3;
  arma::rowvec Sigmajk = Sigma(j, k);
  arma::rowvec SigmaPlus = Sigmajk * CholInv(Sigma(k, k));
  double SigmaStar = arma::as_scalar(Sigma(j, j) - SigmaPlus * arma::trans(Sigmajk));
  arma::vec Phi_k = Phi(arma::vectorise(PhiIndeces.rows(k)));
  arma::vec CondMean = arma::kron(OneM, Delta(j)) + arma::kron(EyeM, SigmaPlus) * (Phi_k - arma::kron(OneM, Delta(k)));
  arma::mat CondRooti = GetRooti(SigmaStar * QInv, EyeM);

  //Initialize objects
  arma::uvec LocInd(1), ColIndeces(2), RowIndeces(Nu);
  arma::colvec YStarWideLoc(Nu), MeanLikelihoodLoc(Nu), MeanLikelihoodLocProposal(Nu);
  arma::mat OmegaLocChol(Nu, Nu), OmegaLocProposalChol(Nu, Nu);
  arma::vec Sigma2LocProposal(Nu), RandU(1), EtaProposal(M);
  arma::mat RootiLikelihood(Nu, Nu), RootiLikelihoodProposal(Nu, Nu);
  arma::mat XThetaProposal(N, 2 * M), XThetaLoc(Nu, 2), XThetaLocProposal(Nu, 2);
  double Component1A, Component1B, Component1;
  double Component2A, Component2B, Component2;
  double TuningSD, LogR, EtaLocProposal, ThetaLocProposal;

  //Indeces
  arma::vec SeqNu(Nu), Seq2(2);
  for (int i = 0; i < Nu; i++) SeqNu(i) = i;
  for (int i = 0; i < 2; i++) Seq2(i) = i;

  //Loop over visits
  for (int Loc = 0; Loc < M; Loc++) {

    //Indeces
    RowIndeces = arma::conv_to<arma::uvec>::from(SeqNu * M + Loc);
    ColIndeces = arma::conv_to<arma::uvec>::from(Seq2 + Loc * 2);

    //Visit specific objects
    LocInd(0) = Loc;
    YStarWideLoc = YStarWide.col(Loc);
    MeanLikelihoodLoc = MuMatrix.col(Loc);
    OmegaLocChol = arma::diagmat(sqrt(Sigma2(RowIndeces)));
    TuningSD = sqrt(MetropEtaVec(Loc));

    //Numerical fix for when a propopsal is rounded to -Inf or Inf
    EtaLocProposal = arma::datum::inf;
    EtaProposal = Eta;
    while ((!arma::is_finite(EtaLocProposal))) {

      //Sample proposal
      EtaLocProposal = arma::as_scalar(rnormRcpp(1, Eta(Loc), TuningSD));
      EtaProposal(Loc) = EtaLocProposal;
      ThetaLocProposal = std::max(std::min(tNu, EtaLocProposal), t1);
      XThetaLoc = GetXThetaLoc(ThetaLocProposal, Time, OneNu, Nu);
      XThetaProposal = XTheta;
      XThetaLocProposal = XThetaProposal.rows(RowIndeces);
      XThetaLocProposal.cols(ColIndeces) = XThetaLoc;
      MeanLikelihoodLocProposal = XThetaLocProposal * Beta;
      Sigma2LocProposal = arma::exp(2 * (XThetaLocProposal * Lambda));
      OmegaLocProposalChol = arma::diagmat(sqrt(Sigma2LocProposal));

    }

    //Likelihood Component
    RootiLikelihoodProposal = arma::solve(arma::trimatu(OmegaLocProposalChol), EyeNu);
    RootiLikelihood = arma::solve(arma::trimatu(OmegaLocChol), EyeNu);
    Component1A = lndMvn(YStarWideLoc, MeanLikelihoodLocProposal, RootiLikelihoodProposal);
    Component1B = lndMvn(YStarWideLoc, MeanLikelihoodLoc, RootiLikelihood);
    Component1 = Component1A - Component1B;

    //Prior components
    Component2A = lndMvn(EtaProposal, CondMean, CondRooti);
    Component2B = lndMvn(Eta, CondMean, CondRooti);
    Component2 = Component2A - Component2B;

    //Log acceptance ratio
    LogR = Component1 + Component2;

    //Metropolis update
    RandU = randuRcpp();
    if (log(RandU(0)) < LogR) {

      //Keep count of acceptances
      AcceptanceEtaVec(Loc)++;

      //Update parameters output
      Eta = EtaProposal;
      Theta(Loc) = ThetaLocProposal;
      XTheta.rows(RowIndeces) = XThetaLocProposal;
      Sigma2.rows(RowIndeces) = Sigma2LocProposal;

    }

    //End loop over locations
  }

  //Update Metropolis object
  MetrObj.AcceptanceEtaVec = AcceptanceEtaVec;

  //Update Para objects
  Para.Eta = Eta;
  Para.Theta = Theta;
  Para.XTheta = XTheta;
  Para.Sigma2 = Sigma2;
  Para.Phi = CreatePhi(Beta, Lambda, Eta, M);
  Para.Omega = arma::diagmat(Sigma2);
  Para.OmegaInv = arma::diagmat(1 / Sigma2);
  Para.Mu = XTheta * Beta;

  //Return final object
  return std::pair<para, metrobj>(Para, MetrObj);

}



//Function to sample latent probit process using Gibbs sampling step------------------------------------------------
arma::vec SampleProbit(datobj DatObj, para Para, dataug DatAug) {

  //Set data objects
  arma::vec YStar = DatObj.YStar;

  //Set parameters
  arma::vec Mu = Para.Mu;
  arma::vec Sigma2 = Para.Sigma2;

  //Set data augmentation objects
  arma::uvec WhichBelow = DatAug.WhichBelow;
  arma::uvec WhichAbove = DatAug.WhichAbove;
  int NBelow = DatAug.NBelow;
  int NAbove = DatAug.NAbove;

  //Moments
  arma::vec Mean = Mu(WhichBelow);
  arma::vec SD = arma::sqrt(Sigma2(WhichBelow));

  //Sample latent Variable from full conditional
  for (int i = 0; i < NBelow; i++) {
    double Temp = rtnormRcpp(Mean(i), SD(i), true);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(Mean(i), SD(i), -100000, 0);
    if (!arma::is_finite(Temp)) Rcpp::stop("infinte value sampled in Tobit sampling step. Most likey cause for this error is that the data being used is inappropriate (i.e., to far from zero) for a Tobit model. Consider scaling towards zero and re-running.");
    YStar(WhichBelow(i)) = Temp;
  }
  for (int i = 0; i < NAbove; i++) {
    double Temp = rtnormRcpp(Mean(i), SD(i), false);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(Mean(i), SD(i), 0, 100000);
    if (!arma::is_finite(Temp)) Rcpp::stop("infinte value sampled in Tobit sampling step. Most likey cause for this error is that the data being used is inappropriate (i.e., to far from zero) for a Tobit model. Consider scaling towards zero and re-running.");
    YStar(WhichAbove(i)) = Temp;
  }
  return YStar;

}




//Function to sample latent tobit process using Gibbs sampling step------------------------------------------------
arma::vec SampleTobit(datobj DatObj, para Para, dataug DatAug) {

  //Set data objects
  arma::vec YStar = DatObj.YStar;

  //Set parameters
  arma::vec Mu = Para.Mu;
  arma::vec Sigma2 = Para.Sigma2;

  //Set data augmentation objects
  int NBelow = DatAug.NBelow;
  arma::uvec WhichBelow = DatAug.WhichBelow;

  //Moments
  arma::vec Mean = Mu(WhichBelow);
  arma::vec SD = arma::sqrt(Sigma2(WhichBelow));

  //Sample latent Variable from full conditional
  for (int i = 0; i < NBelow; i++) {
    double Temp = rtnormRcpp(Mean(i), SD(i), true);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(Mean(i), SD(i), -arma::datum::inf, 0);
    if (!arma::is_finite(Temp)) Rcpp::stop("infinte value sampled in Tobit sampling step. Most likey cause for this error is that the data being used is inappropriate (i.e., to far from zero) for a Tobit model. Consider scaling towards zero and re-running.");
    YStar(WhichBelow(i)) = Temp;
  }
  return YStar;

}




//Function to sample latent process from its full conditional------------------------------------------------------
datobj SampleY(datobj DatObj, para Para, dataug DatAug) {

  //Set data objects
  int FamilyInd = DatObj.FamilyInd;
  int N = DatObj.N;
  int M = DatObj.M;
  int Nu = DatObj.Nu;

  //Sample latent process
  arma::vec YStar(N);
  if (FamilyInd == 1) YStar = SampleProbit(DatObj, Para, DatAug);
  if (FamilyInd == 2) YStar = SampleTobit(DatObj, Para, DatAug);

  //Save output
  DatObj.YStar = YStar;
  DatObj.YStarWide = arma::reshape(YStar, M, Nu);
  return DatObj;

}
