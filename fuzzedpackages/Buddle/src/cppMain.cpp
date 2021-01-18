#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


#include "ActName.h"
#include "Common.h"
#include "Identity.h"
#include "Relu.h"
#include "Dropout.h"
#include "Batchnorm.h"
#include "LeakyRelu.h"
#include "TanH.h"
#include "ArcTan.h"
#include "ArcSinH.h"
#include "ElliotSig.h"
#include "SoftPlus.h"
#include "BentIdentity.h"
#include "Sinusoid.h"
#include "Sigmoid.h"
#include "Gaussian.h"
#include "Sinc.h"
#include "FInv.h"
#include "Link.h"
#include "Affine.h"
#include "SoftmaxLoss.h"
#include "L2loss.h"
#include "Optimization.h"
#include "Layer.h"
#include "Buddle.h"



//'@keywords internal
//[[Rcpp::export]]
List Buddle_Main(arma::mat X_train, arma::mat T_train, arma::mat X_test, arma::mat T_test, int nBatch_Size,
                 int nTotal_Iterations, arma::vec HiddenLayer, int bBatch, int bDrop, double drop_ratio, 
                 double d_learning_rate, double d_init_weight, arma::vec nstrVec, String strOpt, 
                 String strType, int bRand, String strDist, int bDisp){
  
  
  int nTrain_Size = X_train.n_cols;
  int nTest_Size = X_test.n_cols;
  
  int p = X_train.n_rows;
  int r= T_train.n_rows;
  
  arma::mat X(p, nBatch_Size);
  arma::mat T_Mat(r, nBatch_Size);
  
  /////////////////////////////////////// Iteration parameters
  
  int nIter_per_Epoch = nTrain_Size / nBatch_Size;
  
  int nEpoch = nTotal_Iterations/nIter_per_Epoch;
  
  arma::vec train_loss_vec(nEpoch);
  arma::vec train_accuracy_vec(nEpoch);
  arma::vec test_loss_vec(nEpoch);
  arma::vec test_accuracy_vec(nEpoch);
  double dLoss = 0;
  double dAccuracy=0;
  int nRemain=0;
  
  
  ///////////////////////////////////////MLN parameters
  
  
  Layer* Arr_Layer;
  
  int bTest=0;                         /// In test mode for dropout ratio, it should be 1
  int nHiddenLayer = HiddenLayer.n_elem;
  String* strVec = new String[nHiddenLayer];       ///// Activation function name vector
  
  MakeStrVec(nstrVec, strVec);                     ////// Convert to String Vector   
  
  Buddle BClass(p,r,nBatch_Size, nHiddenLayer, HiddenLayer, strType, strVec, strOpt,
          d_learning_rate, d_init_weight, bBatch, bDrop, bTest, drop_ratio, bRand, strDist);
  
  int bTest_predict = 1;
  Buddle BClass_predict(p,r,nTest_Size, nHiddenLayer, HiddenLayer, strType, strVec, strOpt,
           d_learning_rate, d_init_weight, bBatch, bDrop, bTest_predict, drop_ratio, bRand, strDist);
  

  int nInc=0;
  int nInd=0;
  IntegerVector SelVec ; //= RandInts(nBatch_Size, nTrain_Size);
    
  for(int i=1; i<=nTotal_Iterations;i++){
    
    SelVec = RandInts(nBatch_Size, nTrain_Size);
    
    for(int j=1;j<=nBatch_Size;j++){
      nInd = SelVec[j-1];
      X.col(j-1) = X_train.col(nInd-1);
      T_Mat.col(j-1) = T_train.col(nInd-1);
    }

    BClass.forward(X, T_Mat);
    BClass.backward(X, T_Mat);

    BClass.Set_loss();
    dLoss = BClass.Get_loss();
    
    BClass.Set_Accuracy(T_Mat);
    dAccuracy = BClass.Get_Accuracy();
    
    nRemain = GetRemains(i, nIter_per_Epoch);
    
    double d_test_accuracy=0;
    double d_test_loss;
    if(nRemain==0){
      
      Arr_Layer = BClass.Get_Arr_Layer();
      BClass_predict.Set_Arr_Layer(Arr_Layer);
      
      BClass_predict.forward(X_test, T_test);     ///// Use forward not forward_predict
                                                  ///// to get loss.
                                                  ///// Use forward_predict only to predict when no answerkey exists 
      BClass_predict.Set_loss();
      d_test_loss = BClass_predict.Get_loss();
      
      BClass_predict.Set_Accuracy(T_test);
      d_test_accuracy = BClass_predict.Get_Accuracy();
      
      if(bDisp==1){
        Rprintf("Iteration, loss, and accuracy of training are %i, %.3f, %.3f \n", i, dLoss, dAccuracy);
        Rprintf("Epoch, loss and accuracy of test are %i, %.3f, and %.3f \n", nInc, d_test_loss, d_test_accuracy);
        
      }
      
      
      train_loss_vec[nInc] = dLoss;
      train_accuracy_vec[nInc] = dAccuracy;
      
      test_loss_vec[nInc] = d_test_loss;
      test_accuracy_vec[nInc] = d_test_accuracy;
      
      nInc++;

    }
    
  }
  
  
  int nHid = nHiddenLayer+1;
  
  List lst_W(nHid);
  List lst_b(nHid);
  List lst_Result(4);
  
  Arr_Layer = BClass.Get_Arr_Layer();
  
  for(int i=1;i<=nHid;i++){
    lst_W[i-1] = Arr_Layer[i-1].Get_W();
    lst_b[i-1] = Arr_Layer[i-1].Get_b();
  }
  
  lst_Result[0] = train_loss_vec;
  lst_Result[1] = train_accuracy_vec;
  lst_Result[2] = test_loss_vec;
  lst_Result[3] = test_accuracy_vec;
  
  
  List Out(3);
  Out[0] = lst_W;
  Out[1] = lst_b;
  Out[2] = lst_Result;
  
  
  ////////////////////////
  
  /////////////////////   Release dynamic memory
  delete []strVec;
  
  return Out;

}



//'@keywords internal
//[[Rcpp::export]]
List Buddle_Predict(arma::mat X, List lW, List lb, List lParam){

  int r = lParam[0];
  arma::vec HiddenLayer = lParam[1]; 
  int bBatch = lParam[2];  
  int bDrop = lParam[3];  
  double drop_ratio = lParam[4];
  double d_learning_rate = lParam[5]; 
  double d_init_weight = lParam[6]; 
  arma::vec nstrVec = lParam[7];
  String strOpt = lParam[8];
  String strType = lParam[9];
  int bRand =   lParam[10];
  String strDist = lParam[11];
  
  String strClass("Classification");
  String strReg("Regression");
  
  int p = X.n_rows;
  int n= X.n_cols;

  Layer* Arr_Layer;
  
  int nHiddenLayer = HiddenLayer.n_elem;
  String* strVec = new String[nHiddenLayer];       ///// Activation function name vector
  
  MakeStrVec(nstrVec, strVec);                     ////// Convert to String Vector   
  
  int bTest_predict = 1;
  
  Buddle BClass_predict(p,r,n, nHiddenLayer, HiddenLayer, strType, strVec, strOpt,
                  d_learning_rate, d_init_weight, bBatch, bDrop, 
                  bTest_predict, drop_ratio, bRand, strDist);
  
  Arr_Layer = BClass_predict.Get_Arr_Layer();
  
  for(int i=1;i<=(nHiddenLayer+1);i++){
    Arr_Layer[i-1].Set_W( lW[i-1] );
    Arr_Layer[i-1].Set_b( lb[i-1] );
  }
  
  BClass_predict.forward_predict(X);
  arma::mat mOut_predict(r, n);
  arma::mat mOneHot(r, n);
  
  mOut_predict = BClass_predict.Get_mOut();
  
  if(strType == strClass){
    mOneHot = Con2OneHotEncoding(mOut_predict);
  }else{
    mOneHot.zeros();
  }
  
  
  List Out(2);
  
  
  Out[0] = mOut_predict;
  Out[1] = mOneHot;
  
  ////////////////////////
  
  /////////////////////   Release dynamic memory
  delete []strVec;
  
  return Out;
  
  
}






