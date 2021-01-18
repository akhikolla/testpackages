#include "mdr.h"

using namespace Rcpp;
using namespace arma;

RcppExport SEXP mdr(SEXP X, SEXP fold, SEXP status, SEXP t, SEXP cv, SEXP cvp, SEXP top, SEXP na, SEXP fix){

    Rcpp::NumericMatrix XX(X);                 // contains the genotype data
    Rcpp::NumericVector Status(status);        // contains the status
    
    int FOLD = Rcpp::as<int>(fold);            // until what level of interaction should the data be analyzed
    double T = Rcpp::as<int>(t);               // Threshold for low/high risk classification, default: ratio of group sizes
    int CV = Rcpp::as<int>(cv);                // number of cross-validation runs
    double CVP = Rcpp::as<double>(cvp);        // Ratio of training and testing data
    int TOP = Rcpp::as<int>(top);              // Length of output list
    TOP = TOP + 1;
    // int NA = Rcpp::as<int>(na);
    int FIX = Rcpp::as<int>(fix);              // Should one genotyp be forced to be part of the analysis
    
    int k = XX.ncol();
    int p = XX.nrow();
    
    vec index(p);
    for(int i=0;i<p;i++) index(i)=i;
    
    //Create the genotype matrix
    mat XXX(XX.begin(), p, k, false);
    
    mat classifierOne(k,3);
    vec tempCOne(3);
    mat evalClassifierOne(k,7);
    vec tempEvalOne(7);
    vec topOnePrec(TOP);
    vec topOneIndex(TOP);
    topOnePrec.zeros();
    topOneIndex.zeros();
    
    cube evalOneTrain(k,7,CV);
    cube evalOneTest(k,7,CV);
    
    if(FIX==-1){
      if(CV==0){
	for(int i=0;i<k;i++){
	    tempCOne = classOne(XXX.col(i),Status, T);
	    for(int j=0;j<3;j++){
	      classifierOne(i,j) = tempCOne(j);
	    }
	  }
      
	  int posOne=0;

	  for(int i=0;i<k;i++){
	    tempEvalOne = evalClassOne(XXX.col(i),classifierOne.row(i).t(), Status);
	    for(int j=0;j<7;j++){
	      evalClassifierOne(i,j) = tempEvalOne(j);
	    }
	    if(topOnePrec(0)<tempEvalOne(6)){
	      posOne=0;
	      while((topOnePrec(posOne)<tempEvalOne(6))&&posOne<(TOP-1)){
		topOnePrec(posOne) = topOnePrec(posOne+1);
		topOneIndex(posOne) = topOneIndex(posOne+1);
		posOne +=1;
	      }
	      topOnePrec(posOne-1) = tempEvalOne(6);
	      topOneIndex(posOne-1) = i+1;
	    }
	  }    
      } else {
	// This is now the cross-validation for FOLD=1 and FIX=-1
      }
    } else {
        if(CV==0){
	  tempCOne = classOne(XXX.col(FIX),Status, T);
	  for(int j=0;j<3;j++){
	    classifierOne(FIX,j) = tempCOne(j);
	  }
    
	  int posOne=0;

	  tempEvalOne = evalClassOne(XXX.col(FIX),classifierOne.row(FIX).t(), Status);
	  for(int j=0;j<7;j++){
	    evalClassifierOne(FIX,j) = tempEvalOne(j);
	  }
	  if(topOnePrec(0)<tempEvalOne(6)){
	    posOne=0;
	    while((topOnePrec(posOne)<tempEvalOne(6))&&posOne<(TOP-1)){
	      topOnePrec(posOne) = topOnePrec(posOne+1);
	      topOneIndex(posOne) = topOneIndex(posOne+1);
	      posOne +=1;
	    }
	    topOnePrec(posOne-1) = tempEvalOne(6);
	    topOneIndex(posOne-1) = FIX+1;
	  }
	} else {
	 // This is now the cross-validation for FOLD=1 and FIX!=0
	}
    }
      
    cube classifierTwo(3,3,k*(k-1)/2+1);
    mat evalClassifierTwo(k*(k-1)/2+1,7);
    // Toplist requirements:
    vec topTwoPrec(TOP);
    vec topTwoIndex1(TOP);
    vec topTwoIndex2(TOP);
    topTwoPrec.zeros();
    topTwoIndex1.zeros();
    topTwoIndex2.zeros();
    int posTwo=0;
    
    if(FOLD > 1){
      // Calculate the Classifier Two:
      mat tempMatTwo(p,2);
      int runTwo=0;
      if(FIX==-1){
 	if(CV==0){
	  for(int i=0;i<(k-1);i++){
	    for(int j=i+1;j<k;j++){
	      tempMatTwo.col(0)=XXX.col(i);
	      tempMatTwo.col(1)=XXX.col(j);
	      //classifierTwo.slice(runTwo) = mat(3,3).ones();
	      classifierTwo.slice(runTwo) = classTwo(tempMatTwo, Status, T);
	      runTwo +=1;
	    }
	  }
	
	//Evaluate the Classifier Two
	  vec tempEvalTwo(7);
	  runTwo = 0;
	  for(int i=0;i<(k-1);i++){
	    for(int j=i+1;j<k;j++){
	      tempMatTwo.col(0)=XXX.col(i);
	      tempMatTwo.col(1)=XXX.col(j);
	      
	      tempEvalTwo = evalClassTwo(tempMatTwo,classifierTwo.slice(runTwo), Status);
	    
	      for(int a=0;a<7;a++){
		evalClassifierTwo(runTwo,a) = tempEvalTwo(a);
	      }
	      runTwo +=1;
	      if(topTwoPrec(0)<tempEvalTwo(6)){
		posTwo=0;
		while((topTwoPrec(posTwo)<tempEvalTwo(6))&&posTwo<(TOP-1)){
		  topTwoPrec(posTwo) = topTwoPrec(posTwo+1);
		  topTwoIndex1(posTwo) = topTwoIndex1(posTwo+1);
		  topTwoIndex2(posTwo) = topTwoIndex2(posTwo+1);
		  posTwo +=1;
		}
		topTwoPrec(posTwo-1) = tempEvalTwo(6);
		topTwoIndex1(posTwo-1) = i+1;
		topTwoIndex2(posTwo-1) = j+1;
	      }
	    }
	  }
	} else {
	 // This is now the cross-validation part for non-fixed in case of FOLD=2 
	}
      } else {
	// Here is now the fixed part:
	int i=FIX;
	  for(int j=0;j<k;j++){
	    tempMatTwo.col(0)=XXX.col(i);
	    tempMatTwo.col(1)=XXX.col(j);
	    //classifierTwo.slice(runTwo) = mat(3,3).ones();
	    classifierTwo.slice(runTwo) = classTwo(tempMatTwo, Status, T);
	    runTwo +=1;
	  }
      
      //Evaluate the Classifier Two
	vec tempEvalTwo(7);
	runTwo = 0;
	i=FIX;
	  for(int j=0;j<k;j++){
	    tempMatTwo.col(0)=XXX.col(i);
	    tempMatTwo.col(1)=XXX.col(j);
	    
	    tempEvalTwo = evalClassTwo(tempMatTwo,classifierTwo.slice(runTwo), Status);
	  
	    for(int a=0;a<7;a++){
	      evalClassifierTwo(runTwo,a) = tempEvalTwo(a);
	    }
	    runTwo +=1;
	    if(topTwoPrec(0)<tempEvalTwo(6)){
	      posTwo=0;
	      while((topTwoPrec(posTwo)<tempEvalTwo(6))&&posTwo<(TOP-1)){
		topTwoPrec(posTwo) = topTwoPrec(posTwo+1);
		topTwoIndex1(posTwo) = topTwoIndex1(posTwo+1);
		topTwoIndex2(posTwo) = topTwoIndex2(posTwo+1);
		posTwo +=1;
	      }
	      topTwoPrec(posTwo-1) = tempEvalTwo(6);
	      topTwoIndex1(posTwo-1) = i+1;
	      topTwoIndex2(posTwo-1) = j+1;
	    }
	  }
      }
      
      
   } //END OF PART 2

//   // Calculate the Classifier Three
    cube classifierThree(3,3,3);
    mat evalClassifierThree(k*(k-1)*(k-2)/6,7);
    // Toplist requirements:
    vec topThreePrec(TOP);
    vec topThreeIndex1(TOP);
    vec topThreeIndex2(TOP);
    vec topThreeIndex3(TOP);
    topThreePrec.zeros();
    topThreeIndex1.zeros();
    topThreeIndex2.zeros();
    topThreeIndex3.zeros();
    int posThree=0;
 
    
    if(FOLD > 2){
      mat tempMatThree(p,3);
      vec tempEvalThree(7);
      int runThree=0;
      if(FIX==-1){
	for(int i=0;i<(k-2);i++){
	  for(int j=i+1;j<(k-1);j++){
	    for(int m=j+1;m<k;m++){
	      tempMatThree.col(0)=XXX.col(i);
	      tempMatThree.col(1)=XXX.col(j);
	      tempMatThree.col(2)=XXX.col(m);
	      // Just keep the latest classifier;
	      classifierThree = classThree(tempMatThree, Status, T);
	      tempEvalThree = evalClassThree(tempMatThree,classifierThree, Status);
	      for(int a=0;a<7;a++){
		evalClassifierThree(runThree,a) = tempEvalThree(a);
	      }
	      runThree +=1;
	      if(topThreePrec(0)<tempEvalThree(6)){
		posThree=0;
		while((topThreePrec(posThree)<tempEvalThree(6))&&posThree<(TOP-1)){
		  topThreePrec(posThree) = topThreePrec(posThree+1);
		  topThreeIndex1(posThree) = topThreeIndex1(posThree+1);
		  topThreeIndex2(posThree) = topThreeIndex2(posThree+1);
		  topThreeIndex3(posThree) = topThreeIndex3(posThree+1);
		  posThree +=1;
	      }
	      topThreePrec(posThree-1) = tempEvalThree(6);
	      topThreeIndex1(posThree-1) = i+1;
	      topThreeIndex2(posThree-1) = j+1;
	      topThreeIndex3(posThree-1) = m+1;
	      }
	    }
	  }
	}
      } else {
	int i=FIX;
	  for(int j=0;j<(k-1);j++){
	    for(int m=j+1;m<k;m++){
	      tempMatThree.col(0)=XXX.col(i);
	      tempMatThree.col(1)=XXX.col(j);
	      tempMatThree.col(2)=XXX.col(m);
	      // Just keep the latest classifier;
	      classifierThree = classThree(tempMatThree, Status, T);
	      tempEvalThree = evalClassThree(tempMatThree,classifierThree, Status);
	      for(int a=0;a<7;a++){
		evalClassifierThree(runThree,a) = tempEvalThree(a);
	      }
	      runThree +=1;
	      if(topThreePrec(0)<tempEvalThree(6)){
		posThree=0;
		while((topThreePrec(posThree)<tempEvalThree(6))&&posThree<(TOP-1)){
		  topThreePrec(posThree) = topThreePrec(posThree+1);
		  topThreeIndex1(posThree) = topThreeIndex1(posThree+1);
		  topThreeIndex2(posThree) = topThreeIndex2(posThree+1);
		  topThreeIndex3(posThree) = topThreeIndex3(posThree+1);
		  posThree +=1;
	      }
	      topThreePrec(posThree-1) = tempEvalThree(6);
	      topThreeIndex1(posThree-1) = i+1;
	      topThreeIndex2(posThree-1) = j+1;
	      topThreeIndex3(posThree-1) = m+1;
	      }
	    }
	  }
      }
    } // From case 3

    cube classifierFour(3,3,9);
    mat evalClassifierFour(1,7);
    // Toplist requirements:
    vec topFourPrec(TOP);
    vec topFourIndex1(TOP);
    vec topFourIndex2(TOP);
    vec topFourIndex3(TOP);
    vec topFourIndex4(TOP);
    topFourPrec.zeros();
    topFourIndex1.zeros();
    topFourIndex2.zeros();
    topFourIndex3.zeros();
    topFourIndex4.zeros();
    int posFour=0;
    
    if(FOLD > 3){
      // Calculate the Classifier Four:
      mat tempMatFour(p,4);
      vec tempEvalFour(7);
    //  mat evalClassifierFour(k*(k-1)*(k-2)*(k-3)/24,7);
      int runFour=0;
      if(FIX==-1){
      for(int i=0;i<(k-3);i++){
	for(int j=i+1;j<(k-2);j++){
	  for(int m=j+1;m<(k-1);m++){
	    for(int n=m+1;n<k;n++){
	      tempMatFour.col(0)=XXX.col(i);
	      tempMatFour.col(1)=XXX.col(j);
	      tempMatFour.col(2)=XXX.col(m);
	      tempMatFour.col(3)=XXX.col(n);
	      // Just keep the latest classifier;
	      classifierFour = classFour(tempMatFour, Status, T);
	      tempEvalFour = evalClassFour(tempMatFour,classifierFour, Status);
	      for(int a=0;a<7;a++){
		evalClassifierFour(0,a) = tempEvalFour(a);
	      }
	      //while(topVecPrec(pos)<tempEvalFour)
	      runFour +=1;
	      if(topFourPrec(0)<tempEvalFour(6)){
		posFour=0;
		while((topFourPrec(posFour)<tempEvalFour(6))&&posFour<(TOP-1)){
		  topFourPrec(posFour) = topFourPrec(posFour+1);
		  topFourIndex1(posFour) = topFourIndex1(posFour+1);
		  topFourIndex2(posFour) = topFourIndex2(posFour+1);
		  topFourIndex3(posFour) = topFourIndex3(posFour+1);
		  topFourIndex4(posFour) = topFourIndex4(posFour+1);
		  posFour +=1;
		}
		topFourPrec(posFour-1) = tempEvalFour(6);
		topFourIndex1(posFour-1) = i+1;
		topFourIndex2(posFour-1) = j+1;
		topFourIndex3(posFour-1) = m+1;
		topFourIndex4(posFour-1) = n+1;
	      }
	    }
	  }
	}
      }
    } else {
        int i=FIX;
	for(int j=0;j<(k-2);j++){
	  for(int m=j+1;m<(k-1);m++){
	    for(int n=m+1;n<k;n++){
	      tempMatFour.col(0)=XXX.col(i);
	      tempMatFour.col(1)=XXX.col(j);
	      tempMatFour.col(2)=XXX.col(m);
	      tempMatFour.col(3)=XXX.col(n);
	      // Just keep the latest classifier;
	      classifierFour = classFour(tempMatFour, Status, T);
	      tempEvalFour = evalClassFour(tempMatFour,classifierFour, Status);
	      for(int a=0;a<7;a++){
		evalClassifierFour(0,a) = tempEvalFour(a);
	      }
	      //while(topVecPrec(pos)<tempEvalFour)
	      runFour +=1;
	      if(topFourPrec(0)<tempEvalFour(6)){
		posFour=0;
		while((topFourPrec(posFour)<tempEvalFour(6))&&posFour<(TOP-1)){
		  topFourPrec(posFour) = topFourPrec(posFour+1);
		  topFourIndex1(posFour) = topFourIndex1(posFour+1);
		  topFourIndex2(posFour) = topFourIndex2(posFour+1);
		  topFourIndex3(posFour) = topFourIndex3(posFour+1);
		  topFourIndex4(posFour) = topFourIndex4(posFour+1);
		  posFour +=1;
		}
		topFourPrec(posFour-1) = tempEvalFour(6);
		topFourIndex1(posFour-1) = i+1;
		topFourIndex2(posFour-1) = j+1;
		topFourIndex3(posFour-1) = m+1;
		topFourIndex4(posFour-1) = n+1;
	      }
	    }
	  }
	}
      }
    } // From case 4
      
    List res;
    res["classifierOne"] = T;
    res["evalOne"] = evalClassifierOne;
    res["topOnePrec"] = topOnePrec;
    res["topOneIndex"] = topOneIndex;    
    
    if(FOLD>1) res["classifierTwo"] = T;
    if(FOLD>1) res["evalTwo"] = evalClassifierTwo;
    if(FOLD>1) res["topTwoPrec"] = topTwoPrec;
    if(FOLD>1) res["topTwoIndex1"] = topTwoIndex1;    
    if(FOLD>1) res["topTwoIndex2"] = topTwoIndex2;
        
    if(FOLD>2) res["evalThree"] = evalClassifierThree;
    if(FOLD>2) res["topThreePrec"] = topThreePrec;
    if(FOLD>2) res["topThreeIndex1"] = topThreeIndex1;    
    if(FOLD>2) res["topThreeIndex2"] = topThreeIndex2;
    if(FOLD>2) res["topThreeIndex3"] = topThreeIndex3;
    
    
    if(FOLD>3) res["evalFour"] = evalClassifierFour;
    if(FOLD>3) res["topFourPrec"] = topFourPrec;
    if(FOLD>3) res["topFourIndex1"] = topFourIndex1;    
    if(FOLD>3) res["topFourIndex2"] = topFourIndex2;
    if(FOLD>3) res["topFourIndex3"] = topFourIndex3;
    if(FOLD>3) res["topFourIndex4"] = topFourIndex4;
    
    res["evalOneCVTrain"] = evalOneTrain;
    res["evalOneCVTest"] = evalOneTest;
   
    return res;
}

RcppExport SEXP mdrEnsemble(SEXP X, SEXP fold, SEXP status, SEXP t, SEXP top, SEXP oldstatus, SEXP oldX){
  
     Rcpp::NumericMatrix XX(X);
     int kX = XX.ncol();
     int pX = XX.nrow();
     vec indexX(pX);
     for(int i=0;i<pX;i++) indexX(i)=i;
     mat XXX(XX.begin(), pX, kX, false);

     Rcpp::NumericMatrix Top(top);
     int kTop = Top.ncol();
     int pTop = Top.nrow();
     vec indexTop(pTop);     
     for(int i=0;i<pTop;i++) indexTop(i)=i;
     mat TOP(Top.begin(), pTop, kTop, false);

     Rcpp::NumericMatrix OldX(oldX);
     int kOldX = OldX.ncol();
     int pOldX = OldX.nrow();
     vec indexOldX(pOldX);
     for(int i=0;i<pOldX;i++) indexOldX(i)=i;
     mat OLDX(OldX.begin(), pOldX, kOldX, false);

     Rcpp::NumericVector Status(status);
     Rcpp::NumericVector OldStatus(oldstatus);
     int FOLD = Rcpp::as<int>(fold);
     double T = Rcpp::as<int>(t);
      
     vec tempCOne(3);
     mat classLableOne(pX,pTop);
     vec tempEvalOne(7);
     mat evalClassifierOne(pTop,7);

     mat tempCTwo(3,3);
     mat tempMatTwo(pOldX,2);
     mat tempMatTwoNew(pX,2);
     mat classLableTwo(pX,pTop);
     vec tempEvalTwo(7);
     mat evalClassifierTwo(pTop,7);

     cube tempCThree(3,3,3);
     mat tempMatThree(pOldX,3);
     mat tempMatThreeNew(pX,3);
     mat classLableThree(pX,pTop);
     vec tempEvalThree(7);
     mat evalClassifierThree(pTop,7);


     cube tempCFour(3,3,9);
     mat tempMatFour(pOldX,4);
     mat tempMatFourNew(pX,4);
     mat classLableFour(pX,pTop);
     vec tempEvalFour(7);
     mat evalClassifierFour(pTop,7);


    for(int i=0;i<pTop;i++){

       if(FOLD==1){
	tempCOne = classOne(OLDX.col(TOP(i,0)),OldStatus, T);
	classLableOne.col(i) = classifyOne(XXX.col(TOP(i,0)),tempCOne);
	tempEvalOne = evalClassOne(XXX.col(TOP(i,0)), tempCOne, Status);
	for(int j=0;j<7;j++){
	    evalClassifierOne(i,j) = tempEvalOne(j);
	}

       } else if(FOLD==2){
	 tempMatTwo.col(0)=OLDX.col(TOP(i,0));
	 tempMatTwo.col(1)=OLDX.col(TOP(i,1));
	 tempCTwo = classTwo(tempMatTwo,OldStatus, T);
         tempMatTwoNew.col(0)=XXX.col(TOP(i,0));
	 tempMatTwoNew.col(1)=XXX.col(TOP(i,1));
	 classLableTwo.col(i) = classifyTwo(tempMatTwoNew,tempCTwo);
	 tempEvalTwo = evalClassTwo(tempMatTwoNew, tempCTwo, Status);
	 for(int j=0;j<7;j++){
	    evalClassifierTwo(i,j) = tempEvalTwo(j);
	 }

       } else if(FOLD==3){
	 tempMatThree.col(0)=OLDX.col(TOP(i,0));
	 tempMatThree.col(1)=OLDX.col(TOP(i,1));
	 tempMatThree.col(2)=OLDX.col(TOP(i,2));
	 tempCThree = classThree(tempMatThree,OldStatus, T);
         tempMatThreeNew.col(0)=XXX.col(TOP(i,0));
	 tempMatThreeNew.col(1)=XXX.col(TOP(i,1));
	 tempMatThreeNew.col(2)=XXX.col(TOP(i,2));
	 classLableThree.col(i) = classifyThree(tempMatThreeNew,tempCThree);
	 tempEvalThree = evalClassThree(tempMatThreeNew, tempCThree, Status);
	 for(int j=0;j<7;j++){
	    evalClassifierThree(i,j) = tempEvalThree(j);
	 }

       } else if(FOLD==4){
	 tempMatFour.col(0)=OLDX.col(TOP(i,0));
	 tempMatFour.col(1)=OLDX.col(TOP(i,1));
	 tempMatFour.col(2)=OLDX.col(TOP(i,2));
	 tempMatFour.col(3)=OLDX.col(TOP(i,3));
	 tempCFour = classFour(tempMatFour,OldStatus, T);
         tempMatFourNew.col(0)=XXX.col(TOP(i,0));
	 tempMatFourNew.col(1)=XXX.col(TOP(i,1));
	 tempMatFourNew.col(2)=XXX.col(TOP(i,2));
	 tempMatFourNew.col(3)=XXX.col(TOP(i,3));
	 classLableFour.col(i) = classifyFour(tempMatFourNew,tempCFour);
	 tempEvalFour = evalClassFour(tempMatFourNew, tempCFour, Status);
	 for(int j=0;j<7;j++){
	    evalClassifierFour(i,j) = tempEvalFour(j);
	 }
       }
    }
     
     List res;
     res["t"] = T;
     res["top"] = TOP;
     res["classLableOne"] = classLableOne;
     res["classLableTwo"] = classLableTwo;
     res["classLableThree"] = classLableThree;
     res["classLableFour"] = classLableFour;
     res["evalOne"] = evalClassifierOne;
     res["evalTwo"] = evalClassifierTwo;
     res["evalThree"] = evalClassifierThree;
     res["evalFour"] = evalClassifierFour;
     return res;
}

vec classOne(const vec a, const vec status, double T){
  vec tempH(3);
  vec tempC(3);
  tempH.zeros();
  tempC.zeros();
  
  int tempValue;
  for(unsigned int i=0;i<status.n_elem;i++){
    tempValue = a(i);
    if(tempValue<3){
      if(status(i)==0){
	tempH(tempValue) += 1.0;
      } else {
	tempC(tempValue) += 1.0;
      }
    }
  }
  
  vec result(3);
  result = tempC/tempH;
  for (int i=0;i<3;i++){
    if(result(i)>= T){
      result(i) = 1;
    } else {
      result(i) = 0;
    }
  }
  return result;
}

mat classTwo(const mat a, const vec status, double T){
    mat result(3,3);
    mat tempH(3,3);
    mat tempC(3,3);
    tempH.zeros();
    tempC.zeros();
    
    for(unsigned int i=0;i<status.n_elem;i++){
       if((a(i,0)<3) && (a(i,1)<3)){
	   if(status(i)==0){
	     tempH(a(i,0),a(i,1)) += 1.0;
	   } else {
	     tempC(a(i,0),a(i,1)) += 1.0;
	   }
       }
    }
    
    result = tempC/tempH;
    for (int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	if(result(i,j)> T){
	  result(i,j) = 1;
	} else {
	  result(i,j) = 0;
	}
      }
    }
  
    return result;
}

cube classThree(const mat a, const vec status, double T){
    cube result(3,3,3);
    cube tempH(3,3,3);
    cube tempC(3,3,3);
    tempH.zeros();
    tempC.zeros();
    
    for(unsigned int i=0;i<status.n_elem;i++){
       if((a(i,0)<3) && (a(i,1)<3) && (a(i,2)<3)){
	  if(status(i)==0){
	    tempH(a(i,0),a(i,1),a(i,2)) += 1.0;
	  } else {
	    tempC(a(i,0),a(i,1),a(i,2)) += 1.0;
	  }
       }  
    }
    
    result = tempC/tempH;
    for (int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	for(int k=0;k<3;k++){
	  if(result(i,j,k)> T){
	    result(i,j,k) = 1;
	  } else {
	    result(i,j,k) = 0;
	  }
	}
      }
    }
  
    return result;
}

cube classFour(const mat a, const vec status, double T){
    cube result(3,3,9);
    cube tempH(3,3,9);
    cube tempC(3,3,9);
    tempH.zeros();
    tempC.zeros();
    
    for(unsigned int i=0;i<status.n_elem;i++){
       if((a(i,0)<3) && (a(i,1)<3) && (a(i,2)<3) && (a(i,3)<3)){
	  if(status(i)==0){
	    tempH(a(i,0),a(i,1),a(i,2)+3*a(i,3)) += 1.0;
	  } else {
	    tempC(a(i,0),a(i,1),a(i,2)+3*a(i,3)) += 1.0;
	  }
      }
    }
    
    result = tempC/tempH;
    for (int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	for(int k=0;k<3;k++){
	  for(int l=0;l<3;l++){
	    if(result(i,j,k+3*l)> T){
	      result(i,j,k+3*l) = 1;
	    } else {
	      result(i,j,k+3*l) = 0;
	    }
	  }
	}
      }
    }
  
    return result;
}

/* Input of this function:
   vec a      : 
   vec one    :
   vec status :
*/
vec evalClassOne(const vec a, const vec one, const vec status){
 vec result(7);
 result.zeros();
 vec classified(a.n_elem);
 classified.zeros();
 for(unsigned int i=0;i<a.n_elem;i++){
   if(a(i)<3) classified(i)=one(a(i));
 }
 double tp = 0.0;
 double fp = 0.0;
 double tn = 0.0;
 double fn = 0.0;
 for(unsigned int i=0;i<a.n_elem;i++){
   if(a(i)<3){
    if((status(i)==0)&(classified(i)==0)){
      tn += 1.0;
    } else if ((status(i)==1)&(classified(i)==1)){
	tp += 1.0;
    } else if ((status(i)==1)&(classified(i)==0)){
	fn += 1.0;
    } else if ((status(i)==0)&(classified(i)==1)){
	fp += 1.0; 
    }
   }
 }
 result(0)=tp;
 result(1)=fp;
 result(2)=tn;
 result(3)=fn;
 result(4)=tp/(tp+fn);
 result(5)=tn/(tn+fp);
 result(6)=(result(4)+result(5))/2;
 
 return result;
}

vec evalClassTwo(const mat a, const mat two, const vec status){
 vec result(7);
 result.zeros();
 vec classified(status.n_elem);
 for(unsigned int i=0;i<a.n_rows;i++){
    if((a(i,0)<3) && (a(i,1)<3)) classified(i)=two(a(i,0),a(i,1));
 }
 double tp = 0.0;
 double fp = 0.0;
 double tn = 0.0;
 double fn = 0.0;
 for(unsigned int i=0;i<a.n_rows;i++){
   if((a(i,0)<3) && (a(i,1)<3)){
    if((status(i)==0)&(classified(i)==0)){
      tn += 1.0;
    } else if ((status(i)==1)&(classified(i)==1)){
	tp += 1.0;
    } else if ((status(i)==1)&(classified(i)==0)){
	fn += 1.0;
    } else if ((status(i)==0)&(classified(i)==1)){
	fp += 1.0; 
    }
 }
 }
 result(0)=tp;
 result(1)=fp;
 result(2)=tn;
 result(3)=fn;
 result(4)=tp/(tp+fn);
 result(5)=tn/(tn+fp);
 result(6)=(result(4)+result(5))/2;
 
 return result;
}


vec evalClassThree(const mat a, const cube three, const vec status){
 vec result(7);
 result.zeros();
 vec classified(status.n_elem);
 for(unsigned int i=0;i<a.n_rows;i++){
   if((a(i,0)<3) && (a(i,1)<3) && (a(i,2)<3))   classified(i)=three(a(i,0),a(i,1),a(i,2));
 }
 double tp = 0.0;
 double fp = 0.0;
 double tn = 0.0;
 double fn = 0.0;
 for(unsigned int i=0;i<a.n_rows;i++){
   if((a(i,0)<3) && (a(i,1)<3) && (a(i,2)<3)){
      if((status(i)==0)&(classified(i)==0)){
	tn += 1.0;
      } else if ((status(i)==1)&(classified(i)==1)){
	  tp += 1.0;
      } else if ((status(i)==1)&(classified(i)==0)){
	  fn += 1.0;
      } else if ((status(i)==0)&(classified(i)==1)){
	  fp += 1.0; 
      }
   }
 }
 result(0)=tp;
 result(1)=fp;
 result(2)=tn;
 result(3)=fn;
 result(4)=tp/(tp+fn);
 result(5)=tn/(tn+fp);
 result(6)=(result(4)+result(5))/2;
 
 return result;
}


vec evalClassFour(const mat a, const cube four, const vec status){
 vec result(7);
 result.zeros();
 vec classified(status.n_elem);
 for(unsigned int i=0;i<a.n_rows;i++){
      if((a(i,0)<3) && (a(i,1)<3) && (a(i,2)<3) && (a(i,3)<3)) classified(i)=four(a(i,0),a(i,1),a(i,2)+3*a(i,3));
 }
 double tp = 0.0;
 double fp = 0.0;
 double tn = 0.0;
 double fn = 0.0;
 for(unsigned int i=0;i<a.n_rows;i++){
   if((a(i,0)<3) && (a(i,1)<3) && (a(i,2)<3) && (a(i,3)<3)){
      if((status(i)==0)&(classified(i)==0)){
	tn += 1.0;
      } else if ((status(i)==1)&(classified(i)==1)){
	  tp += 1.0;
      } else if ((status(i)==1)&(classified(i)==0)){
	  fn += 1.0;
      } else if ((status(i)==0)&(classified(i)==1)){
	  fp += 1.0; 
      }
   }
 }
 result(0)=tp;
 result(1)=fp;
 result(2)=tn;
 result(3)=fn;
 result(4)=tp/(tp+fn);
 result(5)=tn/(tn+fp);
 result(6)=(result(4)+result(5))/2;
 
 return result;
}

vec classifyOne(const vec a, const vec one){
 vec result(a.n_elem);
 result.zeros();
 result += 2;
 for(unsigned int i=0;i<a.n_elem;i++){
    if(a(i)<3) result(i)=one(a(i));
 }
 return result;
}

vec classifyTwo(const mat a, const mat one){
 vec result(a.n_rows);
 result.zeros();
 result += 2;
 for(unsigned int i=0;i<a.n_rows;i++){
    if((a(i,0)<3) && (a(i,1)<3)) result(i)=one(a(i,0),a(i,1));
 }
 return result;
}


vec classifyThree(const mat a, const cube one){
 vec result(a.n_rows);
 result.zeros();
 result += 2;
 for(unsigned int i=0;i<a.n_rows;i++){
    if((a(i,0)<3) && (a(i,1)<3) && (a(i,2)<3)) result(i)=one(a(i,0),a(i,1),a(i,2));
 }
 return result;
}

vec classifyFour(const mat a, const cube one){
 vec result(a.n_rows);
 result.zeros();
 result += 2;
 for(unsigned int i=0;i<a.n_rows;i++){
    if((a(i,0)<3) && (a(i,1)<3) && (a(i,2)<3) && (a(i,3)<3)) result(i)=one(a(i,0),a(i,1),a(i,2)+3*a(i,3));
 }
 return result;
}
