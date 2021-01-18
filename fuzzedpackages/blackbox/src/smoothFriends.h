#ifndef H_KRIG
#define H_KRIG

#include <iostream>
#include "Matern.h" // Matern
#include "intern_newCSmooth.h"
#include <Rcpp.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

//==============================================================================
// FRIENDs de la class CSmooth
// d['e]clar['e]es friend pour ?tre utilisable comme objective functions
// par bisection search, brent...
// apriori a declareravant les  membres de CSmooth qui y font appel (meme si GCC n'est pas strict la dessus)

template<typename Typeforcov>
Typeforcov Krig_fdf(Typeforcov loglambda) {
        Typeforcov z=0.;
        for (typename std::vector<internalTypeEigen>::iterator ii=test->D_invEigVals.begin();ii!=test->D_invEigVals.end();ii++)
            z+=1./(1.+exp(loglambda)*(*ii));
return (z - test->df);
}



template <typename Typeforcov>
Typeforcov Krig_fgcv(Typeforcov lambda) {
  //Rprintf("debut Krig_fgcv: %f \n",lambda);
    long double MSE=0.,tmp;
    Typeforcov dfOver_n=0.;
    std::vector<internalTypeEigen> lD(0);
    for (typename std::vector<internalTypeEigen>::iterator ii=test->D_invEigVals.begin();ii!=test->D_invEigVals.end();ii++)
       lD.push_back((*ii)*lambda);
//ostream_vector(test->D_invEigVals,std::cout);
    int iii=0;
    for (typename std::vector<internalTypeEigen>::iterator ii=lD.begin();ii!=lD.end();ii++,iii++) {
        tmp=test->u[iii]*(*ii)/(1.+(*ii));
        MSE+=tmp*tmp;
    }
//ostream_vector(test->u,std::cout);
    MSE/=lD.size();
/***************
 Experiments confirm that for FIXED covparam and lambda, the MSE at this point is different whether
 two replicates estimates or only their mean was given as input in pointls.
 Hence the weighted MSE depends on the replicates number, it cannot simply be considered a measure
 of prediction error of the mean values.
***************/
//std::cout<<"MSE: "<<MSE;getchar();
//std::cout<<std::endl<<test->pureSS<<" "<<test->KgPtNbr_with_repl-test->KgPtNbr<<" "<<MSE;getchar();
    if (test->pureSS>0) MSE=MSE+test->pureSS/(test->KgPtNbr_with_repl-test->KgPtNbr); /*fields's method for the pure error*/
//std::cout<<std::endl<<MSE;getchar();
    for (typename std::vector<internalTypeEigen>::iterator ii=lD.begin();ii!=lD.end();ii++)
        dfOver_n+=(1./(1.+(*ii)));  //trA... sum estimates df du predictor
//    Typeforcov denexp=2.*lD.size()/(lD.size()-denexp-1.); // Hurvich-Tsai like correction
    Typeforcov denexp=2.;
    dfOver_n=1.-dfOver_n/lD.size(); // now (residual df)/n
    fnevalcounter++;
    if (dfOver_n>0.) return (Typeforcov(MSE/pow(dfOver_n,denexp)));
    else { //should have been prevented at the level of choice of bounds for df in gcv_Krig()
#ifdef NO_R_CONSOLE
        std::cerr<<"(!) From Krig_fgcv() in DLL: !dfOver_n>0."<<std::endl;
        if (false) {
            std::cerr<<"lambda: "<<lambda<<" MSE: "<<MSE<<" dfOver_n: "<<dfOver_n<<std::endl;
            std::cerr<<"D_invEigVals"<<std::endl;
            ostream_vector(test->D_invEigVals,std::cerr);
            std::cerr<<std::endl<<"u vectors"<<std::endl;
            ostream_vector(test->u,std::cerr);
            std::cerr<<std::endl<<"summands of dfOver_n:"<<std::endl;
            for (typename std::vector<internalTypeEigen>::iterator ii=lD.begin();ii!=lD.end();ii++) {
                std::cerr<<(*ii)<<" "<<(1./(1.+(*ii)))<<std::endl;
            }
        }
        if (batchDebug) std::cin.get();
        exit(-1);
#else
{
        std::stringstream stst;
        stst<<"lambda: "<<double(lambda)<<" MSE: "<<double(MSE)<<" dfOver_n: "<<double(dfOver_n)<<std::endl;
        REprintf(stst.str().c_str());
        REprintf("D_invEigVals");
        stst.str("");
        for (typename std::vector<internalTypeEigen>::iterator ii=test->D_invEigVals.begin();ii!=test->D_invEigVals.end();ii++) {
            stst<<double(*ii)<<" ";
        }
        REprintf(stst.str().c_str());
        stst.str("");
        REprintf("summands of dfOver_n:");
        for (typename std::vector<internalTypeEigen>::iterator ii=lD.begin();ii!=lD.end();ii++) {
            stst<<double(*ii)<<" "<<double(1./(1.+(*ii)))<<"; ";
        }
        REprintf(stst.str().c_str());
        stst.str("");

}
Rf_error("(!) (!) From Krig_fgcv() in DLL: !dfOver_n>0.");
#endif
    }

//return(std::numeric_limits<Typeforcov>::quiet_NaN());
return -1;
}

template <typename Typeforcov>
Typeforcov match_fs2hat_pure_error(Typeforcov lambda) {
    //first, obvious, estimate of MSE
    //second, heuristic, estimate computed as RSS/res_df_with_replicates
    long double RSS=0.,tmp;
    Typeforcov res_df_with_replicates=0.,fs2hat;
    std::vector<internalTypeEigen> lD(0);
    for (typename std::vector<internalTypeEigen>::iterator ii=test->D_invEigVals.begin();ii!=test->D_invEigVals.end();ii++)
       lD.push_back((*ii)*lambda);
//ostream_vector(test->D_invEigVals,std::cout);
    int iii=0;
    for (std::vector<internalTypeEigen>::iterator ii=lD.begin();ii!=lD.end();ii++,iii++) {
        tmp=test->u[iii]*(*ii)/(1.+(*ii));
        RSS+=tmp*tmp;
    }
//ostream_vector(test->u,std::cout);
/***************
 Experiments confirm that for FIXED covparam and lambda, the MSE at this point is different whether
 two replicates estimates or only their mean was given as input in pointls.
 Hence the weighted MSE depends on the replicates number, it cannot simply be considered a measure
 of prediction error of the mean values.
***************/
    RSS+=test->pureSS;// sum over replicated positions of the MSE for each replicated position=SSE
    for (std::vector<internalTypeEigen>::iterator ii=lD.begin();ii!=lD.end();ii++)
        res_df_with_replicates+=(1./(1.+(*ii)));  //trA (unfort. not a function of tr[single matrix for all lambda])
    if (res_df_with_replicates>0.) fs2hat=RSS/res_df_with_replicates;
    else { //should have been prevented at the level of choice of bounds for df in gcv_Krig()
#ifdef NO_R_CONSOLE
        std::cerr<<"(!) From match_fs2hat_pure_error() in DLL: !res_df_with_replicates>0."<<std::endl;
        if (false) {
            std::cerr<<"lambda: "<<lambda<<" RSS: "<<RSS<<" res_df_with_replicates: "<<res_df_with_replicates<<std::endl;
            std::cerr<<"D_invEigVals"<<std::endl;
            ostream_vector(test->D_invEigVals,std::cerr);
            std::cerr<<std::endl<<"u vectors"<<std::endl;
            ostream_vector(test->u,std::cerr);
            std::cerr<<std::endl<<"summands of res_df_with_replicates:"<<std::endl;
            for (typename std::vector<internalTypeEigen>::iterator ii=lD.begin();ii!=lD.end();ii++,iii++) {
                std::cerr<<(*ii)<<" "<<(1./(1.+(*ii)))<<std::endl;
            }
        }
        if (batchDebug) std::cin.get();
        exit(-1);
#else
        REprintf("(!) (!) From match_fs2hat_pure_error() in DLL: !res_df_with_replicates>0.\n");
        std::stringstream stst;
        std::string st="";
        st+="RSS: ";
        stst<<RSS;
        st+=stst.str()+" ";
        stst.str("");
        st+="res_df_with_replicates: ";
        stst<<res_df_with_replicates;
        st+=stst.str()+" ";
        stst.str("");
        st+="test->KgPtNbr_with_repl: ";
        stst<<test->KgPtNbr_with_repl;
        st+=stst.str()+" ";
        stst.str("");
        Rprintf("%s\n",st.c_str());
        Rf_error("");
#endif
    }
    fnevalcounter++;
    if (false) {
        std::stringstream stst;
        std::string st="";
        st+="fs2hat: ";
        stst<<fs2hat;
        st+=stst.str()+" ";
        stst.str("");
        st+="pure MSE: ";
        stst<<test->pure_error;
        st+=stst.str()+" ";
        stst.str("");
        st+="return value: ";
        stst<<pow(fs2hat-test->pure_error,(Typeforcov)(2.));
        st+=stst.str()+" ";
        stst.str("");
        Rprintf("%s\n",st.c_str());
    }
    //comparing the two MSE's
    tmp=fs2hat-test->pure_error;
    return(tmp*tmp);
}

#ifdef TEST_OCV
template <typename Typeforcov>
Typeforcov Krig_ocv(Typeforcov lambda) { // returns ordinary leave-one-out CV error.
// (hopefully... gives results close to GCV)
    Typeforcov MSE=0.;
    Typeforcov tmp;
    Typeforcov OCV=0.,trA=0.;
    std::vector<Typeforcov> lD(0);
    int KgPtNbr=test->KgPtNbr;
    for (std::vector<internalTypeEigen>::iterator ii=test->D_invEigVals.begin();ii!=test->D_invEigVals.end();ii++)
       lD.push_back((*ii)*lambda);
    for (std::vector<internalTypeEigen>::iterator ii=lD.begin();ii!=lD.end();ii++)
        trA+=(1./(1.+(*ii)));
    for (int k=0;k<KgPtNbr-1;k++) {
        MSE=test->u[k]*lD[k]/(1.+lD[k]); // numerator, not yet squared
        tmp=0;
        std::vector<Typeforcov>::iterator j=test->TeigVecs2[k].begin(); // kth eigenvector^2 in row
        for (std::vector<internalTypeEigen>::iterator ii=lD.begin();ii!=lD.end();ii++) {
            tmp+=(*j)*(*j)*1./(1.+(*ii)); //\sum_j e_{kj}^2 /(1+lambda D_j) =A_{kk}
        } // this has to match the denominator summands in GCV, which are (1./(1.+(*ii)))
        OCV+=pow(MSE/(1.-tmp),(Typeforcov)(2.));
    }
    // next term as guessed from non constant term in eq 1.8 of GolubHW79
    if (test->pureSS>0) OCV=OCV-2.*(1.-trA/lD.size())*test->pureSS/(test->KgPtNbr_with_repl-test->KgPtNbr);
    fnevalcounter++;
    return(OCV);
/* full mean square error includes + test->pureSS/(test->KgPtNbr_with_repl-test->KgPtNbr)
// and ideally should equal it
// so the above term should approach zero
// and then we can minimize:
    return(pow(OCV,(Typeforcov)(2.)));  */
}
#endif


template<typename Typeforcov>
int CSmooth::fillcovMat(const Typeforcov& locsmoothness) {
    this->smoothness=locsmoothness; //FR->FR important: that gives the smoothness in later calls to CcovFoval; but poorly placed ?
    covMat.resize(KgPtNbr);
    Typeforcov tmp;
    for (typename std::vector<std::vector<Typeforcov> >::iterator ii=covMat.begin();ii<covMat.end();ii++)
       ii->resize(KgPtNbr);

    for (int ii=0;ii<KgPtNbr;ii++) {
      for (int jj=0;jj<=ii;jj++) {   // symetrie importante pour calculs matriciels
        tmp=Matern<Typeforcov>(euclArray[ii][jj],locsmoothness);
        if ( (! R_FINITE(tmp)) || ISNAN(tmp) ) {
#ifdef NO_R_CONSOLE
          std::cerr<<"(!) From CSmooth::fillcovMat(): something wrong with Matern covariance evaluation" << std::endl;
          std::cerr<<"i,j: "<<ii<<" "<<jj<<" Euclidian distance: "<<euclArray[ii][jj]<<std::endl;
          std::cerr<<"Smoothness: "<<locsmoothness<<std::endl;
          std::cerr<<tmp<<std::endl;
          if (batchDebug) std::cin.get();
          exit(-1);
#else
          Rf_error("(!) From CSmooth::fillcovMat(): something wrong with Matern covariance evaluation");
#endif
        }
        covMat[ii][jj]=covMat[jj][ii]=tmp; //iint<=jint semi matrix
        //if (covMat[ii][jj]<0) {std::stringstream stst;stst<<"ii, jj, covMat[ii][jj]<0 !: "<<ii<<" "<<jj<<" "<<covMat[ii][jj]<<std::endl;REprintf(stst.str().c_str());}
      }
    }

//std::cout<<std::endl<<std::endl<<std::endl<<"[";
//ostream_vec_vector(covMat,std::cout);
//std::cout<<"]";
return 0;
}



template<typename Typeforcov>
int CSmooth::fillcovFocal() {
        covFocal.resize(KgPtNbr);
  for (int jj=0;jj<KgPtNbr;jj++) {
    covFocal[jj]=Matern<Typeforcov>(euclFocal[jj],smoothness);
  }
  return 0;
}

template<typename Typeforcov,typename TypeforES>
int CSmooth::Krig_coef(Typeforcov lambda) { //constructs predictor in compact form
//NOTE default lambda (NaN) value in function prototype
// ne depend de tout le bins en GCV qu'? travers lambda
// reste des calculs uniquement bas['e] sur covariances, QR decomp, et yobs

//DEBUG // still output
//note D_invEigVals[0]=0 always...
{
    std::stringstream stst1;stst1<<"(!) Kriging covariance matrix possibly ill-conditioned"<<std::endl;
    // problem with Mingw64: stst<< long double << displays nonsense... at least through REprintf...
    std::stringstream stst2;stst2<<"       (Eigenvalues max: "<<double(1/D_invEigVals[1])<<"; min: "<<double(1/D_invEigVals.back())<<")."<<std::endl;
#ifdef NO_R_CONSOLE
   if (D_invEigVals.back()/D_invEigVals[1]>100000000000.) std::cout<<stst1.str();
   if (verbosity) std::cout<<stst2.str();
#else
   if (D_invEigVals.back()/D_invEigVals[1]>1e11) REprintf(stst1.str().c_str());
   if(verbosity) REprintf(stst2.str().c_str());
#endif
}

    if (ISNAN(lambda)) {// no explict value provided as fn argument
        lambda=lambdaEst;
        if (ISNAN(lambda)) { //neither explicit argument nor GCV value
#ifdef NO_R_CONSOLE
           std::cerr<<"(!) From Krig_coef(): lambda neither explicitly given nor previously computed by GCV";
           if (batchDebug) std::cin.get();
           exit(-1);
#else
           Rf_error("(!) From Krig_coef(): lambda neither explicitly given nor previously computed by GCV\n");
#endif
        }
    }
    if (lambda<(10.+xy.size())*std::numeric_limits<Typeforcov>::epsilon()) {
#ifdef NO_R_CONSOLE
        std::cerr<<std::endl<<"(!) From Krig_coef: low lambda value likely to lead to numerical artefacts"<<std::endl;
#else
        REprintf("(!) From Krig_coef: low lambda value likely to lead to numerical artefacts\n");
#endif
    }
    int idx;
    std::vector<TypeforES> temp(0);
    std::vector<TypeforES> temp_c(0); // will contain fields' beta then fields' temp.c
// computation of coefs for random part of predictor
    typename std::vector<TypeforES>::iterator it2=D_invEigVals.begin();
    for (typename std::vector<TypeforES>::iterator it1=u.begin();it1<u.end();it1++,it2++) {
      temp.push_back((*it1)/(1.+lambda*(*it2)));
    }
    for (typename std::vector<std::vector<TypeforES> >::iterator it=G_for_Krig_c_coef.begin();it<G_for_Krig_c_coef.end();it++) {
      temp_c.push_back(inner_product((*it).begin(),(*it).end(),temp.begin(),0.));
    }
    temp_c.erase(temp_c.begin(),temp_c.begin()+ncolT);
    temp_c.insert(temp_c.begin(),ncolT,0); // seems to me that beta is unchanged after the last two ops... but this is the R code
    // result goes in coeff_random:
    QR_T->Qy<TypeforES,TypeforES>(temp_c,coefs_random);// now this is temp.c in fields
    for (int ii=0;ii<KgPtNbr;ii++) {
      coefs_random[ii]*=W2[ii];
    }
//std::cout<<coefs_random[0]<<std::endl;

// computation of coefs for fixed part of predictor
    coefs_fixed.resize(0);
    //matrix vector multiplication...
    idx=0;
    for (typename std::vector<std::vector<Typeforcov> >::iterator ii=covMat.begin();ii<covMat.end();ii++,idx++) {
       coefs_fixed.push_back(-inner_product((*ii).begin(),(*ii).end(),coefs_random.begin(),-yobs[idx]));
    }//this has computed -(-y+covMat.tempc); //(size=ptNbt)
    for (int ii=0;ii<KgPtNbr;ii++) {
      coefs_fixed[ii]*=W2[ii];
    }
    // least square fit of input-value to return-value, input-value=T.(return-value)+e
    QR_T->coef<TypeforES>(coefs_fixed); ///output container = input;  (size=nrow mais seuls les ncolT premiers sont non-nuls)
return 0;
}



template <typename Typeforcov>
ioType CSmooth::predict(std::vector<ioType> focal,std::string method) {
   fillaxialFocal(focal);
   filleuclFocal();
   fillcovFocal<Typeforcov>();
//for(int ii=0;ii<10;ii++) {std::cout<<covFocal[ii]<<" "<<coefs_random[ii];getchar();}
   //la fixed part est cod['e]e comme Tmatrix %*% (coefs_fixed). Pour l'instant se ram[`e]ne ? coef_fixed[0]
   ioType sol=std::inner_product(covFocal.begin(),covFocal.end(),coefs_random.begin(),coefs_fixed[0]);
   //fixed part + covvec.temp_c
return sol;
}

template <typename Typeforcov,typename TypeforES>
Typeforcov CSmooth::GCV_lamVar_covFix(std::vector<Typeforcov> covparam) {
  //for(typename std::vector<Typeforcov>::iterator it=covparam.begin();it!=covparam.end();it++) { Rprintf("zo %f ",(*it)); }
  std::vector<Typeforcov> covtheta(covparam);
   covtheta.erase(--covtheta.end());
   Krig_engine_default<Typeforcov,TypeforES>(covtheta,covparam.back());
   return(gcv_Krig<Typeforcov>()); //returns GCV value at lambdaEst for fixed covparam
}


template <typename Typeforcov,typename TypeforES> // second type indeed the second type of CEigensystem
int CSmooth::Krig_engine_default(std::vector<Typeforcov> CovTheta, const Typeforcov& locsmoothness) {
//TypeforES zut;
TypeforES tmp;
bool verbose=false;
std::vector<std::vector<Typeforcov> >tempM(KgPtNbr);
std::vector<TypeforES> temp(0);
#ifdef NO_R_CONSOLE
    if (verbose) std::cout<<"debut Krig_engine_default()\n";
#else
    if (verbose) Rprintf("debut Krig_engine_default()\n");
#endif
    CovFnParam=CovTheta;CovFnParam.push_back(locsmoothness); // this sets a default value if not further estimated by minimiseGCV()

    //for (typename std::vector<Typeforcov>::iterator ii=CovFnParam.begin();ii!=CovFnParam.end();ii++) Rprintf("ici %f ",(*ii));


    ncolT=1;
    for (typename std::vector<std::vector<Typeforcov> >::iterator ii=tempM.begin();ii!=tempM.end();ii++)
        (*ii).resize(KgPtNbr);
    this->CovTheta=CovTheta;
    CovTheta2=CovTheta;
    for(typename std::vector<Typeforcov>::iterator it=CovTheta2.begin();it!=CovTheta2.end();it++) { (*it)=(*it)*(*it); }
    filleuclArray();
    fillcovMat<Typeforcov>(locsmoothness);
    Tmatrix.resize(KgPtNbr); // typically KgPtNbr X 1
    for (unsigned int ii=0;ii<Tmatrix.size();ii++) {
        Tmatrix[ii].resize(0);
        Tmatrix[ii].push_back(sqrt(Typeforcov(W[ii])));
        //Rprintf("la %u", W[ii]);
    }
     if (!QR_Tallocated) {
        QR_T=new CQR<CQRtype>(Tmatrix);    //QR_T is the QR decomp of T
        QR_Tallocated=true; // for safe destruction at destruction of (*this)
     }
//std::cout<<W2.size()<<" "<<;getchar();
     for (int ii=0;ii<KgPtNbr;ii++)
         for (int jj=0;jj<KgPtNbr;jj++)
             tempM[ii][jj]=covMat[ii][jj]*W2[ii]*W2[jj];
//ostream_vec_vector(tempM,std::cout);getchar();
#ifdef NO_R_CONSOLE
    if (verbose) std::cout<<"debut tempM=QR_T->Qty(tempM);\n";
#else
    if (verbose) Rprintf("debut tempM=QR_T->Qty(tempM);\n");
#endif
     QR_T->QtY<Typeforcov>(tempM); /// tempM both input and output; qr.qty(qr.T,t(tempM)) en utilisant la symetrie de tempM
     //Rprintf("tempM avant: %f %f %f %f %f %f\n",double(tempM[0][0]),double(tempM[0][1]),double(tempM[1][0]),
     //       double(tempM[KgPtNbr-2][KgPtNbr-2]),double(tempM[KgPtNbr-3][KgPtNbr-2]),double(tempM[KgPtNbr-2][KgPtNbr-3]));
     //ostream_vec_vector(tempM,std::cout);getchar();
//NB pour y retrouver la sortie de qr.yq2(qr.T, tempM) oublier la 1e ligne et commencer par le dernier pt
     tempM.erase(tempM.begin(),tempM.begin()+ncolT); // erases first ncolT ROWS
#ifdef NO_R_CONSOLE
    if (verbose) std::cout<<"debut tempM=QR_T->Qtyt(tempM);\n";
#else
    if (verbose) Rprintf("debut tempM=QR_T->Qtyt(tempM);\n");
#endif
    MatrixXd for_es;
     for_es=QR_T->Qtyt<Typeforcov,Typeforcov>(tempM); //...Q'.*transpose*(tempM)
     //Rprintf("QR_T->Qtyt: %f %f %f %f %f %f\n",for_es(0,0),for_es(0,1),for_es(1,0),
    //         for_es(KgPtNbr-2,KgPtNbr-2),for_es(KgPtNbr-3,KgPtNbr-2),for_es(KgPtNbr-2,KgPtNbr-3));
     // erases first ncolT ROWS:
     unsigned int rowMax=for_es.rows();
     unsigned int colMax=for_es.cols();
     // last indices describe block of size rowMax-ncolT,colMax
     for_es.block(0,0,rowMax-ncolT,colMax) = for_es.block(ncolT,0,rowMax-ncolT,colMax);
     for_es.conservativeResize(rowMax-ncolT,colMax); // erases lines beyond
     //Rprintf("resize: %f %f %f %f %f %f\n",for_es(0,0),for_es(0,1),for_es(1,0),
     //        for_es(KgPtNbr-2,KgPtNbr-2),for_es(KgPtNbr-3,KgPtNbr-2),for_es(KgPtNbr-2,KgPtNbr-3));
#ifdef NO_R_CONSOLE
    if (verbose) std::cout<<"debut CEigensystem<Typeforcov,TypeforES> es(tempM,true)\n";
#else
    if (verbose) Rprintf("debut CEigensystem<Typeforcov,TypeforES> es(tempM,true)\n");
#endif
    es.compute(for_es);
    VectorXd tempXd=es.eigenvalues();
//ostream_vector(temp,std::cout);
/*{
    stringstream stst;
    stst<<"temp.size(): "<<temp.size()<<std::endl;
#ifdef NO_R_CONSOLE
    if (true) Rprintf(stst.str().c_str());
#else
    if (true) REprintf(stst.str().c_str());
#endif
}*/
     D_invEigVals.resize(ncolT,0);  // ncolT=1, puts initial 0.
     for (int ii=0; ii< int(tempXd.size());ii++) {
        if ((tempXd(ii))>0) { // test added 31/01/2012 .. negative eigenvalues only result from numerical errors
            /* the eigenvalues are those of a bloc from an Hermitian matrix 'H', the latter which eigenvalues are those of the original covmat. Hence the computed eigenvalues are interlaced within
            those of 'H', hence > the smallest of the eigenvalues of the covmat*/
            D_invEigVals.push_back(1./(tempXd(ii))); //D<-c(rep(0,nt),1/tempM$values)
        } else {
            D_invEigVals.push_back(1e+30);
        } //small eigenval => large 1/eigenval !
/* {
    std::stringstream stst;
    stst<<"double(D_invEigVals.back()): "<<double(D_invEigVals.back())<<" 1./tempXd(ii): "<<1./tempXd(ii)<<std::endl;
#ifdef NO_R_CONSOLE
    if (true) Rprintf(stst.str().c_str());
#else
    if (true) REprintf(stst.str().c_str());
#endif
} */
     } /// end loop on ii
     MatrixXd eigVecs=es.eigenvectors(); // (KgPtNbr-1)*(KgPtNbr-1) //,contrary to svd(), the code ensures that for symm matrices A=u.d.t(u)
     //Rprintf("EV: %f %f %f %f %f %f\n",eigVecs(0,0),eigVecs(0,1),eigVecs(1,0),
     //         eigVecs(KgPtNbr-2,KgPtNbr-2),eigVecs(KgPtNbr-3,KgPtNbr-2),eigVecs(KgPtNbr-2,KgPtNbr-3));
     //Rprintf("EV: %f %f %f %f %f %f\n",eigVecs(0,KgPtNbr-2),eigVecs(0,KgPtNbr-3),eigVecs(KgPtNbr-3,0),
     //         eigVecs(KgPtNbr-2,0),eigVecs(KgPtNbr-3,0),eigVecs(0,KgPtNbr-3));
#ifdef TEST_OCV
     TeigVecs2.resize(KgPtNbr-1);
     for (int jj=0;jj<KgPtNbr-1;jj++) {
        TeigVecs2[jj].resize(KgPtNbr-1);
        for (int kk=0;kk<KgPtNbr-1;kk++) {
            tmp=eigVecs(kk,jj);
            TeigVecs2[jj][kk]=tmp*tmp; // note transposition
        }
     }
#endif
//     ostream_vec_vector(eigVecs,std::cout);
     G_for_Krig_c_coef.resize(KgPtNbr);
     G_for_Krig_c_coef[0].resize(KgPtNbr,0); // pour le cas ncolT=1 !!!
     for (int ii=1;ii<KgPtNbr;ii++) { // 1 !!!
       G_for_Krig_c_coef[ii].resize(KgPtNbr,0);
         for (int jj=1;jj<KgPtNbr;jj++) { // -1 !!!
           G_for_Krig_c_coef[ii][jj]=eigVecs(ii-1, jj-1)*D_invEigVals[jj];
         } //G[(nt+1):np,(nt+1):np]<-tempM$vectors
     } //G[(nt+1):np,(nt+1):np]<-tempM$vectors
     temp.resize(0);
     for (int ii=0;ii<KgPtNbr;ii++) {
       temp.push_back(W2[ii]*yobs[ii]); // out$W2 %*d% out$yM
       //Rprintf("| %f %f %f|",W2[ii],yobs[ii],temp[ii]);
     }
#ifdef NO_R_CONSOLE
    if (verbose) std::cout<<"debut appel QR_T->Qy(temp,temp)\n";
#else
    if (verbose) Rprintf("debut appel QR_T->Qy(temp,temp)\n");
#endif
    //for (typename std::vector<TypeforES>::iterator ii=temp.begin();ii!=temp.end();ii++) Rprintf("la %f ",(*ii));

    QR_T->Qy<TypeforES,TypeforES>(temp,temp);  // <..> = the types of temp; note that the second argument is the (address of) the result.

     //for (typename std::vector<TypeforES>::iterator ii=temp.begin();ii!=temp.end();ii++) Rprintf("ici %f ",(*ii));
     //Rf_error("FR stop debug");

     temp.erase(temp.begin(),temp.begin()+1); //cf def qr.q2ty:  qr.qty(qr, y)[(rank + 1):dy[1], ]
//    ostream_vector(temp,std::cout);
     // where rank of T must be one... qy is the # of rows of (column vector) y
     // well, we could count nonzero eigenvalues since they directly play in the next computation
     int eigdim=eigVecs.cols();    // not .size(); !!!!!!!
     //Rprintf("ICI %u \n",eigdim);
     u.resize(1,0); // 1!
     for (int ii=0;ii<eigdim;ii++) {
         tmp=0;
         for (int jj=0;jj<eigdim;jj++) {
             tmp+=eigVecs(jj,ii)*temp[jj];
         }
         u.push_back(tmp);
     }
/*    std::cout<<std::endl;
    ostream_vector(u,std::cout);
    std::cout<<std::endl;
    ostream_vector(D_invEigVals,std::cout);
    std::cout<<std::endl;
    ostream_vec_vector(G,std::cout);
    std::cout<<std::endl;
*/
#ifdef NO_R_CONSOLE
    if (verbose) std::cout<<"fin Krig_engine_default()\n";
#else
    if (verbose) Rprintf("fin Krig_engine_default()\n");
#endif
return 0; //D_invEigVals, eigVecs, u et G sont conserv['e]s comme variables membres.
}




template<typename Typeforcov>
Typeforcov CSmooth::Krig_df_to_lambda(Typeforcov explicitdf) { // default argument value is -1.
    if (explicitdf>0) df=explicitdf; // else df (a member of CSmooth!) keeps it current value
    Typeforcov (*fdfPtr)(Typeforcov);
    fdfPtr=&Krig_fdf<Typeforcov>;
    Typeforcov l1=1.;
    Typeforcov l2=1.;
    Typeforcov tr;
    for (int k=1;k<26;k++) {
        tr=0.;
        for (typename std::vector<internalTypeEigen>::iterator ii=D_invEigVals.begin();ii!=D_invEigVals.end();ii++)
            tr+=1./(1.+l1*(*ii));
        if (tr <= df) break;
        l1*=4;
    }
    for (int k=1;k<26;k++) {
        tr=0.;
        for (typename std::vector<internalTypeEigen>::iterator ii=D_invEigVals.begin();ii!=D_invEigVals.end();ii++)
            tr+=1./(1.+l2*(*ii));
        if (tr >= df) break;
        l2/=4;
    }
    return(bisection_search(fdfPtr, log(l1), log(l2))); // fdfPtr(log(l1))*fdfPtr(log(l2)) must be <0
}

#ifdef TEST_OCV
template<typename Typeforcov>
Typeforcov CSmooth::ocv_Krig() {//estimates lambda by OCV, returns OCV criterion
//    Typeforcov (*fgcvPtr)(Typeforcov) = &Krig_ocv;
    fgcvPtr=&Krig_ocv<Typeforcov>;
    Typeforcov objective;
    int maxit=100;
    int it=0;
    Typeforcov tmp;
    Typeforcov hi,lo=1; // lambda value here
    Typeforcov previous=(*fgcvPtr)(lo);
//starts from some value, seek a low lambda value at which objective fn increases
    lo/=10.;
    while (it<maxit && (tmp=(*fgcvPtr)(lo))<previous) {previous=tmp;lo/=10.;it++;}
    if (it==maxit) {
#ifdef NO_R_CONSOLE
        std::cout<<"(!) From CSmooth::ocv_Krig(): CV search gives a minimum at the lower endpoint of the search.";
#else
        REprintf("(!) From CSmooth::ocv_Krig(): CV search gives a minimum at the lower endpoint of the search.\n");
#endif
        lambdaEst=lo;
        objective=(*fgcvPtr)(lo);
        return(objective);
    }
// then go backwards; previous is at the minimum and lo lower than minimum
    it=0;
    hi=lo*100;
    while (it<maxit && (tmp=(*fgcvPtr)(hi))<previous) {previous=tmp;hi*=10.;it++;}
    if (it==maxit) {
#ifdef NO_R_CONSOLE
        std::cout<<"(!) From CSmooth::ocv_Krig(): CV search gives a minimum at the upper endpoint of the search.";
#else
        REprintf("(!) From CSmooth::ocv_Krig(): CV search gives a minimum at the upper endpoint of the search.\n");
#endif
        lambdaEst=hi;
        objective=(*fgcvPtr)(hi);
        return(objective);
    }
//ELSE
//DEBUG
//{ std::stringstream stst;stst<<lo<<" "<<hi<<std::endl;REprintf(stst.str().c_str());}
    objective=brent<brentTypedef>(fgcvPtr,lo,hi/10.,hi,&lambdaEst); //fills lambdaEst
return objective;
}
#endif



template<typename Typeforcov>
Typeforcov CSmooth::gcv_Krig() {//estimates lambda by GCV, returns GCV criterion
#ifndef NO_R_CONSOLE
    R_CheckUserInterrupt();
#endif
    Typeforcov objective;
    grid_df.resize(0);
    gcv_grid.resize(0);
    lambda_grid.resize(0);
    Typeforcov hi=KgPtNbr*0.95;
    Typeforcov lo=1; //# col T dans fields ??
    if (hi<lo) {
#ifdef NO_R_CONSOLE
       std::cerr<<"(!) From CSmooth::gcv_Krig(): problem with bounds. Seek this message in source and compare to Migraine code";
       if (batchDebug) std::cin.get();
       exit(-1);
#else
       Rf_error("(!) From CSmooth::gcv_Krig(): problem with bounds. Seek this message in source and compare to Migraine code\n");
#endif
    }
    lo+=(hi-lo)/80000.;
/** Krig_df_to_lambda(hi) calls bisection_search that finds a lambda value that 'fits' df=hi;
    exp(Krig_df_to_lambda(hi) is the lambda value matching the df value.
    The aim of the first loop is to find a df value low enough that the matching lambda is
    small enough according to the (retrospectively dubious) loop test criterion
    The aim of the second loop is to find a df value high enough that the matching lambda is small enough... */
    //{ std::stringstream stst;stst<<lo<<" "<<hi<<std::endl;REprintf(stst.str().c_str());}
    Typeforcov step=(hi-lo)/30.;
    df=lo;
    for (int jj=0;jj<8;jj++,df+=step*pow(double(2),double(std::min(jj,7-jj)))) { //adds steps* 1 2 4 8 8 4 2 1 (sum=30)
        grid_df.push_back(df);//df made a class member...
        Typeforcov tmp=exp(Krig_df_to_lambda<Typeforcov>()); //donne lambda(df) (bisection search)
        lambda_grid.push_back(tmp);
        gcv_grid.push_back(tmp);
    }
    //next few lines ~ Krig.find.gcvmin
    for (typename std::vector<Typeforcov>::iterator ii=gcv_grid.begin();ii!=gcv_grid.end();ii++)
        (*ii)=(*fgcvPtr)((*ii)); //gcv values for the grid of lambda values
    typename std::vector<Typeforcov>::iterator min_ptr=min_element(gcv_grid.begin(),gcv_grid.end());
    int minidx=int(min_ptr-gcv_grid.begin());
    if ( min_ptr!=gcv_grid.begin() && min_ptr!=(--gcv_grid.end()) ) {
        objective=brent<brentTypedef>(fgcvPtr,lambda_grid[minidx-1],lambda_grid[minidx],lambda_grid[minidx+1],&lambdaEst); //fills lambdaEst
    } else {
#ifdef NO_R_CONSOLE
      if (verbosity>1) std::cout<<"(*) From CSmooth::gcv_Krig(): GCV search gives a minimum at the endpoints of the grid search.";
#else
     if (verbosity>1) REprintf("(*) From CSmooth::gcv_Krig(): GCV search gives a minimum at the endpoints of the grid search.\n");
#endif
        lambdaEst=lambda_grid[minidx]; //!! lambda_GCV <=> phi_HGLM/lambda_HGLM
        objective=(*fgcvPtr)(*min_ptr);
    }
//std::cout<<"lambdaEst="<<lambdaEst;//getchar();
//std::cout<<objective;
return objective;
}


#endif
