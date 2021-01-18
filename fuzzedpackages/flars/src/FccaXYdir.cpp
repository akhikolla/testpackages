#include <Rmath.h>
#include <RcppEigen.h>


// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                 
using Eigen::MatrixXd;            
using Eigen::VectorXd;            
using Eigen::SelfAdjointEigenSolver;
using Eigen::GeneralizedSelfAdjointEigenSolver;
using Eigen::MatrixXi;
using Eigen::Lower;
using Eigen::ArrayXd;
using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::export]]
List FccaXYdir0(List Lx, List Lp, List LPhi, NumericVector resp, 
      VectorXd l1, VectorXd l2, int cv=1){
  int nx=Lx.size();
  const Map<MatrixXd> M=Lx[0];
  const int nrowx(M.rows());
  NumericVector y(resp);
  
  VectorXd xidx(VectorXd(nx).setZero());
  //find the ncol for each of x(t) in th list
  for(int i=0;i<nx;i++){
    const Map<MatrixXd> M=Lx[i];
    const double ncolx(M.cols());
    xidx[i]=ncolx;
  }
  VectorXd PHIidx(VectorXd(nx).setZero());
  //find the ncol for each of phi in th list
  for(int i=0;i<nx;i++){
    const Map<MatrixXd> Phi=LPhi[i];
    const double ncolPhi(Phi.cols());
    PHIidx[i]=ncolPhi;
  }
  
  MatrixXd x = MatrixXd::Ones(nrowx,xidx.sum());
  
  VectorXd cumsumxidx(VectorXd(nx+1).setZero());
  //find the cummetive ncol for each of x(t) in th list
  int acc = 0;
  cumsumxidx[0]=acc;
  for(int i = 1; i < cumsumxidx.size(); i++){
     acc += xidx[i-1];
     cumsumxidx[i] = acc;
   }
  
  VectorXd cumsumPHIidx(VectorXd(nx+1).setZero());
  //find the cummetive ncol for each of phi in th list
  int accP = 0;
  cumsumPHIidx[0]=accP;
  for(int i = 1; i < cumsumPHIidx.size(); i++){
     accP += PHIidx[i-1];
   	cumsumPHIidx[i] = accP;
 	}
  
  
  //centre both x and y
  y=y-mean(y);
  MatrixXd eigenY(x.rows(),1);
  for(int i=0;i<x.rows();i++){
    eigenY(i,0)=y[i];
  }
  
  List Z(nx), zy(nx);
  for(int i=0;i<nx;i++){
    MatrixXd M=Lx[i];
    MatrixXd Phi=LPhi[i];
    MatrixXd meanMi(M.colwise().mean());
    M=M.rowwise() - M.colwise().mean();
    
    MatrixXd ZM=M*Phi;
    Z[i]=ZM;
    zy[i]=ZM.adjoint()*eigenY;
  }
  MatrixXd cz = MatrixXd::Ones(nrowx,PHIidx.sum());
  MatrixXd czy = MatrixXd::Ones(PHIidx.sum(),1);
  //cz is cbind z
  //czy is connect zy
  for(int p=0;p<nx;p++){
    int tmp1=cumsumPHIidx[p];
    int tmp2=cumsumPHIidx[p+1];
    MatrixXd ZM0=Z[p];
    MatrixXd zy0=zy[p];
    for(int i=tmp1;i<tmp2;i++){
      czy(i,0)=zy0(i-tmp1,0);
      for(int j=0;j<nrowx;j++){
        cz(j,i)=ZM0(j,i-tmp1);
      }
    }
  }
  MatrixXd Vx=cz.adjoint()*cz;
  
  List L2L1(nx);
  List L2L2(nx);
  
  double GCV=-1; double l_1=0.1; double TrHat1=0; double TrHat2=0;
  double GCV2=-1; double l_2cv=l2(0); double l_1cv=l1(0);
  MatrixXd idI=MatrixXd::Identity(Vx.rows(), Vx.cols());
  MatrixXd H1_idI_H1=MatrixXd::Identity(eigenY.size(), eigenY.size());
  MatrixXd H1_idI_H1cv=MatrixXd::Identity(eigenY.size(), eigenY.size());
  MatrixXd H_gcv=MatrixXd::Identity(eigenY.size(), eigenY.size());
  MatrixXd H_gcv2=MatrixXd::Identity(eigenY.size(), eigenY.size());  
  MatrixXd VxCV=Vx;
  
  MatrixXd GCV_mat=MatrixXd::Identity(l2.size(),l1.size());
  
  //build penalty matrix list
  for(int j2=0;j2<l2.size();j2++){
    double l_2=l2(j2);
    GCV=-1; l_1=0.1;
    
    for(int i=0;i<nx;i++){
      MatrixXd l0=Lp[i];
      MatrixXd id0 = MatrixXd::Identity(PHIidx[i], PHIidx[i]);
      
      L2L1[i]=l0.adjoint() * l0;
      L2L2[i]=l_2*id0.adjoint() * id0;
    }
    
    MatrixXd K1=MatrixXd::Identity(Vx.rows(), Vx.cols());
    MatrixXd K2=MatrixXd::Identity(Vx.rows(), Vx.cols());
    
    
    
    //joint penalty matrices on the diagonal block matrices
    for(int p=0;p<nx;p++){
      int tmp1=cumsumPHIidx[p];
      int tmp2=cumsumPHIidx[p+1];
      const Map<MatrixXd> ML2L1=L2L1[p];
      const Map<MatrixXd> ML2L2=L2L2[p];
      for(int i=tmp1;i<tmp2;i++){
        for(int j=tmp1;j<tmp2;j++){
          K1(j,i)=ML2L1(j-tmp1,i-tmp1);
          K2(j,i)=ML2L2(j-tmp1,i-tmp1);
        }
      }
    }
    
    Vx=Vx+K2;
    
    if(cv>0){
      MatrixXd I=MatrixXd::Identity(eigenY.size(), eigenY.size());
      GeneralizedSelfAdjointEigenSolver<MatrixXd> es(K1, Vx);
      MatrixXd d(es.eigenvalues());
      MatrixXd dm(es.eigenvectors());
      
      MatrixXd H1=dm.adjoint()*cz.adjoint();
      
      for(int i=0; i<l1.size();i++){
        double lam1=l1(i);
        for(int j=0; j<idI.rows();j++){
          idI(j,j)=1/(1+lam1*d(j,0));
        }
        MatrixXd H=H1.adjoint()*idI*H1;
        VectorXd trHV(VectorXd(H.rows()).setZero());
        for(int k=0; k<eigenY.size(); k++){
          trHV(k)=H(k,k);
        }
        double TrHat=trHV.sum();
        double C=(1-TrHat/nrowx);
        C=nrowx*C*C;
        double GCV_tmp=(eigenY.adjoint()*(I-H)*(I-H)*eigenY/C).sum();
        GCV_mat(j2,i)=GCV_tmp;
        if(GCV<0){
          GCV=GCV_tmp;
          l_1=lam1;
          H_gcv=H;
          H1_idI_H1=dm*idI*dm.adjoint();
          TrHat1=TrHat;
        }
        if(GCV>0){
          if(GCV_tmp<GCV){
            GCV=GCV_tmp;
            l_1=lam1;
            H_gcv=H;
            H1_idI_H1=dm*idI*dm.adjoint();
            TrHat1=TrHat;
          }
        }
      }
    }
    
    Vx=Vx+l_1*K1;
    
    if(GCV2<0){
      GCV2=GCV;
      l_1cv=l_1;
      H_gcv2=H_gcv;
      H1_idI_H1cv=H1_idI_H1;
      l_2cv=l_2;
      VxCV=Vx;
      TrHat2=TrHat1;
    }
    if(GCV2>0){
      if(GCV<GCV2){
        GCV2=GCV;
        l_1cv=l_1;
        H_gcv2=H_gcv;
        H1_idI_H1cv=H1_idI_H1;
        l_2cv=l_2;
        VxCV=Vx;
        TrHat2=TrHat1;
      }
    }
  }  
  VectorXd betahat(H1_idI_H1cv*czy);
  
  return List::create(Named("betahat")=betahat,
                      Named("K")=VxCV,
                      Named("S")=H_gcv2,
                      Named("lambda1")=l_1cv,
                      Named("lambda2")=l_2cv,
                      Named("M")=M,
                      Named("GCV")=GCV2,
                      Named("GCV_mat")=GCV_mat,
                      Named("HatMatTrace")=TrHat2
                      );
}

// [[Rcpp::export]]
VectorXd FccaXYdir(List Lx, List Lp, List LPhi, NumericVector resp, 
      VectorXd l1, VectorXd l2, int cv=1){
  int nx=Lx.size();
  const Map<MatrixXd> M=Lx[0];
  const int nrowx(M.rows());
  NumericVector y(resp);
  
  VectorXd xidx(VectorXd(nx).setZero());
  //find the ncol for each of x(t) in th list
  for(int i=0;i<nx;i++){
    const Map<MatrixXd> M=Lx[i];
    const double ncolx(M.cols());
    xidx[i]=ncolx;
  }
  VectorXd PHIidx(VectorXd(nx).setZero());
  //find the ncol for each of phi in th list
  for(int i=0;i<nx;i++){
    const Map<MatrixXd> Phi=LPhi[i];
    const double ncolPhi(Phi.cols());
    PHIidx[i]=ncolPhi;
  }
  
  MatrixXd x = MatrixXd::Ones(nrowx,xidx.sum());
  
  VectorXd cumsumxidx(VectorXd(nx+1).setZero());
  //find the cummetive ncol for each of x(t) in th list
  int acc = 0;
  cumsumxidx[0]=acc;
  for(int i = 1; i < cumsumxidx.size(); i++){
     acc += xidx[i-1];
   	cumsumxidx[i] = acc;
 	}
  
  VectorXd cumsumPHIidx(VectorXd(nx+1).setZero());
  //find the cummetive ncol for each of phi in th list
  int accP = 0;
  cumsumPHIidx[0]=accP;
  for(int i = 1; i < cumsumPHIidx.size(); i++){
     accP += PHIidx[i-1];
   	cumsumPHIidx[i] = accP;
 	}
  
  
  //centre both x and y
  y=y-mean(y);
  MatrixXd eigenY(x.rows(),1);
  for(int i=0;i<x.rows();i++){
    eigenY(i,0)=y[i];
  }
  
  List Z(nx), zy(nx);
  for(int i=0;i<nx;i++){
    MatrixXd M=Lx[i];
    MatrixXd Phi=LPhi[i];
    MatrixXd meanMi(M.colwise().mean());
    for(int iCol=0;iCol<M.cols();iCol++){
      for(int iRow=0;iRow<M.rows();iRow++){
        M(iRow,iCol)=M(iRow,iCol)-meanMi(0,iCol);
      }
    }
    
    MatrixXd ZM=M*Phi;
    Z[i]=ZM;
    zy[i]=ZM.adjoint()*eigenY;
  }
  
  MatrixXd cz = MatrixXd::Ones(nrowx,PHIidx.sum());
  MatrixXd czy = MatrixXd::Ones(PHIidx.sum(),1);
  //cz is cbind z
  //czy is connect zy
  for(int p=0;p<nx;p++){
    int tmp1=cumsumPHIidx[p];
    int tmp2=cumsumPHIidx[p+1];
    MatrixXd ZM0=Z[p];
    MatrixXd zy0=zy[p];
    for(int i=tmp1;i<tmp2;i++){
      czy(i,0)=zy0(i-tmp1,0);
      for(int j=0;j<nrowx;j++){
        cz(j,i)=ZM0(j,i-tmp1);
      }
    }
  }
  MatrixXd Vx=cz.adjoint()*cz;
  

  List L2L1(nx);
  List L2L2(nx);
  
  double GCV=-1; double l_1=0.1;
  double GCV2=-1; 
  MatrixXd idI=MatrixXd::Identity(Vx.rows(), Vx.cols());
  MatrixXd H1_idI_H1=MatrixXd::Identity(eigenY.size(), eigenY.size());
  MatrixXd H1_idI_H1cv=MatrixXd::Identity(eigenY.size(), eigenY.size());
  MatrixXd H_gcv=MatrixXd::Identity(eigenY.size(), eigenY.size());
  MatrixXd H_gcv2=MatrixXd::Identity(eigenY.size(), eigenY.size());  
  MatrixXd VxCV=Vx;
  
  MatrixXd GCV_mat=MatrixXd::Identity(l2.size(),l1.size());
  
  //build penalty matrix list
  
  for(int j2=0;j2<l2.size();j2++){
    double l_2=l2(j2);
    GCV=-1; l_1=0.1;
    
    for(int i=0;i<nx;i++){
      MatrixXd l0=Lp[i];
      MatrixXd id0 = MatrixXd::Identity(PHIidx[i], PHIidx[i]);
      
      L2L1[i]=l0.adjoint() * l0;
      L2L2[i]=l_2*id0.adjoint() * id0;
    }
    
    MatrixXd K1=MatrixXd::Identity(Vx.rows(), Vx.cols());
    MatrixXd K2=MatrixXd::Identity(Vx.rows(), Vx.cols());
    //joint penalty matrices on the diagonal block matrices
    for(int p=0;p<nx;p++){
      int tmp1=cumsumPHIidx[p];
      int tmp2=cumsumPHIidx[p+1];
      const Map<MatrixXd> ML2L1=L2L1[p];
      const Map<MatrixXd> ML2L2=L2L2[p];
      for(int i=tmp1;i<tmp2;i++){
        for(int j=tmp1;j<tmp2;j++){
          K1(j,i)=ML2L1(j-tmp1,i-tmp1);
          K2(j,i)=ML2L2(j-tmp1,i-tmp1);
        }
      }
    }
    
    Vx=Vx+K2;
    
    if(cv>0){
      MatrixXd I=MatrixXd::Identity(eigenY.size(), eigenY.size());
      GeneralizedSelfAdjointEigenSolver<MatrixXd> es(K1, Vx);
      MatrixXd d(es.eigenvalues());
      MatrixXd dm(es.eigenvectors());
      
      MatrixXd H1=dm.adjoint()*cz.adjoint();
      
      for(int i=0; i<l1.size();i++){
        double lam1=l1(i);
        for(int j=0; j<idI.rows();j++){
          idI(j,j)=1/(1+lam1*d(j,0));
        }
        MatrixXd H=H1.adjoint()*idI*H1;
        VectorXd trHV(VectorXd(H.rows()).setZero());
        for(int k=0; k<eigenY.size(); k++){
          trHV(k)=H(k,k);
        }
        double C=(1-trHV.sum()/nrowx);
        C=nrowx*C*C;
        double GCV_tmp=(eigenY.adjoint()*(I-H)*(I-H)*eigenY/C).sum();
        GCV_mat(j2,i)=GCV_tmp;
        if(GCV<0){
          GCV=GCV_tmp;
          l_1=lam1;
          H_gcv=H;
          H1_idI_H1=dm*idI*dm.adjoint();
        }
        if(GCV>0){
          if(GCV_tmp<GCV){
            GCV=GCV_tmp;
            l_1=lam1;
            H_gcv=H;
            H1_idI_H1=dm*idI*dm.adjoint();
          }
        }
      }
    }
    
    Vx=Vx+l_1*K1;
    
    if(GCV2<0){
      GCV2=GCV;
      
      H_gcv2=H_gcv;
      H1_idI_H1cv=H1_idI_H1;
      
      VxCV=Vx;
    }
    if(GCV2>0){
      if(GCV<GCV2){
        GCV2=GCV;
        
        H_gcv2=H_gcv;
        H1_idI_H1cv=H1_idI_H1;
        
        VxCV=Vx;
      }
    }
  }  
  
  VectorXd betahat(H1_idI_H1cv*czy);
  
  
  return(betahat);
}

