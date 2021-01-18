#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <random>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Core>
#include <Eigen/Cholesky>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXf;
using Eigen::VectorXf;
using Eigen::VectorXi;
std::mt19937 mt;

struct IdLess {					//internal function.
    template <typename T>
    IdLess(T iter) : values(&*iter) {}
    bool operator()(int left,int right){
        return values[left]<values[right];
    }
    float const* values;
};
float GetUniform(){
    static std::uniform_real_distribution<float> Dist(0,1);
    return Dist(mt);
}
void GetSmallest(const VectorXf& r,const int& h,const MatrixXf& x,VectorXf& y,MatrixXf& xSub,VectorXf& ySub,VectorXi& RIndex){
	const int n=x.rows();
	VectorXi SIndx2(n);
	SIndx2.setLinSpaced(n,0,n-1);
	std::nth_element(SIndx2.data(),SIndx2.data()+h,SIndx2.data()+SIndx2.size(),IdLess(r.data()));
	for (int i=0;i<h;i++){
	 	xSub.row(i)=x.row(SIndx2(i));
		ySub(i)=y(SIndx2(i));
	}
	RIndex.head(h)=SIndx2.head(h);	
}
VectorXi SampleR(const int m,const int p){
	int i,j,nn=m;
	VectorXi ind(nn);
	VectorXi y(p);
	ind.setLinSpaced(nn,0,nn-1);
    	for(i=0;i<p;i++){
		j=GetUniform()*nn;
		y(i)=ind(j);
		ind(j)=ind(--nn);
    	}
	return y;		
}
VectorXf FindLine(const MatrixXf& xSub,const VectorXf& ySub,const int h){
	const int p=xSub.cols();
	VectorXi  QIndexp(p);
	VectorXf  bt=VectorXf::Ones(p);
	QIndexp=SampleR(h,p);
	MatrixXf  A(p,p);
	for(int i=0;i<p;i++){
		A.row(i)=xSub.row(QIndexp(i));
		bt(i)=ySub(QIndexp(i));
	}
	return(A.lu().solve(bt));
}
VectorXf OneProj(const MatrixXf& x,const VectorXf& y,const MatrixXf& xSub,VectorXf& ySub,const int h,const VectorXi& RIndex,const int h_m){
	const int n=x.rows();
	VectorXf praj(n);
	praj=((x*FindLine(xSub,ySub,h)).array()-y.array()).array().abs2();
	float prem=0.0f,tol=1e-8;
	for(int i=0;i<h;i++)	prem+=praj(RIndex(i));
	prem=prem/(float)h;
	if(prem<tol){	
		const int n=praj.size();
		VectorXf d_resd=VectorXf::Zero(n);
		d_resd=(praj.array()<tol).select(1.0,d_resd);
		if((d_resd.sum())>=h_m){
			prem=1.0;
		} else {
			float maxin=praj.maxCoeff();
			d_resd=(praj.array()<tol).select(maxin,praj);
			prem=d_resd.minCoeff();
		}
	}
	return praj/=prem;
}
float SubsetRankFun(const MatrixXf& x,const VectorXf& y,const MatrixXf& xSub,const VectorXf& ySub,const int h,const VectorXi& RIndex){
	const int n=x.rows();
	VectorXf praj(n);
	VectorXf prej(h);
	praj=((x*FindLine(xSub,ySub,h)).array()-y.array()).array().abs2();
	for(int i=0;i<h;i++)	prej(i)=praj(RIndex(i));
	nth_element(praj.data(),praj.data()+h,praj.data()+praj.size());	
	float prem=praj.head(h).mean(),fin=(prem>1e-7)?(prej.head(h).mean()/prem):(1.0);
	return fin;
}
float Main(MatrixXf& x,VectorXf& y,const int k0,const int J,const int k1,VectorXf& dP,const int h_m,VectorXi& samset,const VectorXi& hl,VectorXi& RIndex){
	int p=x.cols(),n=x.rows(),h=p+1,ni=samset.size();
	MatrixXf xSub(h_m,p);
	VectorXf ySub(h_m);
	VectorXf fin(k1);
	RIndex.head(h)=SampleR(ni,h);			//draws random p-subset
	for(int i=0;i<h;i++){
		xSub.row(i)=x.row(samset(RIndex(i)));			
		ySub(i)=y(samset(RIndex(i)));
	}
	for(int j=0;j<J;j++){					//growing step
		dP=VectorXf::Zero(n);
		for(int i=0;i<k0;i++) dP+=OneProj(x,y,xSub,ySub,hl(j),RIndex,h_m);
		h=hl(j+1);
		GetSmallest(dP,hl(j+1),x,y,xSub,ySub,RIndex);
	}
	for(int i=0;i<k1;i++) fin(i)=SubsetRankFun(x,y,xSub,ySub,hl(J),RIndex);
	return fin.array().log().mean(); 
}
void CStep(VectorXi& dI,const MatrixXf& x,const VectorXf& y,const int h,const int h0){
	const int n=x.rows(),p=x.cols();
	float w1=0,w0=0;
	int w2=1,i,j=0;
	MatrixXf xSub(h,p);
	VectorXf ySub(h);
	VectorXf m_cf(p);
	MatrixXf b=MatrixXf::Identity(p,p);
	MatrixXf Sig(p,p);

	for(i=0;i<h0;i++) 	xSub.row(i)=x.row(dI(i));
	for(i=0;i<h0;i++) 	ySub(i)=y(dI(i));
	ColPivHouseholderQR<MatrixXf> QR(xSub.topRows(h0));
	m_cf=QR.solve(ySub.head(h0));
	VectorXf dP(n);
	dP=((x*m_cf).array()-y.array()).abs2();
	w1=std::numeric_limits<float>::max();
	while(w2){	
		dI.setLinSpaced(n,0,n-1);
		std::nth_element(dI.data(),dI.data()+h,dI.data()+dI.size(),IdLess(dP.data()));
		for(i=0;i<h;i++) 	xSub.row(i)=x.row(dI(i));
		for(i=0;i<h;i++) 	ySub(i)=y(dI(i));
		HouseholderQR<MatrixXf> QR(xSub);
		m_cf=QR.solve(ySub);
		dP=((x*m_cf).array()-y.array()).abs2();
		w0=w1;
		j++;
		w1=log(((xSub*m_cf).array()-ySub.array()).abs2().sum()/(float)(h-1));
		((w0-w1<1e-3) | (j>10))?(w2=0):(w2=1);
	}
} 
extern "C"{
	void fastrcs(
			int* n,		//1
			int* p,		//2
			int* k0,	//3
			float* xi,	//4
			float* yi,	//5
			int* k1,	//6
			float* DpC,	//7
			int* nsamp,	//8
			int* J,		//9
			float* objfunC,	//10
			int* seed,	//11
			int* ck,	//12
			int* ni,	//13
			int* n1,	//14
			int* n2,	//15
			int* h,		//16
			int* hf		//17
		){
		const int ik0=*k0,iJ=*J,ik1=*k1,h_m=*h,h_f=*hf;
		float objfunA,objfunB=*objfunC;
		mt.seed(*seed);
		MatrixXf x=Map<MatrixXf>(xi,*n,*p);	
		VectorXi icK=Map<VectorXi>(ck,*ni);
		VectorXf y=Map<VectorXf>(yi,*n);	
		VectorXf DpA=VectorXf::Zero(*n);
		VectorXf DpB=VectorXf::Zero(*n);
		VectorXi RIndexi(*n);
		VectorXi RIndexf(*n);
		VectorXi hl(*J+1);
		hl.setLinSpaced(*J+1,*p+1,h_m);
		hl(*J)=h_m;
		for(int i=0;i<*nsamp;i++){			//for i=0 to i<#p-subsets.
			objfunA=Main(x,y,ik0,iJ,ik1,DpA,h_m,icK,hl,RIndexi);
			if(objfunA<objfunB){
				objfunB=objfunA;
				DpB=DpA;
				RIndexf.head(h_m)=RIndexi.head(h_m);
			}
		}
		Map<VectorXi>(n1,h_m)=RIndexf.head(h_m).array()+1;
		Map<VectorXf>(DpC,*n)=DpB.array();		
		CStep(RIndexf,x,y,h_f,h_m);
		Map<VectorXi>(n2,h_f)=RIndexf.head(h_f).array()+1;
		*objfunC=objfunB;
	}
}
