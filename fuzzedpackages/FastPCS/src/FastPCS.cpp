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

#include <Eigen/QR>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Cholesky>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXf;
using Eigen::VectorXf;
using Eigen::VectorXi;
using Eigen::RowVectorXf;
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
void GetSmallest(const VectorXf& r,const int& h,const MatrixXf& x,MatrixXf& xSub,VectorXi& RIndex){
	const int n=x.rows();
	VectorXi SIndx2(n);
	SIndx2.setLinSpaced(n,0,n-1);
	std::nth_element(SIndx2.data(),SIndx2.data()+h,SIndx2.data()+SIndx2.size(),IdLess(r.data()));
	for (int i=0;i<h;i++) 	xSub.row(i)=x.row(SIndx2(i));
	RIndex.head(h)=SIndx2.head(h);	
}
VectorXi SampleR(const int m,const int p){
	int i,j,nn=m;
	VectorXi ind(nn);
	ind.setLinSpaced(nn,0,nn-1);
	VectorXi y(p);
    	for(i=0;i<p;i++){
		j=GetUniform()*nn;
		y(i)=ind(j);
		ind(j)=ind(--nn);
    	}
	return y;		
}
VectorXf FindLine(const MatrixXf& xSub,const int h){
	const int p=xSub.cols();
	VectorXi  QIndexp(p);
	QIndexp=SampleR(h,p);
	MatrixXf A(p,p);
	for(int i=0;i<p;i++)	A.row(i)=xSub.row(QIndexp(i));
	VectorXf bt=VectorXf::Ones(p);
	return(A.lu().solve(bt));
}
VectorXf OneProj(const MatrixXf& x,const MatrixXf& xSub,const int h,const VectorXi& RIndex,const int h_m){
	const int p=x.cols(),n=x.rows();
	VectorXf beta(p);
	beta=FindLine(xSub,h);
	VectorXf praj(n);
	praj=((x*beta).array()-1.0f).array().abs2();
	praj/=beta.squaredNorm();
	float prem=0.0f,tol=1e-8;
	for(int i=0;i<h;i++)	prem+=praj(RIndex(i));
	prem=prem/(float)h;
	if(prem<tol){	
		const int n=praj.size();
		VectorXf d_resd=VectorXf::Zero(n);
		d_resd=(praj.array()<tol).select(1.0f,d_resd);
		if((d_resd.sum())>=h_m){
			prem=1.0f;
		} else {
			float maxin=praj.maxCoeff();
			d_resd=(praj.array()<tol).select(maxin,praj);
			prem=d_resd.minCoeff();
		}
	}
	return praj/=prem;
}
float SubsetRankFun(const MatrixXf& x,const MatrixXf& xSub,const int h,const VectorXi& RIndex){
	const int p=x.cols(),n=x.rows();
	VectorXf beta(p);
	beta=FindLine(xSub,h);
	VectorXf praj(n);
	praj=((x*beta).array()-1.0f).array().abs2();
	praj/=beta.squaredNorm();
	VectorXf prej(h);
	for(int i=0;i<h;i++)	prej(i)=praj(RIndex(i));
	std::nth_element(praj.data(),praj.data()+h,praj.data()+praj.size());	
	float prem=praj.head(h).mean(),fin=(prem>1e-7)?(prej.head(h).mean()/prem):(1.0f);
	return fin;
}
float Main(const MatrixXf& x,const int k0,const int J,const int k1,VectorXf& dP,const int h_m,VectorXi& samset,const VectorXi& hl,VectorXi& RIndex){
	int p=x.cols(),n=x.rows(),h=p+1,ni=samset.size();
	RIndex.head(h)=SampleR(ni,h);
	MatrixXf xSub(h_m,p);
	for(int i=0;i<h;i++) xSub.row(i)=x.row(samset(RIndex(i)));			
	for(int j=0;j<J;j++){					//growing step
		dP=VectorXf::Zero(n);
		for(int i=0;i<k0;i++) dP+=OneProj(x,xSub,hl(j),RIndex,h_m);
		GetSmallest(dP,hl(j+1),x,xSub,RIndex);
	}
	VectorXf fin(k1);
	for(int i=0;i<k1;i++) fin(i)=SubsetRankFun(x,xSub,hl(J),RIndex);
	return fin.array().log().mean();
}
void CStep(VectorXi& dIn,MatrixXf& x,const int h,const int h0){
	const int n=x.rows(),p=x.cols();
	float w1,w0;
	int w2=1,i;
	MatrixXf xSub(h,p);
	for(i=0;i<h0;i++) 	xSub.row(i)=x.row(dIn(i));
	RowVectorXf xSub_mean(p);
	xSub_mean=xSub.topRows(h0).colwise().mean();	
	xSub.topRows(h0).rowwise()-=xSub_mean;
	x.rowwise()-=xSub_mean;
	MatrixXf Sig(p,p);
	Sig=xSub.topRows(h0).adjoint()*xSub.topRows(h0);
	Sig.array()/=(float)(h0-1);
	LDLT<MatrixXf> chol=Sig.ldlt();
	MatrixXf b=MatrixXf::Identity(p,p);
	chol.solveInPlace(b);
	w1=chol.vectorD().array().minCoeff();
	VectorXf dP(n);
	if(w1>1e-6){
		w1=std::numeric_limits<float>::max();
		dP=((x*b).cwiseProduct(x)).rowwise().sum();
	} else {
		w2=0;
	}
	while(w2){	
		dIn.setLinSpaced(n,0,n-1);
		std::nth_element(dIn.data(),dIn.data()+h,dIn.data()+dIn.size(),IdLess(dP.data()));
		for(i=0;i<h;i++) 	xSub.row(i)=x.row(dIn(i));
		xSub_mean=xSub.colwise().mean();	
		xSub.rowwise()-=xSub_mean;
		x.rowwise()-=xSub_mean;
		Sig=xSub.adjoint()*xSub;
		Sig.array()/=(float)(h-1);
		LDLT<MatrixXf> chol=Sig.ldlt();
		b=MatrixXf::Identity(p,p);
		chol.solveInPlace(b);
		if(chol.vectorD().array().minCoeff()>1e-6){
			w0=w1;
			w1=chol.vectorD().array().log().sum()*2.0f;
			dP=((x*b).cwiseProduct(x)).rowwise().sum();
			(w0-w1<1e-3)?(w2=0):(w2=1);
		} else {
			w2=0;
		}
	}
} 
extern "C"{
	void fastpcs(
			int* n,		//1
			int* p,		//2
			int* k0,	//3
			float* xi,	//4
			int* k1,	//5
			float* DpC,	//6
			int* nsamp,	//7
			int* J,		//8
			float* objfunC,	//9
			int* seed,	//10
			int* ck,	//11
			int* ni,	//12
			int* n1,	//13
			int* n2,	//14
			int* hm,	//15
			int* hf		//16
		){
		const int ik0=*k0,iJ=*J,ik1=*k1,ih_m=*hm,ih_f=*hf;
		float objfunA,objfunB=*objfunC;
		mt.seed(*seed);
		MatrixXf x=Map<MatrixXf>(xi,*n,*p);	
		VectorXi icK=Map<VectorXi>(ck,*ni);
		VectorXf DpA=VectorXf::Zero(*n);
		VectorXf DpB=VectorXf::Zero(*n);
		VectorXi RIndexi(*n);
		VectorXi RIndexf(*n);
		VectorXi hl(*J+1);
		hl.setLinSpaced(*J+1,*p+1,ih_m);
		hl(*J)=ih_m;

		for(int i=0;i<*nsamp;i++){			//for i=0 to i<#p-subsets.
			objfunA=Main(x,ik0,iJ,ik1,DpA,ih_m,icK,hl,RIndexi);
			if(objfunA<objfunB){
				objfunB=objfunA;
				DpB=DpA;
				RIndexf.head(ih_m)=RIndexi.head(ih_m);
			}
		}
		Map<VectorXi>(n1,*hm)=RIndexf.head(ih_m).array()+1;		
 		Map<VectorXf>(DpC,*n)=DpB.array();
		CStep(RIndexf,x,ih_f,ih_m);
		Map<VectorXi>(n2,*hf)=RIndexf.head(ih_f).array()+1;
		*objfunC=objfunB;
	}
}
