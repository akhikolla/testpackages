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
#include <R.h>
#include <Rmath.h>
#include "mystruc.h"

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXf;
using Eigen::VectorXf;
using Eigen::VectorXi;
using Eigen::RowVectorXf;
std::mt19937 mt;

float GetUniform(){
    static std::uniform_real_distribution<float> Dist(0,1);
    return Dist(mt);
}
void GetSmallest(
			const VectorXf& r,
			const int& h,
			const MatrixXf& x,
			MatrixXf& xSub,
			VectorXi& RIndex
		){
	const int n=x.rows();
	VectorXi SIndx2(n);
	SIndx2.setLinSpaced(n,0,n-1);
	std::nth_element(SIndx2.data(),SIndx2.data()+h,SIndx2.data()+SIndx2.size(),IdLess(r.data()));
	for (int i=0;i<h;i++) 	xSub.row(i)=x.row(SIndx2(i));
	RIndex.head(h)=SIndx2.head(h);	
}
VectorXi SampleR(
			const int& m,
			const int& p
		){
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
VectorXf FindLine(
			const MatrixXf& xSub,
			const int& h
		){
	const int p=xSub.cols();
	VectorXi QIndexp=SampleR(h,p);
	MatrixXf A(p,p);
	for(int i=0;i<p;i++)	A.row(i)=xSub.row(QIndexp(i));
	VectorXf bt=VectorXf::Ones(p);
	return(A.lu().solve(bt));
}
void OneProj(
		const MatrixXf& x,
		const MatrixXf& xSub,
		const int& h,
		const VectorXi& RIndex,
		const int& h_m,
		int& EF,
 		int& EF2,
		VectorXf& Dep
	){
	const int p=x.cols(),n=x.rows();
	VectorXf beta(p);
	beta=FindLine(xSub,h);
	float betan=beta.norm();
	if((isnan(betan)) | (betan>std::numeric_limits<float>::max())){
		EF2++;
		return(void());
	}
	VectorXf praj=((x*beta).array()-1.0f).array().abs();
	float prem=0.0f,tol=1e-8;
	for(int i=0;i<h;i++)	prem+=praj(RIndex(i));
	prem=prem/(float)h;
	if(prem<tol){	
		const int n=praj.size();
		VectorXf d_resd=VectorXf::Zero(n);
		d_resd=(praj.array()<tol).select(1.0f,d_resd);
		if((d_resd.sum())>=h_m){
			Dep=praj;
			EF=1;
		} else {
			float maxin=praj.maxCoeff();
			d_resd=(praj.array()<tol).select(maxin,praj);
			prem=d_resd.minCoeff();
			praj/=prem;
			Dep+=praj;
		}
	} else {
		praj/=prem;
		Dep+=praj;
	}
}
float SubsetRankFun(
			const MatrixXf& x,
			const MatrixXf& xSub,
			const int& h,
			const VectorXi& RIndex,
			int& EF,
			int& EF2
		){
	const int p=x.cols(),n=x.rows();
	VectorXf beta(p);
	beta=FindLine(xSub,h);
	float betan=beta.norm();
	if((isnan(betan)) | (betan>std::numeric_limits<float>::max())){
		EF2++;
		return(std::numeric_limits<float>::max());
	}
	VectorXf praj=((x*beta).array()-1.0f).array().abs();
	VectorXf prej(h);
	for(int i=0;i<h;i++)	prej(i)=praj(RIndex(i));
	std::nth_element(praj.data(),praj.data()+h,praj.data()+praj.size());	
	float prem=praj.head(h).mean();
	float fin=(prem>1e-8)?prej.head(h).mean()/prem:1.0f;
	return fin;
}
float Main(
		MatrixXf& x,
		const int& k0,
		const int& J,
		const int& k1,
		const int& K,
		VectorXf& dP,
		const int& h_m,
		const VectorXi& hl,
		VectorXi& hm,
		VectorXi& RIndex,
		int& EF,
		int& EF2
	){
	const int p=x.cols(),n=x.rows(),h=K+1;
	RIndex.head(h)=SampleR(n,h);
	hm=RIndex.head(h);
	MatrixXf xSub(h_m,p);
	for(int i=0;i<h;i++) xSub.row(i)=x.row(RIndex(i));	
	RowVectorXf xSub_mean(p);
	xSub_mean=xSub.topRows(h).colwise().mean();	
	x.rowwise()-=xSub_mean;					
	xSub.rowwise()-=xSub_mean;
	JacobiSVD<MatrixXf> svd(xSub.topRows(h),ComputeThinV);
	MatrixXf w(n,K);
	w=(x*svd.matrixV().topLeftCorner(p,K));
	MatrixXf wSub(h_m,K);
	for(int i=0;i<h;i++) wSub.row(i)=w.row(RIndex(i));		
	for(int j=0;j<J;j++){			//growing step
		dP=VectorXf::Zero(n);
		EF2==0;
		for(int i=0;i<k0;i++) OneProj(w,wSub,hl(j),RIndex,h_m,EF,EF2,dP);	//const w!
		if(EF2==k0)	break;
		if(EF==1){
			GetSmallest(dP,h_m,w,wSub,RIndex);				//const w!
			break;
		}
		if(EF==0)	GetSmallest(dP,hl(j+1),w,wSub,RIndex);				//const w!
	}
	VectorXf fin=VectorXf::Ones(k1);
	if((EF==0) & (EF2<k0)){
		EF2=0;
		for(int i=0;i<k1;i++) fin(i)=SubsetRankFun(w,wSub,hl(J),RIndex,EF,EF2);	//const w!
	}
	return(fin.array().log().mean());
}
extern "C"{
	void FastHCS(
			int* rn,	//1
			int* p,		//2
			int* k0,	//3
			float* xi,	//4
			int* k1,	//5
			int* k2,	//6
			int* nsmp,	//7
			int* J,		//8
			float* objfunC,	//9
			int* seed,	//10
			int* rraw,	//13
			int* rh,	//14
			float* pco,	//15
			int* hf,	//16
			int* rra2,	//17	
			int* EF1out,	//18
			int* EF2out	//19
		){
		const int n=*rn,ik0=*k0,iJ=*J,ik1=*k1,K=*k2,ih_m=*rh;
		float objfunA,objfunB=std::numeric_limits<float>::max();
		mt.seed(*seed);
		MatrixXf x=Map<MatrixXf>(xi,n,*p);	
		VectorXf DpA=VectorXf::Zero(n);
		VectorXf DpB=VectorXf::Zero(n);
		VectorXi hl(*J+1);
		VectorXi hk(K+1);
		VectorXi hm(K+1);
		VectorXi RIndex0(ih_m);
		VectorXi RIndex1(ih_m);
		hl.setLinSpaced(*J+1,K+1,ih_m);
		hl(*J)=ih_m;
		int EF=0,EF2=0;
		for(int i=0;i<*nsmp;i++){			//for i=0 to i<#p-subsets.
			objfunA=Main(x,ik0,iJ,ik1,K,DpA,ih_m,hl,hk,RIndex0,EF,EF2);
			if(objfunA<objfunB){
				objfunB=objfunA;
				RIndex1=RIndex0;
				DpB=DpA;
				hm=hk;
			}
			if((EF) | (EF2>min(ik0,ik1)))	break;
		}
		*objfunC=objfunB;
		*EF1out=EF;
		*EF2out=EF2;
		Map<VectorXi>(rra2,K+1)=hm.array()+1;
		Map<VectorXf>(pco,n)=DpB;
		Map<VectorXi>(rraw,ih_m)=RIndex1.array()+1;
	}
}
