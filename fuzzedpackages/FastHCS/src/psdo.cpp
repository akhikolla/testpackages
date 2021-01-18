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


using namespace Eigen;
using Eigen::VectorXf;
std::mt19937 mt2;

extern VectorXi SampleR(
			const int& m,
			const int& p
			);
void unimcd(
		const VectorXf& m_resd,
		const int& h,
		VectorXf& m_medn,
		const int& rw=1
	){
	const int n1=m_resd.size(),len=n1-h+1;
	float initmean=0.0,initcov=0.0,sumw=0.0;
	int minone;

	VectorXf y=m_resd;
	VectorXf ay(len);
	VectorXf ay2(len);
	VectorXf sq(len);
	VectorXf y2(n1);

	std::sort(y.data(),y.data()+y.size());
	ay(0)=y.head(h).sum();	
	for(int samp=1;samp<len;samp++) ay(samp)=ay(samp-1)-y(samp-1)+y(samp+h-1);
	ay2=ay.array().square()/(float)h;
	y2=y.array().square();
	sq(0)=y2.head(h).sum()-ay2(0);
	for(int samp=1;samp<len;samp++) sq(samp)=sq(samp-1)-y2(samp-1)+y2(samp+h-1)-ay2(samp)+ay2(samp-1);
	initcov=sq.minCoeff(&minone);
	initcov/=(float)(h-1);
	initmean=ay(minone)/(float)h;
	if(rw){
		y2=(m_resd.array()-initmean).array().abs2()/initcov;
		std::nth_element(y2.data(),y2.data()+h-1,y2.data()+y2.size());	
		initcov*=y2(h-1)/qchisq(0.975,h/(float)n1,1,0);
		y=(m_resd.array()-initmean).array().abs2()/initcov;
		y2.setOnes();
		y2=(y.array()>5.023886f).select(0.0f,y2);
		sumw=y2.array().sum();
		initmean=(y2.array()*m_resd.array()).sum()/sumw;
		y=(m_resd.array()-initmean).array().abs2();
		initcov=(y2.array()*y.array()).sum()/(sumw-1.0f);
	}
	m_medn<<initmean,sqrt(initcov);
}
float quantiles(
			VectorXf& x,
			const float& quant
		){
	const int n=x.size();
	float lq,uq,fq;
	const float q1=n*(float)quant+0.5;
	const int index1=floor(q1);
	const int index2=ceil(q1);
	const float index3=(float)index2-q1;
	std::nth_element(x.data(),x.data()+index1-1,x.data()+x.size());
	lq=x(index1-1);
	if(index1==index2){
		fq=lq;
	} else {
		uq=x.segment(index1,x.size()-index1-1).minCoeff();
		fq=lq*index3+uq*(1.0-index3);
	}
	return(fq);
}
void mad(
		const VectorXf& m_resd,
		const int& h,
		VectorXf& m_medn
	){
	const int n=m_resd.size();
	float temp0=0.5,hf=h/(float)n;
	VectorXf i_resd=m_resd;
	float temp1=quantiles(i_resd,temp0);
	i_resd.array()-=temp1;
	i_resd=i_resd.array().abs();
	m_medn<<temp1,quantiles(i_resd,temp0);
}
void psdo_dir(
			const MatrixXf& x,
			int& m_r,
			VectorXf& m_coef
		){
	const int n=x.rows();
	VectorXi QIndexpin(2);
	QIndexpin=SampleR(n,2);
	m_coef=x.row(QIndexpin(0))-x.row(QIndexpin(1));
	float m_medn=m_coef.norm();
	if(m_medn>1e-8){
		m_coef.array()/=m_medn;	
		m_r=1;
	}
}
void main_psdo(
			const MatrixXf& x,
			int& ndir,
			int& sdr,
			VectorXf& outlyingness,
			const int& h
		){
	const int p=x.cols(),n=x.rows();
	int m_r=0,j=0;
	VectorXf m_medn(2);
	VectorXf m_coef=VectorXf::Ones(p);
	VectorXf m_resd=VectorXf::Zero(n);
	while(j<ndir){
		m_r=0;
		psdo_dir(x,m_r,m_coef);
		j++;
		if(m_r==1){
	    		m_resd=x*m_coef;
			mad(m_resd,h,m_medn);
			m_resd.array()-=m_medn(0);
			m_resd=m_resd.array().abs();
			if(m_medn(1)>1e-12){
				m_resd.array()/=m_medn(1);
				outlyingness=outlyingness.cwiseMax(m_resd);
			} else {
				j=ndir;
				VectorXf w_resd=VectorXf::Zero(n);
				outlyingness=(m_resd.array()>1e-12).select(1.0,w_resd);
			}
		} else {
			sdr++;
		}
	}
}
extern "C"{
	void r_psdo(
			int* n,		//1
			int* p,		//2
			int* ndir,	//3
			float* xi,	//4
			float* pout,	//5
			int* sdr,	//6
			int* ih,	//7
			int* seed,	//8
			int* poutind	//9
		){
		const int h=*ih;
		int indir=*ndir,isdr=0;
		mt2.seed(*seed);
		MatrixXf x=Map<MatrixXf>(xi,*n,*p);	
		VectorXf outin=VectorXf::Zero(*n);
		main_psdo(x,indir,isdr,outin,h);
 		Map<VectorXf>(pout,*n)=outin;
		VectorXi SIndx2(*n);
		SIndx2.setLinSpaced(*n,0,*n-1);
		std::nth_element(SIndx2.data(),SIndx2.data()+h,SIndx2.data()+SIndx2.size(),IdLess(outin.data()));
		Map<VectorXi>(poutind,h)=SIndx2.head(h).array()+1;
		*sdr=isdr;
	}
}
