#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <iomanip> 

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>

using namespace std;
using namespace Eigen;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::RowVectorXd;

struct IdLess {
    template <typename T>
    IdLess(T iter) : values(&*iter) {}
    bool operator()(int left,int right){
        return values[left]<values[right];
    }
    double const* values;
};
double whimed_i(
			VectorXd& a,
			VectorXi& w,
			int& n,
			VectorXd& a_cand,
			VectorXd& a_psrt,
			VectorXi& w_cand
		){
//Eigen-ized version of whimed_i. See citation("pcaPP")
	int n2,i,k_cand,nn=n;/* sum of weights: `int' do overflow when  n ~>= 1e5 */
	int64_t wleft,wmid,wright,w_tot,wrest=0;
	double trial;
	w_tot=w.head(n).sum();
	do{
		a_psrt.head(nn)=a.head(nn);
		n2=nn/2;
		std::nth_element(a_psrt.data(),a_psrt.data()+n2,a_psrt.data()+nn);
		trial=a_psrt(n2);
		wleft=0;wmid=0;wright=0;
		for (i=0;i<nn;++i){
		    if(a(i)<trial)				wleft+=w(i);
		    else if(a(i)>trial)				wright+=w(i);
		    else wmid+=w(i);
		}
		k_cand=0;
		if((2*(wrest+wleft))>w_tot){
			for (i=0;i<nn;++i){
				if(a(i)<trial){
				    	a_cand(k_cand)=a(i);
				    	w_cand(k_cand)=w(i);
					++k_cand;
				}
			}
		}
		else if((2*(wrest+wleft+wmid))<=w_tot){
			for (i=0;i<nn;++i){
				if (a(i)>trial){
					a_cand(k_cand)=a(i);
					w_cand(k_cand)=w(i);
					++k_cand;
				}
		    	}
			wrest+=wleft+wmid;
		} else {
		        return(trial);
		}
		nn=k_cand;
		a.head(nn)=a_cand.head(nn);
		w.head(nn)=w_cand.head(nn);
	} while(1); 
}
double qn(VectorXd& x){
//Eigen-ized version of qn. See citation("pcaPP")
		const int n=x.rows();
		VectorXd y=x;
		VectorXd work=VectorXd::Zero(n);
    		VectorXi p=VectorXi::Zero(n);
    		VectorXi q=VectorXi::Zero(n);
		VectorXi left=VectorXi::Zero(n);
    		VectorXi right=VectorXi::Zero(n);
    		VectorXi weight=VectorXi::Zero(n);
		VectorXd a_cand(n);
		VectorXd a_psrt(n);
		
		bool found;
		int h,i,j,jj,jh;
		double ri,trial,dn=1.0;
	    	int64_t k,knew,nl,nr,sump,sumq;

	    	h=n/2+1;
    		k=(int64_t)h*(h-1)/2;
    		for(i=0;i<n;++i){
			left(i)=n-i+1;
			right(i)=(i<=h)?n:n-(i-h);
		}
		std::sort(y.data(),y.data()+y.size());

		nl=(int64_t)n*(n+1)/2;
		nr=(int64_t)n*n;
		knew=k+nl;/*=k +(n+1 \over 2) */
		found=false;	
		while(!found && nr-nl>n){
			j=0;
			for(i=1;i<n;++i){
				if(left(i)<=right(i)){
					weight(j)=right(i)-left(i)+1;
					jh=left(i)+weight(j)/2;
					work(j)=(y(i)-y(n-jh));
					++j;
				}
			}
			trial=whimed_i(work,weight,j,a_cand,a_psrt,p);
			j=0;
			for(i=n-1;i>=0;--i){
				while(j<n&&((y(i)-y(n-j-1)))<trial) ++j;
				p(i)=j;
			}
			j=n+1;
			for(i=0;i<n;++i){
				while((y(i)-y(n-j+1))>trial)	--j;
				q(i)=j;
			}
			sump=p.sum();
			sumq=q.sum()-n;
			if(knew<=sump){
				right=p;
				nr=sump;
			} else if(knew>sumq){
				left=q;
				nl=sumq;
			} else { /* sump < knew <= sumq */
			    	found=true;
			}
    		} /* while */
		if(found) ri=trial;
		else{
			j=0;
			for(i=1;i<n;++i){
		    		for(jj=left(i);jj<=right(i);++jj){
					work(j)=y(i)-y(n-jj);
					j++;
		    		}/* j will be=sum_{i=2}^n(right[i] - left[i] + 1)_{+}  */
			}
			std::nth_element(work.data(),work.data()+knew-nl-1,work.data()+j);
			ri=work(knew-nl-1);
		}
    		ri*=2.21914446598508;//ri*=2.2219;
		if(n<=9){
			if(n==2)	dn=.399;
	    		else if(n==3)	dn=.994;
	    		else if(n==4)	dn=.512;
	    		else if(n==5)	dn=.844;
	    		else if(n==6)	dn=.611;
	    		else if(n==7)	dn=.857;
	    		else if(n==8)	dn=.669;
	    		else if(n==9)	dn=.872;
		} else {	
	    		if(n%2==1) dn=n/(n+1.4);
	    		else dn=n/(n+3.8);
		}
		ri*=dn;
		return(ri);
}
double Fmedian(VectorXd& x){
	int n=x.rows();
	int half=(n+1)/2-1;				
	double med;
	std::nth_element(x.data(),x.data()+half,x.data()+x.size());	
	if((n%2)==1){
		med=x(half);
	} else {
		double tmp0=x(half);
		double tmp1=x.segment(half+1,half-1).minCoeff(); 
		med=0.5*(tmp0+tmp1);
	}
	return med;
}
double pCalc(VectorXd& x,double (*pCalcMethod)(VectorXd&)){
        return(pCalcMethod(x));
}
double scaleTau2(VectorXd& x){
	int n=x.rows();
	const double C1=4.5;
  	const double Es2c=0.9247153921761315;
	double med=0.0,sigma0=0.0;
	VectorXd dwork1(n);
	VectorXd dwork2(n);

	dwork1=x;
	med=Fmedian(dwork1);
	dwork1=(x.array()-med).array().abs();
	dwork2=dwork1;
	sigma0=Fmedian(dwork2);

	dwork1.array()/=(sigma0*C1);
	dwork2=1-dwork1.array().square();
	dwork2=((dwork2.array().abs()+dwork2.array())*0.5).array().square();
	med=(dwork2.array()*x.array()).sum()/dwork2.sum();

	x=(x.array()-med).array()/sigma0;
	dwork1=x.array().square();
	dwork1=(dwork1.array()>9.0).select(9.0,dwork1);	
	return(sigma0*sqrt(dwork1.sum()/(n*Es2c)));
}
VectorXi Finitset(
			const MatrixXd& x,
			const MatrixXd& P,
			const int& calcM,
			const int& h0,
			int& Svol
		){
	const int n=x.rows(),p=x.cols();	
	double (*pFo[])(VectorXd&)={&qn,&scaleTau2}; 
	SelfAdjointEigenSolver<Eigen::MatrixXd> solved(P);
	MatrixXd Data2=x*solved.eigenvectors();
	RowVectorXd lamba(p);
	VectorXd lambb(n);
	for(int i=0;i<p;++i){
		lambb=Data2.col(i);
		lamba(i)=pCalc(lambb,pFo[calcM]);
	}
	Data2=x*(solved.eigenvectors()*lamba.asDiagonal().inverse()*solved.eigenvectors().transpose());
	RowVectorXd lambc(p);
	for(int i=0;i<p;++i){
		lambb=Data2.col(i);
		lambc(i)=Fmedian(lambb);
	}
	lambc*=(solved.eigenvectors()*lamba.asDiagonal()*solved.eigenvectors().transpose());
	Data2=x;
	Data2.rowwise()-=lambc;
	Data2*=solved.eigenvectors();
	Data2*=lamba.asDiagonal().inverse();
	lambb=Data2.cwiseAbs2().rowwise().sum();
	VectorXi SIndx1(n);
	SIndx1.setLinSpaced(n,0,n-1);
	std::nth_element(SIndx1.data(),SIndx1.data()+h0,SIndx1.data()+SIndx1.size(),IdLess(lambb.data()));

	MatrixXd xSub(h0,p);
	for(int i=0;i<h0;i++) 	xSub.row(i)=x.row(SIndx1(i));
	RowVectorXd xSub_mean(p);
	xSub_mean=xSub.colwise().mean();	
	xSub.rowwise()-=xSub_mean;
	MatrixXd Sig(p,p);
	Sig.setZero().selfadjointView<Lower>().rankUpdate(xSub.transpose());
	Sig.array()/=(double)(h0-1);
	SelfAdjointEigenSolver<Eigen::MatrixXd> solvedSig(Sig);
	if(solvedSig.eigenvalues().minCoeff()>1e-8){
		Data2=x;
		Data2.rowwise()-=xSub_mean;	
		Data2*=solvedSig.eigenvectors();
		Data2*=solvedSig.eigenvalues().cwiseSqrt().asDiagonal().inverse();
		lambb=Data2.cwiseAbs2().rowwise().sum();
		SIndx1.setLinSpaced(n,0,n-1);
		std::nth_element(SIndx1.data(),SIndx1.data()+h0,SIndx1.data()+SIndx1.size(),IdLess(lambb.data()));
	} else {
		Svol=1;
	}
	return(SIndx1.head(h0));
}
void cov_CStep(
			VectorXi& dIn,
			const MatrixXd& x,
			const int& h,
			const int& h0,
			int& w3,
			double& w1
		){
	const int n=x.rows(),p=x.cols();
	double w0;
	int w2=1,i;
	MatrixXd xSub(h,p);
	MatrixXd Data2(n,p);
	for(i=0;i<h0;i++) 	xSub.row(i)=x.row(dIn(i));
	RowVectorXd xSub_mean(p);
	xSub_mean=xSub.topRows(h0).colwise().mean();	
	xSub.topRows(h0).rowwise()-=xSub_mean;
	MatrixXd Sig(p,p);
	Sig.setZero().selfadjointView<Lower>().rankUpdate(xSub.topRows(h0).transpose());
	Sig.array()/=(double)(h0-1);
	SelfAdjointEigenSolver<Eigen::MatrixXd> solved(Sig);
	w1=solved.eigenvalues().minCoeff();
	VectorXd dP(n);	
	w3=0;
	w1=std::numeric_limits<double>::max();
	Data2=x;
	Data2.rowwise()-=xSub_mean;	
	Data2*=solved.eigenvectors();
	Data2*=solved.eigenvalues().cwiseSqrt().asDiagonal().inverse();
	dP=Data2.cwiseAbs2().rowwise().sum();
	while(w2){	
		dIn.setLinSpaced(n,0,n-1);
		std::nth_element(dIn.data(),dIn.data()+h,dIn.data()+dIn.size(),IdLess(dP.data()));
		for(i=0;i<h;i++) 	xSub.row(i)=x.row(dIn(i));
		xSub_mean=xSub.colwise().mean();	
		xSub.rowwise()-=xSub_mean;
		Sig.setZero().selfadjointView<Lower>().rankUpdate(xSub.transpose());
		Sig.array()/=(double)(h-1);
		SelfAdjointEigenSolver<Eigen::MatrixXd> solved(Sig);
		w3++;
		if(solved.eigenvalues().minCoeff()>1e-8){
			w0=w1;
			w1=solved.eigenvalues().array().log().sum();
			Data2=x;
			Data2.rowwise()-=xSub_mean;	
			Data2*=solved.eigenvectors();
			Data2*=solved.eigenvalues().cwiseSqrt().asDiagonal().inverse();
			dP=Data2.cwiseAbs2().rowwise().sum();
			(w0-w1<1e-6)?(w2=0):(w2=1);
		} else {
			w2=0;
			w1=solved.eigenvalues().array().log().sum();
		}
	}
}  
VectorXi CStepAll(
				const MatrixXd& xi,
				const int& hf,
				int& wo,
				int& nit,
				MatrixXi& hsup,
				double& bestobjfun
			){
	const int n=xi.rows(),p=xi.cols(),h0=hsup.rows(),S=hsup.cols();	
	int w3;
	double objfun;
	bestobjfun=std::numeric_limits<double>::max();
	VectorXi dIn(n);
	VectorXi dOut(hf);
	for(int i=0;i<S;i++){
		dIn.head(h0)=hsup.col(i);
		cov_CStep(dIn,xi,hf,h0,w3,objfun);
		if(objfun<bestobjfun){
			wo=i;
			nit=w3;
			bestobjfun=objfun;
			dOut=dIn.head(hf);
		}
	}
	return(dOut);
}  
extern "C"{
	void R_FastR(
			int* n,		//1
			int* p,		//2
			double* X,	//3
			double* Q,	//4
			int* cMet,	//5
			int* h,		//6
			int* hL,	//7
			int* wM,	//8
			int* BestD,	//9
			int* BestI,	//10
			int* pQ,	//11
			int* svol,	//12
			double* obj,	//13
			int* wMef,	//14
			int* Hmin,	//15
			int* exft	//16				//int* mW		//14
		){
		const int CalcMet=*cMet,q=(*p)*(*p),hlt=*hL,hmin=*Hmin;
		VectorXd objy(hlt);
		VectorXi hi=Map<VectorXi>(h,hlt);
		MatrixXd xi=Map<MatrixXd>(X,*n,*p);	
		MatrixXd qi=Map<MatrixXd>(Q,q,*pQ);
		MatrixXi hsub=MatrixXi::Zero(*n,hlt);
		MatrixXi hsup(hmin,*pQ);
		VectorXi bit(hlt);
		VectorXi bdt(hlt);
		VectorXi Svol=VectorXi::Zero(*pQ);
		for(int i=0;i<*pQ;++i){
			Map<const MatrixXd>q2(qi.col(i).data(),*p,*p);
			hsup.col(i)=Finitset(xi,q2,CalcMet,hmin,Svol(i));	
		}
		//Map<MatrixXi>(mW,hi(0),*pQ)=hsup.array()+1;	//all the indexes.
		Map<VectorXi>(svol,*pQ)=Svol;			
		int wef,ismax=Svol.maxCoeff(&wef);
		if(ismax<1){
			for(int i=0;i<hlt;i++)	hsub.col(i).head(hi(i))=CStepAll(xi,hi(i),bdt(i),bit(i),hsup,objy(i));
			Map<VectorXi>(BestD,hlt)=bdt.array()+1;		//which one was chosen
			Map<VectorXi>(BestI,hlt)=bit.array()+1;		//nit
			Map<MatrixXi>(wM,*n,hlt)=hsub.array()+1;	//all the indexes.
			Map<VectorXd>(obj,hlt)=objy;
		} else {
			Map<VectorXi>(wMef,*Hmin)=hsup.col(wef).array()+1;
			*exft=1;
			bdt(0)=wef;
			Map<VectorXi>(BestD,hlt)=bdt.array()+1;		//which one was chosen
		}
	}
}
