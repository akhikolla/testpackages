#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
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

struct IdLess {					//internal function.
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
			int n,
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
		    if(a(i)<trial)					wleft+=w(i);
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
double qn(
		Ref<VectorXd> x,
		int retMed
	){//Eigen-ized version of qn. See citation("pcaPP")
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
void cov_CStep(
		VectorXi& dIn,
		MatrixXd& x,
		const int h,
		const int h0,
		int w3
	){
	const int n=x.rows(),p=x.cols();
	double w1,w0;
	int w2=1,i;
	MatrixXd xSub(h,p);
	for(i=0;i<h0;i++) 	xSub.row(i)=x.row(dIn(i));
	RowVectorXd xSub_mean(p);
	xSub_mean=xSub.topRows(h0).colwise().mean();	
	xSub.topRows(h0).rowwise()-=xSub_mean;
	x.rowwise()-=xSub_mean;
	MatrixXd Sig(p,p);
	//Sig=xSub.topRows(h0).adjoint()*xSub.topRows(h0);
	Sig.setZero().selfadjointView<Lower>().rankUpdate(xSub.topRows(h0).transpose());
	Sig.array()/=(double)(h0-1);
	LDLT<MatrixXd> chol=Sig.ldlt();
	MatrixXd b=MatrixXd::Identity(p,p);
	chol.solveInPlace(b);
	w1=chol.vectorD().array().minCoeff();
	VectorXd dP(n);	
	if(w1>1e-6){
		w1=std::numeric_limits<double>::max();
		dP=((x*b).cwiseProduct(x)).rowwise().sum();
	} else {
		w2=0;
		w3=0;
		w1=chol.vectorD().array().log().sum()*2.00;
	}
	while(w2){	
		dIn.setLinSpaced(n,0,n-1);
		std::nth_element(dIn.data(),dIn.data()+h,dIn.data()+dIn.size(),IdLess(dP.data()));
		for(i=0;i<h;i++) 	xSub.row(i)=x.row(dIn(i));
		xSub_mean=xSub.colwise().mean();	
		xSub.rowwise()-=xSub_mean;
		x.rowwise()-=xSub_mean;
		//Sig=xSub.adjoint()*xSub;
		Sig.setZero().selfadjointView<Lower>().rankUpdate(xSub.transpose());
		Sig.array()/=(double)(h-1);
		chol=Sig.ldlt();
		b=MatrixXd::Identity(p,p);
		chol.solveInPlace(b);
		if(chol.vectorD().array().minCoeff()>1e-6){
			w0=w1;
			w1=chol.vectorD().array().log().sum()*2.00;
			dP=((x*b).cwiseProduct(x)).rowwise().sum();
			(w0-w1<1e-3)?(w2=0):(w2=1);
		} else {
			w2=0;
			w3=0;
		}
	}
} 
int unimcd(
		VectorXd& y,
		const int& n,
		const int& h,
		const int& len,
		double& initmean,
		double& initcov
	){
	//does NOT accept ties in the y's!
	std::sort(y.data(),y.data()+y.size());
	VectorXd ay(len);
	ay(0)=y.head(h).sum();	
	for(int samp=1;samp<len;samp++) ay(samp)=ay(samp-1)-y(samp-1)+y(samp+h-1);
	VectorXd ay2(len);
	ay2=ay.array().square()/(double)h;
	VectorXd y2(n);
	y2=y.array().square();
	VectorXd sq(len);
	sq(0)=y2.head(h).sum()-ay2(0);
	for(int samp=1;samp<len;samp++) sq(samp)=sq(samp-1)-y2(samp-1)+y2(samp+h-1)-ay2(samp)+ay2(samp-1);
	int minone;
	initcov=sq.minCoeff(&minone);
	initcov/=(double)(h-1);
	initmean=ay(minone)/(double)h;
	return(minone);
}
double Fmedian(
		Ref<VectorXd> x,
		int retMed
	){
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
double pCalc(
		Ref<VectorXd> x,
		int retMed,
		double (*pCalcMethod)(Ref<VectorXd>,int)
	){
		return(pCalcMethod(x,retMed));
}
double qCalc(
		Ref<VectorXd> x,
		int retMed,
		double (*qCalcMethod)(Ref<VectorXd>,int)
	){
		return(qCalcMethod(x,retMed));
}
double scaleTau2(
		Ref<VectorXd> x,
		int retMed
	){
	int n=x.rows();
	const double C1=4.5,tol=1e-8;
  	const double Es2c=0.9247153921761315;
	double med=0.0,sigma0=0.0;
	VectorXd dwork1(n);
	dwork1=x;
	med=Fmedian(dwork1,0);
	dwork1=(x.array()-med).array().abs();
	VectorXd dwork2=dwork1;
	sigma0=Fmedian(dwork2,0);
	if(sigma0<tol)	return(0.0);
	double sigma1=1/(C1*sigma0);
	dwork1.array()*=(sigma1);
	dwork2=1-dwork1.array().square();
	dwork2=((dwork2.array().abs()+dwork2.array())*0.5).array().square();
	med=(dwork2.array()*x.array()).sum()/dwork2.sum();
	if(retMed)	return(med);
	dwork1=(x.array()-med);
	sigma1=1/sigma0;
	dwork1.array()*=(sigma1);
	dwork1=dwork1.array().abs2();
	dwork1=(dwork1.array()>9.0).select(9.0,dwork1);
	sigma1=dwork1.sum();
	if(sigma1<tol)	return(0.0);
	sigma0*=sqrt(sigma1/(n*Es2c));
	return(sigma0);
}
void FFOgkBasis(
		const MatrixXd& xi,
		const int& calcM,
		const int& intercept,
		VectorXi& warn,
		const int& h,
		VectorXi& dIn,
		int w3
	){
	double (*pFo[])(Ref<VectorXd>,int)={&qn,&scaleTau2}; 
	double (*qFo[])(Ref<VectorXd>,int)={&Fmedian,&scaleTau2}; 
	const int p=xi.cols(),n=xi.rows(),h0=(n+1)/2;	
	double b1=0.0,b2=0.0;
	const double tol=1e-8;
	int i,j;
	MatrixXd x=xi;
	RowVectorXd lamba(p);
	MatrixXd x2=x;
	if(intercept){
		for(i=0;i<p;i++)	lamba(i)=qCalc(x2.col(i),1,qFo[calcM]);
		x.rowwise()-=lamba;	
	}
	for(i=0;i<p;i++)		lamba(i)=pCalc(x2.col(i),0,pFo[calcM]);
	for(i=0;i<p;i++)		warn(i)=(lamba(i)<tol)?1:0;
	i=warn.sum();
	if(i>0)				return;
	for(i=0;i<p;i++)		x.col(i).array()/=lamba(i);
	VectorXd dvec1=VectorXd::Ones(p);
	MatrixXd U=dvec1.asDiagonal();
	VectorXd sYi(n);
	VectorXd sYj(n);
	VectorXd dY(n);
	for(i=0;i<p;++i){
		sYi=x.col(i);
		for(j=0;j<i;++j){
			sYj=x.col(j);
			dY=sYi+sYj;
			b1=pCalc(dY,0,pFo[calcM]);
			b1*=b1;
			dY=sYi-sYj;
			b2=pCalc(dY,0,pFo[calcM]);
			b2*=b2;
			U(i,j)=0.25*(b1-b2);
			U(j,i)=U(i,j);	
		}		
	}
	JacobiSVD<MatrixXd> svd(U,ComputeThinV);
	x2=x*svd.matrixV();
	for(i=0;i<p;i++)		lamba(i)=pCalc(x2.col(i),0,pFo[calcM]);
	for(i=0;i<p;i++)		warn(i)=(lamba(i)<tol)?1:0;
	i=warn.sum();
	if(i>0)				return;
	for(i=0;i<p;i++)		x2.col(i).array()/=lamba(i);
	dY=x2.array().abs2().rowwise().sum();
	dIn.setLinSpaced(n,0,n-1);
	std::nth_element(dIn.data(),dIn.data()+h,dIn.data()+dIn.size(),IdLess(dY.data()));
	cov_CStep(dIn,x,h,h,w3);
	return;
}
double unimcd_in(
		const VectorXd& m_resd,
		const int& h
	){
	const int n1=m_resd.size(),len=n1-h+1;
	double initmean=0.0,initcov=0.0,sumw=0.0;
	int minone;
	if(h==n1){
		initmean=m_resd.sum()/(double)h;
		initcov=(m_resd.array()-initmean).abs2().sum()/(double)(h-1);
		return(sqrt(initcov));
	}
	VectorXd y=m_resd;
	VectorXd ay(len);
	VectorXd ay2(len);
	VectorXd sq(len);
	VectorXd y2(n1);

	std::sort(y.data(),y.data()+y.size());
	ay(0)=y.head(h).sum();	
	for(int samp=1;samp<len;samp++) ay(samp)=ay(samp-1)-y(samp-1)+y(samp+h-1);
	ay2=ay.array().square()/(double)h;
	y2=y.array().square();
	sq(0)=y2.head(h).sum()-ay2(0);
	for(int samp=1;samp<len;samp++) sq(samp)=sq(samp-1)-y2(samp-1)+y2(samp+h-1)-ay2(samp)+ay2(samp-1);
	initcov=sq.minCoeff(&minone);
	initcov/=(double)(h-1);
	initmean=ay(minone)/(double)h;
	return(initmean);
}
void CstepPrep(
		const MatrixXd& x,
		const VectorXd& y,
		const VectorXi& dIn,
		VectorXd& dP,
		const int& h0
	){
	const int n=x.rows(),p=x.cols();
	double w1,w0;
	int w2=1,i,j=0;
	MatrixXd xSub(h0,p);
	VectorXd ySub(h0);
	for(i=0;i<h0;i++) 	xSub.row(i)=x.row(dIn(i));
	for(i=0;i<h0;i++) 	ySub(i)=y(dIn(i));
	MatrixXd Sig(p,p);
	//Sig=xSub.adjoint()*xSub;
	//LDLT<MatrixXd> chol=Sig.ldlt();
	Sig.setZero().selfadjointView<Lower>().rankUpdate(xSub.transpose());
	LLT<MatrixXd> chol=Sig.llt();
	VectorXd m_cf(p);
	m_cf=chol.solve(xSub.adjoint()*ySub);
	dP=((x*m_cf).array()-y.array()).abs();
}
VectorXi CStep(
		const VectorXd& dPin,
		const MatrixXd& x,
		const VectorXd& y,
		const VectorXi& h,
		VectorXd& objfun,
		const int& I,
		VectorXi& citer
	){
	const int n=x.rows(),p=x.cols();
	double w1,w0;
	int w2=1,i,j=0;
	MatrixXd xSub(h(I),p);
	VectorXd ySub(h(I));
	VectorXd m_cf(p);
	MatrixXd Sig(p,p);
	VectorXd dP(n);
	VectorXi dIn(n);
	dP=dPin;
	w1=std::numeric_limits<double>::max();
	while(w2){	
		dIn.setLinSpaced(n,0,n-1);
		std::nth_element(dIn.data(),dIn.data()+h(I),dIn.data()+n,IdLess(dP.data()));
		for(i=0;i<h(I);i++) 	xSub.row(i)=x.row(dIn(i));
		for(i=0;i<h(I);i++) 	ySub(i)=y(dIn(i));
		//Sig=xSub.adjoint()*xSub;
		//LDLT<MatrixXd> chol=Sig.ldlt();
		Sig.setZero().selfadjointView<Lower>().rankUpdate(xSub.transpose());
		LLT<MatrixXd> chol=Sig.llt();
		m_cf=chol.solve(xSub.adjoint()*ySub);
		dP=y.array()-(x.leftCols(p-1)*m_cf.head(p-1)).array();
		m_cf(p-1)=unimcd_in(dP,h(I));
		dP=((x*m_cf).array()-y.array()).abs();
		w0=w1;
		j++;
		w1=log(((xSub*m_cf).array()-ySub.array()).abs().mean());
		(w0-w1<1e-3)?(w2=0):(w2=1);
	}
	objfun(I)=w1;
	citer(I)=j;
	return(dIn.head(h(I)).array());
} 
extern "C"{
	void R_FastOGK(
			int* rn,	//01
			int* rp,	//02
			double* X,	//03
			int* P,		//04
			int* cMet,	//05
			int* rint,	//06
			int* rwarn,	//07
			int* h,		//08
			int* lh,	//09
			int* rraw,	//10
			double* Ofuns,	//11
			int* doCstep,	//12
			int* Citer,	//13
			int* wout3	//14
		){
		const int CalcMet=*cMet,intercept=*rint,n=*rn,p=*rp,lhi=*lh,DoCsteps=*doCstep;
		int w3=1;
		MatrixXd xi=Map<MatrixXd>(X,n,p);
		VectorXi hi=Map<VectorXi>(h,lhi);
		VectorXd of(lhi);
		VectorXi ct(lhi);
		VectorXi warn(p);
		VectorXi dIn(n);
		MatrixXi hsub(n,lhi);
		warn.setZero();
		
		FFOgkBasis(xi,CalcMet,intercept,warn,hi(0),dIn,w3);
		Map<VectorXi>(rraw,hi(0))=dIn.head(hi(0)).array()+1;
		Map<VectorXi>(rwarn,p)=warn.array();
		int ws=warn.sum();
		*wout3=w3;
		if(ws==0 && DoCsteps==1 && w3==1){
			VectorXd dP(n);
			VectorXd yi=xi.col(p-1);
			xi.col(p-1).setOnes();
			CstepPrep(xi,yi,dIn,dP,hi(0));
			for(int i=0;i<lhi;i++)	hsub.col(i).head(hi(i))=CStep(dP,xi,yi,hi,of,i,ct);
			Map<MatrixXi>(P,n,lhi)=hsub.array()+1;
			Map<VectorXd>(Ofuns,lhi)=of.array().exp();
			Map<VectorXi>(Citer,lhi)=ct;
		} 
	}
}
