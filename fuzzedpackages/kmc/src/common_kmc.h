/* 
 * File:   kmc_comm.h
 * Author: Yifan Yang
 * Newton Raphson - root finding
 * Jacob Updating: oordinate descent
 * Created on November 2, 2014, 5:39 PM
 */

#ifndef KMC_COMM_H
#define	KMC_COMM_H

#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;


inline double perturb_c(double value,double pert) {
 /* A small, positive value, not too small*/
    double re;
    value= value>0 ? value : -value;
    re = value>1 ? value*pert:pert;
    return(re);
}



void rootfiding_kmc(
        void f(double* ,double *,vector<double>&),// para, data, re
        double * init_value,
        double *para,
        double *data,
        int *n, // sample size
        int * p,// number of parameters
        int * maxiter

){
    int pp=p[0];
    Eigen::MatrixXd Jacb(pp,pp);
    vector<double> delt(pp);
    Eigen::MatrixXd relchange(pp,1);
    Eigen::MatrixXd refxEigen(pp,1);
    double * x;
    double * refx;
    refx=new double[pp];
    x=new double[pp];
    
    vector<double> fx(pp);
    vector<double> fx1(pp);
    //Initialization:
    for (int i=0;i<pp;i++){
        x[i]=init_value[i];
    }
    
    int MAXITER=maxiter[0];
    
    for(int iter=0;iter<MAXITER;iter++){
        
    for (int i=0;i<pp;i++) {
        delt[i]= perturb_c(x[i],1e-8);
        refx[i]=x[i];
    }

    f(x,data,fx1);//reffx
    
    for (int j=0;j<pp;j++){
        x[j]+=delt[j];
        f(x,data,fx);
         for (int i=0;i<pp;i++){        
             Jacb(i,j)=(fx[i]-fx1[i])/delt[j];
        }
        x[j]=refx[j];//restore
     }
     
     //pp>1
     if(pp>1){
        for (int i=0;i<pp;i++) refxEigen(i,0)=fx1[i];
        relchange=(Jacb.ldlt().solve(-1.*refxEigen));
     }else{
         relchange(0,0)=-refxEigen(0,0)/Jacb(0,0);
     }
     for(int i=0;i<pp;i++){
         x[i]+=relchange(i,0);
     }
    }
    for (int i=0;i<pp;i++) para[i]=x[i];
     delete[] refx;
     delete[] x;
}

#endif	/* KMC_COMM_H */

