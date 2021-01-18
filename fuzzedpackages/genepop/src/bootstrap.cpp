/***************************************************************************
@ F. Rousset 2005-2007

francois.rousset@umontpellier.fr

This file is part of Genepop'007
This software is a computer program whose purpose is to perform statistical analyses.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

 ***************************************************************************/

#include <iostream>
#include <limits>
#include "GenepopS.h"
#include "F_est.h"
#include "bootstrap.h"
#include "settings.h"
#include "proba.h"
#include "myutils.h"
#include "tools.h"
#ifdef COMPATIBILITYRCPP
#include <Rcpp.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#endif

using namespace std;

/** in no way specific to bootstrap**/
vector<double> bisection_search(double (*func)(double d),double x1,double x2,bool verbose) {
//R: function (x1, x2, f, tol = 1e-07, niter = 25, f.extra = NA, upcross.level = 0)
// rtbis, Numerical recipes
//Using bisection, find the root of a function func known to lie between x1 and x2. The root,
//returned as rtbis, will be re...ned until its accuracy is ...xacc.
    vector<double>v(0); // [0] : indicateur probleme eventuel; [1]: actual P value
    const double xacc=numeric_limits<double>::epsilon()*(abs(x1)+abs(x2))/2; //suggestion Num. Rec.
    const int maxit=-2*int(log(numeric_limits<double>::epsilon())/log(2.));
    int j;
    double dx,f,fmid,xmid,rtb;
    f=(*func)(x1);
    fmid=(*func)(x2);
    if (f*fmid >= 0.0) {
        if (verbose) {
                noR_cout<<"(!) From CKrig::bisection_search() : Root must be bracketed for bisection. "<<endl;
                noR_cout<<"   x1, f(x1), x2, f(x2) were "<<x1<<" "<<f<<" "<<x2<<" "<<fmid<<endl;
        }
        //if (cinGetOnError) cin.get();
        v.push_back(-1);
        return(v); // length 1
    }
    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2); //Orient the search so that f>0
/*std::cout<<x1<<" ";
std::cout<<x2<<" ";
std::cout<<dx<<" ";
std::cout<<rtb;getchar();*/
    for (j=1;j<=maxit;j++) { //lies at x+dx.
        fmid=(*func)(xmid=rtb+(dx *= 0.5)); //Bisection loop.
        if (fmid <= 0.0) rtb=xmid;
        if (fabs(dx) < xacc || fmid == 0.0) {
            v.push_back(0);
            v.push_back(rtb);
            return(v);
        }
    }
    if (f*fmid >= 0.0) { //mainly means !=0...
        if (verbose) {
                noR_cout<<"(!) From CKrig::bisection_search() : Too many bisections. "<<endl;
        }
        //if (cinGetOnError) cin.get();
        v.push_back(-1);
        v.push_back(numeric_limits<double>::quiet_NaN());  // visible problem; otherwise perhps xmid ?
        return(v);
    }
// should never reach:
v.push_back(numeric_limits<double>::signaling_NaN());
return(v);
}


static CABCbootstrap* ABCptr;

CABCbootstrap::CABCbootstrap(size_t units) {
    nb_units=units;
    z=std::numeric_limits<double>::quiet_NaN();
}

int CABCbootstrap::bootstrapOverLoci(double (*estimatingFnPtr)(vector<double> d),string legend,
                                     string bootOutfile,bool clear_screen) {
    ofstream bootOut;
    estimFnPtr=estimatingFnPtr;
    testLegend=legend; // copy to namespace variable used by other function
    //nb_units is a numbere of relevant loci (in particular with the relevant ploidy)
    // it is the function called through estimatingFnPtr that must map the weights to the chosen subset of loci
    // using a set of indices for the mapping (see varForBootstrapGenepop::idxPloid in wrapper() )
    if (nb_units==0) {
        #ifdef COMPATIBILITYRCPP
            Rcpp::Rcerr<<"(!) 0 replicates to bootstrap over";
        #else
            cerr<<"(!) 0 replicates to bootstrap over";
        #endif
        return(0);
    }
	delta.resize(nb_units);
	if (clear_screen) effacer_ecran();
	vector<double> locABCweight(nb_units);
	vector<double> dt(nb_units);
	vector<double> ddt(nb_units);
	vector<double>tminput(nb_units);
	vector<double>tpinput(nb_units);
    double epsn,tm,tp,cq,bhat,curv,machin;//ABC bootstrap variables
	double sigmahat=0.0;
    double epsn_value=0.00002;
    seuil_inf=(1.0-widthCI)/2.0;
    seuil_sup=1.0-seuil_inf;


     //point estimate for original data (ABCweigth constant)
	for(size_t loc=0;loc<nb_units;loc++) {
           locABCweight[loc]=1.0/nb_units;
    }
	t0=(*estimatingFnPtr)(locABCweight);  //first call to isoldeproc sets _first_of_repl=false (and if it was =true, also clears screen !)
    // does not delete the .MIG file here
    if (std::isnan(t0)) {
        tinf=numeric_limits<double>::quiet_NaN();
        tsup=numeric_limits<double>::quiet_NaN();
        noR_cout<<"\n No valid estimate for original data. Skipping further bootstrap computation.\n";
        if (pauseGP) {
            noR_cout<<"(Return) to continue";cin.get();
        }
        genepop_exit(-1, "No valid estimate for original data. Skipping further bootstrap computation.");
    }
    if (nb_units<3) {
        bootOut.open(bootOutfile.c_str(),ios::app);
        if ( ! bootOut.is_open()) {
            #ifdef COMPATIBILITYRCPP
                // Rcpp::Rcerr<<"(!) From bootstrapOverLoci(): error while opening file "<<bootOutfile<<endl;
            #else
                cerr<<"(!) From bootstrapOverLoci(): error while opening file "<<bootOutfile<<endl;
                if (cinGetOnError) cin.get();
            #endif

	       genepop_exit(-1, "(!) From bootstrapOverLoci(): error while opening file ");
        }
        noR_cout<<"\nNot enough loci to bootstrap over\n";
        bootOut<<"\nNot enough loci to bootstrap over\n";
        bootOut.close();
        if (pauseGP) {
            noR_cout<<"(Return) to continue"<<endl; getchar();
        }
        genepop_exit(-1, "Not enough loci to bootstrap over");
    }
    /**ELSE**/
	int consy=wherey();
    _gotoxy(0,consy+2);
    noR_cout<<" Computing [ "<<seuil_inf<<"- "<<seuil_sup<<" ] confidence interval"<<legend<<":";
	//pour borne inf
	epsn=-epsn_value;
#ifdef COMPATIBILITYRCPP
	Rcpp::Rcerr<<"Computing confidence interval"<<legend<<":"<<endl;
	Progress progbar(2*nb_units, true);
#endif
	for(ABCloc=0;ABCloc<nb_units;ABCloc++) {
        for(size_t loc=0;loc<nb_units;loc++)
            locABCweight[loc]= - epsn/nb_units+1.0/nb_units;
        locABCweight[ABCloc]+= epsn;
//ostream_vector(locABCweight,cout);getchar();
        tminput[ABCloc]=(*estimatingFnPtr)(locABCweight);
//cout<<tminput[ABCloc];getchar();
#ifdef COMPATIBILITYRCPP
progbar.increment();
#else
_gotoxy(0,consy+3);
noR_cout<<" Computing confidence interval... about      % done             ";
_gotoxy(41,consy+3);
noR_cout<< int(100*(double(ABCloc)+1.0)/(2*double(nb_units)+5.0));
#endif
    }

	//pour borne sup
	epsn=epsn_value;
	for(ABCloc=0;ABCloc<nb_units;ABCloc++) {
        for(size_t loc=0;loc<nb_units;loc++)
            locABCweight[loc]= - epsn/nb_units+1.0/nb_units;
        locABCweight[ABCloc]+= epsn;
        tpinput[ABCloc]=(*estimatingFnPtr)(locABCweight);
#ifdef COMPATIBILITYRCPP
        progbar.increment();
#else
        _gotoxy(0,consy+3);
        noR_cout<<" Computing confidence interval... about      % done            ";
        _gotoxy(41,consy+3);
        noR_cout<< int(100*(double(ABCloc)+double(nb_units)+1.0)/(2*double(nb_units)+5.0));
#endif
    }

	for(size_t loc=0;loc<nb_units;loc++) {
	  dt[loc]=(tpinput[loc]-tminput[loc])/(2*epsn);
		ddt[loc]=(tpinput[loc]+tminput[loc]-2*t0)/pow(epsn,2);
//cout<<tpinput[loc]<<" "<<dt[loc]<<" "<<ddt[loc]<<endl;
		}
//getchar();


	for(size_t loc=0;loc<nb_units;loc++) sigmahat+=pow(dt[loc],2);
	sigmahat=sqrt(sigmahat)/nb_units;
    if (sigmahat==0) {
        bootOut.open(bootOutfile.c_str(),ios::app);
        if ( ! bootOut.is_open()) {
            #ifdef COMPATIBILITYRCPP
                // Rcpp::Rcerr<<"(!) From bootstrapOverLoci(): error while opening file "<<bootOutfile<<endl;
            #else
                cerr<<"(!) From bootstrapOverLoci(): error while opening file "<<bootOutfile<<endl;
                if (cinGetOnError) cin.get();
            #endif

	       genepop_exit(-1, "(!) From bootstrapOverLoci(): error while opening file ");
        }
        noR_cout<<"\nNot enough data for constructing confidence interval\n";
       bootOut<<"\nNot enough data for constructing confidence interval\n";
       bootOut.close();
       if (pauseGP) {
            noR_cout<<"(Return) to continue"<<endl; getchar();
       }
       genepop_exit(-1, "Not enough data for constructing confidence interval");
    }

    ahat=0;
	for(size_t loc=0;loc<nb_units;loc++) ahat+=pow(dt[loc],3);
	ahat=ahat/(6*pow(double(nb_units),3)*pow(sigmahat,3));
	for(size_t loc=0;loc<nb_units;loc++) delta[loc]=dt[loc]/(pow(double(nb_units),2)*sigmahat);

	epsn=-epsn_value;
	for(size_t loc=0;loc<nb_units;loc++)
		locABCweight[loc]=delta[loc]*epsn+(1.0/nb_units);
	tp=(*estimatingFnPtr)(locABCweight);
	for(size_t loc=0;loc<nb_units;loc++)
		locABCweight[loc]=-delta[loc]*epsn+(1.0/nb_units);
	tm=(*estimatingFnPtr)(locABCweight);
    _gotoxy(0,consy+3);
    noR_cout<<" Computing confidence interval... about      % done              ";
    _gotoxy(41,consy+3);
    noR_cout<< int(100*(2*double(nb_units)+2)/(2*double(nb_units)+5));
	cq=(tp+tm-2*t0)/(2*sigmahat*pow(epsn,2));
	bhat=0.0;
	for(size_t loc=0;loc<nb_units;loc++) bhat+=ddt[loc];
	bhat/=(2*pow(double(nb_units),2));
	curv=(bhat/sigmahat)-cq;
	machin=2*ndtr(ahat)*ndtr(-curv);

	if((machin>=1)||(machin<=0.0)) {//z=INFINITY;
        bootOut.open(bootOutfile.c_str(),ios::app);
        if ( ! bootOut.is_open()) {
            #ifdef COMPATIBILITYRCPP
                // Rcpp::Rcerr<<"(!) From bootstrapOverLoci(): error while opening file "<<bootOutfile<<endl;
            #else
                cerr<<"(!) From bootstrapOverLoci(): error while opening file "<<bootOutfile<<endl;
                if (cinGetOnError) cin.get();
            #endif

	       genepop_exit(-1, "(!) From bootstrapOverLoci(): error while opening file ");
        }

         noR_cout<<"Putative confidence interval"<<legend<<" has infinite range.\n Bootstrap computation aborted.";
	     bootOut<<"Putative confidence interval"<<legend<<" has infinite range.\n Bootstrap computation aborted.";
	     bootOut.close();
	   	 if (pauseGP) {
            noR_cout<<"\n(Return) to continue"<<endl; getchar();
         }
	   	 return 1; // exit() !
    } else z=ndtri(machin);

	bidullevel=(z+ndtri(seuil_inf))/(pow(1-ahat*(z+ndtri(seuil_inf)),2));
	for(size_t loc=0;loc<nb_units;loc++) //pour borne inf
		locABCweight[loc]=delta[loc]*bidullevel+1.0/nb_units;
	tinf=(*estimatingFnPtr)(locABCweight);
    _gotoxy(0,consy+3);
    noR_cout<<" Computing confidence interval... about      % done             ";
    _gotoxy(41,consy+3);
    noR_cout<< int(100*(2*double(nb_units)+4)/(2*double(nb_units)+5));

	bidullevel=(z+ndtri(seuil_sup))/(pow(1-ahat*(z+ndtri(seuil_sup)),2));
	for(size_t loc=0;loc<nb_units;loc++) //pour borne sup
		locABCweight[loc]=delta[loc]*bidullevel+1.0/nb_units;
	tsup=(*estimatingFnPtr)(locABCweight);
	if (tsup<tinf) {double stat=tsup;tsup=tinf;tinf=stat;}

    _gotoxy(0,consy+3);
    noR_cout<<" Computing confidence interval... 100 % done                      ";
    noR_cout<<"\n ABC bootstrap results"<<":";
    noR_cout<<"\n Point estimate and "<<100*widthCI<<"% confidence interval"<<legend<<":"<<endl;
    noR_cout<<t0<<" [ "<<tinf<<" , "<<tsup<<" ]\n";

    bootOut.open(bootOutfile.c_str(),ios::app);
    if ( ! bootOut.is_open()) {
        #ifdef COMPATIBILITYRCPP
            // Rcpp::Rcerr<<"(!) From bootstrapOverLoci(): error while opening file "<<bootOutfile<<endl;
        #else
            cerr<<"(!) From bootstrapOverLoci(): error while opening file "<<bootOutfile<<endl;
            if (cinGetOnError) cin.get();
        #endif

        genepop_exit(-1, "(!) From bootstrapOverLoci(): error while opening file ");
    }
	if ( ! perf) {	 //ecriture des resultats
        bootOut<<"\n ABC bootstrap results"<<legend<<":";
	    bootOut<<"\n Point estimate and "<<100*widthCI<<"% confidence interval:\n"<<t0<<" [ "<<tinf<<" , "<<tsup<<" ]\n";
 	}
 	if ( std::strcmp(legend.c_str(), " for SLOPE")==0 && ! std::isnan(testPointslope)) {
 	    //FR->FR comme on n'est pas surde travailler sur des slopes ici, il faut un mecanisme genre migraine sur les testPoint
        testPointPvalue=Pvalue(testPointslope,true,false);
        if (std::isnan(testPointPvalue)) Pvalue(testPointslope,true,true); // verbose version

        noR_cout<<"\nUnidirectional P-value for [ slope="<<testPointslope<<" ] : "<<testPointPvalue<<endl;

        bootOut<<"\nUnidirectional P-value for [ slope="<<testPointslope<<" ] : "<<testPointPvalue<<endl;
    }
    bootOut.close();
    ///// bug corrected 2015/12/8 : added lines:
	for(size_t loc=0;loc<nb_units;loc++) {locABCweight[loc]=1.0/nb_units;}
	t0=(*estimatingFnPtr)(locABCweight);  //need to reset any variabes modified by estimatingFnPtr (slope() -> datamatrix::data used later by mantel test)
return 0;
}

double CABCbootstrap::cancelland(double unidirPvalue) { //wrapper from (double) to (vector<double>)
    vector<double>weights(nb_units);
	bidullevel=(z+ndtri(unidirPvalue))/(pow(1-ahat*(z+ndtri(unidirPvalue)),2));
	for(size_t loc=0;loc<nb_units;loc++) //pour borne inf
		weights[loc]=delta[loc]*bidullevel+1.0/nb_units;
	double t=(*estimFnPtr)(weights);
	return(t-testPoint); // bisection_search must seek the zero of this. testPoint is global within the file...
}

double cancellandWrapper(double unidirPvalue) {
    //allows a member function to be called through an ordinary function pointer ie as first argument of bootstrapOverLoci
    return(ABCptr->cancelland(unidirPvalue)); // of course one needs a global pointer to the appropriate class
}

double CABCbootstrap::Pvalue(double testPt,bool unidir=false,bool verbose=true) { // gives (by default unidirectional) P-value associated with (currently) slope value
    /** returns a valid P value, or else
        NaN -> invalid call, no valid CI previouly computed
        -1 -> failure of bisection search
    **/

    vector<double> resu(1,-1);
    double guess=1;
    ABCptr=this;
    testPoint=testPt;  // for later uses outside Pvalue()
    double sig,level1=0.0,level2=0.0,pvalue;
    if(std::isnan(z)) {
        noR_cout<<"Attempt to compute P value by ABC bootstrap\n    while confidence interval computation was not called, or failed.";
	     //(need to have bootOut open... )
	     //bootOut<<"Attempt to compute P value by ABC bootstrap\n    while confidence interval computation was not called, or failed.";
	     //bootOut.close();
	   	 if (pauseGP) {

            noR_cout<<"\n(Return) to continue"<<endl; getchar();

        }
	   	 return(numeric_limits<double>::quiet_NaN()); // not exit() !
    } /**ELSE **/

        noR_cout<<" Computing test"<<testLegend<<"= "<<testPoint<<"; beginning..";


    int it=1;
    while(resu.size()==1 && it <50) {  //resu.size()=1means bad starting values
        // note that resu.size()=2 either when solution is found or for bisection failure to converge (visible NaN)
        if (verbose && it==2) {
                noR_cout<<"(*) From Pvalue(): Problem finding starting values for bisection search"<<endl;
                noR_cout<<"tinf, t0, tsup were "<<tinf<<" "<<t0<<" "<<tsup<<endl;
                noR_cout<<"Initial levels were "<<level1<<" "<<level2<<endl;
        }
        guess*=10;
        /** just initial values... **/
        if (testPt<t0) {
            if (tinf<testPt) { // testtPt in ]tinf,t0[
                level1=seuil_inf;level2=0.5+0.01*it;
            } else { // testtPt in ]-inf,tinf]
                sig=(t0-tinf)/(-ndtri(seuil_inf)); // 1 if N(t0,1) CI => estimation sigma in N(t0,sigma2) CI
                //seuil_inf = 0.025 du CI would result in unidirPvalue=0.05
                level2=1-(1-0.01*it)*(1-2*seuil_inf);
                //level1 better not =0...
                level1=max(ndtr(guess*(testPt-t0)/sig),min(level2/2,pow(0.01,double(it)/5+1)));
            }
        } else {
            if (testPt<tsup) { // testtPt in [t0,tsup[
                level1=0.5-0.01*it;level2=seuil_sup;
            } else { // testtPt in [tsup,inf[
                //seuil_sup = 0.975 du CI would result in unidirPvalue=0.95
                level1=(1-0.01*it)*(1-2*(1-seuil_sup));
                sig=(tsup-t0)/(ndtri(seuil_sup)); // 1 if N(t0,1) CI => estimation sigma in N(t0,sigma2) CI
                //level2 better not =1...
                level2=min(ndtr(guess*(testPt-t0)/sig),max(1.-(1-level1)/2,1-pow(0.01,double(it)/5+1)));
            }
        }
        if (it>1 && verbose) {

            noR_cout<<"New initial levels "<<level1<<" "<<level2<<endl;

        }
        it++;
        /** main stuff **/
        resu=bisection_search(cancellandWrapper,level1,level2,verbose);
    }
    if (resu.size()==1) { // Exited while() without good starting values
        if (verbose) {
                noR_cout<<"(!) From Pvalue(): Failed to find starting values for bisection search";
                noR_cout<<"tinf, t0, tsup were "<<tinf<<" "<<t0<<" "<<tsup<<endl;
        }
        // now a reasonable guess as to why good starting values were not found...
        if (testPt<tinf) resu.push_back(0);
        else if (tsup<testPt) resu.push_back(1);
        else resu.push_back(numeric_limits<double>::quiet_NaN()); // no clear reason to reach this point => make it visible
    }
    if ( ! unidir) pvalue=2*(min(resu[1],1.-resu[1])); else pvalue=resu[1];
    return(pvalue);
}
