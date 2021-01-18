#include "runTest.h"
#include "runFunctions.h"

using namespace Rcpp ;
using namespace std ;

vector<vector<double> > convertToVVd(SEXP const& rMatrix)
{
  NumericMatrix data(rMatrix);
  int n = data.nrow(), p=data.ncol();
  vector<vector<double> > output(n,vector<double> (p));
  for(int i=0;i<n;i++)
    for(int j=0;j<p;j++)
      output[i][j]=data[j*n+i];

  return output;
}

vector<vector<int> > convertToVVi(SEXP const& rMatrix)
{
  NumericMatrix data(rMatrix);
  int n = data.nrow(), p=data.ncol();
  vector<vector<int> > output(n,vector<int> (p));
  for(int i=0;i<n;i++)
    for(int j=0;j<p;j++)
      output[i][j]=data[j*n+i];

  return output;
}

vector<Rank> downUniVariateRank(NumericMatrix XR)
{
    int const n(XR.nrow());
	int const m(XR.ncol());
    vector<Rank> donnees(n);
    set<int> element;

    for(int i(1);i<m+1;i++)
		element.insert(i);

    for(int j(0);j<n;j++)
		donnees[j].rank.resize(m);

    for(int j(0);j<n;j++)
    {
		donnees[j].missingNumber=element;

        donnees[j].isPartial=false;

        for(int i(0);i<m;i++)
        {
			donnees[j].rank[i]=XR[j+i*n];
			if(donnees[j].rank[i]==0)
			{
				donnees[j].isPartial=true;
				donnees[j].missingIndex.push_back(i);
			}
			else
				donnees[j].missingNumber.erase(donnees[j].rank[i]);
		}
	}

    return donnees;
}


RcppExport SEXP kullback(SEXP m,SEXP mu1,SEXP mu2,SEXP p1, SEXP p2,SEXP proportion1,SEXP proportion2)
{
	//conversion
    NumericVector proportion1R(proportion1),proportion2R(proportion2);
	NumericVector mR(m);
    vector<int> mC=as<vector<int> > (mR);
	vector<double> prop1=as<vector<double> > (proportion1R);
	vector<double> prop2=as<vector<double> > (proportion2R);
	vector<vector<double> > p1C,p2C;
	p1C=convertToVVd(p1);
	p2C=convertToVVd(p2);
	NumericMatrix mu1R(mu1),mu2R(mu2);
	vector<vector<vector<int> > > mu1C,mu2C;
	mu1C=numMat2vvvInt(mu1,mC);

	mu2C=numMat2vvvInt(mu2,mC);

	//
	double dkl(0);
	dkl=divKL(mC,mu1C,mu2C,p1C,p2C,prop1,prop2);
	return wrap(dkl);
}


RcppExport SEXP adkhi2(SEXP donnees, SEXP p, SEXP proportion, SEXP mu, SEXP nBootstrap)
{
	//conversion

	int nBoot=as<int>(nBootstrap);
	vector<double> prop=as<vector<double> >(proportion);
	vector<double> pC=as<vector<double> >(p);
	vector<vector<int> > muC=convertToVVi(mu);
	vector<vector<int> > data=convertToVVi(donnees);


	//
	double pval(0);
	pval=khi2(data,pC,prop,muC,nBoot);
	return wrap(pval);
}

RcppExport SEXP adkhi2partial(SEXP donnees, SEXP p, SEXP proportion, SEXP mu, SEXP nBootstrap)
{
	//conversion
	int nBoot=as<int>(nBootstrap);
	vector<double> prop=as<vector<double> >(proportion);
	vector<double> pC=as<vector<double> >(p);
	vector<vector<int> > muC=convertToVVi(mu);
	NumericMatrix donneesR(donnees);
	vector<Rank > data=downUniVariateRank(donneesR);

	//
	double pval(0);
	pval=khi2partial(data,pC,prop,muC,nBoot);
	return wrap(pval);
}



