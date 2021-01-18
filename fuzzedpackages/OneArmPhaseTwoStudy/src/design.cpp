#include "design.h"

Design::Design()
{
    logFacLookupTable = new std::map<int, long double>();
    binLookupTable = new std::map<BinomKey, long double>();
}

Design::~Design()
{
    delete logFacLookupTable;
    delete binLookupTable;
}


long double Design::logFact(int n)
{
    std::map<int,long double>::iterator it = logFacLookupTable->find(n);
    if(it != logFacLookupTable->end())
    {
        return it->second;
    }
    else
    {
        if(n <= 0)
        {
            logFacLookupTable->insert(std::pair<int,long double>(n,log((float)1)));
        }
        long double sum = 0;
        int i;
        for(i=n; i>1;i--)
            sum = sum + log((float)i);
        logFacLookupTable->insert(std::pair<int,long double>(n,sum));
        return sum;
    }
}

long double Design::binomKoeff(int n, int k)
{
    long double koeff = exp(logFact(n) - logFact(n - k) - logFact(k));
    return koeff;
}


long double Design::bin(int n, int r, long double p)
{
    if(r>n)
        return 0;
    else
    {
        BinomKey nrp(n,r, p);
        std::map<BinomKey, long double>::iterator it = binLookupTable->find(nrp);
        if(it != binLookupTable->end())
        {
            return it->second;
        }
        else
        {
            long double test = exp(logFact(n) - logFact(n-r) -logFact(r));
            test = test *pow(p,r);
            test = test *pow(1-p,n-r);
            binLookupTable->insert(std::pair<BinomKey, long double>(nrp, test));
            return test;
        }
    }
}

long double Design::binsum(int n, int r, long double p)
{
    long double sum =0;
    int i;
    for(i=0;i<=r;i++)
        sum += bin(n,i,p);
    return sum;
}


void Design::setAlpha(double a){ alpha = a;}

double Design::getAlpha(){return alpha;}

void Design::setBeta(double b){ beta = b;}

void Design::setP0(double p0){ this->p0 = p0;}

void Design::setP1(double p1){ this->p1 = p1;}

long double Design::calcAlpha(int n1, int r1, int n, int r, long double p0)
{
    return 0;
}

long double Design::calcBeta(int n1, int r1, int n, int r, long double p1)
{
    return 0;
}

long double Design::calcAlpha(int n1, int r1, int n, int r, int s, long double pc0, long double pt0)
{
  return 0;
}

long double Design::calcBeta(int n1, int r1, int n, int r, int s, long double pc1, long double pt1)
{
  return 0;
}
