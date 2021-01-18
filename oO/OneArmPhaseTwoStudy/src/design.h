#ifndef DESIGN_H
#define DESIGN_H

#include <map>
#include <vector>
#include <math.h>


//Base class for the implemented design types (simon's two-stage design, subset design)
class Design
{
public:
    //Used as key in the "binLookupTable".
    struct BinomKey
    {
    public:
        BinomKey() :  n_int(0), r_int(0), p_int(0) {}
        BinomKey(int n, int r, long double p) : n_int(n), r_int(r), p_int(p) {}
        int n_int,r_int;
        long double p_int;

        friend bool operator <(const Design::BinomKey &lhs, const Design::BinomKey &rhs)
        {
            if(lhs.n_int != rhs.n_int){
                return lhs.n_int < rhs.n_int;
            }else if(lhs.r_int != rhs.r_int){
                return lhs.r_int < rhs.r_int;
            }else{
                return lhs.p_int < rhs.p_int;
            }
        }
    };

    Design();
    ~Design();
    // Sets the maximal type I error rate.
    virtual void setAlpha(double a);
    // Returns the maximal type I error rate.
    virtual double getAlpha();
    // Sets the maximal type II error rate.
    virtual void setBeta(double b);
    // Sets the response probability under the null hypothesis.
    virtual void setP0(double p0);
    // Sets the response probability under the alternative hypothesis.
    virtual void setP1(double p1);

protected:
    int n1, n, maxn;
    int miniMaxPos, optimalPos;
    double alpha, beta, p0, p1;
    std::map<int,long double> *logFacLookupTable;
    std::map<BinomKey, long double> *binLookupTable;
    
    // Calculates the factorial of n.
    virtual long double logFact(int n);
    virtual long double binomKoeff(int n, int k);
    virtual long double bin(int n, int r, long double p);
    virtual long double binsum(int n, int r, long double p);
    virtual long double calcAlpha(int n1, int r1, int n, int r, long double p0);
    virtual long double calcBeta(int n1, int r1, int n, int r, long double p1);
    virtual long double calcAlpha(int n1, int r1, int n, int r, int s, long double pc0, long double pt0);
    virtual long double calcBeta(int n1, int r1, int n, int r, int s, long double pc1, long double pt1);

};

#endif // DESIGN_H
