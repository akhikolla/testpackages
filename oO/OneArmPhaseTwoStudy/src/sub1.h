#ifndef SUB1DESIGN_H
#define SUB1DESIGN_H

#include "simon.h"
#include <Rcpp.h>
#include <map>
#include "design.h"
#include "resultsub1.h"
#include <math.h>

// An object of this class implements all methods needed for planning, monitoring and analysing subset designs. 
// Also all methods needed to apply (non-)stochastic curtailment are implemented.
class Sub1Design: public Design
{
public:
    //Used as key in the "multinomialLookupTable".
    struct MultiKey
    {
    public:
        MultiKey() : x1(0), x2(0), n_int(0), p1(0), p2(0){}
        MultiKey(int x_1, int x_2, int n, long double p_1, long double p_2) : x1(x_1), x2(x_2), n_int(n), p1(p_1), p2(p_2){}
        int x1, x2, n_int;
        long double p1, p2;

        friend bool operator <(const Sub1Design::MultiKey &lhs, const Sub1Design::MultiKey &rhs)
        {
            if(lhs.x1 != rhs.x1)
            {
                return lhs.x1 < rhs.x1;
            }else if(lhs.x2 != rhs.x2){
                return lhs.x2 < rhs.x2;
            }else if(lhs.n_int != rhs.n_int){
                return lhs.n_int < rhs.n_int;
            }else if(lhs.p1 != rhs.p1){
                return lhs.p1 < rhs.p1;
            }else{
                return lhs.p2 < rhs.p2;
            }
        }
    };
    
    //Used as key in the "alphaLookupTable" and "betaLookupTable".
    struct AlphaBetaKey
    {
    public:
      AlphaBetaKey() : sumboarder3(0), sumboarder4(0), n2(0){}
      AlphaBetaKey(int sumboarder3, int sumboarder4, int n2) : sumboarder3(sumboarder3), sumboarder4(sumboarder4), n2(n2){}
      int sumboarder3, sumboarder4, n2;

      friend bool operator <(const Sub1Design::AlphaBetaKey &lhs, const Sub1Design::AlphaBetaKey &rhs)
      {
        if(lhs.sumboarder3 != rhs.sumboarder3)
        {
          return lhs.sumboarder3 < rhs.sumboarder3;
        }
        else if (lhs.sumboarder4 != rhs.sumboarder4){
          return lhs.sumboarder4 < rhs.sumboarder4;
        }
        else
        {
          return lhs.n2 < rhs.n2;
        }
      }

      void setValues(int sumboarder3, int sumboarder4, int n2)
      {
        this->sumboarder3 = sumboarder3;
        this->sumboarder4 = sumboarder4;
        this->n2 = n2;
      }
    };


    Sub1Design();
    ~Sub1Design();
    
    // Sets the maximal type I error rate.
    void setAlpha(double a);
    // Returns the maximal type I error rate.
    double getAlpha();
    // Set the maximal type II error rate
    void setBeta(double b);
    // Sets the probability of an event under the null hypothesis.
    void setP0(double p0);
    // Sets the probability of an event under the alternative hypothesis.
    void setP1(double p1);
    // Sets the probability of an event in the subset endpoint under the null hypothesis.
    void setPc0(double pc0);
    // Sets the probability of an event in the superset endpoint under the null hypothesis.
    void setPt0(double pt0);
    // Sets the probability of an event in the subset endpoint under the alternative hypothesis.
    void setPc1(double pc1);
    // Sets the probability of an event in the superset endpoint under the alternative hypothesis.
    void setPt1(double pt1);
    // Approximates "maxN" (maximal number of patients to be recruited) in a way that it is most likely to identify the optimal design.
    // The value approximated is internally set for the member "maxn" and also returned by this function.
    long double aproximateMaxN();
    // Calculates possible designs for the set values "alpha", "beta", "pc0", "pc1", "pt0", "pt1".
    // skipS: Skips the iteration over "s" at certian points to improve calculation speed (finds less designs)
    // skipR: Skips the iteration over "r" at certian points to improve calculation speed (finds less designs)
    // skipN1: Skips the iteration over "n1" at certian points to improve calculation speed (finds less designs and it is impossible to determine the optimalization criteria of the found designs)
    // lowerBorder: Sets a minimal value for "n" (number of patients to be recruited)
    // upperBorder: Sets a maximal value for "n" (number of patients to be recruited) if set to "0" "maxn" is used instead
    void calculateStudySolutions(bool skipS, bool skipR, bool skipN1, int lowerBorder, int upperBorder);
    // Returns the number of found designs/solutions.
    int getSolutionCount();
    // Returns an R-representation of the identified designs.
    SEXP getResultsForR();    
    // Implements the binomial distribution.
    long double binsum(int n, int r, long double p);
    // Returns the probability of early termination for the given parameter set if the null hypothesis is true.
    // r1: Critical value for the first stage (more than "r1" responses in the subset endpoint needed to proceed to the second stage).
    // n1: Sample size for the first stage.
    // pc0: The response probability for the subset endpoint under the null hypothesis.
    long double calcPet_p0(int r1, int n1, long double pc0);
    // Returns the expected sample size for the given parameter set if the null hypothesis is true.
    // r1: Critical value for the first stage (more than "r1" responses in the subset endpoint needed to proceed to the second stage).
    // n1: Sample size for the first stage.
    // pet_p0: Probability of early termination under the null hypothesis.
    long double calcEN_p0(int n1, int n, long double pet_p0);    
    // Calculates the probability for a type I error.    
    // n1: Sample size for the first stage.
    // r1: Critical value for the first stage.
    // n: Overall sample size.
    // r: Critical value for the subset endpoint.
    // s: Critical value for the superset endpoint.
    // pc0: The response probability under the null hypothesis for the subset endpoint.
    // pt0: The response probability under the null hypothesis for the superset endpoint.
    long double calcAlpha(int n1, int r1, int n, int r, int s, long double pc0, long double pt0);
    // Calculates the probability for a type II error.
    // n1: Sample size for the first stage.
    // r1: Critical value for the first stage.
    // n: Overall sample size.
    // r: Critical value for the subset endpoint.
    // s: Critical value for the superset endpoint.
    // pc1: The response probability under the alternative hypothesis for the subset endpoint.
    // pt1: The response probability under the alternative hypothesis for the superset endpoint.
    long double calcBeta(int n1, int r1, int n, int r, int s, long double pc1, long double pt1);
    // Returns the exact p-value.
    // t: Observed responses in the subset endpoint.
    // u: Observed responses in the superset endpoint.
    // r1: Critical value for the first stage.
    // n1: Sample size for the first stage.
    // n: Overall sample size.
    // pc0: The response probability under the null hypothesis for the subset endpoint.
    // pt0: The response probability under the null hypothesis for the superset endpoint.
    long double get_p_exact(int t, int u, int r1, int n1, int n, long double pc0, long double pt0);
    // Returns the conditional power given that there were "t" responses observed in the subset endpoint and "u" responses in the superset endpoint.
    // t: Observed responses in the subset endpoint.
    // u: Observed responses in the superset endpoint.
    // enrolled: Number of patients enrolled so far.
    // r1: Critical value for the first stage.
    // n1: Sample size for the first stage.
    // r: Critical value for the subset endpoint.
    // s: Critical value for the superset endpoint.
    // n: Overall sample size.
    // pc1: The response probability under the alternative hypothesis for the subset endpoint.
    // pt1: The response probability under the alternative hypothesis for the superset endpoint.
    long double get_conditionalPower(int t, int u, int enrolled, int r1, int n1, int r, int s, int n, long double pc1, long double pt1);
    // Estimates the effect of (non-)stochastic curtailment for the design with the ID "resID".
    // cut: Sets the "cut point" used to calculate the effect of (non-)stochastic curtailment. A study is stopped if the conditional power falls below the value of "cut".
    // reps: Number of replications used to estimate the effect of (non-)stochastic curtailment.
    // all: If true the effect of (non-)stochastic curtailment will be calculated for different cut points in 0.05 steps starting with the value of the parameter "cut".
    void calculateSC(int resID, double cut, int reps , bool all);
    // Returns the ID of the minimax design under the found designs. (only used by the GUI for the R-package)
    int getMinimaxPos();
    // Returns the ID of the optimal design under the found designs. (only used by the GUI for the R-package)
    int getOptimalPos();
    
protected:
    long double logFact(int n);
    long double bin(int n, int r, long double p);
    
    
private:
    long double multinomial(int x1, int x2, int n, long double p1, long double p2);
    long double multinomialTest(int r1, int r, int s, int n1, int n, long double p1, long double p2);
    long double aproximateMaxNInternal(long double alpha, long double beta, long double pc0, long double pt0, long double pc1, long double pt1, int n1, int startValue);
    // Returns true if a design was found and saves the identified design internally.
    bool getDesign(int r1, int n1, int r, int s, int n, bool skipN1);
    ResultSub1::Curtailment_SubD1 calcSCintern(int r1, int r, int s, int n1, int n, long double pc1, long double pt1, long double pc0, long double pt0, long double cut, int rep);
    double calculateIntersection(double slope1, double yIntercept1, double slope2, double yIntercept2);
    // Identifies admissible designs based on the aproach deskribed in
    // "Kunz C (2011) To-Stage Designs for Phase II Trials with One or Two Endpoints. 
    // PH.D. thesis, Medical Faculty of Heidelberg, Rupecht-Karls-Universitaet" on page 16.
    void setAdmissible(std::vector<ResultSub1*>* results);
        

    std::map<MultiKey, long double> *multinomialLookupTable;
    std::map<AlphaBetaKey, long double> *alphaLookupTable;
    std::map<AlphaBetaKey, long double> *betaLookupTable;
    
    SimonDesign *simon;
    int r1, r, s;
    double pc0, pt0, pc1, pt1;
    long double enCurrent;
    int miniN;
    bool minMaxFound;
    bool firstDesignFound;
    
    std::vector<ResultSub1*> *allSub1Results;

};

#endif // SUB1DESIGN_H
