#ifndef SIMONDESIGN_H
#define SIMONDESIGN_H

#include <Rcpp.h>
#include <map>
#include <math.h>
#include "result.h"
#include "design.h"


// An object of this class implements all methods needed for planning, monitoring and analysing simon's two-stage designs. 
// Also all methods needed to apply (non-)stochastic curtailment are implemented.
class SimonDesign: public Design
{

public:


  SimonDesign();
  ~SimonDesign();
  
  // Sets the maximal type I error rate.
  void setAlpha(double a);
  // Returns the maximal type I error rate.
  double getAlpha();
  // Sets the maximal type II error rate.
  void setBeta(double b);
  // Returns the maximal type II error rate.
  double getBeta();
  // Sets the response probability under the null hypothesis.
  void setP0(double p0);
  // Returns the response probability under the null hypothesis.
  double getP0();
  // Sets the response probability under the alternative hypothesis.
  void setP1(double p1);
  // Returns the response probability under the alternative hypothesis.
  double getP1();
  // Sets the maximal possible value for "n".
  void setMaxN(int maxn);
  // Returns the maximal value for "n" which was set  through "setMaxN".
  int getMaxN();
  // Returns the number of found Solutions (which were identified through "calculateStudySolutions").
  int getSolutionCount();
  // Returns the ID of the minimax design under the found designs. (only used by the GUI for the R-package).
  int getMinimaxPos();
  // Returns the ID of the optimal design under the found designs. (only used by the GUI for the R-package).
  int getOptimalPos();
  // Approximate the maximal value of "n" for the previously set values of "alpha", "beta", "p1" and "p0" in a way that the optimal design should be found through "calculateStudySolutions".
  long double aproximateMaxN();
  // Calculates the type I error rate for the given design ("r1", "n1", "r", "n"), where "p0" represents the probability of an event under the null hypothesis.
  // r1: Critical value for the first stage (more than "r1" responses needed to proceed to the second stage).
  // n1: Number of patients to be enrolled in the first stage.
  // r: Critical value for the whole trial (more than "r" responses needed at the end of the study to reject the null hypothesis).
  // n: Number of patients to be enrolled in the whole trial.
  // p0: Response probability under the null hypothesis.
  long double calcAlpha(int n1, int r1, int n, int r, long double p0);
  // Calculates the type II error rate for the given design ("r1", "n1", "r", "n").
  // r1: Critical value for the first stage (more than "r1" responses needed to proceed to the second stage).
  // n1: Number of patients to be enrolled in the first stage.
  // r: Critical value for the whole trial (more than "r" responses needed at the end of the study to reject the null hypothesis).
  // n: Number of patients to be enrolled in the whole trial.
  // p1: Response probability under the alternative hypothesis.
  long double calcBeta(int n1, int r1, int n, int r, long double p1);
  // Identifies and internally saves all designs which fit to the previously set values of "alpha", "beta", "p1" and "p0". 
  void calculateStudySolutions();
  // Returns the conditional power for the given design ("r1", "n1", "r", "n") if "rk" out of "k" enrolled patients had an event.
  // r1: Critical value for the first stage (more than "r1" responses needed to proceed to the second stage).
  // n1: Number of patients to be enrolled in the first stage.
  // r: Critical value for the whole trial (more than "r" responses needed at the end of the study to reject the null hypothesis).
  // n: Number of patients to be enrolled in the whole trial.
  // p1: Response probability under the alternative hypothesis.
  double getConditionalPower(int rk, int k, int r1, int n1, int r, int n, double p1);
  // Estimates the effect of (non-)stochastic curtailment for the design with the ID "resID".
  // cut: Sets the "cut point" used to calculate the effect of (non-)stochastic curtailment. A study is stopped if the conditional power falls below the value of "cut".
  // reps: Number of replications used to estimate the effect of (non-)stochastic curtailment.
  // all: If true the effect of (non-)stochastic curtailment will be calculated for different cut points in 0.05 steps starting with the value of the parameter "cut".
  void calculateSC(int resID, double cut, int reps , bool all);
  // Returns an R-representation of the identified designs.
  SEXP getResultsForR();
  // If present returns the estimated effect of (non-)stochastic curtailment for the design with the ID "id".
  SEXP getCurResultForR(int id);
  
  
private:
  long double aproximateMaxNInternal(long double alpha, long double beta, long double p0, long double p1, int n1);
  Result::Curtailment calcSCIntern(int resID, double cut, int reps);//, int seed);
  
  double calculateIntersection(double slope1, double yIntercept1, double slope2, double yIntercept2);
  // Identifies admissible designs based on the approach described in
  // "Kunz C (2011) To-Stage Designs for Phase II Trials with One or Two Endpoints. 
  // PH.D. thesis, Medical Faculty of Heidelberg, Rupecht-Karls-Universitaet" on page 16.
  void setAdmissible(std::vector<Result*>* results);
    
  long double logFact(int n);
  long double bin(int n, int r, long double p);
  long double binsum(int n, int r, long double p);
  
  int miniMaxPos, optimalPos, maxn;
  std::vector<Result*> *allResults;
};

#endif // SIMONDESIGN_H
