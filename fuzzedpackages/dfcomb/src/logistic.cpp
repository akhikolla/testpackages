#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <algorithm>

#include <R.h>
#include <Rcpp.h>

#include <progress.hpp>

#define BOOST_DISABLE_ASSERTS

#include <boost/random/exponential_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/linear_congruential.hpp>

extern "C" {
  #include "arms.h"
  #include "dfcomb.h"
}

using namespace std;

namespace dfcomb { namespace logistic {

enum para_ty {
  PARA_d0,
  PARA_a0,
  PARA_b0,
  PARA_c0
} PARA;
int NDOSE1, NDOSE2;
double TARGET, TARGET_MIN, TARGET_MAX, WEEK_incl, TIMEFULL;
int NCOHORT, COHORT;
bool TITE;
double ESCP, DESP, ARRET, COH_MIN;
int NBURN, NITER;
double CTARG, COVER;
int COH_STOP_EARLY;
int COH_FIN;

boost::random::minstd_rand r;

struct datastru{
  // Doses scales from priors
  vector<double> dose_scale_1, dose_scale_2;

  // Description of included patients
  int pat_incl;
  vector<vector<int> > y, n;
  vector<bool> delta;
  vector<int> dose_adm1, dose_adm2;
  // For TITE only:
  vector<double> time_ev, time_follow, time_min;

  // Current dose for inclusion
  int cdose1, cdose2;

  // Estimated parameters
  double d0, a0, b0, c0;
  vector<vector<double> > pi, ptox, ptox_inf_targ, ptox_targ, ptox_sup_targ;

  datastru(int NDOSE1, int NDOSE2, double prior1[], double prior2[]) :
    dose_scale_1(NDOSE1), dose_scale_2(NDOSE2),
    y(NDOSE1, vector<int>(NDOSE2)), n(NDOSE1, vector<int>(NDOSE2)),
    pi(NDOSE1, vector<double>(NDOSE2)), ptox(NDOSE1, vector<double>(NDOSE2)),
    ptox_inf_targ(NDOSE1, vector<double>(NDOSE2)), ptox_targ(NDOSE1, vector<double>(NDOSE2)), ptox_sup_targ(NDOSE1, vector<double>(NDOSE2))
  {
    for(int i=0; i<NDOSE1; i++) this->dose_scale_1[i] = log(prior1[i]/(1-prior1[i]));
    for(int i=0; i<NDOSE2; i++) this->dose_scale_2[i] = log(prior2[i]/(1-prior2[i]));
  }

  void new_trial() {
    this->pat_incl = 0;
    for(int i=0; i<NDOSE1; i++)
      for(int j=0; j<NDOSE2; j++){
        this->y[i][j]=0;
        this->n[i][j]=0;
      }
    this->delta.clear();
    this->dose_adm1.clear();
    this->dose_adm2.clear();
    if(TITE) {
      this->time_ev.clear();
      this->time_follow.clear();
      this->time_min.clear();
    }
    this->cdose1 = 0;
    this->cdose2 = 0;
  }

  void update_tite() {
    this->delta.resize(this->pat_incl);
    this->time_min.resize(this->pat_incl);

    for(int i=0; i<NDOSE1; i++)
      for(int j=0; j<NDOSE2; j++)
        this->y[i][j] = 0;

    for(int i=0; i<this->pat_incl; i++){
      this->delta[i] = this->time_ev[i] <= min(this->time_follow[i], TIMEFULL);
      this->time_min[i] = min(this->time_ev[i], min(this->time_follow[i], TIMEFULL));
      this->y[this->dose_adm1[i]][this->dose_adm2[i]] += (int)this->delta[i];
    }
  }
};

// Cohort inclusion
void genpopn(datastru * datapt, const vector<vector<double> >& piV){
  if(TITE) {
    double lambda = -log(1-piV[datapt->cdose1][datapt->cdose2])/TIMEFULL;
    boost::random::exponential_distribution<double> exp_rng;
    for(int i=0; i<COHORT; i++){
      datapt->dose_adm1.push_back(datapt->cdose1);
      datapt->dose_adm2.push_back(datapt->cdose2);
      datapt->time_ev.push_back(exp_rng(r, lambda));

      double time_btw_incl = exp_rng(r, WEEK_incl);
      datapt->time_follow.push_back(0);
      for(int j = 0; j <= datapt->pat_incl+i; j++)
        datapt->time_follow[j] += time_btw_incl;
    }

    datapt->n[datapt->cdose1][datapt->cdose2] += COHORT;
    datapt->pat_incl += COHORT;
    datapt->update_tite();
  } else {
    double pi_ij = piV[datapt->cdose1][datapt->cdose2];

    boost::random::uniform_real_distribution<double> uni_rng;
    for(int i=0; i<COHORT; i++){
      datapt->dose_adm1.push_back(datapt->cdose1);
      datapt->dose_adm2.push_back(datapt->cdose2);
      bool tox = uni_rng(r)<pi_ij;
      datapt->y[datapt->cdose1][datapt->cdose2] += (int)tox;
      datapt->delta.push_back(tox);
    }

    datapt->n[datapt->cdose1][datapt->cdose2] += COHORT;
    datapt->pat_incl += COHORT;
  }
}

// Start-up phase
// Pas de startup
void startup0(datastru *datapt, const vector<vector<double> >& piV){
}


// On augmente chaque dose en même temps
void startup1(datastru *datapt, const vector<vector<double> >& piV){
  while(true) {
    genpopn(datapt, piV);

    if(datapt->cdose1==NDOSE1-1 && datapt->cdose2==NDOSE2-1)
      break;

    if(datapt->y[datapt->cdose1][datapt->cdose2])
      break;

    if(datapt->cdose1<NDOSE1-1) datapt->cdose1++;
    if(datapt->cdose2<NDOSE2-1) datapt->cdose2++;
  }
}

// On fait les (0, x) puis les (x, 0)
void startup2(datastru *datapt, const vector<vector<double> >& piV){
  while(true) {
    genpopn(datapt, piV);
    if(datapt->cdose1 == NDOSE1-1)
      break;
    if(datapt->y[datapt->cdose1][datapt->cdose2])
      break;
    datapt->cdose1++;
  }
  if(NDOSE2 > 1 && !datapt->y[0][0]) {
    datapt->cdose1 = 0;
    datapt->cdose2 = 1;
    while(true) {
      genpopn(datapt, piV);
      if(datapt->cdose2 == NDOSE2-1)
        break;
      if(datapt->y[datapt->cdose1][datapt->cdose2])
        break;
      datapt->cdose2++;
    }
  }
}

// On augmente une dose, puis l'autre alternativement
void startup3(datastru *datapt, const vector<vector<double> >& piV){
  bool stat = true;
  while(true){
    genpopn(datapt, piV);

    if(datapt->cdose1==NDOSE1-1 && datapt->cdose2==NDOSE2-1)
      break;

    if(datapt->y[datapt->cdose1][datapt->cdose2])
      break;

    if(stat) {
      if(datapt->cdose1<NDOSE1-1) datapt->cdose1++;
      else datapt->cdose2++;
    } else {
      if(datapt->cdose2<NDOSE2-1) datapt->cdose2++;
      else datapt->cdose1++;
    }
    stat = !stat;
  }
}

// Toxicity probability model
double proba_tox(double x_p, double y_q, double d0, double a0, double b0, double c0){
  double x = d0+a0*x_p+b0*y_q+c0*x_p*y_q;
  return( exp(x)/(1+exp(x)) );
}

double density(double x,  void *datapt){
  struct datastru *dp;
  dp = (datastru *) datapt;
  double d0, a0, b0, c0;
  d0=dp->d0; a0=dp->a0; b0=dp->b0; c0=dp->c0;

  double logprior;
  switch(PARA) {
    case PARA_d0:
      d0=x;
      logprior=-0.05*d0*d0;
      break;
    case PARA_a0:
      a0=x;
      logprior=-a0;
      break;
    case PARA_b0:
      b0=x;
      logprior=-b0;
      break;
    case PARA_c0:
      c0=x;
      logprior=-0.05*c0*c0;
      break;
    default:
      throw std::logic_error("Internal error: invalid PARA.");
  }

  double loglike=0;

  if(TITE) {
    int card2=0;
    for(int j=0; j<dp->pat_incl; j++)
      if(dp->time_follow[j] >= TIMEFULL && dp->delta[j])
        card2 += 1;

    for(int p=0; p<dp->pat_incl; p++){
      double p_tox = proba_tox(dp->dose_scale_1[dp->dose_adm1[p]],
                               dp->dose_scale_2[dp->dose_adm2[p]], d0, a0, b0, c0);
      double weight;
      if(dp->delta[p])
        weight = 1;
      else{
        int card1=0;
        for(int j=0; j<dp->pat_incl; j++)
          if(dp->time_follow[j] >= TIMEFULL &&
             dp->time_ev[j] <= min(dp->time_follow[p],TIMEFULL))
            card1 += 1;

        double w_ini = dp->time_min[p] / TIMEFULL;
        weight = (double)(card1+w_ini)/(card2+1);
        if(weight < 0 || weight > 1)
          throw std::logic_error("Internal error: invalid weight.");
      }
      loglike += dp->delta[p] ? log(weight*p_tox) : log(1-weight*p_tox);
    }
  } else
    for(int i=0; i<NDOSE1; i++)
      for(int j=0; j<NDOSE2; j++)
        if(dp->n[i][j] != 0){
          double pi_ij = proba_tox(dp->dose_scale_1[i], dp->dose_scale_2[j], d0, a0, b0, c0);
          loglike += dp->y[i][j]*log(pi_ij) + (dp->n[i][j]-dp->y[i][j])*log(1-pi_ij);
        }

  return(loglike+logprior);
}

boost::random::minstd_rand r_estim;
boost::random::uniform_real_distribution<double> uni_rng_estim;

double u_random() {
  return uni_rng_estim(r_estim);
}

// Estimations
void estimation(datastru * datapt){
  r_estim.seed(53425);
  uni_rng_estim.reset();

  datapt->d0=1, datapt->a0=1, datapt->b0=1, datapt->c0=0;

  // paramters for ARMS
  int ninit=4, dometrop=1;
  double d0L=-8, d0R=8, d0prev=1, d0prev2=1;
  double a0L=0.01, a0R=8, a0prev=1, a0prev2=1;
  double b0L=0.01, b0R=8, b0prev=1, b0prev2=1;
  double c0L=-8, c0R=8, c0prev=0, c0prev2=0;

  for(int i=0; i<NDOSE1; i++){
    for(int j=0; j<NDOSE2; j++){
      datapt->pi[i][j]=0;
      datapt->ptox[i][j]=0;
      datapt->ptox_inf_targ[i][j]=0;
      datapt->ptox_targ[i][j]=0;
      datapt->ptox_sup_targ[i][j]=0;
    }
  }

  for(int iter=0; iter<NBURN+NITER; iter++){
    int err=0;
    double d0samp;
    double a0samp;
    double b0samp;
    double c0samp;

    while(true){
      PARA = PARA_d0;
      err = arms_simple(ninit, &d0L, &d0R, density, datapt, dometrop, &d0prev, &d0samp, u_random);
      if(err>0) {
        char buf[100]; sprintf(buf, "%d", err);
        throw std::logic_error(string("arms internal error (d0): ") + string(buf));
      }
      datapt->d0 = d0samp;
      d0prev2 = d0prev;
      d0prev = d0samp;

      PARA = PARA_a0;
      err = arms_simple(ninit, &a0L, &a0R, density, datapt, dometrop, &a0prev, &a0samp, u_random);
      if(err>0) {
        char buf[100]; sprintf(buf, "%d", err);
        throw std::logic_error(string("arms internal error (a0): ") + string(buf));
      }
      datapt->a0 = a0samp;
      a0prev2 = a0prev;
      a0prev =  a0samp;

      PARA = PARA_b0;
      err = arms_simple(ninit, &b0L, &b0R, density, datapt, dometrop, &b0prev, &b0samp, u_random);
      if(err>0) {
        char buf[100]; sprintf(buf, "%d", err);
        throw std::logic_error(string("arms internal error (b0): ") + string(buf));
      }
      datapt->b0 = b0samp;
      b0prev2 = b0prev;
      b0prev = b0samp;

      PARA = PARA_c0;
      err = arms_simple(ninit, &c0L, &c0R, density, datapt, dometrop, &c0prev, &c0samp, u_random);
      if(err>0) {
        char buf[100]; sprintf(buf, "%d", err);
        throw std::logic_error(string("arms internal error (c0): ") + string(buf));
      }
      datapt->c0 = c0samp;
      c0prev2 = c0prev;
      c0prev = c0samp;

      if(a0samp+c0samp*(datapt->dose_scale_2[0]) >= 0 &&
         a0samp+c0samp*(datapt->dose_scale_2[NDOSE2-1]) >= 0 &&
         b0samp+c0samp*(datapt->dose_scale_1[0]) >= 0 &&
         b0samp+c0samp*(datapt->dose_scale_1[NDOSE1-1]) >= 0)
        break;

      datapt->d0 = d0prev2;
      d0prev = d0prev2;
      datapt->a0 = a0prev2;
      a0prev = a0prev2;
      datapt->b0 = b0prev2;
      b0prev = b0prev2;
      datapt->c0 = c0prev2;
      c0prev = c0prev2;
    }

    if(iter>=NBURN){
      for(int i=0; i<NDOSE1; i++){
        for(int j=0; j<NDOSE2; j++){
          double pi_ij;
          pi_ij=proba_tox(datapt->dose_scale_1[i], datapt->dose_scale_2[j], d0samp, a0samp, b0samp, c0samp);
          datapt->pi[i][j] += pi_ij;
          datapt->ptox[i][j] += (pi_ij<TARGET)?1:0;
          datapt->ptox_inf_targ[i][j] += (pi_ij<TARGET_MIN)?1:0;
          datapt->ptox_targ[i][j] += (pi_ij>=TARGET_MIN && pi_ij<=TARGET_MAX)?1:0;
          datapt->ptox_sup_targ[i][j] += (pi_ij>TARGET_MAX)?1:0;
        }
      }
    }
  }
  for(int i=0; i<NDOSE1; i++){
    for(int j=0; j<NDOSE2; j++){
      datapt->pi[i][j] /= NITER;
      datapt->ptox[i][j] /= NITER;
      datapt->ptox_inf_targ[i][j] /= NITER;
      datapt->ptox_targ[i][j] /= NITER;
      datapt->ptox_sup_targ[i][j] /= NITER;
    }
  }
}


void take_if_better(datastru * datapt, int& recomdose1, int& recomdose2, int candidate_dose1, int candidate_dose2) {
  if(recomdose1 == -1 && recomdose2 == -1) {
    recomdose1 = candidate_dose1;
    recomdose2 = candidate_dose2;
    return;
  }

  double candidateTarg = fabs(datapt->pi[candidate_dose1][candidate_dose2]-TARGET);
  double bestTarg = fabs(datapt->pi[recomdose1][recomdose2]-TARGET);

  if(candidateTarg < bestTarg) {
    recomdose1 = candidate_dose1;
    recomdose2 = candidate_dose2;
  }
}

// Stop for over/underdosing
bool over_under_stop(datastru * datapt) {
  if(datapt->cdose1 == 0 && datapt->cdose2 == 0 &&
     datapt->n[datapt->cdose1][datapt->cdose2] >= COH_MIN*COHORT &&
     1-datapt->ptox[datapt->cdose1][datapt->cdose2] >= ARRET)
    return true;

  if(datapt->cdose1 == NDOSE1-1 && datapt->cdose2 == NDOSE2-1 &&
     datapt->n[datapt->cdose1][datapt->cdose2] >= COH_MIN*COHORT &&
     datapt->ptox[datapt->cdose1][datapt->cdose2] >= ARRET)
    return true;

  return false;
}

// Determination of the next combination
void alloc_rule1(datastru * datapt){
  int recomdose1=-1, recomdose2=-1;
  int cdose1 = datapt->cdose1;
  int cdose2 = datapt->cdose2;

  if(datapt->ptox[cdose1][cdose2] > ESCP){
    int d1[4] = { +1, +1, 0, -1 };
    int d2[4] = { -1, 0, +1, +1 };
    for(int i = 0; i < 4; i++) {
      int candidate_dose1 = cdose1 + d1[i];
      int candidate_dose2 = cdose2 + d2[i];
      if(candidate_dose1 < 0 || candidate_dose2 < 0 ||
         candidate_dose1 >= NDOSE1 || candidate_dose2 >= NDOSE2)
        continue;
      if(datapt->pi[candidate_dose1][candidate_dose2] <= datapt->pi[cdose1][cdose2])
        continue;
      take_if_better(datapt, recomdose1, recomdose2, candidate_dose1, candidate_dose2);
    }
    if(recomdose1==-1 && recomdose2==-1){
      recomdose1 = cdose1;
      recomdose2 = cdose2;
    }
  }
  else if(datapt->ptox[cdose1][cdose2] < DESP){
    int d1[4] = { +1, 0, -1, -1 };
    int d2[4] = { -1, -1, 0, +1 };
    for(int i = 0; i < 4; i++) {
      int candidate_dose1 = cdose1 + d1[i];
      int candidate_dose2 = cdose2 + d2[i];
      if(candidate_dose1 < 0 || candidate_dose2 < 0 ||
         candidate_dose1 >= NDOSE1 || candidate_dose2 >= NDOSE2)
        continue;
      if(datapt->pi[candidate_dose1][candidate_dose2] >= datapt->pi[cdose1][cdose2])
        continue;
      take_if_better(datapt, recomdose1, recomdose2, candidate_dose1, candidate_dose2);
    }
    if(recomdose1==-1 && recomdose2==-1) {
      recomdose1 = cdose1;
      recomdose2 = cdose2;
    }
  }
  else{
    recomdose1 = cdose1;
    recomdose2 = cdose2;
  }

  datapt->cdose1 = recomdose1;
  datapt->cdose2 = recomdose2;
}

void alloc_rule2(datastru * datapt){
  int recom1 = -1, recom2 = -1;
  double max_ptox_targ=0.0;
  int d1[8] = { -1, -1,  0,  0,  0, +1, +1, +1 };
  int d2[8] = { 0, +1,  -1,  0, +1, -1,  0, +1 };
  for(int dose1=0; dose1<NDOSE1; dose1++)
    for(int dose2=0; dose2<NDOSE2; dose2++) {
      bool admissible = false;
      for(int i = 0; i < 8 && !admissible; i++)
        if(dose1+d1[i] >= 0 && dose2+d2[i] >= 0 &&
           dose1+d1[i] < NDOSE1 && dose2+d2[i] < NDOSE2 &&
           datapt->n[dose1+d1[i]][dose2+d2[i]] > 0)
          admissible = true;
      if(!admissible || datapt->ptox_sup_targ[dose1][dose2] >= COVER)
        continue;
      double x = datapt->ptox_targ[dose1][dose2];
      if(x >= max_ptox_targ){
        max_ptox_targ = x;
        recom1=dose1;
        recom2=dose2;
      }
    }

  if(recom1 == -1 || recom2 == -1)
    recom1 = recom2 = 0;
  datapt->cdose1 = recom1;
  datapt->cdose2 = recom2;
}

void alloc_rule3(datastru * datapt){
  int recom1 = -1, recom2 = -1;
  double max_ptox_targ=0.0;
  int d1[8] = { -1, -1, -1,  0,  0,  0, +1, +1 };
  int d2[8] = { -1,  0, +1, -1,  0, +1, -1,  0 };
  for(int i = 0; i < 8; i++) {
    int candidate_dose1 = datapt->cdose1 + d1[i];
    int candidate_dose2 = datapt->cdose2 + d2[i];
    if(candidate_dose1 < 0 || candidate_dose2 < 0 ||
       candidate_dose1 >= NDOSE1 || candidate_dose2 >= NDOSE2)
      continue;
    if(datapt->ptox_sup_targ[candidate_dose1][candidate_dose2] >= COVER)
      continue;
    double x = datapt->ptox_targ[candidate_dose1][candidate_dose2];
    if(x >= max_ptox_targ){
      max_ptox_targ = x;
      recom1=candidate_dose1;
      recom2=candidate_dose2;
    }
  }

  if(recom1 == -1 || recom2 == -1) {
    recom1 = max(0, datapt->cdose1 - 1);
    recom2 = max(0, datapt->cdose2 - 1);
  }
  datapt->cdose1 = recom1;
  datapt->cdose2 = recom2;
}

void alloc_rule(datastru * datapt, int id) {
  switch(id) {
    case 1: alloc_rule1(datapt); break;
    case 2: alloc_rule2(datapt); break;
    case 3: alloc_rule3(datapt); break;
    default: throw std::logic_error("Unknown alloc rule ID.");
  }
}

// No stop rule : we go to NCOHORT cohorts (max sample size)
bool early_finding_rule1(datastru * datapt){
  return false;
}

// At least 2 cohorts at recommended dose and ptox_targ >= CTARG and ptox_sup_targ < COVER
bool early_finding_rule2(datastru * datapt){
  return datapt->n[datapt->cdose1][datapt->cdose2] >= COH_STOP_EARLY * COHORT &&
         datapt->ptox_targ[datapt->cdose1][datapt->cdose2] >= CTARG &&
         datapt->ptox_sup_targ[datapt->cdose1][datapt->cdose2] < COVER;
}

// At least 3 cohorts at recommended dose
bool early_finding_rule3(datastru * datapt){
  return datapt->n[datapt->cdose1][datapt->cdose2] >= COH_STOP_EARLY * COHORT;
}

bool early_finding_rule(datastru * datapt, int id){
  switch(id) {
    case 1: return early_finding_rule1(datapt);
    case 2: return early_finding_rule2(datapt);
    case 3: return early_finding_rule3(datapt);
    default: throw std::logic_error("Unknown early finding rule ID.");
  }
}

void final_recom(datastru * datapt) {
  int recom1 = -1, recom2 = -1;
  double max_ptox_targ=0.0;
  for(int i=0; i<NDOSE1; i++)
    for(int j=0; j<NDOSE2; j++)
      if(datapt->n[i][j] >= COHORT * COH_FIN){
        double x = datapt->ptox_targ[i][j];
        if(x >= max_ptox_targ){
          max_ptox_targ = x;
          recom1=i;
          recom2=j;
        }
      }
  if(recom1 == -1 || recom2 == -1)
    throw std::logic_error("Internal error: no recommended dose.");
  datapt->cdose1 = recom1;
  datapt->cdose2 = recom2;
}

}}

using namespace dfcomb::logistic;

R_NativePrimitiveArgType logistic_sim_args[] =
  {LGLSXP,
   INTSXP, INTSXP,
   REALSXP, REALSXP,
   REALSXP,
   REALSXP, REALSXP, REALSXP,
   REALSXP, REALSXP,
   INTSXP, INTSXP, INTSXP,
   REALSXP, REALSXP, REALSXP,
   REALSXP, REALSXP,
   INTSXP, INTSXP, INTSXP,
   INTSXP,
   INTSXP, INTSXP, INTSXP,
   INTSXP, INTSXP,

   REALSXP, REALSXP, REALSXP,
   REALSXP, REALSXP,
   REALSXP};
const int logistic_sim_nargs = 34;

void logistic_sim(int* tite,
                  int* ndose1, int* ndose2,
                  double* timefull, double* week_incl,
                  double* piv,
                  double* target, double* target_max, double* target_min,
                  double* prior_1, double* prior_2,
                  int* ncohort, int* cohort, int* ntrial,
                  double* escp, double* desp, double* arret,
                  double* ctarg, double* cover,
                  int* coh_min, int* coh_stop_early, int* coh_fin,
                  int* seed,
                  int* startup, int* alloc_rule_id, int* early_finding_rule_id,
                  int* nburn, int* niter,

                  double* nrecdoserat, double* nptsdoserat, double* ntoxrat,
                  double* inconc_rat, double* early_finding_rat,
                  double* n_treated_tab)
{
  try {

  TITE = *tite;
  if(TITE) {
    TIMEFULL = *timefull;
    WEEK_incl = *week_incl;
  }

  ESCP = *escp;
  DESP = *desp;
  ARRET = *arret;
  COH_MIN = *coh_min;
  r.seed(*seed);
  NBURN = *nburn;
  NITER = *niter;
  CTARG = *ctarg;
  COVER = *cover;
  COH_STOP_EARLY = *coh_stop_early;
  COH_FIN = *coh_fin;

  NDOSE1 = *ndose1;
  NDOSE2 = *ndose2;

  struct datastru data(NDOSE1, NDOSE2, prior_1, prior_2);

  vector<vector<double> > piV(NDOSE1, vector<double>(NDOSE2));

  for(int i=0; i<NDOSE1; i++)
    for(int j=0; j<NDOSE2; j++)
      piV[i][j] = piv[i + j*NDOSE1];

  TARGET = *target;
  TARGET_MIN = *target_min;
  TARGET_MAX = *target_max;
  NCOHORT = *ncohort;
  COHORT = *cohort;

  int inconc=0;
  int early_finding=0;
  vector<vector<int> > nptsdose(NDOSE1, vector<int>(NDOSE2, 0)),
    ntox(NDOSE1, vector<int>(NDOSE2, 0)), nrecdose(NDOSE1, vector<int>(NDOSE2, 0));

  // Trials simulations
  Progress prog(*ntrial);
  for(int trial=0; trial<*ntrial; trial++) {
    data.new_trial();

    // Start-up phase
    switch(*startup) {
      case 0:
        startup0(&data, piV);
        break;
      case 1:
        startup1(&data, piV);
        break;
      case 2:
        startup2(&data, piV);
        break;
      case 3:
        startup3(&data, piV);
        break;
      default:
        throw std::logic_error("Unknown startup ID.");
    }

    // Cohort inclusions and estimations, based on the model
    while(true){
      if(prog.check_abort()) {
        *ntrial = trial;
        goto aborted;
      }

      estimation(&data);

      if(over_under_stop(&data)) {
        inconc++;
        break;
      }

      alloc_rule(&data, *alloc_rule_id);

      if(early_finding_rule(&data, *early_finding_rule_id)) {
        early_finding++;
        nrecdose[data.cdose1][data.cdose2]++;
        break;
      }

      if(data.pat_incl >= COHORT * NCOHORT) {
        // FIXME : we should perhaps do that for early finding or inconclusive trials
        if(TITE) {
          for(int i=0; i<data.pat_incl; i++)
            data.time_follow[i] = INFINITY;
          data.update_tite();
          estimation(&data);
        }

        final_recom(&data);
        nrecdose[data.cdose1][data.cdose2]++;
        break;
      }

      genpopn(&data, piV);
    }

    n_treated_tab[trial] = data.pat_incl;
    for(int i=0; i<NDOSE1; i++)
      for(int j=0; j<NDOSE2; j++){
        ntox[i][j] += data.y[i][j];
        nptsdose[i][j] += data.n[i][j];
      }


    prog.increment();
  }

  aborted:

  for(int i=0; i<NDOSE1; i++)
    for(int j=0; j<NDOSE2; j++){
      nrecdoserat[i + j*NDOSE1] = (double)nrecdose[i][j] / *ntrial;
      nptsdoserat[i + j*NDOSE1] = (double)nptsdose[i][j] / *ntrial;
      ntoxrat[i + j*NDOSE1] = (double)ntox[i][j] / *ntrial;
    }

  *inconc_rat = (double)inconc / *ntrial;
  *early_finding_rat = (double)early_finding / *ntrial;
  }
  catch (std::logic_error &e) { error("Internal error in dfcomb (details: %s)", e.what()); }
  catch (...) { error("Internal error in dfcomb"); }

  return;
}

R_NativePrimitiveArgType logistic_next_args[] =
  {LGLSXP,
   INTSXP, INTSXP,
   REALSXP,
   REALSXP, REALSXP, REALSXP,
   REALSXP, REALSXP,
   INTSXP, LGLSXP,
   REALSXP, REALSXP, REALSXP,
   REALSXP, REALSXP,
   INTSXP, INTSXP, INTSXP,
   INTSXP, INTSXP,
   INTSXP, INTSXP,

   INTSXP,
   INTSXP, INTSXP,
   INTSXP, INTSXP,
   REALSXP, REALSXP,
   LGLSXP,

   LGLSXP, LGLSXP,
   REALSXP, REALSXP, REALSXP,
   REALSXP, REALSXP };
const int logistic_next_nargs = 38;

void logistic_next(int* tite,
                   int* ndose1, int* ndose2,
                   double* timefull,
                   double* target, double* target_max, double* target_min,
                   double* prior_1, double* prior_2,
                   int* cohort, int* trial_end,
                   double* escp, double* desp, double* arret,
                   double* ctarg, double* cover,
                   int* coh_min, int* coh_stop_early,  int* coh_fin,
                   int* early_finding_rule_id, int* alloc_rule_id,
                   int* nburn, int* niter,

                   int* pat_incl,
                   int* cdose1, int* cdose2,
                   int* dose_adm1, int* dose_adm2,
                   double* time_ev, double* time_follow /* Only for tite */,
                   int* delta /* Only for non-tite */,

                   int* inconc, int* early_finding,
                   double* pi, double* ptox, double* ptox_inf_targ,
                   double* ptox_targ, double* ptox_sup_targ)
{
  try {

  TITE = *tite;
  if(TITE) TIMEFULL = *timefull;

  ESCP = *escp;
  DESP = *desp;
  ARRET = *arret;
  COH_MIN = *coh_min;
  NBURN = *nburn;
  NITER = *niter;
  CTARG = *ctarg;
  COVER = *cover;
  COH_STOP_EARLY = *coh_stop_early;
  COH_FIN = *coh_fin;

  NDOSE1 = *ndose1;
  NDOSE2 = *ndose2;

  *inconc = 0;
  *early_finding = 0;

  struct datastru data(NDOSE1, NDOSE2, prior_1, prior_2);

  TARGET = *target;
  TARGET_MIN = *target_min;
  TARGET_MAX = *target_max;
  COHORT = *cohort;
  NCOHORT = *trial_end ? *pat_incl : *pat_incl + 1;

  data.pat_incl = *pat_incl;
  data.cdose1 = *cdose1;
  data.cdose2 = *cdose2;
  for(int i=0; i < *pat_incl; i++) {
    data.dose_adm1.push_back(dose_adm1[i]);
    data.dose_adm2.push_back(dose_adm2[i]);
    data.n[dose_adm1[i]][dose_adm2[i]]++;
    if(TITE) {
      data.time_ev.push_back(time_ev[i]);
      if(*trial_end && time_follow[i] < TIMEFULL)
        error("dfcomb : the final recommendation cannot be computed when "
              "all the patients have not been fully followed");
      data.time_follow.push_back(time_follow[i]);
    } else {
      data.delta.push_back(delta[i]);
      data.y[dose_adm1[i]][dose_adm2[i]] += (int)data.delta[i];
    }
  }
  if(TITE) data.update_tite();

  estimation(&data);

  for(int i=0; i<NDOSE1; i++)
    for(int j=0; j<NDOSE2; j++){
      pi[i + j*NDOSE1] = data.pi[i][j];
      ptox[i + j*NDOSE1] = data.ptox[i][j];
      ptox_inf_targ[i + j*NDOSE1] = data.ptox_inf_targ[i][j];
      ptox_targ[i + j*NDOSE1] = data.ptox_targ[i][j];
      ptox_sup_targ[i + j*NDOSE1] = data.ptox_sup_targ[i][j];
    }

  if(over_under_stop(&data)) {
    *inconc = 1;
    return;
  }

  alloc_rule(&data, *alloc_rule_id);

  if(early_finding_rule(&data, *early_finding_rule_id))
    *early_finding = 1;
  else if(*trial_end)
    final_recom(&data);

  *cdose1 = data.cdose1;
  *cdose2 = data.cdose2;
  }
  catch (std::logic_error &e) { error("Internal error in dfcomb (details: %s)", e.what()); }
  catch (...) { error("Internal error in dfcomb"); }
}
