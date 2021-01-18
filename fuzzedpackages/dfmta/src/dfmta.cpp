#include<R.h>

#define BOOST_DISABLE_ASSERTS
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <RcppArmadillo.h>

#undef dnorm

#include <cmath>
#include <string>
#include <functional>

#include <cppbugs/cppbugs.hpp>

#include <boost/random/exponential_distribution.hpp>
#include <cstdlib>

#ifdef _OPENMP
#include<omp.h>
#endif

#include <progress.hpp>

using namespace std;
using namespace arma;
using namespace cppbugs;

namespace dfmta {

double TARG_SUP, EFF_MIN;
bool HAS_TIME;
double TIMEFULL, CYCLE;
int COHORT_START, COHORT;

/* True toxicity and response rates from input */
struct true_data {
  vector<double> piV;
  vector<vector<double>> respV;
  vector<double> incl_per_week;
};

/* A struct representing the parameters of the toxicity model */
struct toxicity_parameters {
  double beta0, beta1;

  toxicity_parameters(double beta0, double beta1): beta0(beta0), beta1(beta1) {}
  toxicity_parameters(): beta0(0), beta1(0) {}

  template<typename U>
  void proba_tox(const U& dose_tox, U& out) {
    out = 1-1/(1+cppbugs::exp_approx(beta0 + beta1*dose_tox));
  }
};

/* A struct representing the parameters of the efficacy model */
struct efficacy_parameters {
  double gamma0, gamma1;
  int tau;

  efficacy_parameters(double gamma0, double gamma1, int tau):
    gamma0(gamma0), gamma1(gamma1), tau(tau) {}
  efficacy_parameters():
    gamma0(0), gamma1(0), tau(0) {}

  template<typename U>
  void responseRate(const vector<U>& dose_eff_tau, U& out) const {
    unsigned tau = min((int)dose_eff_tau.size()-1, max(0, this->tau));
    out = 1-1/(1+cppbugs::exp_approx(gamma0 + gamma1*dose_eff_tau[tau]));
  }
};

struct estimations {
  /* Estimated parameters */
  toxicity_parameters tox_params;
  efficacy_parameters eff_params;

  /* Estimated information on doses */
  vector<double> pi, ptox_inf;
  vector<double> resp, qeff_inf, resp2;
  vector<double> proba_tau;

  estimations(int ndose):
    pi(ndose), ptox_inf(ndose), resp(ndose),
    qeff_inf(ndose), resp2(ndose), proba_tau(ndose)
  { }
};


struct trial_data;

typedef estimations (*estimator)(const trial_data& trial_data, uint_fast64_t seed, int group, bool final);

/* All the data needed dring a trial */
struct trial_data{
  estimator est;

  /* Transformation of dose levels to a more appropriate scale for estimation */
  vector<double> doseT;
  vector<vector<double>> doseE;

  /* State information for inclusion process */
  vector<int> cdose; /* Current dose. -1 if stopped because of too high toxivity */
  vector<int> startup_end;
  double time_cur;
  int pat_incl;
  vector<int> pat_incl_group;
  vector<unsigned> dose_adm;
  vector<int> group;
  vector<double> time_eff, time_incl;
  vector<int> efficacy, toxicity;

  double s_2;

  vector<double> s_1;

  std::mt19937_64 r;

  trial_data(estimator est, const vector<double>& doseTV, const vector<vector<double>>& doseEV,
             int ngroups, double s_2, vector<double> s_1, uint_fast64_t seed):
    est(est),
    doseT(doseTV.size()), doseE(doseEV.size(), vector<double>(ngroups)),
    cdose(ngroups, 0), startup_end(ngroups, -1),
    time_cur(0), pat_incl(0),
    pat_incl_group(ngroups, 0),
    s_2(s_2),
    s_1(s_1),
    r(seed)
  {
    for(int i=0; i < (int)doseT.size(); i++) {
      doseT[i] = log(doseTV[i]/(1-doseTV[i]));
      for(int j = 0; j < ngroups; j++)
        doseE[i][j] = log(doseEV[i][j]/(1-doseEV[i][j]));
    }
  }
};

/* Results of a set of trials */
struct results {
  vector<int> inconc;
  vector<vector<int>> n_pat_dose, rec_dose;
  int n_pat_tot;
  vector<vector<int>> n_tox, n_eff;
  int tox_tot, eff_tot;
  vector<int> n_pat_mtd;
  double duration;
  int nb_trials;

  results(int ndose, int ngroups) :
    inconc(ngroups, 0),
    n_pat_dose(ndose, vector<int>(ngroups, 0)), rec_dose(ndose, vector<int>(ngroups, 0)),
    n_pat_tot(0),
    n_tox(ndose, vector<int>(ngroups, 0)),      n_eff(ndose, vector<int>(ngroups, 0)),
    tox_tot(0), eff_tot(0),
    n_pat_mtd(ngroups, 0),
    duration(0), nb_trials(0)
  { }

  void accumulate(const trial_data& trial_data, const vector<int>& recom) {
    for(int group=0; group<recom.size(); group++)
      if(recom[group] == -1)
        inconc[group]++;
      else
        rec_dose[recom[group]][group]++;

    for(int pat=0; pat<trial_data.pat_incl; pat++) {
      int dose_pat = trial_data.dose_adm[pat];
      int group_pat = trial_data.group[pat];
      int tox = trial_data.toxicity[pat];
      int eff = HAS_TIME ? trial_data.time_eff[pat] < TIMEFULL
                         : trial_data.efficacy[pat];
      n_pat_dose[dose_pat][group_pat]++;
      n_pat_tot++;
      n_tox[dose_pat][group_pat] += tox;
      tox_tot += tox;
      n_eff[dose_pat][group_pat] += eff;
      eff_tot += eff;
      if(dose_pat == recom[group_pat])
        n_pat_mtd[group_pat]++;
    }
    if(HAS_TIME)
      duration += trial_data.time_cur;
    nb_trials++;
  }
};

// Wait for the next patient to arrive.
// Rejects it if the trial has been stopped for its group, or if we did not
// spend CYCLE week since the last inclusion in this cohort.
// Returns true if the trial ends.
// Puts in group the group of the next patient.
bool wait_patient(trial_data& trial_data, const true_data& true_data, int& group) {
  boost::exponential_distribution<double> exp_rng;
  double rate = 0;
  const int ngroups = trial_data.cdose.size();
  vector<double> hist(ngroups);
  for(int g = 0; g < ngroups; g++)
    if(trial_data.cdose[g] < 0)
      // We do not chose this group, either because it is too toxic,
      // or because we have included all the patients of this group.
      hist[g] = 0;
    else {
      hist[g] = true_data.incl_per_week[g];
      rate += hist[g];
    }

  // There is no longer any admissible group. The trial is finished.
  if(rate == 0) {
    if(HAS_TIME &&
       trial_data.pat_incl > 0 &&
       trial_data.time_cur -
       trial_data.time_incl[trial_data.pat_incl-1] <= TIMEFULL)
      trial_data.time_cur = TIMEFULL + 0.01 + trial_data.time_incl[trial_data.pat_incl-1];
    return true;
  } else {
    std::discrete_distribution<int> group_rng(hist.begin(), hist.end());

    if(!HAS_TIME)
      group = group_rng(trial_data.r);
    else {
      reject_patient:

      /* How much time do we wait for this patient ? */
      double time_incl = exp_rng(trial_data.r, rate);

      /* We shift everything by this amount of time. */
      trial_data.time_cur += time_incl;

      /* What is the group of the new patient ? */
      group = group_rng(trial_data.r);

      /* Should we reject this patient ? */
      int startup_end = trial_data.startup_end[group];
      if( (startup_end == -1 && trial_data.pat_incl_group[group] % COHORT_START == 0) ||
          (startup_end != -1 && (trial_data.pat_incl_group[group] - startup_end) % COHORT == 0) )
        for(int pat = trial_data.pat_incl-1; pat >= 0; pat--)
          if(trial_data.group[pat] == group) {
            if(trial_data.time_cur - trial_data.time_incl[pat] < CYCLE)
              /* Because we do not already know whether it is toxic for the latest
                 cohort ? */
              goto reject_patient;
            break;
          }
    }
    return false;
  }
}

void include_patient(trial_data& trial_data, const true_data& true_data, int group) {
  trial_data.dose_adm.push_back(trial_data.cdose[group]);
  trial_data.group.push_back(group);

  /* Is there a toxicity ? */
  double pi_ij = true_data.piV[trial_data.cdose[group]];
  std::uniform_real_distribution<double> uni_rng;
  trial_data.toxicity.push_back(uni_rng(trial_data.r) < pi_ij);

  /* When will there be an efficacy ? */
  double resp_ij = true_data.respV[trial_data.cdose[group]][group];
  if(HAS_TIME) {
    boost::exponential_distribution<double> exp_rng;
    double param_exp = -log(1-resp_ij)/TIMEFULL;
    trial_data.time_eff.push_back(exp_rng(trial_data.r, param_exp));

    trial_data.time_incl.push_back(trial_data.time_cur);
  } else
    trial_data.efficacy.push_back(uni_rng(trial_data.r) < resp_ij);

  trial_data.pat_incl++;
  trial_data.pat_incl_group[group]++;
}

template<typename T>
typename T::value_type median(const T& x) {
  if(x.size() == 0)
    return 0;
  std::vector<typename T::value_type> v(x.begin(), x.end());
  std::nth_element(v.begin(), v.begin()+v.size()/2, v.end());
  if(v.size() % 2 == 1)
    return v[v.size()/2];
  else
    return (v[v.size()/2] +
            *std::max_element(v.begin(), v.begin()+v.size()/2))/2.;
}

/* Estimate toxicity and efficacy parameters.
   The efficacy is estimated only for the given group. */
estimations estimate_ra(const trial_data& trial_data, uint_fast64_t seed, int group, bool final) {
  const int ndose = trial_data.doseT.size();

  vector<double> doseE(ndose);
  for(int d = 0; d < ndose; d++)
    doseE[d] = trial_data.doseE[d][group];

  /* Read usefull data from trial_data */
  int NPatientsCycle;
  if(!HAS_TIME)
    NPatientsCycle = trial_data.pat_incl;
  else
    for(NPatientsCycle = 0; NPatientsCycle < trial_data.pat_incl; NPatientsCycle++)
      if(trial_data.time_cur - trial_data.time_incl[NPatientsCycle] < CYCLE)
        break;

  const vec toxicity = conv_to<vec>::from(trial_data.toxicity).head(NPatientsCycle);
  const vec dose_tox = vec(trial_data.doseT)
    .elem(conv_to<uvec>::from(trial_data.dose_adm).head(NPatientsCycle));

  uvec in_group = find(conv_to<ivec>::from(trial_data.group).head(NPatientsCycle) == group);
  vec efficacy, weights;
  if(HAS_TIME) {
    efficacy = vec(in_group.size());
    vec time_min_eff(in_group.size());

    int card2 = 0;
    for(unsigned i = 0; i < in_group.size(); i++) {
      double time_follow = min(TIMEFULL, trial_data.time_cur - trial_data.time_incl[in_group(i)]);
      if(trial_data.time_eff[in_group(i)] < time_follow) {
        efficacy(i) = 1;
        time_min_eff(i) = trial_data.time_eff[in_group(i)];
        if(time_follow == TIMEFULL)
          card2++;
      } else {
        efficacy(i) = 0;
        time_min_eff(i) = time_follow;
      }
    }

    weights = vec(in_group.size());
    for(unsigned i = 0; i < in_group.size(); i++) {
      if(efficacy(i)) {
        weights(i) = 1;
      } else {
        double time_follow = min(TIMEFULL, trial_data.time_cur - trial_data.time_incl[in_group(i)]);
        int card1 = 0;
        for(unsigned j = 0; j < in_group.size(); j++)
          if(trial_data.time_cur - trial_data.time_incl[in_group(j)] >= TIMEFULL &&
             trial_data.time_eff[in_group(j)] <= time_follow)
            card1++;
        double w_ini = time_min_eff(i) / TIMEFULL;
        weights(i) = (card1+w_ini)/(card2+1);
      }
    }
  } else
    efficacy = conv_to<vec>::from(trial_data.efficacy).elem(in_group);

  const uvec dose_adm_eff = conv_to<uvec>::from(trial_data.dose_adm).elem(in_group);

  vector<vec> dose_eff_tau;
  for(int tau = 0; tau < ndose; tau++)
    dose_eff_tau.push_back(vec(doseE).elem(clamp(dose_adm_eff, 0, tau)));

  /* MCMC initalization and sampling */
  efficacy_parameters eff_params(0, 0.001, ndose-1);
  toxicity_parameters tox_params(0, 0.001);

  vec pi_jk, resp_jk;

  std::function<void ()> model = [&]() {
    if(HAS_TIME) {
      eff_params.responseRate(dose_eff_tau, resp_jk);
      resp_jk = weights % resp_jk;
    } else
      eff_params.responseRate(dose_eff_tau, resp_jk);
    tox_params.proba_tox(dose_tox, pi_jk);
  };
  MCModel<std::mt19937> m(model);

  m.track<Normal>(eff_params.gamma0).dnorm(0, 0.01);
  m.track<Exponential>(eff_params.gamma1).dexp(1);
  vec distr(ndose);
  for(int i = 0; i < ndose; i++)
    distr[i] = 1.0;
  m.track<Discrete>(eff_params.tau).ddiscr(distr);
  m.track<ObservedBernoulli>(efficacy).dbern(resp_jk);

  m.track<ObservedBernoulli>(toxicity).dbern(pi_jk);
  m.track<Normal>(tox_params.beta0).dnorm(0, 0.01);
  m.track<Exponential>(tox_params.beta1).dexp(1);

  m.sample(1e6, 1e5, 1e4, 10);

  /* Get list of eff_params and tox_params */
  vector<pair<efficacy_parameters, toxicity_parameters>> params_draws_list;
  auto itbeta0 = m.getNode(tox_params.beta0).history.begin(),
    itbeta0end = m.getNode(tox_params.beta0).history.end();
  auto itbeta1 = m.getNode(tox_params.beta1).history.begin();
  auto itgamma0 = m.getNode(eff_params.gamma0).history.begin();
  auto itgamma1 = m.getNode(eff_params.gamma1).history.begin();
  auto ittau = m.getNode(eff_params.tau).history.begin();

  while(itbeta0 != itbeta0end) {
    params_draws_list.emplace_back(
      efficacy_parameters(*itgamma0, *itgamma1, *ittau),
      toxicity_parameters(*itbeta0, *itbeta1));
    itbeta0++; itbeta1++;
    itgamma0++; itgamma1++; ittau++;
  }

  int n_samp = params_draws_list.size();

  /* Estimations */
  estimations result(ndose);
  result.tox_params.beta0 = median(m.getNode(tox_params.beta0).history);
  result.tox_params.beta1 = median(m.getNode(tox_params.beta1).history);

  for(int d = 0; d < ndose; d++) {
    int count_tox = 0;
    std::vector <double> pi_median;
    for(auto draw: params_draws_list){
      double proba_tox;
      draw.second.proba_tox(trial_data.doseT[d], proba_tox);
      count_tox += proba_tox < TARG_SUP;
      pi_median.push_back(proba_tox);
    }
    result.ptox_inf[d] = double(count_tox)/n_samp;
    result.pi[d] = median(pi_median);
  }

  vector<int> counts(ndose, 0);
  for(int tau: m.getNode(eff_params.tau).history) {
    counts[tau]++;
  }
  for(int d = 0; d < ndose; d++) {
    result.proba_tau[d] = (double)counts[d]/n_samp;
  }

  if(final) {
    result.eff_params.tau = max_element(counts.begin(), counts.end())-counts.begin();
  } else {
    vector<int> tau_high_prob;
    vector<double> prob_tau_high;
    double max_val = *max_element(result.proba_tau.begin(), result.proba_tau.end());
    for(int d=0; d<ndose; d++){
      if(max_val-result.proba_tau[d] <= trial_data.s_1[trial_data.pat_incl_group[group]]){
        tau_high_prob.push_back(d);
        prob_tau_high.push_back(result.proba_tau[d]);
      }
    }
    std::discrete_distribution<int> tau_rng(prob_tau_high.begin(), prob_tau_high.end());
    mt19937_64 r(seed);
    int sample = tau_rng(r);
    result.eff_params.tau = tau_high_prob[sample];
  }

  result.eff_params.gamma0 = median(m.getNode(eff_params.gamma0).history);
  result.eff_params.gamma1 = median(m.getNode(eff_params.gamma1).history);

  for(int d = 0; d < ndose; d++) {
    int count_eff = 0;
    std::vector <double> q_median;
    vector<double> dose_eff_tau;
    for(int tau = 0; tau < ndose; tau++)
      dose_eff_tau.push_back(doseE[min(tau, d)]);
    for(auto draw: params_draws_list){
      double resp_rate;
      draw.first.responseRate(dose_eff_tau, resp_rate);
      count_eff += resp_rate < EFF_MIN;
      q_median.push_back(resp_rate);
    }
    result.qeff_inf[d] = double(count_eff)/n_samp;
    result.resp[d] = median(q_median);
    result.resp2[d] = median(q_median);
  }

  for(int d = 0; d < ndose; d++) {
    if(d > result.eff_params.tau)
      result.resp[d] = result.resp[d-1];
  }

  return result;
}

/* Estimate toxicity and efficacy parameters.
   The efficacy is estimated only for the given group. */
estimations estimate_pm(const trial_data& trial_data, uint_fast64_t, int group, bool) {
  const int ndose = trial_data.doseT.size();

  if(!HAS_TIME)
    throw std::logic_error("Internal error: HAS_TIME in PM.");

  vector<double> doseE(ndose);
  for(int d = 0; d < ndose; d++)
    doseE[d] = trial_data.doseE[d][group];

  /* Read usefull data from trial_data */
  int NPatientsCycle;
  for(NPatientsCycle = 0; NPatientsCycle < trial_data.pat_incl; NPatientsCycle++)
    if(trial_data.time_cur - trial_data.time_incl[NPatientsCycle] < CYCLE)
      break;

  const vec toxicity = conv_to<vec>::from(trial_data.toxicity).head(NPatientsCycle);
  const vec dose_tox = vec(trial_data.doseT)
    .elem(conv_to<uvec>::from(trial_data.dose_adm).head(NPatientsCycle));

  uvec in_group =
    find(conv_to<ivec>::from(trial_data.group).head(NPatientsCycle) == group);
  vec efficacy(in_group.size()), time_min_eff(in_group.size());

  int card2 = 0;
  for(unsigned i = 0; i < in_group.size(); i++) {
    double time_follow = min(TIMEFULL, trial_data.time_cur - trial_data.time_incl[in_group(i)]);
    if(trial_data.time_eff[in_group(i)] < time_follow) {
      efficacy(i) = 1;
      time_min_eff(i) = trial_data.time_eff[in_group(i)];
      if(time_follow == TIMEFULL)
	card2++;
    } else {
      efficacy(i) = 0;
      time_min_eff(i) = time_follow;
    }
  }

  vec weights(in_group.size());
  for(unsigned i = 0; i < in_group.size(); i++) {
    if(efficacy(i)) {
      weights(i) = 1;
    } else {
      double time_follow = min(TIMEFULL, trial_data.time_cur - trial_data.time_incl[in_group(i)]);
      int card1 = 0;
      for(unsigned j = 0; j < in_group.size(); j++)
        if(trial_data.time_cur - trial_data.time_incl[in_group(j)] >= TIMEFULL &&
           trial_data.time_eff[in_group(j)] <= time_follow)
          card1++;
      double w_ini = time_min_eff(i) / TIMEFULL;
      weights(i) = (card1+w_ini)/(card2+1);
    }
  }
  const uvec dose_adm_eff = conv_to<uvec>::from(trial_data.dose_adm).elem(in_group);

  vector<vec> dose_eff_tau;
  for(int tau = 0; tau < ndose; tau++)
    dose_eff_tau.push_back(vec(doseE).elem(clamp(dose_adm_eff, 0, tau)));

  /* MCMC initalization and sampling */
  efficacy_parameters eff_params(0, 0.001, ndose-1);
  toxicity_parameters tox_params(0, 0.001);

  vec pi_jk, resp_jk, resp_jk_weight;

  std::function<void ()> model = [&]() {
    eff_params.responseRate(dose_eff_tau, resp_jk);
    resp_jk_weight = weights % resp_jk;
    tox_params.proba_tox(dose_tox, pi_jk);
  };
  MCModel<std::mt19937> m(model);

  m.track<Normal>(eff_params.gamma0).dnorm(0, 0.01);
  m.track<Exponential>(eff_params.gamma1).dexp(1);
  vec distr(ndose);
  for(int i = 0; i < ndose; i++)
    distr[i] = 1.0;
  m.track<Discrete>(eff_params.tau).ddiscr(distr);
  m.track<ObservedBernoulli>(efficacy).dbern(resp_jk_weight);

  m.track<ObservedBernoulli>(toxicity).dbern(pi_jk);
  m.track<Normal>(tox_params.beta0).dnorm(0, 0.01);
  m.track<Exponential>(tox_params.beta1).dexp(1);

  m.sample(1e6, 1e5, 1e4, 10);

  /* Get list of eff_params and tox_params */
  vector<pair<efficacy_parameters, toxicity_parameters>> params_draws_list;
  auto itbeta0 = m.getNode(tox_params.beta0).history.begin(),
    itbeta0end = m.getNode(tox_params.beta0).history.end();
  auto itbeta1 = m.getNode(tox_params.beta1).history.begin();
  auto itgamma0 = m.getNode(eff_params.gamma0).history.begin();
  auto itgamma1 = m.getNode(eff_params.gamma1).history.begin();
  auto ittau = m.getNode(eff_params.tau).history.begin();

  while(itbeta0 != itbeta0end) {
    params_draws_list.emplace_back(
      efficacy_parameters(*itgamma0, *itgamma1, *ittau),
      toxicity_parameters(*itbeta0, *itbeta1));
    itbeta0++; itbeta1++;
    itgamma0++; itgamma1++; ittau++;
  }

  int n_samp = params_draws_list.size();

  /* Estimations */
  /* Estimation of toxicity */
  estimations result(ndose);
  result.tox_params.beta0 = median(m.getNode(tox_params.beta0).history);
  result.tox_params.beta1 = median(m.getNode(tox_params.beta1).history);

  for(int d = 0; d < ndose; d++) {
    int count_tox = 0;
    std::vector <double> pi_median;
    for(auto draw: params_draws_list){
      double proba_tox;
      draw.second.proba_tox(trial_data.doseT[d], proba_tox);
      count_tox += proba_tox < TARG_SUP;
      pi_median.push_back(proba_tox);
    }
    result.ptox_inf[d] = double(count_tox)/n_samp;
    result.pi[d] = median(pi_median);
  }

  /* Estimation of efficacy */
  vector<vector<double>> qeff_inf_tau(ndose, vector<double>(ndose));
  vector<vector<double>> resp_tau(ndose, vector<double>(ndose));
  for(int tau = 0; tau < ndose; tau++) {
    int count = 0;
    vector<double> gammas0, gammas1;
    for(auto draw: params_draws_list) {
      if(tau == draw.first.tau) {
        count++;
        gammas0.push_back(draw.first.gamma0);
        gammas1.push_back(draw.first.gamma1);
      }
    }
    if(count > 0) {
      result.proba_tau[tau] = (double)count/n_samp;

      double prevresp = 0/*Dummy*/;
      for(int d = 0; d < ndose; d++) {
        int count_eff = 0;
        std::vector <double> q_median;
        vector<double> dose_eff_tau;
        for(int tau = 0; tau < ndose; tau++)
          dose_eff_tau.push_back(doseE[min(tau, d)]);
        for(auto draw: params_draws_list){
          if(tau == draw.first.tau) {
            double resp_rate;
            draw.first.responseRate(dose_eff_tau, resp_rate);
            count_eff += resp_rate < EFF_MIN;
            q_median.push_back(resp_rate);
          }
        }
        qeff_inf_tau[d][tau] = double(count_eff)/count;
        if(d <= tau)
          resp_tau[d][tau] = prevresp = median(q_median);
        else
          resp_tau[d][tau] = prevresp;
      }
    }
    else{
      result.proba_tau[tau] = 0;
    }
  }
  vector<double> proba_resp_BMA(ndose, 0);
  for(int d = 0; d < ndose; d++) {
    for(int tau = 0; tau < ndose; tau++) {
      proba_resp_BMA[d] += resp_tau[d][tau]*result.proba_tau[tau];
    }
    result.resp2[d] = proba_resp_BMA[d];
  }
  result.eff_params.tau = 0;
  for(int tau=ndose-1; tau>0; tau--){
    if(abs(proba_resp_BMA[tau]-proba_resp_BMA[tau-1]) > trial_data.s_2){
        result.eff_params.tau = tau;
        break;
    }
  }

  for(int d = 0; d < ndose; d++) {
    result.qeff_inf[d] = qeff_inf_tau[d][result.eff_params.tau];
    result.resp[d] = resp_tau[d][result.eff_params.tau];
  }

  return result;
}

void take_if_better(const estimations& estim, int& nextdose, int candidate_dose) {
  const int ndose = estim.pi.size();

  if(nextdose == -1) {
    nextdose = candidate_dose;
    return;
  }

  if(nextdose < 0 || candidate_dose < 0 || nextdose >= ndose || candidate_dose >= ndose)
    throw std::logic_error("Internal error: invalid nextdose or candidate_dose.");

  double candidateResp = estim.resp[candidate_dose];
  double bestResp = estim.resp[nextdose];

  if(candidateResp > bestResp ||
     (candidateResp == bestResp && estim.pi[candidate_dose] < estim.pi[nextdose]) ) {
    nextdose = candidate_dose;
  }
}

// Choose the next dose for the given group.
// Returns the dose, or -1 if the trial must be stopped for this group.
int find_next_dose(trial_data& trial_data, int group, double c_tox, double c_eff,
                   bool final, estimations* estimRet = NULL){
  const int ndose = trial_data.doseT.size();

  if(trial_data.startup_end[group] == -1 && final)
    trial_data.startup_end[group] = trial_data.pat_incl;

  // Warning : this code should detect whether the startup has ended, even
  // if it was a long time ago.
  if(trial_data.startup_end[group] == -1) {
    if(trial_data.pat_incl_group[group] == 0)
      return 0;

    int dlt = 0;
    for(int pat = 0; pat < trial_data.pat_incl; pat++)
      if(trial_data.group[pat] == group && trial_data.toxicity[pat])
        dlt++;

    if(dlt == 0 && trial_data.cdose[group] != ndose-1)
      return trial_data.cdose[group] + 1;

    trial_data.startup_end[group] = trial_data.pat_incl_group[group];
  }

  estimations estim = trial_data.est(trial_data, trial_data.r(), group, final);
  if(estimRet != NULL)
    *estimRet = estim;

  int nextdose = -1;
  const int cdose = trial_data.cdose[group];

  vector<int> pat_explo_dose(ndose, 0);
  int highest_tested=0;
  for(int pat = 0; pat < trial_data.pat_incl; pat++){
    if(trial_data.group[pat] == group){
      pat_explo_dose[trial_data.dose_adm[pat]]++;
      if((int)trial_data.dose_adm[pat] > highest_tested){
        highest_tested = trial_data.dose_adm[pat];
      }
    }
  }

  vector<int> pats;
  for(int pat = 0; pat < trial_data.pat_incl; pat++) {
    if(trial_data.group[pat] == group)
      pats.push_back(pat);
  }
  int dlt = 0;
  for(int pat = 0; pat < min(2*max(COHORT,COHORT_START),trial_data.pat_incl_group[group]); pat++){
    if(trial_data.toxicity[pats[pat]]){
        dlt++;
    }
  }

  int max_candidate_dose;
  if(final) max_candidate_dose = highest_tested;
  else max_candidate_dose = max(min(ndose-1, cdose+1), highest_tested);

  for(int candidate_dose = 0; candidate_dose <= max_candidate_dose; candidate_dose++) {
    // Determination of the next dose or stopping rule
    if( (1-estim.ptox_inf[candidate_dose] < c_tox &&
         (1-estim.qeff_inf[candidate_dose] >= c_eff ||
          (!final && pat_explo_dose[candidate_dose] < 2*max(COHORT,COHORT_START)))) ||
        (!final &&
         candidate_dose==0 &&
         trial_data.pat_incl_group[group] < 2*max(COHORT,COHORT_START) &&
         dlt==1) ) {
      take_if_better(estim, nextdose, candidate_dose);
    }
  }

  return nextdose;
}

}

using namespace dfmta;

extern "C" {

static R_NativePrimitiveArgType dfmta_next_args[] =
  {LGLSXP, LGLSXP,
   INTSXP, INTSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   INTSXP, INTSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   INTSXP,

   INTSXP,
   REALSXP,
   INTSXP,
   INTSXP, INTSXP,
   REALSXP, REALSXP,
   LGLSXP, LGLSXP,
   LGLSXP,

   LGLSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   REALSXP };
static const int dfmta_next_nargs = 33;

static void dfmta_next(int* has_time, int* ra_pm /* true for RA, false for PM */,
                       int* ngroups, int* ndose,
                       double* targ_sup, double* eff_min,
                       double* doseTV0, double* doseEV0,
                       double* time_full, double* cycle,
                       int* cohort_start, int* cohort,
                       double* c_tox, double* c_eff,
                       double* s_1, double* s_2,
                       int* seed,

                       int* cdose /* -1 if group stopped */,
                       double* time_cur,
                       int* pat_incl_group,
                       int* dose_adm, int* group,
                       double* time_eff, double* time_incl,
                       int* efficacy, int* toxicity,
                       int* final,

                       int* in_startup,
                       double* pi, double* ptox_inf,
                       double* resp2, double* qeff_inf,
                       double* proba_tau) {
  try {

  CYCLE = *cycle;
  TARG_SUP = *targ_sup;
  EFF_MIN = *eff_min;
  HAS_TIME = *has_time;
  if(HAS_TIME) TIMEFULL = *time_full;

  int pat_incl = 0;
  for(int g = 0; g < *ngroups; g++)
    pat_incl += pat_incl_group[g];
  int group_cur = group[pat_incl];

  vector<double> doseTV(doseTV0, doseTV0+*ndose);
  vector<vector<double>> doseEV(*ndose, vector<double>(*ngroups, 0));
  for(int i = 0; i < *ndose; i++)
    doseEV[i][group_cur] = doseEV0[i];

  COHORT_START = *cohort_start;
  COHORT = *cohort;

  estimator est = *ra_pm ? estimate_ra : estimate_pm;
  vector<double> _s_1;
  if(!*final) {
    _s_1.resize(pat_incl_group[group_cur]+1);
    _s_1[pat_incl_group[group_cur]] = *s_1;
  }
  trial_data trial_data(est, doseTV, doseEV, *ngroups, *s_2, _s_1, *seed);
  trial_data.pat_incl = pat_incl;
  trial_data.cdose = vector<int>(*ngroups, 0);
  trial_data.cdose[group_cur] = *cdose;
  trial_data.startup_end = vector<int>(*ngroups, 0);
  // Assume the startup has not ended : if it has, it will detect it.
  trial_data.startup_end[group_cur] = -1;
  trial_data.time_cur = *time_cur;
  trial_data.pat_incl_group = vector<int>(pat_incl_group, pat_incl_group+*ngroups);
  trial_data.dose_adm = vector<unsigned int>(dose_adm, dose_adm+pat_incl);
  trial_data.group = vector<int>(group, group+pat_incl);
  if(HAS_TIME) {
    trial_data.time_eff = vector<double>(time_eff, time_eff+pat_incl);
    trial_data.time_incl = vector<double>(time_incl, time_incl+pat_incl);
  } else {
    trial_data.efficacy = vector<int>(efficacy, efficacy+pat_incl);
  }
  trial_data.toxicity = vector<int>(toxicity, toxicity+pat_incl);

  estimations estim(*ndose);
  *cdose = find_next_dose(trial_data, group_cur, *c_tox, *c_eff, *final, &estim);
  *in_startup = trial_data.startup_end[group_cur] == -1;

  if(!*in_startup)
    for(int d = 0; d < *ndose; d++) {
      pi[d] = estim.pi[d];
      ptox_inf[d] = estim.ptox_inf[d];
      resp2[d] = estim.resp2[d];
      qeff_inf[d] = estim.qeff_inf[d];
      proba_tau[d] = estim.proba_tau[d];
    }
  }
  catch (std::logic_error &e) {
    error("Internal error in dfmta (details: %s)", e.what());
  }
  catch (...) {
    error("Internal error in dfmta");
  }
}

static R_NativePrimitiveArgType dfmta_simu_args[] =
  {LGLSXP, LGLSXP,
   INTSXP, INTSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   REALSXP, INTSXP,
   REALSXP, REALSXP,
   INTSXP, INTSXP, INTSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   INTSXP, INTSXP,

   INTSXP, INTSXP, INTSXP,
   INTSXP, INTSXP, INTSXP,
   INTSXP, INTSXP, INTSXP,
   REALSXP};
static const int dfmta_simu_nargs = 33;

static void dfmta_simu(int* has_time, int* ra_pm /* true for RA, false for PM */,
                       int* ngroups, int* ndose,
                       double* piV, double* respV,
                       double* targ_sup, double* eff_min,
                       double* doseTV0, double* doseEV0,
                       double* incl_per_week, int* npatientspergroup,
                       double* time_full, double* cycle,
                       int* cohort_start, int* cohort, int* ntrial,
                       double* c_tox, double* c_eff,
                       double* s_1, double* s_2,
                       int* seed, int* nthreads,

                       int* inconc, int* n_pat_dose, int* rec_dose,
                       int* n_pat_tot, int* n_tox, int* n_eff,
                       int* tox_tot, int* eff_tot, int* n_pat_mtd,
                       double* duration) {
  string errstr;

  {
    struct true_data true_data;

    CYCLE = *cycle;
    HAS_TIME = *has_time;
    if(HAS_TIME) TIMEFULL = *time_full;

    vector<double> _s_1(s_1, s_1+*npatientspergroup);
    true_data.piV = vector<double>(piV, piV+*ndose);

    true_data.respV = vector<vector<double> >(*ndose, vector<double>(*ngroups));
    for(int g = 0; g < *ngroups; g++)
      for(int i = 0; i < *ndose; i++)
        true_data.respV[i][g] = respV[i**ngroups + g];

    TARG_SUP = *targ_sup;
    EFF_MIN = *eff_min;

    vector<double> doseTV(doseTV0, doseTV0+*ndose);
    vector<vector<double> > doseEV(*ndose, vector<double>(*ngroups, 0));
    for(int g = 0; g < *ngroups; g++)
      for(int i = 0; i < *ndose; i++)
        doseEV[i][g] = doseEV0[i**ngroups + g];
    true_data.incl_per_week = vector<double>(incl_per_week, incl_per_week+*ngroups);

    COHORT_START = *cohort_start;
    COHORT = *cohort;

    std::mt19937_64 global_rng(*seed);
    vector<uint_fast64_t> seeds;
    for(int trial=0; trial<*ntrial; trial++)
      seeds.push_back(global_rng());

    struct results results(*ndose, *ngroups);

#ifdef _OPENMP
    if(*nthreads > 0)
      omp_set_num_threads(*nthreads);
    else
      omp_set_num_threads(omp_get_num_procs());
#endif

    bool err = false;
    Progress prog(*ntrial);
#pragma omp parallel for schedule(dynamic, 1)
    for(int trial=0; trial<*ntrial; trial++) {
      try {
        bool err2;
#pragma omp critical
        err2 = err;
        if(prog.check_abort() || err2)
          continue;

        estimator est = *ra_pm ? estimate_ra : estimate_pm;
        struct trial_data trial_data(est, doseTV, doseEV, *ngroups, *s_2, _s_1, seeds[trial]);

        while(true) {
#pragma omp critical
          err2 = err;
          if(prog.check_abort() || err2)
            goto aborted;

          int group;
          if(wait_patient(trial_data, true_data, group))
            break;

          int startup_end = trial_data.startup_end[group];
          if( (startup_end == -1 && trial_data.pat_incl_group[group] % COHORT_START == 0) ||
              (startup_end != -1 && (trial_data.pat_incl_group[group] - startup_end) % COHORT == 0) )
            trial_data.cdose[group] = find_next_dose(trial_data, group, *c_tox, *c_eff, false);

          if(trial_data.cdose[group] >= 0) {
            include_patient(trial_data, true_data, group);
            if(trial_data.pat_incl_group[group] >= *npatientspergroup)
              trial_data.cdose[group] = -2; /* No more patiens for this group */
          }
        }

        vector<int> recom;
        for(int group = 0; group < *ngroups; group++)
          if(trial_data.cdose[group] == -2)
            recom.push_back(find_next_dose(trial_data, group, *c_tox, *c_eff, true));
          else
            recom.push_back(-1);

#pragma omp critical
        results.accumulate(trial_data, recom);

      } catch (std::logic_error &e) {
#pragma omp critical
        { err = true;
          errstr = e.what(); }
      } catch(...) {
#pragma omp critical
        err = true;
      }

#pragma omp critical
      prog.increment();
      aborted: ;
    }

    if(err) goto errlbl;

    copy(results.inconc.begin(), results.inconc.end(), inconc);
    for(int i = 0; i < *ndose; i++)
      for(int j = 0; j < *ngroups; j++) {
        n_pat_dose[i**ngroups + j] = results.n_pat_dose[i][j];
        rec_dose[i**ngroups + j]   = results.rec_dose[i][j];
        n_tox[i**ngroups + j]      = results.n_tox[i][j];
        n_eff[i**ngroups + j]      = results.n_eff[i][j];
      }
    *n_pat_tot = results.n_pat_tot;
    *tox_tot   = results.tox_tot;
    *eff_tot   = results.eff_tot;
    *duration  = results.duration;
    *ntrial    = results.nb_trials;
    copy(results.n_pat_mtd.begin(), results.n_pat_mtd.end(), n_pat_mtd);

  }

  if(false) {
  errlbl:
    error("Internal error in dfmta (details: %s)", errstr.c_str());
  }
}

static const R_CMethodDef cMethods[] = {
  {"dfmta_next", (DL_FUNC) &dfmta_next,
   dfmta_next_nargs, dfmta_next_args},
  {"dfmta_simu", (DL_FUNC) &dfmta_simu,
   dfmta_simu_nargs, dfmta_simu_args},
  {NULL, NULL, 0}
};

void R_init_dfmta(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

}
