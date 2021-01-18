#ifndef BSA_BINARY_WALKER_H
#define BSA_BINARY_WALKER_H

#include "bsalite.h"
#include "log_binary_functor.h"

class McmcReparametrizingBinarySampler {
public:
  McmcReparametrizingBinarySampler(const LogBinaryLik& loglik_, 
                                   const LogBinaryPri& logpri_)
    : loglik(loglik_), 
      logpri(logpri_)
  {}
  
  bool sample(ThetaBinary& theta, double& log_lik, double& log_pri) const;
protected:
  virtual ThetaBinary reparametrize(const ThetaBinary& theta_cur) const = 0;
private:
  const LogBinaryLik& loglik;
  const LogBinaryPri& logpri;
};

class ReparametrizeBinaryAlpha : public McmcReparametrizingBinarySampler {
public:
  ReparametrizeBinaryAlpha(const LogBinaryLik& loglik_, 
                           const LogBinaryPri& logpri_,
                           const arma::vec& sampler_jump)
    : McmcReparametrizingBinarySampler::McmcReparametrizingBinarySampler(loglik_, logpri_),
      sampler_jump(sampler_jump)
  {}
protected:
  virtual ThetaBinary reparametrize(const ThetaBinary& theta_cur) const;
private:
  const arma::vec sampler_jump;
};

class ReparametrizeBinaryBetaZ : public McmcReparametrizingBinarySampler {
public:
  ReparametrizeBinaryBetaZ(const LogBinaryLik& loglik_, 
                           const LogBinaryPri& logpri_,
                           const int p_,
                           const arma::vec& sampler_jump)
    : McmcReparametrizingBinarySampler::McmcReparametrizingBinarySampler(loglik_, logpri_),
      p(p_),
      sampler_jump(sampler_jump)
  {}
protected:
  virtual ThetaBinary reparametrize(const ThetaBinary& theta_cur) const;
private:
  const int p;
  const arma::vec sampler_jump;
};

class ReparametrizeBinaryTauSq : public McmcReparametrizingBinarySampler {
public:
  ReparametrizeBinaryTauSq(const LogBinaryLik& loglik_, 
                           const LogBinaryPri& logpri_,
                           const int p_,
                           const arma::vec& sampler_jump)
    : McmcReparametrizingBinarySampler::McmcReparametrizingBinarySampler(loglik_, logpri_),
      p(p_),
      sampler_jump(sampler_jump)
  {}
protected:
  virtual ThetaBinary reparametrize(const ThetaBinary& theta_cur) const;
private:
  const int p;
  const arma::vec sampler_jump;
};

class ReparametrizeBinaryBetaUGammaX : public McmcReparametrizingBinarySampler {
public:
  ReparametrizeBinaryBetaUGammaX(const LogBinaryLik& loglik_, 
                                 const LogBinaryPri& logpri_,
                                 const double sampler_jump,
                                 const double el2_)
    : McmcReparametrizingBinarySampler::McmcReparametrizingBinarySampler(loglik_, logpri_),
      sampler_jump(sampler_jump),
      el2(el2_)
  {}
protected:
  virtual ThetaBinary reparametrize(const ThetaBinary& theta_cur) const;
private:
  const double sampler_jump;
  const double el2;
};

class ReparametrizeBinaryGammaZ : public McmcReparametrizingBinarySampler {
public:
  ReparametrizeBinaryGammaZ(const LogBinaryLik& loglik_, 
                            const LogBinaryPri& logpri_,
                            const int p_,
                            const arma::vec& sampler_jump)
    : McmcReparametrizingBinarySampler::McmcReparametrizingBinarySampler(loglik_, logpri_),
      p(p_),
      sampler_jump(sampler_jump)
  {}
protected:
  virtual ThetaBinary reparametrize(const ThetaBinary& theta_cur) const;
private:
  const int p;
  const arma::vec sampler_jump;
};

#endif
