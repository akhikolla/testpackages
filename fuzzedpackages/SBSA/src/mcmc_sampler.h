#ifndef BSA_WALKER_H
#define BSA_WALKER_H

#include "bsalite.h"
#include "log_functor.h"

class McmcReparametrizingSampler {
public:
  McmcReparametrizingSampler(const LogLik& loglik_, 
                             const LogPri& logpri_)
    : loglik(loglik_), 
      logpri(logpri_) {}
  
  bool sample(Theta& theta, double& log_lik, double& log_pri) const;
protected:
  virtual Theta reparametrize(const Theta& theta_cur) const = 0;
private:
  const LogLik& loglik;
  const LogPri& logpri;
};

class ReparametrizeAlpha : public McmcReparametrizingSampler {
public:
  ReparametrizeAlpha(const LogLik& loglik_, 
                     const LogPri& logpri_,
                     const arma::vec& sigma_)
    : McmcReparametrizingSampler::McmcReparametrizingSampler(loglik_, logpri_),
      sigma(sigma_) {}
protected:
  virtual Theta reparametrize(const Theta& theta_cur) const;
private:
  const arma::vec sigma;
};

class ReparametrizeBetaZ : public McmcReparametrizingSampler {
public:
  ReparametrizeBetaZ(const LogLik& loglik_, 
                     const LogPri& logpri_,
                     const int p_,
                     const arma::vec& sigma_)
    : McmcReparametrizingSampler::McmcReparametrizingSampler(loglik_, logpri_),
      p(p_),
      sigma(sigma_) {}
protected:
  virtual Theta reparametrize(const Theta& theta_cur) const;
private:
  const int p;
  const arma::vec sigma;
};

class ReparametrizeSigmaSq : public McmcReparametrizingSampler {
public:
  ReparametrizeSigmaSq(const LogLik& loglik_, 
                       const LogPri& logpri_,
                       const double sigma_)
    : McmcReparametrizingSampler::McmcReparametrizingSampler(loglik_, logpri_),
      sigma(sigma_) {}
protected:
  virtual Theta reparametrize(const Theta& theta_cur) const;
private:
  const double sigma;
};

class ReparametrizeTauSq : public McmcReparametrizingSampler {
public:
  ReparametrizeTauSq(const LogLik& loglik_, 
                     const LogPri& logpri_,
                     const int p,
                     const arma::vec& sigma_)
    : McmcReparametrizingSampler::McmcReparametrizingSampler(loglik_, logpri_),
      p(p),
      sigma(sigma_) {}
protected:
  virtual Theta reparametrize(const Theta& theta_cur) const;
private:
  const int p;
  const arma::vec sigma;
};

class ReparametrizeBetaUGammaX : public McmcReparametrizingSampler {
public:
  ReparametrizeBetaUGammaX(const LogLik& loglik_, 
                           const LogPri& logpri_,
                           const double sigma_,
                           const double el2_)
    : McmcReparametrizingSampler::McmcReparametrizingSampler(loglik_, logpri_),
      sigma(sigma_), el2(el2_) {}
protected:
  virtual Theta reparametrize(const Theta& theta_cur) const;
private:
  const double sigma;
  const double el2;
};

class ReparametrizeGammaZ : public McmcReparametrizingSampler {
public:
  ReparametrizeGammaZ(const LogLik& loglik_, 
                      const LogPri& logpri_,
                      const int p_,
                      const arma::vec& sigma_)
    : McmcReparametrizingSampler::McmcReparametrizingSampler(loglik_, logpri_),
      p(p_),
      sigma(sigma_) {}
protected:
  virtual Theta reparametrize(const Theta& theta_cur) const;
private:
  const int p;
  const arma::vec sigma;
};

#endif
