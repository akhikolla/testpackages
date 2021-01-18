#include "fssemR.h"

RcppExport SEXP RidgeReg(SEXP X_, SEXP Y_, SEXP S_, SEXP gamma_)
{
BEGIN_RCPP
  Rcpp::RObject res;
  Rcpp::RNGScope rcpp_rngScope;
  MatrixXf X(Rcpp::as<MatrixXf>(X_));  // n x k
  MatrixXf Y(Rcpp::as<MatrixXf>(Y_));  // n x p
  MatrixXf xm, ym;
  const double gamma = Rcpp::as<double>(gamma_);
  const int n = X.rows();
  const int p = Y.cols();
  const int k = X.cols();
  std::vector<ArrayXd> S = Rcpp::as<std::vector<ArrayXd> >(S_);
  MatrixXf B(p, p);
  std::vector<MatrixXf> fi(p);
  MatrixXf mu;
  center(X, xm);
  center(Y, ym);
  B.setZero();
  MatrixXf Xi, Yi, yi, bi;
  double error = 0, sigma2;
#ifdef _OPENMP
  #pragma omp parallel for private(Xi, Yi, yi, bi) shared(B, fi) reduction(+:error)
#endif
  for(int i = 1; i <= p; i++)
  {
    Xi = get_Cols(X, S[i-1]);
    Yi = rm_Col(Y, i);
    yi = Y.col(i-1);
    bi = LR_Shrink(Xi, Yi, yi, gamma, n, p, k);
    set_Row(B, i, bi);
    fi[i-1] = LR_Fi(Xi, Yi, yi, bi);
    error += (yi - Yi * bi - Xi * fi[i-1]).norm();
  }
  sigma2 = error / (n * p - 1);
  mu     = (MatrixXf::Identity(p, p) - B) * ym;
#ifdef _OPENMP
  #pragma omp parallel for shared(mu)
#endif
  for (int i = 1; i <= p; i++) {
    mu(i-1, 0) -= (xm.transpose() * fi[i-1])(0,0);
  }
  MatrixXf F = get_Fs(fi, S, k);
  double err = LR_Objerr(X, Y, B, F, mu);
  res = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("B") = B,
                                      Rcpp::Named("f") = fi,
                                      Rcpp::Named("F") = F,
                                      Rcpp::Named("mu") = mu)
  );
  return res;
END_RCPP
}

RcppExport SEXP L2lamax(SEXP Xs_, SEXP Ys_, SEXP S_, SEXP n_, SEXP p_, SEXP k_)
{
BEGIN_RCPP
  std::vector<MatrixXf> Xs = Rcpp::as<std::vector<MatrixXf> >(Xs_);
  std::vector<MatrixXf> Ys = Rcpp::as<std::vector<MatrixXf> >(Ys_);
  std::vector<ArrayXd>  S  = Rcpp::as<std::vector<ArrayXd> >(S_);
  // const int n = Rcpp::as<int>(n_);
  std::vector<int> n = Rcpp::as<std::vector<int> >(n_);
  const int p = Rcpp::as<int>(p_);
  const int k = Rcpp::as<int>(k_);
  const int m = Xs.size();
  double lambda = 0;
  for(int i = 0; i < m; i++)
  {
    lambda = max(lambda, L2lamax(Xs[i], Ys[i], S, n[i], p, k));
  }
  return Rcpp::wrap(lambda);
END_RCPP
}

RcppExport SEXP MultiReg(SEXP Xs_, SEXP Ys_, SEXP S_, SEXP gamma_, SEXP n_, SEXP p_, SEXP k_)
{
BEGIN_RCPP
  std::vector<MatrixXf> Xs = Rcpp::as<std::vector<MatrixXf> >(Xs_);
  std::vector<MatrixXf> Ys = Rcpp::as<std::vector<MatrixXf> >(Ys_);
  const double gamma       = Rcpp::as<double>(gamma_);
  std::vector<ArrayXd>  S  = Rcpp::as<std::vector<ArrayXd> >(S_);
  // const int n = Rcpp::as<int>(n_);
  std::vector<int> n = Rcpp::as<std::vector<int> >(n_);
  const int p = Rcpp::as<int>(p_);
  const int k = Rcpp::as<int>(k_);
  const int m = Xs.size();
  std::vector<MatrixXf> Bs(m);
  std::vector<std::vector<MatrixXf> > fs(m);
  std::vector<MatrixXf> Fs(m);
  std::vector<MatrixXf> mu(m);
  double error2 = 0, df = 0;
  for(int i = 0; i < m; i++)
  {
    Bs[i].resize(p, p);
    Bs[i].setZero();
    fs[i].resize(p);
    error2 += L2lr(Xs[i], Ys[i], S, Bs[i], fs[i], mu[i], gamma, n[i], p, k);
    df += n[i] * p - 1;
  }
  double sigma2 = error2 / df;
  for(int i = 0; i < m; i++)
  {
    Fs[i] = get_Fs(fs[i], S, k).transpose();
  }
  return Rcpp::List::create(
    Rcpp::Named("Bs") = Bs,
    Rcpp::Named("fs") = fs,
    Rcpp::Named("Fs") = Fs,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("sigma2") = sigma2
  );
END_RCPP
}

RcppExport SEXP ObjMultiReg(SEXP Xs_, SEXP Ys_, SEXP fit_)
{
BEGIN_RCPP
  std::vector<MatrixXf> Xs = Rcpp::as<std::vector<MatrixXf> >(Xs_);
  std::vector<MatrixXf> Ys = Rcpp::as<std::vector<MatrixXf> >(Ys_);
  Rcpp::List fit(fit_);
  std::vector<MatrixXf> Bs = Rcpp::as<std::vector<MatrixXf> >(fit["Bs"]);
  std::vector<MatrixXf> Fs = Rcpp::as<std::vector<MatrixXf> >(fit["Fs"]);
  std::vector<MatrixXf> mu = Rcpp::as<std::vector<MatrixXf> >(fit["mu"]);
  double error2 = 0;
  int m = Xs.size();
#ifdef _OPENMP
  #pragma omp parallel for shared(Xs, Ys, Bs, Fs, mu) reduction(+:error2)
#endif
  for(int i = 0; i < m; i++) {
    error2 += LR_Objerr(Xs[i], Ys[i], Bs[i], Fs[i], mu[i]);
  }
  double err = std::sqrt(error2);
  return Rcpp::wrap(err);
END_RCPP
}

RcppExport SEXP MultiFSSEMiPALM(SEXP Xs_, SEXP Ys_, SEXP Bs_, SEXP Fs_, SEXP S_, SEXP sigma2_,
                                SEXP lambda_, SEXP rho_, SEXP options_,
                                SEXP n_, SEXP p_, SEXP k_, SEXP m_)
{
BEGIN_RCPP
  std::vector<MatrixXf> Xs = Rcpp::as<std::vector<MatrixXf> >(Xs_);
  std::vector<MatrixXf> Ys = Rcpp::as<std::vector<MatrixXf> >(Ys_);
  std::vector<ArrayXd>  S  = Rcpp::as<std::vector<ArrayXd> >(S_);
  double sigma2 = Rcpp::as<double>(sigma2_);
  double lambda = Rcpp::as<double>(lambda_);
  double rho    = Rcpp::as<double>(rho_);
  const std::vector<int> n = Rcpp::as<std::vector<int> >(n_);
  const int p   = Rcpp::as<int>(p_);
  const int k   = Rcpp::as<int>(k_);
  const int m   = Rcpp::as<int>(m_);
  // option vars
  Rcpp::List options(options_);
  int maxit   = Rcpp::as<int>(options["maxit"]);
  double threshold = Rcpp::as<double>(options["threshold"]);
  int i, j;
  std::vector<MatrixXf> Wl(m);
  MatrixXf Wf;
  // center
  std::vector<MatrixXf> Xm(m);
  std::vector<MatrixXf> Ym(m);
#ifdef _OPENMP
  #pragma omp parallel for private(i) shared(Xs, Ys, Xm, Ym)
#endif
  for(i = 0; i < m; i++) {
    center(Xs[i], Xm[i], true);                                                 // X => k x n
    center(Ys[i], Ym[i], true);                                                 // Y => p x n
  }
  // Bs & Fs
  std::vector<MatrixXf> Bs = Rcpp::as<std::vector<MatrixXf> >(Bs_);             // B => p x p
  std::vector<MatrixXf> Fs = Rcpp::as<std::vector<MatrixXf> >(Fs_);             // F => k x p
  std::vector<MatrixXf> mu(m);
  std::vector<MatrixXf> Fx(m);
  std::vector<VectorXf> Y2norm(m);
  std::vector<VectorXf> B2norm(m);
    MatrixXf Xi, Yi, yi, P;
  std::vector<std::vector<MatrixXf> > f0(m);
  std::vector<std::vector<MatrixXf> > f1(m);
  for(i = 0; i < m; i++)
  {
    Fx[i] = Fs[i].transpose() * Xs[i];
    Wl[i] = (Bs[i].cwiseInverse()).cwiseAbs();                                  // lasso weight
    f0[i].resize(p);
    f1[i].resize(p);
    Y2norm[i].resize(p);
    B2norm[i].resize(p);
#ifdef _OPENMP
    #pragma omp parallel for private(j, Xi, P, Yi, yi) shared(Y2norm, f0, f1)
#endif
    for(j = 1; j <= p; j++)
    {
      Xi = get_Rows(Xs[i], S[j-1]);
      P  = ProjectMat(Xi);
      Yi = rm_Row(Ys[i], j);
      yi = Ys[i].row(j-1);
      Y2norm[i][j-1] = yi.squaredNorm();
      f0[i][j-1] = P * yi.transpose();
      f1[i][j-1] = P * Yi.transpose();
      B2norm[i][j-1] = Bs[i].col(j-1).squaredNorm();
    }
  }
  Wf = (((Bs[1].array() - Bs[0].array()).matrix()).cwiseInverse()).cwiseAbs();  // fused lasso weight
  int niter = 0;
  std::vector<MatrixXf> ImBs(m);
  std::vector<double> Dets(m);
  std::vector<MatrixXf> IBinv(m);
  for(i = 0; i < m; i++)
  {
    Bs[i].setZero();
    ImBs[i]  = MatrixXf::Identity(p, p) - Bs[i];
    IBinv[i] = inverse(ImBs[i], Dets[i]);
  }
  double L = loglikelihoodFSSEM(Bs, Wl, Wf, Dets, lambda, rho, n, p, m, sigma2);
  // aux vars
  std::vector<std::vector<MatrixXf> > Bs_last{Bs, Bs};
  std::vector<MatrixXf> Fs_last;
  double L_last;
  // status check
  bool zeroDet, nonincreased = true, converged;
  // temp vars
  double t, minit = 1, c2, Li, tau, error, lerror;
  std::vector<MatrixXf> Bt(m), ci(m), bi(m), Ri(m), gi(m);
  std::vector<MatrixXf> ui(m), w(m), xi(m), dB(m);
  std::vector<double> Detm(m);
  MatrixXf grad, r, fi;
  VectorXf li(m);
  while(niter <= maxit)
  {
    t = inert(niter);
    Fs_last = Fs;
    L_last  = L;
    if(!nonincreased)
    {
      if(minit <= 1e6)
        minit *= 1.2;
      else
        break;
    }
    for(i = 0; i < m; i++)
    {
      Bt[i] = Bs_last[1][i] + t * (Bs_last[1][i] - Bs_last[0][i]);
    }
    // block coordinate descent
    for(i = 1; i <= p; i++)
    {
      for(j = 0; j < m; j++)
      {
        ci[j] = IBinv[j].row(i-1);
        bi[j] = Bt[j].col(i-1);
        Ri[j] = Ys[j] - rm_Col(Bs[j], i) * rm_Row(Ys[j], i) - Fx[j];
        grad  = n[j] * ci[j].transpose() + (Y2norm[j][i-1] * bi[j] - Ri[j] * Ys[j].row(i-1).transpose()) / sigma2;
        gi[j] = rm_Row(grad, i);
        c2    = (ci[j] * Dets[j]).squaredNorm();
        li[j] = n[j] * c2 / std::pow(1 - std::sqrt(c2 * B2norm[j][i-1]), 2) + Y2norm[j][i-1] / sigma2;
      }
      Li = li.norm();
      // Armoji-scheme
      tau = minit;
      while(true) {
        for(j = 0; j < m; j++)
        {
          ui[j] = rm_Row(bi[j], i) - 1 / (tau * Li) * gi[j];
          w[j]  = rm_Row(Wl[j], i).col(i-1);
        }
        r = rm_Row(Wf, i).col(i-1);
        prox_flsa(ui, xi, lambda, rho, tau * Li, w, r);
        zeroDet = false;
        for(j = 0; j < m; j++)
        {
          Detm[j] = (IBinv[j](i-1, i-1) - (rm_Col(IBinv[j], i).row(i-1) * xi[j])(0, 0)) * Dets[j];
          if (Detm[j] == 0)
          {
            zeroDet = true;
          }
        }
        if (!zeroDet) {
          break;
        }
        tau = tau * 1.2;
      }
      // update B
      for(j = 0; j < m; j++)
      {
        set_Col(Bs[j], i, xi[j]);
        ImBs[j] = MatrixXf::Identity(p, p) - Bs[j];
        IBinv[j] = inverse(ImBs[j], Dets[j]);
      }
    }
    // update Fs
    for(j = 0; j < m; j++)
    {
      for(i = 1; i <= p; i++)
      {
        fi = f0[j][i-1] - f1[j][i-1] * rm_Col(Bs[j], i).row(i-1);
        set_Row(Fs[j], fi, S[i-1], i-1);
      }
      Fx[j] = Fs[j].transpose() * Xs[j];
    }
    // update sigma2
    sigma2 = get_sigma2(Xs, Ys, Bs, Fs, m);
    // converge
    error = 0;
    for(j = 0; j < m; j++)
    {
      error += get_err(Bs[j], Bs_last[1][j]) + get_err(Fs[j], Fs_last[j]);
    }
    L = loglikelihoodFSSEM(Bs, Wl, Wf, Dets, lambda, rho, n, p, m, sigma2);
    lerror = std::abs(L_last - L) / (1 + std::abs(L_last));
    nonincreased = (L_last > L);
    converged = (error < threshold) && (lerror < threshold);
    // record Bs history
    Bs_last[0] = Bs_last[1];
    Bs_last[1] = Bs;
    niter++;
    // debug printf
    Rprintf("FSSEM: iter = %d, error = %.5f, loglik = %.5f, sigma2 = %.5f\n", niter, error, L, sigma2);
    if(converged || niter > maxit)
    {
      break;
    }
  }
  return Rcpp::List::create(Rcpp::Named("Bs") = Rcpp::wrap(Bs),
                            Rcpp::Named("Fs") = Rcpp::wrap(Fs),
                            Rcpp::Named("sigma2") = Rcpp::wrap(sigma2));
END_RCPP
}

void R_init_fssemR(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
