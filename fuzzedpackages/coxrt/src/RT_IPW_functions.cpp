#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat getGamma(const arma::vec expbZ, const arma::vec X, const arma::vec T,
                   const arma::mat Z, const arma::vec wh)
{
  // expbZ - exp(beta %*% Z)
  // X - target lifetime
  // T - right truncation variable
  // Z - matrix with covariates in columns mXp
  // wh - weights=P(R>T_i)

  int p = Z.n_cols, m = Z.n_rows;
  arma::mat Gamma(p,p), S2_i(p,p), S1S1(p,p);
  double S0_i=0;
  arma::vec S1_i(p), tmp(m);

  Gamma.fill(0);
  for (int ii = 0; ii < m; ii++)
  {
    tmp = wh(ii)/wh % expbZ % as<arma::vec>(wrap(X>=X(ii)));
    S0_i = sum(tmp);
    for (int hr=0; hr<p; hr++) // number of covariates
    {
      S1_i(hr) = sum(tmp % Z.col(hr));

      for (int hc=0; hc<=hr; hc++) // number of covariates
      {
        S2_i(hr,hc)=sum(tmp % Z.col(hr) % Z.col(hc));
        S2_i(hc,hr)=S2_i(hr,hc);
        S1S1(hr,hc)=S1_i(hr) * S1_i(hc);
        S1S1(hc,hr)= S1S1(hr,hc);
        Gamma(hr,hc)=Gamma(hr,hc)+S2_i(hr,hc)/S0_i - S1S1(hr,hc)/std::pow(S0_i,2);
        Gamma(hc,hr)=Gamma(hr,hc);
      }
    }
  }
  return(Gamma);
}


// [[Rcpp::export]]
arma::mat getSigma_cpp(const arma::mat Z, const arma::mat TMP, const arma::vec S0,
                       const arma::mat xi, const arma::mat dMi_Tl)
{
  int p = Z.n_cols, m = Z.n_rows;
  arma::mat Sigma(p,p), S1(m,p), S1_S0(m,p);
  arma::rowvec e(m), iid(p), a(p);

  Sigma.fill(0);
  for (int ii = 0; ii < m; ii++)
  {
    for (int g=0; g<p; g++)
    {
      e = TMP.row(ii);
      S1(ii,g) = mean(e.t() % Z.col(g) );
    }
    S1_S0.row(ii) = S1.row(ii)/S0(ii);
  }


  for (int ii = 0; ii < m; ii++)
  {
    iid.fill(0);
    a.fill(0);

    for (int ll = 0; ll < m; ll++) //=j
    {
      if (ll==ii)
        iid = iid +(Z.row(ii) - S1_S0.row(ll)) * ( 1 - dMi_Tl(ll,ii));
      else
        iid = iid - (Z.row(ii) - S1_S0.row(ll)) * dMi_Tl(ll,ii);

      for (int k = 0; k < m; k++)
      {
        a=a + TMP(ll,k) * Z.row(k) *(xi(ii,ll)-xi(ii,k))/S0(ll);
      }
    }

    iid=iid+a/std::pow((double)m,2);
    Sigma = Sigma + iid.t() * iid;
  } // ii loop

  return(Sigma);
}


// [[Rcpp::export]]
arma::mat getVar(const arma::vec exp_bZ, const arma::vec X, const arma::vec T, const arma::mat Z,
                        const arma::vec wh)
{
  // X - event time, T - right truncation time
  int m =Z.n_rows;

  arma::mat Tk(m,m), Ti(m,m), exp_bZk(m,m), S_R_k(m,m),
  S_R_i(m,m), TMP(m,m), Rk(m,m), Ri(m,m), dMi_Tl(m,m), xi(m,m), B(m,m),
  C(m,m), S0_j(m,m), Qk(m,m), Gamma(m,m);
  arma::vec S0(m), Q_atR_v(m);
  arma::mat M1, M2, Sigma, Var;
  arma::rowvec Q_atT_v(m);

  B.fill(0);
  C.fill(0);
  Tk.each_row() = X.t();
  Ti.each_col() = X;
  exp_bZk.each_row() = exp_bZ.t();
  S_R_k.each_row() = wh.t();
  S_R_i.each_col() = wh;
  Rk.each_row() = T.t();
  Ri.each_col() = T;

  TMP = S_R_i/S_R_k % exp_bZk % (Tk>=Ti);

  for(int i = 0; i < m; i++)
  {
    S0(i) = mean(TMP.row(i));
  }

  S0_j.each_col()=S0;
  S0_j = m*S0_j;
  dMi_Tl = TMP/S0_j;
  M1 =  as<arma::mat>(wrap( (Ti <= Tk) % (Tk <= Ri) ));
  M2=as<arma::mat>(wrap((Ti <= Rk) % (Rk <= Ri)));
  for(int i = 0; i < m; i++)
  {
    Q_atT_v(i) = mean(M1.col(i));
    Q_atR_v(i) = mean(M2.col(i));
  }
  Qk.each_row() = Q_atT_v;
  xi=(Ri<=Tk)/Qk;
  for (int q = 0; q < m; q++)
  {
    C = C + as<arma::mat>(wrap((Ti <= T(q)) % (T(q)<= Ri) % (T(q)<=Tk)))/std::pow(Q_atR_v(q),2);
    B = B + as<arma::mat>(wrap((Ri <= X(q)) % (X(q)<= Tk)))/std::pow(Q_atT_v(q),2) - as<arma::mat>(wrap((Ri <= T(q)) % (T(q)<= Tk)))/std::pow(Q_atR_v(q),2);
  }
  xi = xi + (B-C)/m;

  Sigma = getSigma_cpp(Z, TMP, S0, xi, dMi_Tl);
  Gamma = getGamma(exp_bZ, X, T, Z, wh);
  Var = inv_sympd(Gamma) * Sigma * inv_sympd(Gamma);
  return(Var);
}



