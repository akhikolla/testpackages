#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame method1_cpp(const int init_pop_size,
                      const IntegerVector& vaccinations,
                      const IntegerVector& cases,
                      const NumericVector& ve) {
  int nt = cases.size();
  IntegerVector pops(nt), popn(nt), cases_novac(nt), avert(nt);
  NumericVector pvac(nt), pflu(nt), vc_lag(nt);
  for (int i = 0; i < nt; i++)
    pvac[i] = vaccinations[i] / double(init_pop_size);
  vc_lag[0] = pvac[0] / 2;
  pops[0] = R::fround(init_pop_size * (1 - vc_lag[0] * ve[0]), 0);
  pflu[0] = double(cases[0]) / pops[0];
  popn[0] = init_pop_size - cases[0];
  cases_novac[0] = R::fround(pflu[0] * popn[0], 0);
  avert[0] = cases_novac[0] - cases[0];
  for (int i = 1; i < nt; i++) {
    vc_lag[i] = (pvac[i] + pvac[i - 1]) / 2;
    pops[i] = R::fround(
      (pops[i - 1] - cases[i - 1]) * (1 - vc_lag[i] * ve[i]), 0
    );
    pflu[i] = double(cases[i]) / pops[i];
    popn[i] = popn[i - 1] - cases_novac[i - 1];
    cases_novac[i] = R::fround(pflu[i] * popn[i], 0);
    avert[i] = cases_novac[i] - cases[i];
  }
  DataFrame method1 = DataFrame::create(
    _["cases"] = cases,
    _["vaccinations"] = vaccinations,
    _["ve"] = ve,
    _["pvac"] = pvac,
    _["vc_lag"] = vc_lag,
    _["pops"] = pops,
    _["pflu"] = pflu,
    _["popn"] = popn,
    _["cases_novac"] = cases_novac,
    _["avert"] = avert
  );
  return method1;
}

// [[Rcpp::export]]
DataFrame method3_cpp(const int init_pop_size,
                      const IntegerVector& vaccinations,
                      const IntegerVector& cases,
                      const NumericVector& ve) {
  int nt = cases.size();
  IntegerVector b(nt), A(nt), C(nt), D(nt), E(nt), F(nt),
  popn(nt), cases_novac(nt), avert(nt);
  NumericVector pvac(nt), pflu(nt);
  pflu[0] = cases[0] / double(init_pop_size);
  pvac[0] = vaccinations[0] / double(init_pop_size);
  b[0] = R::fround(init_pop_size * pvac[0], 0);
  A[0] = R::fround(init_pop_size * (1 - pflu[0]) - b[0], 0);
  C[0] = R::fround(b[0] * (1 - ve[0]), 0);
  D[0] = R::fround(b[0] * ve[0], 0);
  E[0] = R::fround(init_pop_size * pflu[0], 0);
  F[0] = 0;
  cases_novac[0] = R::fround(init_pop_size * pflu[0], 0);
  popn[0] = init_pop_size - cases_novac[0];
  avert[0] = cases_novac[0] - cases[0];
  for (int i = 1; i < nt; i++) {
    pflu[i] = double(cases[i]) / (A[i - 1] + C[i - 1]);
    pvac[i] = double(vaccinations[i]) / (A[i - 1] + E[i - 1]);
    b[i] = R::fround(A[i - 1] * pvac[i], 0);
    A[i] = R::fround(A[i - 1] * (1 - pflu[i]) - b[i], 0);
    C[i] = R::fround(C[i - 1] * (1 - pflu[i]) + b[i] * (1 - ve[i]), 0);
    D[i] = R::fround(D[i - 1] + b[i] * ve[i], 0);
    E[i] = R::fround(E[i - 1] * (1 - pvac[i]) + A[i - 1] * pflu[i], 0);
    F[i] = R::fround(F[i - 1] + C[i - 1] * pflu[i] + E[i - 1] * pvac[i], 0);
    cases_novac[i] = R::fround(pflu[i] * popn[i - 1], 0);
    popn[i] = popn[i - 1] - cases_novac[i];
    avert[i] = cases_novac[i] - cases[i];
  }
  DataFrame method3 = DataFrame::create(
    _["cases"] = cases,
    _["vaccinations"] = vaccinations,
    _["ve"] = ve,
    _["b"] = b,
    _["A"] = A,
    _["C"] = C,
    _["D"] = D,
    _["E"] = E,
    _["F"] = F,
    _["pvac"] = pvac,
    _["pflu"] = pflu,
    _["popn"] = popn,
    _["cases_novac"] = cases_novac,
    _["avert"] = avert
  );
  return method3;
}
