#include "expsbm.h"

void expsbm::Print()
{
  std::ostringstream strs;
  strs << "\n\nclass expbm\n";
  strs << "\nN\t=\t" << N << "\n";
  strs << "\ndirected\t=\t" << directed << "\n";
  strs << "\ntrunc\t=\t" << trunc << "\n";
  // strs << "\nEdgelist:\n";
  // edgelist.print(strs);
  // strs << "\nA:\n";
  // A.print(strs);
  // strs << "\nX:\n";
  // X.print(strs);
  strs << "\nW:\n";
  W.print(strs);
  strs << "\nK\t=\t" << K << "\n";
  strs << "\nZ:\n";
  Z.print(strs);
  strs << "\nlambda:\n";
  lambda.t().print(strs);
  strs << "\nmu:\n";
  mu.print(strs);
  strs << "\nnu:\n";
  nu.print(strs);
  strs << "\nA1:\n";
  A1.print(strs);
  strs << "\nA0:\n";
  A0.print(strs);
  strs << "\nX1:\n";
  X1.print(strs);
  strs << "\nX0:\n";
  X0.print(strs);
  strs << "\nL_mu:\n";
  L_mu.print(strs);
  strs << "\nL_nu:\n";
  L_nu.print(strs);
  strs << "\neta:\n";
  eta.print(strs);
  strs << "\nzeta:\n";
  zeta.print(strs);
  strs << "\nelbo_value\t=\t" << elbo_value << "\n";
  Rcpp::Rcout << strs.str() << std::endl << std::endl;
}

