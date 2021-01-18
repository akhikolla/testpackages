#include "dblpm.h"

void dblpm::Print()
{
  if (debug) Rcpp::Rcout << "dblpm::Print has been called" << std::endl;
  Rcpp::Rcout << "\n\nClass dblpm\n\nDimensions:\n" << std::endl;
  Rcpp::Rcout << "T\t=\t" << T << std::endl;
  Rcpp::Rcout << "N\t=\t" << N << std::endl;
  Rcpp::Rcout << "M\t=\t" << M << std::endl;
  Rcpp::Rcout << "L\t=\t" << L << std::endl;
  Rcpp::Rcout << "D\t=\t" << D << std::endl << std::endl << std::endl;
  Rcpp::Rcout << "Hyperparameters:\n" << std::endl;
  Rcpp::Rcout << "taux\t=\t" << taux << std::endl;
  Rcpp::Rcout << "delta\t=\t" << delta << std::endl;
  Rcpp::Rcout << "aw\t=\t" << aw << "\t\t\t\tbw\t=\t" << bw << std::endl;
  Rcpp::Rcout << "agamma\t=\t" << agamma << "\t\t\t\tbgamma\t=\t" << bgamma << std::endl;
  Rcpp::Rcout << "abeta\t=\t" << abeta << "\t\t\t\tbbeta\t=\t" << bbeta << std::endl;
  Rcpp::Rcout << "\n\nPrecisions:\n" << std::endl;
  Rcpp::Rcout << "tauw\t=\t" << tauw << "\t\t\ttauw0\t=\t" << tauw0 << std::endl;
  Rcpp::Rcout << "tauga\t=\t" << taugamma << "\t\t\t\ttauga0\t=\t" << taugamma0 << std::endl;
  Rcpp::Rcout << "taube\t=\t" << taubeta << "\t\t\t\ttaube0\t=\t" << taubeta0 << std::endl << std::endl << std::endl;
  Rcpp::Rcout << "Likelihood parameters:" << std::endl;
  Rcpp::Rcout << "\nx\t=\n\n";
  x.print();
  Rcpp::Rcout << "\nw\t=\n\n";
  w.print();
  Rcpp::Rcout << "\nw0_ss = " << w0_ss << "\t\t\tw_innovation_ss = " << w_innovation_ss << "\n\n";
  Rcpp::Rcout << "\ngamma\t=\n\n";
  gamma.t().print();
  Rcpp::Rcout << "\ngamma_innovation_ss = " << gamma_innovation_ss << "\n\n";
  Rcpp::Rcout << "\nbeta\t=\n\n";
  beta.t().print();
  Rcpp::Rcout << "\nbeta_innovation_ss = " << beta_innovation_ss << "\n\n";
  Rcpp::Rcout << "\n\nAdjacency cube:\n\n";
  y.print();
  Rcpp::Rcout << "\nEdgelist:\n\n";
  edgelist.print();
  Rcpp::Rcout << "\n\nOut-degrees:\n\n";
  out_degrees.print();
  Rcpp::Rcout << "\n\nTotal out-degrees:\n\n";
  out_tot_degrees.t().print();
  Rcpp::Rcout << "\n\nIn-degrees:\n\n";
  in_degrees.print();
  Rcpp::Rcout << "\n\nTotal in-degrees:\n\n";
  in_tot_degrees.t().print();
  Rcpp::Rcout << "\n\ni_activity_table (directors' activity table):\n\n";
  i_activity_table.print();
  Rcpp::Rcout << "\n\nLists of active directors per time frame\n" << std::endl;
  unsigned int t;
  for (t=0; t<T; ++t)
  {
    Rcpp::Rcout << "i_activity_list_" << t << ":\t";
    i_activity_list.at(t).t().print();
  }
  Rcpp::Rcout << "\n\nj_activity_table (boards' activity table):\n\n";
  j_activity_table.print();
  Rcpp::Rcout << "\n\nj_first_activity (boards' first time frame of activity)\t\n\n";
  j_first_activity.t().print();
  Rcpp::Rcout << "\n\nj_last_activity (boards' last time frame of activity)\t\n\n";
  j_last_activity.t().print();
  Rcpp::Rcout << "\n\n" << N_active << " active directors:\t";
  i_active.t().print();
  Rcpp::Rcout << M_active << " active boards:\t";
  j_active.t().print();
  Rcpp::Rcout << "\n\n\n\n";
  if (debug) Rcpp::Rcout << "dblpm::Print has terminated" << std::endl;
  
}
