// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix nn_matchC(const IntegerVector& treat,
                        const IntegerVector& ord,
                        const IntegerVector& ratio,
                        const bool& replace,
                        const LogicalVector& discarded,
                        const Nullable<NumericVector>& distance_ = R_NilValue,
                        const Nullable<IntegerVector>& exact_ = R_NilValue,
                        const Nullable<double>& caliper_dist_ = R_NilValue,
                        const Nullable<NumericVector>& caliper_covs_ = R_NilValue,
                        const Nullable<NumericMatrix>& calcovs_covs_mat_ = R_NilValue,
                        const Nullable<NumericMatrix>& mah_covs_ = R_NilValue,
                        const Nullable<NumericMatrix>& mahSigma_inv_ = R_NilValue,
                        const bool& disl_prog = false)
  {

  // Initialize

  NumericVector distance, caliper_covs;
  double caliper_dist;
  NumericMatrix calcovs_covs_mat, mah_covs, mahSigma_inv, mah_covs_c;
  IntegerVector exact;

  Environment pkg = Environment::namespace_env("stats");
  Function mah = pkg["mahalanobis"];

  bool use_exact = false;
  bool use_caliper_dist = false;
  bool use_caliper_covs = false;
  bool use_mah_covs = false;

  int n = treat.size();
  int n1 = sum(treat);
  int n0 = n - n1;

  IntegerVector ind = Range(0, n - 1);
  IntegerVector ind1 = ind[treat == 1];

  int max_rat = max(ratio);

  // Output matrix with sample indices of C units
  IntegerMatrix mm(n1, max_rat);
  mm.fill(NA_INTEGER);
  rownames(mm) = as<CharacterVector>(treat.names())[ind1];

  // Store who has been matched; exclude discarded
  LogicalVector matched = clone(discarded);

  // Indices of units to form control pool
  IntegerVector c_pool = ind;

  int t, t_ind, min_ind, c_chosen, num_eligible, cal_len, t_rat;
  double dt, cal_var_t;

  NumericVector cal_var, cal_diff, ps_diff, diff, mah_covs_t, mah_covs_col,
                match_distance(n0);

  LogicalVector c_eligible(n);

  IntegerVector c_available(n0), indices(n0);

  if (distance_.isNotNull()) {
    distance = distance_;
  }
  if (exact_.isNotNull()) {
    exact = exact_;
    use_exact = true;
  }
  if (caliper_dist_.isNotNull()) {
    caliper_dist = as<double>(caliper_dist_);
    use_caliper_dist = true;
    ps_diff = NumericVector(n);
  }
  if (caliper_covs_.isNotNull()) {
    caliper_covs = caliper_covs_;
    use_caliper_covs = true;
    cal_len = caliper_covs.size();
    cal_diff = NumericVector(n);
  }
  if (calcovs_covs_mat_.isNotNull()) {
    calcovs_covs_mat = as<NumericMatrix>(calcovs_covs_mat_);
  }
  if (mah_covs_.isNotNull()) {
    mah_covs = as<NumericMatrix>(mah_covs_);
    use_mah_covs = true;
  } else {
    ps_diff = NumericVector(n);
  }
  if (mahSigma_inv_.isNotNull()) {
    mahSigma_inv = as<NumericMatrix>(mahSigma_inv_);
  }

  bool ps_diff_assigned = false;

  //progress bar
  int prog_length;
  if (replace) prog_length = n1 + 1;
  else prog_length = max_rat*n1 + 1;
  Progress p(prog_length, disl_prog);

  //Counters
  int rat, i, x, j;

  //Matching
  for (rat = 0; rat < max_rat; ++rat) {
    for (i = 0; i < n1; ++i) {

      p.increment();

      if (all(matched).is_true()) {
        break;
      }

      t = ord[i] - 1;     // index among treated
      t_rat = ratio[t];

      if (t_rat < rat + 1) {
        continue;
      }

      t_ind = ind1[t]; // index among sample

      if (discarded[t_ind]) {
        continue;
      }

      c_eligible = treat == 0; // index among sample

      c_eligible[matched] = false;

      if (use_exact) {
        c_eligible[exact != exact[t_ind]] = false;
      }

      if (any(c_eligible).is_false()) {
        continue;
      }

      if (use_caliper_dist) {
        dt = distance[t_ind];
        diff = Rcpp::abs(as<NumericVector>(distance[c_eligible]) - dt);

        ps_diff[c_eligible] = diff;
        ps_diff_assigned = true;

        c_eligible[ps_diff > caliper_dist] = false;

        if (any(c_eligible).is_false()) {
          continue;
        }
      }

      if (use_caliper_covs) {
        for (x = 0; (x < cal_len) && any(c_eligible).is_true(); ++x) {
          cal_var = calcovs_covs_mat( _ , x );

          cal_var_t = cal_var[t_ind];

          diff = Rcpp::abs(as<NumericVector>(cal_var[c_eligible]) - cal_var_t);

          cal_diff[c_eligible] = diff;

          c_eligible[cal_diff > caliper_covs[x]] = false;
        }

        if (any(c_eligible).is_false()) {
          continue;
        }
      }

      //Compute distances among eligible
      c_available = as<IntegerVector>(c_pool[c_eligible]);
      num_eligible = c_available.size();

      if (use_mah_covs) {
        mah_covs_c = NumericMatrix(num_eligible, mah_covs.ncol());
        for (j = 0; j < mah_covs.ncol(); ++j) {
          mah_covs_col = mah_covs.column(j);
          mah_covs_c(_,j) = as<NumericVector>(mah_covs_col[c_available]);
        }
        mah_covs_t = mah_covs( t_ind , _ );
        match_distance = sqrt(as<NumericVector>(mah(mah_covs_c, mah_covs_t, mahSigma_inv, true))); //mahalanobis in R

      } else {
        if (ps_diff_assigned) {
          match_distance = ps_diff[c_eligible];
          ps_diff_assigned = false; // reset for next iter
        } else {
          dt = distance[t_ind];
          match_distance = Rcpp::abs(as<NumericVector>(distance[c_eligible]) - dt);
        }
      }

      if (!replace) {
        min_ind = which_min(match_distance);
        c_chosen = c_available[min_ind];

        mm( t , rat ) = c_chosen + 1; // + 1 because C indexing starts at 0 but mm is sent to R

        matched[c_chosen] = true;
      }
      else {
        if (num_eligible <= t_rat) {
          for (j = 0; j < num_eligible; ++j) {
            mm( t , j ) = c_available[j] + 1;
          }
        }
        else {
          //When matching w/ replacement, get t_rat closest control units
          indices = Range(0, num_eligible - 1);

          std::partial_sort(indices.begin(), indices.begin() + t_rat, indices.end(),
                            [&match_distance](int k, int j) {return match_distance[k] < match_distance[j];});

          for (j = 0; j < t_rat; ++j) {
            min_ind = indices[j];
            mm( t , j ) = c_available[min_ind] + 1;
          }
        }
      }
    }

    if (replace) break;
  }

  p.update(prog_length);

  return mm;
}
