#include "outlier_utils.h"

/*****************************************
 * In:
 * - x: Vector of strings
 * Out: A namedvector of counts of all unique
 *      elements in x
 ***************************************/
// [[Rcpp::export]]
RIV count_unique(VS  x) { // std::unordered_map<std::string, int>
  std::map<std::string, int> tab;
  int n = x.size();
  for (int i = 0; i < n; i++) {
    auto s = x[i];
    tab[s]++;
  }
  return Rcpp::wrap(tab);
}

// [[Rcpp::export]]
VS matpr(Rcpp::CharacterMatrix A) {
  // Paste rows in a character matrix
  int n = A.nrow();
  VS  x(n);
  for( int i = 0; i < n; i++ ) {
    auto row = A(i, Rcpp::_);
    std::string s;
    s = std::accumulate(row.begin(),row.end(), s);
    x[i] = s;
  }
  return x;
}

/*****************************************
 * In:
 * - A: A character matrix with all a-variables (and only them)
 * -    A needs to have dimnames = list(NULL, colnames)
 * Out: The a-marginal table with attribute = variable names
 ***************************************/
// [[Rcpp::export]]
RIV n_a(RCM & A) {
  VS x = matpr(A);
  auto na = count_unique(x);
  // Rcpp::List Delta_A = A.attr("dimnames"); // Use colnames(A) ?
  na.attr("vars") = Rcpp::colnames(A); //Delta_A[1];
  return na;
}

// [[Rcpp::export]]
int na_ya(RIV & na, std::string ya) {
  RCV na_names = na.names();
  bool in_na = std::find(na_names.begin(), na_names.end(), ya) != na_names.end();
  if (in_na ) return na[ya];
  else return 0;
}

/*****************************************
 * In:
 * - na: The a-marginal table
 * -  b: Named vector with positions of the b's
 * Out:  The b'th slice of the a-marginal table
 ***************************************/
// [[Rcpp::export]]
RIV n_b(RIV & na, RIV & b) {

  /**
   *  - Assert that max(b) <= nv. Otherwise an error is produced.
   */
  RIV n_b_out;
  VS  cells     = na.names();
  VS  vars      = na.attr("vars");
  int nc        = cells.size();
  int nv        = vars.size();
  int nb        = b.size();
  VS  names_b   = b.names();
  VI  b_sorted  = Rcpp::as<VI>(b);
  VS  fix_b(nv);
  
  // ensures correct deletion of b index
  std::sort(b_sorted.rbegin(), b_sorted.rend()); 
  std::fill(fix_b.begin(), fix_b.end(), ".");

  // store the conditional information
  VS  vars_cond_names;
  RCV vars_cond = b.names();
  
  for (int i = 0; i < b.size(); i++) {
    int bi_idx = b[i] - 1;
    fix_b[bi_idx] = names_b[i];
    vars_cond_names.push_back(vars[bi_idx]);
  }

  vars_cond.names() = vars_cond_names;
  
  std::string fb;
  fb = std::accumulate(fix_b.begin(),fix_b.end(), fb);
  std::regex look(fb);

  for (int c = 0; c < nc; c++) {
    if( std::regex_match(cells[c], look) ) {
      std::string new_cell = cells[c];
      int new_val = na[c];
      for (int r = 0 ; r < nb; r++) {
    	new_cell.erase(new_cell.begin() + b_sorted[r] - 1); 
      }
      n_b_out[new_cell] = new_val;
    }    
  }

  VS new_vars = vars;
  for (int j = 0; j < nb; j++) {
    new_vars.erase(new_vars.begin() + b_sorted[j] - 1); 
  }
  n_b_out.attr("vars") = new_vars;
  n_b_out.attr("vars_cond") = vars_cond;
  return n_b_out;
}


// [[Rcpp::export]]
VD subtract_one(VD x) {
  std::transform(x.begin(), x.end(), x.begin(), [](double i) {
    return i - 1.0;
  });
  return x;
}

// [[Rcpp::export]]
VD Gx_(VD x) {
  std::transform(x.begin(), x.end(), x.begin(), [](double j) {
    if ( j > 0 ) {
      return j*std::log(j);
    } else {
      return 0.0;
    }
  });
  return x;
}

// [[Rcpp::export]]
VD Hx_(VD x) {
  VD v1 = Gx_(subtract_one(x)), v2 = Gx_(x), v;
  std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(v),
    [](double i, double j) { return i-j; });
  return v;
}

// [[Rcpp::export]]
RCM subM( RCM & A, RCV & x ) {
  RCV A_Delta = Rcpp::colnames(A);
  RIV sub_idx = Rcpp::match(x, A_Delta);
  RCM A_sub(A.nrow(), x.size());
  for (int i = 0; i < sub_idx.size(); i++) {
    A_sub(Rcpp::_ , i ) = A(Rcpp::_, sub_idx[i] - 1);
  }
  Rcpp::colnames(A_sub) = x;
  return A_sub;
}

/*****************************************
 * In:
 * - A  : Character matrix (the data)
 * - am : A list of variables (typically cliques or separators with RIP ordering)
 * Out: A list with sparse marginal tables corresponding to variables in list am
 ***************************************/
// [[Rcpp::export]]
RL a_marginals( RCM A, RL & am ) {
  int n = am.size();
  RL  out(n);
  for (int i = 0; i < n; i++) {
    if ( am[i] == R_NilValue ) {
      out[i] = (am[i]);      
    }
    else {
      RCV z = Rcpp::as<RCV>(am[i]);
      RCM A_sub = subM(A, z);
      out[i] = n_a(A_sub);
    }
  }
  return out;
}

/*****************************************
 * In:
 * - y : A named vector (named according to data)
 * - C_marginals : Clique marginal tables (use a_marginals function)
 * - S_marginals : Separator marginal tables (use a_marginals function)
 * Out: The affine value T(y) of -2 log likelihood-ratio
 ***************************************/
// [[Rcpp::export]]
double TY(RCV y, RL & C_marginals, RL & S_marginals) {
  int nC = C_marginals.size();
  RCV y_names = y.names();
  VD  CS(nC), SS(nC); // SS[0] initialize as zero according to NULL

  // The first clique
  RIV nC0       = C_marginals[0];
  RCV nC0_Delta = nC0.attr("vars");
  RIV yC_idx    = Rcpp::match(nC0_Delta, y_names);
  RCV yC0       = y[yC_idx - 1];  // -1 to account for the R side
  std::string yC0_;
  yC0_ = std::accumulate(yC0.begin(), yC0.end(), yC0_);
  CS[0] = na_ya(nC0, yC0_); //nC0[yC0_];

  for (int i = 1; i < nC; i++) {
    // Cliques
    RIV nCi       = C_marginals[i];
    RCV nCi_Delta = nCi.attr("vars");
    RIV yC_idx    = Rcpp::match(nCi_Delta, y_names);
    RCV yCi       = y[yC_idx - 1];
    std::string yCi_;
    yCi_  = std::accumulate(yCi.begin(), yCi.end(), yCi_);
    CS[i] = na_ya(nCi, yCi_); // nCi[yCi_];

    // Separators
    RIV nSi = S_marginals[i];
    // Handling the empty separator = M = |n| with no vars attribute
    if ( Rf_isNull(nSi.attr("vars")) ) {
      SS[i] = nSi[0]; // 0.0;
    }
    else {
      RCV nSi_Delta = nSi.attr("vars");
      RIV yS_idx    = Rcpp::match(nSi_Delta, y_names);
      RCV ySi       = y[yS_idx - 1];
      std::string ySi_;
      ySi_  = std::accumulate(ySi.begin(), ySi.end(), ySi_);
      SS[i] = na_ya(nSi, ySi_); // nSi[ySi_];
    }
  }
  VD H_CS = Hx_(CS), H_SS = Hx_(SS);
  double sum_HCS = std::accumulate(H_CS.begin(), H_CS.end(), 0.0);
  double sum_HSS = std::accumulate(H_SS.begin(), H_SS.end(), 0.0);
  // Rcpp::Rcout << "HC - HS = " << sum_HCS << " - " << sum_HSS << "\n";
  return sum_HCS - sum_HSS;
}
