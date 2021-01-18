// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include <Rcpp.h>
using namespace Rcpp;


static const double neg_inf = -std::numeric_limits<double>::infinity();


inline double logsum(const double &a, const double &b) {
  const double mx = std::max(a, b);
  if (mx==neg_inf) return a;
  return mx+log(exp(a-mx)+exp(b-mx));
}


/* Precompute joint probabilities jointp[i-1][n1] = P(b^i) */
std::vector<std::vector<double> > precompute_joint_p(NumericVector logprior, Progress &pbar, bool divideByBinom = true) {

  const int n = logprior.length()-1;
  
  std::vector<std::vector<double> > jointp;
  jointp.resize(n);
  for (int i = 1; i <= n; i++) {
    jointp[i-1].resize(i+1, neg_inf);
  }
  
  // Initialise for i=n
  
  if (logprior[0] != logprior[0] || logprior[n] != logprior[n])
    stop("logprior contains NaNs!\n");
  // NOTE: log(n choose 0) = 0, so we don't have to subtract in case of divideByBinom
  jointp[n-1][0] = logprior[0];
  jointp[n-1][n] = logprior[n];
  double logcnk = 0;  // log choose(n,k) starting at k=0
  
  // Use choose(n,n-k) = choose(n,k)
  for (int k = 1; k <= n/2.0; k++) {
    if (logprior[k] != logprior[k]) stop("logprior contains NaNs!\n");
    if (divideByBinom) { 
      logcnk = logcnk + log(n-k+1.0) - log((double)k);
      jointp[n-1][k] = logprior[k] - logcnk;
      jointp[n-1][n-k] = logprior[n-k] - logcnk;
    } else {
      jointp[n-1][k] = logprior[k];
      jointp[n-1][n-k] = logprior[n-k];
    }
  }
  
  // Recursively complete jointp
  for (int i = n-1; i >= 1; i--) {
    for (int n1 = 0; n1 <= i; n1++) {
      jointp[i-1][n1] = logsum(jointp[i][n1], jointp[i][n1+1]);
    }
    if (Progress::check_abort()) stop("User abort");
    pbar.increment();
  }
  
  return jointp;
}  





// [[Rcpp::export]]
NumericVector HierarchicalPriorC(NumericVector logphi, NumericVector logpsi, NumericVector logprior,
                                 bool showProgress = true, bool divideByBinom = true) {
  
  int n = logphi.length();
  if (logpsi.length()   != n)   stop("Lengths of logpsi and logphi disagree!");
  if (logprior.length() != n+1) stop("Lengths of data and prior disagree!");
  
    /* Progress bar:
   * First phase:  calculate jointp, takes n-1 steps
   * Second phase: forward pass,     takes n-1 steps
   * Third phase:  backward pass,    takes n-1 steps
   * Fourth phase: marginals,        takes n steps
   */
  Progress pbar(4*n-3, showProgress);
  
  std::vector<std::vector<double> > jointp = precompute_joint_p(logprior, pbar, divideByBinom);
  
  /* Forward pass */
  
  // pf[i][j][n1] = P(S_{i+1}, X^{i+1}) for S_{i+1} = (j,n1), 
  // where j in {0,1} signals if there is a slab at time i+1 and n1 in
  // {0,...,i+1} counts the number of slabs up to and including time i+1
  std::vector<std::vector<std::vector<double> > > pf;
  pf.resize(n);
  for (int i = 0; i < n; i++) {
    pf[i].resize(2);
    for (int j = 0; j < 2; j++) {
      pf[i][j].resize(i+2, neg_inf);
    }
    
  }
  
  // Initialize for i = 0
  pf[0][0][0] = logphi[0] + jointp[0][0];
  pf[0][1][1] = logpsi[0] + jointp[0][1];
  pf[0][0][1] = neg_inf;
  pf[0][1][0] = neg_inf;
  
  for (int i = 1; i < n; i++) {
    for (int n1 = 0; n1 <= i+1; n1++) {
      // Case that j=0 in round i+1 (can only happen if n1 < i+1)
      if (n1 < i+1) {
        pf[i][0][n1] = logphi[i] + jointp[i][n1] - jointp[i-1][n1] + (logsum(pf[i-1][0][n1],pf[i-1][1][n1]));
      } else {
        pf[i][0][n1] = neg_inf; // fill out unreachable case to avoid case-distinction in the other case
      }
      
      // Case that j=1 in round i+1 (can only happen if n1 > 0)
      if (n1 > 0) {
        pf[i][1][n1] = logpsi[i] + jointp[i][n1] - jointp[i-1][n1-1] + (logsum(pf[i-1][0][n1-1],pf[i-1][1][n1-1]));
      } else {
        pf[i][1][n1] = neg_inf; // fill out unreachable case to avoid case-distinction in the other case
      }
    }
    if (Progress::check_abort()) stop("User abort");
    pbar.increment();
  }
  
  /* Backward pass */
  
  // logarithmic representation:
  // pb[i][j][n1] = P(X_{i+2}^N | S_{i+1}) for S_{i+1} = (j,n1),
  // where j in {0,1} signals if there is a slab at time i+1 and n1 in
  // {0,...,i+1} counts the number of slabs up to and including time i+1
  std::vector<std::vector<std::vector<double> > > pb;
  pb.resize(n);
  for (int i = 0; i < n; i++) {
    pb[i].resize(2);
    for (int j = 0; j < 2; j++) {
      pb[i][j].resize(i+2, neg_inf);
    }
  }
  
  // Initialize for i = N-1
  for (int j = 0; j < 2; j++) {
    for (int n1 = 0; n1 <= n; n1++) {
      pb[n-1][j][n1] = 0;
    }
  }
  
  for (int i = n-2; i >= 0; i--) {
    for (int n1 = 0; n1 <= i+1; n1++) {
      // Cases j=0 and j=1 get the same value v
      double v = logsum(logphi[i+1] + jointp[i+1][n1]   - jointp[i][n1] + pb[i+1][0][n1],
                        logpsi[i+1] + jointp[i+1][n1+1] - jointp[i][n1] + pb[i+1][1][n1+1]);
      
      // Case that j=0 in round i+1 (can only happen if n1 < i+1)
      // fill out unreachable case to avoid case-distinction in computation of v
      pb[i][0][n1] = n1<i+1 ? v : neg_inf;
      
      // Case that j=1 in round i+1 (can only happen if n1 > 0)
      // fill out unreachable case to avoid case-distinction in computation of v
      pb[i][1][n1] = n1>0 ? v : neg_inf;
    }
    if (Progress::check_abort()) stop("User abort");
    pbar.increment();
  }
  
  /* Combine the forward and backward steps */
  
  // Compute P(b_{i+1}=1 | X^n) for all i
  NumericVector postb1(n);
  for (int i = 0; i < n; i++) {
    // Compute joint P(X^N, b_{i+1}) for b_{i+1}=0 and b_{i+1}=1
    double pjoint_spikes = neg_inf;
    double pjoint_slabs  = neg_inf;
    for (int n1 = 0; n1 <= i+1; n1++) {
      // NB Pr(X^n,S_{i+1}) = pf[i][j][n1] * pb[i][j][n1]
      pjoint_spikes = logsum(pjoint_spikes,pf[i][0][n1] + pb[i][0][n1]);
      pjoint_slabs  = logsum(pjoint_slabs,pf[i][1][n1]  + pb[i][1][n1]);
    }

    postb1[i] = exp(pjoint_slabs-logsum(pjoint_spikes, pjoint_slabs));
    if (Progress::check_abort()) stop("User abort");
    pbar.increment();
  }
  
  return postb1;
}


// [[Rcpp::export]]
NumericVector DiscreteSpikeSlabPriorC(NumericVector logphi, NumericVector logpsi, NumericVector grid, NumericVector logGridPrior, 
                                      bool showProgress = true) {
  int n = logphi.length();
  int ndisc = grid.length(); // nr of discretization points
  
  if (logpsi.length() != n) stop("Lengths of logpsi and logphi disagree!");
  if (logGridPrior.length() != ndisc) stop("Lengths of grid and logGridPrior disagree!");
  
  /* Progress bar:
   * Phase 1: posterior on alpha grid points, n steps
   * Phase 2: marginal probabilities,         n steps
   */
  
  Progress pbar(2*n, showProgress);
  
  
  /* Compute posterior pf on alpha */
  
  // We will recursively maintain
  // pf[j] = P(alpha_j, X^i) for i = 0, ..., n
  std::vector<double> pf;
  pf.resize(ndisc, neg_inf);
  
  // Initialize to prior
  for (int j = 0; j < ndisc; j++) {
    pf[j] = logGridPrior[j];
  }
  
  // Multiply pf by likelihood of the data
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < ndisc; j++) {
      double logspike = log(1.0 - grid[j]) + logphi[i];
      double logslab = log(grid[j]) + logpsi[i];
      pf[j] = pf[j] + logsum(logspike, logslab);
    }
    if (Progress::check_abort()) stop("User abort");
    pbar.increment();
  }
  
  // Compute normalization: log sum_j exp(pf[j]) = P(X^n)
  double lognorm = neg_inf;
  for (int j=0; j < ndisc; j++) {
    lognorm = logsum(lognorm, pf[j]);
  }
  
  /* Compute marginal posterior probabilities */
  
  NumericVector postb1(n);
  for (int i = 0; i < n; i++) {
    double acc = neg_inf;
    for (int j = 0; j < ndisc; j++) {
      // Compute condb1 = P(b_i = 1 | alpha_j, X_i)
      //  = -log(1 + (1-alpha_j)phi_i/(alpha_j psi_i) )
      double logspike = log(1.0 - grid[j]) + logphi[i];
      double logslab = log(grid[j]) + logpsi[i];
      double condb1 = -logsum(0, logspike - logslab);
      
      acc = logsum(acc, condb1 + pf[j] - lognorm);
    }
    postb1[i] = exp(acc);
    if (Progress::check_abort()) stop("User abort");
    pbar.increment();
  }
  
  return postb1;
}
