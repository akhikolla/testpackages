// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

using namespace arma;
using namespace Rcpp;

//' Simulate a sample path from an endpoint conditioned CTMC by modified
//' rejection sampling.
//'
//' @param a,b States at the interval endpoints, provided as integers
//'    corresponding to rows of the CTMC rate matrix.
//' @param t0,t1 times of the interval endpoints
//' @param Q CTMC rate matrix
//'
//' @return matrix whose first column is the sequence of transition times
//' bookended by interval endpoints, and whose second column is the sequence of
//' states
// [[Rcpp::export]]
arma::mat sample_path_mr(const int a, const int b, const double t0, const double t1, const Rcpp::NumericMatrix& Q) {

        // Get the number of states and initialize vector of states
        int n_states = Q.nrow();
        Rcpp::IntegerVector states = Rcpp::seq_len(n_states);

        // Initialize booleans for whether to keep simulating and whether a
        // valid path has been obtained.
        bool valid_path = false;

        // Initialize objects for storing the sequences of times and states
        std::vector<double> time_vec;
        std::vector<int> state_vec;

        // Sample paths until a valid path has been obtained
        while(valid_path == false) {

                // Set boolean to initiate forward sampling
                bool keep_going = true;

                // insert the initial time and state
                time_vec.push_back(t0);
                state_vec.push_back(a);

                // set the current time and state
                Rcpp::NumericVector cur_time(1, t0);
                Rcpp::IntegerVector cur_state(1, a);
                double cur_rate = -Q(cur_state[0] - 1, cur_state[0] - 1);

                // get the state transition probabilities
                Rcpp::NumericVector state_probs = pmax(Q(cur_state[0] - 1, _ ), 0);

                // If the beginning and end states don't match, sample first transition
                if(a != b) {

                        // sample the first transition time
                        cur_time += -log(1 - Rcpp::runif(1, 0, 1) * (1 - exp(-(t1-t0) * cur_rate))) / cur_rate;

                        // sample the next state
                        cur_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs);

                        // update the rate of transition out of the new state
                        // and update the state transition probabilities
                        cur_rate  = -Q(cur_state[0] - 1, cur_state[0] - 1);
                        state_probs = pmax(Q(cur_state[0] - 1, _ ), 0);

                        // Insert the next state and transition time into the
                        // appropriate vectors
                        time_vec.push_back(cur_time[0]);
                        state_vec.push_back(cur_state[0]);
                }

                // Proceed with forward sampling algorithm
                while(keep_going == true) {

                        // check if the state is an absorbing state
                        if(is_true(all(state_probs == 0))) {

                                // stop sampling forward
                                keep_going = false;

                                if(cur_state[0] == b) {
                                        valid_path = true;
                                } else {
                                        valid_path = false;
                                        time_vec.clear();
                                        state_vec.clear();
                                }

                                break;
                        }

                        // Sample the next transition time
                        cur_time += Rcpp::rexp(1, cur_rate);

                        // If the next time is after the right endpoint, stop
                        // sampling and determine if the path is valid
                        if(cur_time[0] > t1) {

                                // Stop forward sampling
                                keep_going = false;

                                // Determine if the path is valid
                                if(cur_state[0] == b) {
                                        valid_path = true;

                                } else {
                                        valid_path = false;
                                        time_vec.clear();
                                        state_vec.clear();
                                }

                                // If the next transition occurs before the right
                                // endpoint, sample the next state
                        } else {
                                // sample the next state
                                cur_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs);

                                // update the rate of transition out of the new state
                                // and update the state transition probabilities
                                cur_rate  = -Q(cur_state[0] - 1, cur_state[0] - 1);
                                state_probs = pmax(Q(cur_state[0] - 1, _ ), 0);

                                // update the state and time vectors
                                time_vec.push_back(cur_time[0]);
                                state_vec.push_back(cur_state[0]);
                        }
                }
        }

        // Add the time and state at the right endpoint
        time_vec.push_back(t1);
        state_vec.push_back(b);

        // Collect the sequences of states and times into a matrix
        arma::mat path(time_vec.size(), 2);;
        path.col(0) = arma::conv_to<arma::colvec>::from(time_vec);
        path.col(1) = arma::conv_to<arma::colvec>::from(state_vec);

        return path;
}

//' Simulate a sample path from an endpoint conditioned CTMC by uniformization.
//'
//' @param a,b States at the interval endpoints, provided as integers
//'    corresponding to rows of the CTMC rate matrix.
//' @param t0,t1 times of the interval endpoints
//' @param Q CTMC rate matrix
//'
//' @return matrix whose first column is the sequence of transition times
//' bookended by interval endpoints, and whose second column is the sequence of
//' states
// [[Rcpp::export]]
arma::mat sample_path_unif(const int a, const int b, const double t0, const double t1, const arma::mat& Q) {

        // Get the number of states and initialize vector of states
        int n_states = Q.n_rows;
        Rcpp::IntegerVector states = Rcpp::seq_len(n_states);

        // Get the length of the interval and the largest diagonal element of Q
        double T = t1 - t0;
        double m = max(abs(Q.diag()));

        // Construct the transition probability matrix and extract the a,b elem.
        arma::mat tpm = arma::expmat(Q * T);
        double p_ab = tpm(a-1, b-1);

        // Generate the auxilliary transition matrix
        arma::mat R = arma::eye(n_states, n_states) + Q/m;

        // Sample threshold for determining the number of states
        Rcpp::NumericVector n_thresh = Rcpp::runif(1, 0, 1);

        // Initialize number of jumps and conditional probability of n jumps
        int n_jumps = 0;
        double c_prob = exp(-m*T) * (a == b) / p_ab;

        // proceed with sampling by uniformization
        // first the cast when there are no jumps
        if(c_prob > n_thresh[0]) {

                // initialize matrix
                arma::mat path(2,2);

                // fill matrix
                path(0,0) = t0; path(0,1) = t1;
                path(1,0) = a;  path(1,1) = b;

                return path;

        } else {

                // increment the number of jumps and compute c_prob
                n_jumps += 1;
                c_prob += exp(-m*T) * pow(m*T, n_jumps) / Rcpp::internal::factorial(n_jumps) * R(a-1, b-1) / p_ab;

                // if there is exactly one jump
                if(c_prob > n_thresh[0]) {

                        // if the endpoints match, the only jump is a virtual one
                        if(a == b) {

                                // initialize matrix
                                arma::mat path(2,2);

                                // fill matrix
                                path(0,0) = t0; path(0,1) = t1;
                                path(1,0) = a;  path(1,1) = b;

                                return path;


                        // if the endpoints don't match, the jump is real
                        } else {

                                // initialize matrix
                                arma::mat path(3,2);

                                // fill matrix
                                path(0,0) = t0; path(0,1) = a;
                                path(1,0) = Rcpp::runif(1,t0,t1)[0]; path(1,1) = b;
                                path(2,0) = t1; path(2,1) = b;

                                return path;
                        }

                // Else, there are at least two jumps
                } else {

                        // Initialize a cube for storing powers of R
                        arma::cube R_pow(n_states, n_states, 8);
                        int R_pow_size = R_pow.n_slices;
                        R_pow.slice(0) = arma::eye(size(R));
                        R_pow.slice(1) = R;

                        // Initialize a vector for storing the transition probabilities
                        Rcpp::NumericVector state_probs(n_states);

                        // keep calculating the conditional probability of n jumps until
                        // the threshold is exceeded. store powers of R accordingly.
                        while(c_prob < n_thresh[0]) {

                                // increment the number of jumps
                                n_jumps += 1;

                                // check whether to add additional slices to the R_pow cube
                                if(n_jumps == R_pow_size) {
                                        R_pow.insert_slices(R_pow.n_slices, 8);
                                        R_pow_size = R_pow.n_slices;
                                }

                                // Add the new power of R to the cube and calculate c_prob
                                R_pow.slice(n_jumps) = R_pow.slice(n_jumps - 1) * R;
                                c_prob += exp(-m*T) * pow(m*T, n_jumps) / Rcpp::internal::factorial(n_jumps) * R_pow.slice(n_jumps)(a-1, b-1) / p_ab;
                        }

                        // initialize the path matrix
                        int path_nrows = n_jumps + 2;
                        arma::mat path(path_nrows, 2);
                        path(0,0) = t0;
                        path(0,1) = a;
                        path(path_nrows - 1, 0) = t1;
                        path(path_nrows - 1, 1) = b;

                        // transition times are uniformly distributed in the
                        // interval. Sample them, sort them, and place in path.
                        arma::colvec transitions = Rcpp::runif(n_jumps, t0, t1);
                        std::sort(transitions.begin(), transitions.end());
                        path(arma::span(1,n_jumps), 0) = transitions;

                        // Sample the states at the transition times
                        for(int j = 1; j < n_jumps + 1; ++j) {
                                state_probs = arma::trans(R(path(j-1, 1) - 1, span::all)) % R_pow.slice(n_jumps-j)(span::all, b-1) / R_pow(path(j-1, 1)-1, b-1, n_jumps-j+1);
                                path(j, 1) = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs)[0];
                        }

                        // Determine which transitions are virtual transitions
                        arma::vec keep_inds(path_nrows, arma::fill::ones);
                        for(int j = 1; j < n_jumps + 1; ++j) {
                                if(path(j, 1) == path(j-1, 1)) {
                                        keep_inds[j] = 0;
                                }
                        }

                        // create a matrix for the complete path without virtual jumps
                        arma::mat path_comp = path.rows(arma::find(keep_inds == 1));

                        return path_comp;
                }
        }
}

//' Simulate a sample path from an endpoint conditioned CTMC by uniformization
//' using pre-computed eigen-values (assumes that all eigenvalues are real).
//'
//' @param a,b States at the interval endpoints, provided as integers
//'    corresponding to rows of the CTMC rate matrix.
//' @param t0,t1 times of the interval endpoints
//' @param Q CTMC rate matrix
//' @param eigen_vals vector of eigen values of Q.
//' @param eigen_vecs matrix of eigen vectors of Q.
//' @param inverse_vecs inverse of the eigen vector matrix.
//'
//' @return matrix whose first column is the sequence of transition times
//' bookended by interval endpoints, and whose second column is the sequence of
//' states
// [[Rcpp::export]]

arma::mat sample_path_unif2(const int a, const int b, const double t0, const double t1, const arma::mat& Q, const arma::vec& eigen_vals, const arma::mat& eigen_vecs, const arma::mat& inverse_vecs) {

        // Get the number of states and initialize vector of states
        int n_states = Q.n_rows;
        Rcpp::IntegerVector states = Rcpp::seq_len(n_states);

        // Get the length of the interval and the largest diagonal element of Q
        double T = t1 - t0;
        double m = max(abs(Q.diag()));

        // Construct the transition probability matrix and extract the a,b elem.
        arma::mat tpm = eigen_vecs * arma::diagmat(exp(eigen_vals * T)) * inverse_vecs;
        double p_ab = tpm(a-1, b-1);

        // Generate the auxilliary transition matrix
        arma::mat R = arma::eye(n_states, n_states) + Q/m;

        // Sample threshold for determining the number of states
        Rcpp::NumericVector n_thresh = Rcpp::runif(1, 0, 1);

        // Initialize number of jumps and conditional probability of n jumps
        int n_jumps = 0;
        double c_prob = exp(-m*T) * (a == b) / p_ab;

        // proceed with sampling by uniformization
        // first the case when there are no jumps
        if(c_prob > n_thresh[0]) {

                // initialize matrix
                arma::mat path(2,2);

                // fill matrix
                path(0,0) = t0; path(0,1) = t1;
                path(1,0) = a;  path(1,1) = b;

                return path;

        } else {

                // increment the number of jumps and compute c_prob
                n_jumps += 1;
                c_prob += exp(-m*T) * pow(m*T, n_jumps) / Rcpp::internal::factorial(n_jumps) * R(a-1, b-1) / p_ab;

                // if there is exactly one jump
                if(c_prob > n_thresh[0]) {

                        // if the endpoints match, the only jump is a virtual one
                        if(a == b) {

                                // initialize matrix
                                arma::mat path(2,2);

                                // fill matrix
                                path(0,0) = t0; path(0,1) = t1;
                                path(1,0) = a;  path(1,1) = b;

                                return path;


                                // if the endpoints don't match, the jump is real
                        } else {

                                // initialize matrix
                                arma::mat path(3,2);

                                // fill matrix
                                path(0,0) = t0; path(0,1) = a;
                                path(1,0) = Rcpp::runif(1,t0,t1)[0]; path(1,1) = b;
                                path(2,0) = t1; path(2,1) = b;

                                return path;
                        }

                        // Else, there are at least two jumps
                } else {

                        // Initialize a cube for storing powers of R
                        arma::cube R_pow(n_states, n_states, 8);
                        int R_pow_size = R_pow.n_slices;
                        R_pow.slice(0) = arma::eye(size(R));
                        R_pow.slice(1) = R;

                        // Initialize a vector for storing the transition probabilities
                        Rcpp::NumericVector state_probs(n_states);

                        // keep calculating the conditional probability of n jumps until
                        // the threshold is exceeded. store powers of R accordingly.
                        while(c_prob < n_thresh[0]) {

                                // increment the number of jumps
                                n_jumps += 1;

                                // check whether to add additional slices to the R_pow cube
                                if(n_jumps == R_pow_size) {
                                        R_pow.insert_slices(R_pow.n_slices, 8);
                                        R_pow_size = R_pow.n_slices;
                                }

                                // Add the new power of R to the cube and calculate c_prob
                                R_pow.slice(n_jumps) = R_pow.slice(n_jumps - 1) * R;
                                c_prob += exp(-m*T) * pow(m*T, n_jumps) / Rcpp::internal::factorial(n_jumps) * R_pow.slice(n_jumps)(a-1, b-1) / p_ab;
                        }

                        // initialize the path matrix
                        int path_nrows = n_jumps + 2;
                        arma::mat path(path_nrows, 2);
                        path(0,0) = t0;
                        path(0,1) = a;
                        path(path_nrows - 1, 0) = t1;
                        path(path_nrows - 1, 1) = b;

                        // transition times are uniformly distributed in the
                        // interval. Sample them, sort them, and place in path.
                        arma::colvec transitions = Rcpp::runif(n_jumps, t0, t1);
                        std::sort(transitions.begin(), transitions.end());
                        path(arma::span(1,n_jumps), 0) = transitions;

                        // Sample the states at the transition times
                        for(int j = 1; j < n_jumps + 1; ++j) {
                                state_probs = arma::trans(R(path(j-1, 1) - 1, span::all)) % R_pow.slice(n_jumps-j)(span::all, b-1) / R_pow(path(j-1, 1)-1, b-1, n_jumps-j+1);
                                path(j, 1) = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs)[0];
                        }

                        // Determine which transitions are virtual transitions
                        arma::vec keep_inds(path_nrows, arma::fill::ones);
                        for(int j = 1; j < n_jumps + 1; ++j) {
                                if(path(j, 1) == path(j-1, 1)) {
                                        keep_inds[j] = 0;
                                }
                        }

                        // create a matrix for the complete path without virtual jumps
                        arma::mat path_comp = path.rows(arma::find(keep_inds == 1));

                        return path_comp;
                }
          }
}

//' Simulate a sample path from an endpoint conditioned CTMC by uniformization
//' using a pre-computed transition probability matrix.
//'
//' @param a,b States at the interval endpoints, provided as integers
//'    corresponding to rows of the CTMC rate matrix.
//' @param t0,t1 times of the interval endpoints
//' @param Q CTMC rate matrix
//' @param P CTMC transition probability matrix over the interval.
//'
//' @return matrix whose first column is the sequence of transition times
//' bookended by interval endpoints, and whose second column is the sequence of
//' states
// [[Rcpp::export]]

arma::mat sample_path_unif3(const int a, const int b, const double t0, const double t1, const arma::mat& Q, const arma::mat& P) {

          // Get the number of states and initialize vector of states
          int n_states = Q.n_rows;
          Rcpp::IntegerVector states = Rcpp::seq_len(n_states);

          // Get the length of the interval and the largest diagonal element of Q
          double T = t1 - t0;
          double m = max(abs(Q.diag()));

          // Construct the transition probability matrix and extract the a,b elem.
          double p_ab = P(a-1, b-1);

          // Generate the auxilliary transition matrix
          arma::mat R = arma::eye(n_states, n_states) + Q/m;

          // Sample threshold for determining the number of states
          Rcpp::NumericVector n_thresh = Rcpp::runif(1, 0, 1);

          // Initialize number of jumps and conditional probability of n jumps
          int n_jumps = 0;
          double c_prob = exp(-m*T) * (a == b) / p_ab;

          // proceed with sampling by uniformization
          // first the case when there are no jumps
          if(c_prob > n_thresh[0]) {

                    // initialize matrix
                    arma::mat path(2,2);

                    // fill matrix
                    path(0,0) = t0; path(0,1) = t1;
                    path(1,0) = a;  path(1,1) = b;

                    return path;

          } else {

                    // increment the number of jumps and compute c_prob
                    n_jumps += 1;
                    c_prob += exp(-m*T) * pow(m*T, n_jumps) / Rcpp::internal::factorial(n_jumps) * R(a-1, b-1) / p_ab;

                    // if there is exactly one jump
                    if(c_prob > n_thresh[0]) {

                              // if the endpoints match, the only jump is a virtual one
                              if(a == b) {

                                        // initialize matrix
                                        arma::mat path(2,2);

                                        // fill matrix
                                        path(0,0) = t0; path(0,1) = t1;
                                        path(1,0) = a;  path(1,1) = b;

                                        return path;


                                        // if the endpoints don't match, the jump is real
                              } else {

                                        // initialize matrix
                                        arma::mat path(3,2);

                                        // fill matrix
                                        path(0,0) = t0; path(0,1) = a;
                                        path(1,0) = Rcpp::runif(1,t0,t1)[0]; path(1,1) = b;
                                        path(2,0) = t1; path(2,1) = b;

                                        return path;
                              }

                              // Else, there are at least two jumps
                    } else {

                              // Initialize a cube for storing powers of R
                              arma::cube R_pow(n_states, n_states, 8);
                              int R_pow_size = R_pow.n_slices;
                              R_pow.slice(0) = arma::eye(size(R));
                              R_pow.slice(1) = R;

                              // Initialize a vector for storing the transition probabilities
                              Rcpp::NumericVector state_probs(n_states);

                              // keep calculating the conditional probability of n jumps until
                              // the threshold is exceeded. store powers of R accordingly.
                              while(c_prob < n_thresh[0]) {

                                        // increment the number of jumps
                                        n_jumps += 1;

                                        // check whether to add additional slices to the R_pow cube
                                        if(n_jumps == R_pow_size) {
                                                  R_pow.insert_slices(R_pow.n_slices, 8);
                                                  R_pow_size = R_pow.n_slices;
                                        }

                                        // Add the new power of R to the cube and calculate c_prob
                                        R_pow.slice(n_jumps) = R_pow.slice(n_jumps - 1) * R;
                                        c_prob += exp(-m*T) * pow(m*T, n_jumps) / Rcpp::internal::factorial(n_jumps) * R_pow.slice(n_jumps)(a-1, b-1) / p_ab;
                              }

                              // initialize the path matrix
                              int path_nrows = n_jumps + 2;
                              arma::mat path(path_nrows, 2);
                              path(0,0) = t0;
                              path(0,1) = a;
                              path(path_nrows - 1, 0) = t1;
                              path(path_nrows - 1, 1) = b;

                              // transition times are uniformly distributed in the
                              // interval. Sample them, sort them, and place in path.
                              arma::colvec transitions = Rcpp::runif(n_jumps, t0, t1);
                              std::sort(transitions.begin(), transitions.end());
                              path(arma::span(1,n_jumps), 0) = transitions;

                              // Sample the states at the transition times
                              for(int j = 1; j < n_jumps + 1; ++j) {
                                        state_probs = arma::trans(R(path(j-1, 1) - 1, span::all)) % R_pow.slice(n_jumps-j)(span::all, b-1) / R_pow(path(j-1, 1)-1, b-1, n_jumps-j+1);
                                        path(j, 1) = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs)[0];
                              }

                              // Determine which transitions are virtual transitions
                              arma::vec keep_inds(path_nrows, arma::fill::ones);
                              for(int j = 1; j < n_jumps + 1; ++j) {
                                        if(path(j, 1) == path(j-1, 1)) {
                                                  keep_inds[j] = 0;
                                        }
                              }

                              // create a matrix for the complete path without virtual jumps
                              arma::mat path_comp = path.rows(arma::find(keep_inds == 1));

                              return path_comp;
                    }
          }
}

//' Compute the matrix exponential.
//'
//' @param Q matrix
//'
//' @return Matrix exponential of Q
// [[Rcpp::export]]
arma::mat comp_expmat(const arma::mat& Q) {
          return arma::expmat(Q);
}
