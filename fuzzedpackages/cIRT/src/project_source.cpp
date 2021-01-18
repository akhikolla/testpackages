#include <RcppArmadillo.h>
#include "armadillo_manipulations.h"
#include "rgens.h"
#include "total_cpp.h"

//' Two Parameter Choice IRT Model MCMC
//'
//' Performs an MCMC routine for a two parameter IRT Model using Choice Data
//'
//' @param unique_subject_ids A `vector` with length \eqn{N \times 1}{N x 1}
//'                           containing unique subject IDs.
//' @param subject_ids        A `vector` with length \eqn{NK \times 1}{N*K x 1}
//'                           containing subject IDs.
//' @param choices_nk         A `vector` with length \eqn{NK \times 1}{N*K x 1}
//'                           containing subject choices.
//' @param fixed_effects      A `matrix` with dimensions
//'                           \eqn{NK \times P_1}{N*K x P_1} containing
//'                           fixed effect design matrix without theta.
//' @param B                  A \eqn{V} dimensional column `vector`
//'                           relating \eqn{\theta_i}{theta_i} and
//'                           \eqn{\zeta_i}{zeta_i}.
//' @param rv_effects_design  A `matrix` with dimensions
//'                           \eqn{NK \times V}{N*K x V} containing
//'                           random effect variables.
//' @param gamma              A `vector` with dimensions \eqn{P \times 1}{P x 1}
//'                           containing fixed parameter estimates,
//'                           where \eqn{P = P_1 + P_2}
//' @param beta               A `vector` with dimensions \eqn{P_2}
//'                           containing random parameter estimates.
//' @param zeta_rv            A `matrix` with dimensions \eqn{N \times V}{N x V}
//'                           containing random parameter estimates.
//' @param Sigma_zeta_inv     A `matrix` with dimensions
//'                           \eqn{P_2 \times P_2}{P_2 x P_2}.
//' @param Y                  A `matrix` of dimensions \eqn{N \times J}{N x J}
//'                           for Dichotomous item responses
//' @param theta0             A `vector` of length \eqn{N \times 1}{N x 1}
//'                           for latent theta.
//' @param a0                 A `vector` of length \eqn{J}
//'                           for item discriminations.
//' @param b0                 A `vector` of length \eqn{J}
//'                           for item locations.
//' @param mu_xi0             A `vector` of dimension 2 (i.e. c(0,1)) that is a
//'                           prior for item parameter means.
//' @param Sig_xi0            A `matrix` of dimension 2x2 (i.e. diag(2)) that
//'                           is a prior for item parameter vc matrix.
//'
//' @return
//' A `list` that contains:
//' \describe{
//'   \item{\code{ai1}}{A `vector` of length J}
//'   \item{\code{bi1}}{A `vector` of length J}
//'   \item{\code{theta1}}{A `vector` of length N}
//'   \item{\code{Z_c}}{A `matrix` of length NK}
//'   \item{\code{Wzeta_0}}{A `matrix` of length NK}
//' }
//'
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @seealso [cIRT()], [rmvnorm()], and [riwishart()]
//' @export
//' @examples \dontrun{
//' # Call with the following data:
//' TwoPLChoicemcmc(cogDAT, theta0, a0, b0, mu_xi0, Sig_xi0)
//' }
// [[Rcpp::export]]
Rcpp::List TwoPLChoicemcmc(/* New Parameters */
                           const arma::vec &unique_subject_ids,
                           const arma::vec &subject_ids,
                           const arma::vec &choices_nk,
                           const arma::mat &fixed_effects, const arma::vec &B,
                           const arma::mat &rv_effects_design,
                           const arma::vec &gamma, const arma::vec &beta,
                           const arma::mat &zeta_rv,
                           const arma::mat &Sigma_zeta_inv,
                           /* Old Parameters */
                           const arma::mat &Y, const arma::vec &theta0,
                           const arma::vec &a0, const arma::vec &b0,
                           const arma::vec &mu_xi0, const arma::mat &Sig_xi0)
{
    // Get dimensions of the item matrix
    // Number of participants
    unsigned int N = Y.n_rows;
    // Number of items
    unsigned int J = Y.n_cols;

    unsigned int NK = choices_nk.n_elem;

    // Number of fixed effects no theta
    //  unsigned int P = fixed_effects_n.n_cols;
    //  unsigned int P2 = fixed_effects_theta.n_cols;

    // Support for eq. 35
    // Variance of theta (1 x 1)
    double sigma2_theta = 1;
    // Mean of theta (1 x 1)
    double mu_theta = 0;

    /* design matrix */
    arma::mat X0 = theta0; // need to make copy
    X0.insert_cols(0, arma::ones<arma::vec>(N) * -1.0);

    /* matrix of item parameters */
    arma::mat xi0 = arma::join_rows(b0, a0);

    /* mean for augmented data  */
    arma::mat mZ0 = X0 * (xi0.t());

    /* sampling augmented data  */
    arma::mat u = arma::randu<arma::mat>(N, J);

    /* Handles how to convert (-1)^Y to rcpp
     * The operation is defined as:
     * If element = 1, -1 is put in place.
     * If element = 0, 1 is put in its place.
     */
    arma::mat p = Y * -1;
    p.elem(arma::find(p == 0)).fill(1);

    /* Brief history lesson before loop:
     * calling q/p/d/rnorm functions requires first parameter to be vectorized
     * calling the function with Rf_ before it should require only a double in
     * first parameter. However, Rf_qnorm seems to be broken. Hence,since 0.10
     * there are R::qnorm and R::pnorm, which still only access the c-level
     * functions
     */
    arma::mat Z0(N, J);
    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < J; j++) {
            Z0(i, j) = R::qnorm(
                Y(i, j) * u(i, j) +
                    u(i, j) * p(i, j) * R::pnorm(0.0, mZ0(i, j), 1.0, 1, 0) +
                    Y(i, j) * R::pnorm(0.0, mZ0(i, j), 1.0, 1, 0),
                mZ0(i, j), 1.0, 1, 0);
        }
    }

    /* Inject new updation code here! */
    /* Sampling thetas from truncated normal */
    arma::vec theta1(N);

    /* Z_ci = (Z_ci1, ... Z_ciK), which is K * 1 */
    /* Z_c = (Z_c1, ... Z_cN), which is N*K x 1 */

    arma::vec Z_c(NK);
    arma::vec Wzeta_0(NK);
    /* The loop is meant to identify and extract candidates */
    for (unsigned int i = 0; i < N; i++) {

        /* Subsetting for the current subject being considered */
        arma::uvec current_subject = find(subject_ids == unique_subject_ids(i));

        /* Isolate all the data we need */

        /* Fixed effects without theta for subject i (k x p) */
        //    arma::mat Xi_ntheta = fixed_effects_n.rows(current_subject);

        /* Fixed effects with theta for subject i (k x p_2) */
        //    arma::mat Xi_theta = fixed_effects_theta.rows(current_subject);

        //    double ptheta = theta0(i);
        //    Xi_theta = ptheta*Xi_theta;

        /* Fixed effects combined */
        //    arma::mat Xi = arma::join_rows(Xi_ntheta, Xi_theta);
        arma::mat Xi = fixed_effects.rows(current_subject);

        /* Random Effects for subject i (k x V) */
        arma::mat Wi = rv_effects_design.rows(current_subject);

        /* Random effect zeta_rv for subject i (1 x v) */
        arma::mat zeta_i = zeta_rv.row(i);

        /* Choices give 1 if hard or 0 if easy (k x 1). */
        arma::vec c_ik = choices_nk.elem(current_subject);

        // Cache W*Zeta_0
        Wzeta_0.elem(current_subject) = Wi * trans(zeta_rv.row(i)); // nk x 1

        /* Augmented Data for subject i (J x 1) */

        /* Obtain the mean by summing k x 1 elements of fixed
         * and random effects for the normal distribution.
         */
        arma::vec mu_ik = Xi * gamma + Wzeta_0.elem(current_subject);

        /* Handles how to convert (-1)^Y to rcpp
         * The operation is defined as:
         * If element = 1, -1 is put in place.
         * If element = 0, 1 is put in its place.
         */
        arma::vec p = c_ik * -1;
        p.elem(arma::find(p == 0)).fill(1);

        /* Creates Z_ci (aka new mean)*/
        unsigned int temp_length = Xi.n_rows;

        /* Update z_c with k x 1 */
        arma::vec Z_ci(temp_length);

        /* Brief history lesson before loop:
         * calling q/p/d/rnorm functions requires first parameter to be
         * vectorized Since 0.10 there are R::qnorm and R::pnorm which give
         * point returns, which still only access the c-level functions
         */
        for (unsigned int ii = 0; ii < temp_length; ii++) {
            double u = R::runif(0.0, 1.0);

            Z_ci(ii) = R::qnorm(
                c_ik(ii) * u + u * p(ii) * R::pnorm(0.0, mu_ik(ii), 1.0, 1, 0) +
                    c_ik(ii) * R::pnorm(0.0, mu_ik(ii), 1.0, 1, 0),
                mu_ik(ii), 1.0, 1, 0);
        }

        // Cache Z_c
        Z_c.elem(current_subject) = Z_ci;

        // Gamma no theta (P1 x 1), we have an offset since the index starts at
        // 0.
        //     arma::vec gamma_ntheta = gamma_fixed.rows(0, P1-1);

        // Gamma theta (P2 x 1), we have an offset since the index starts at 0.
        //     arma::vec gamma_theta = gamma_fixed.rows(P1, P1+P2-1);

        // End Extraction

        // This satisfies eq.32
        // Z_tilda_yi gives a (1 x J) vector
        arma::mat Z_tilda_yi = Z0.row(i) + trans(b0);
        // End eq 32

        // Start eq. 37 for full conditional for zeta_i
        // z_ci is given by k x 1, (k x 1) - (k x p1)*(p1 x 1) - (k X v)*(V x k)
        //     arma::vec Z_tilda_ci = Z_ci - Xi_ntheta*gamma_ntheta -
        //     Wi*trans(zeta_i);
        arma::vec Z_tilda_ci = Z_ci - Xi * gamma;

        // End eq. 33

        // Start eq. 35

        // Store results
        // (k x p2) * (p2 x 1) = (k x 1)
        arma::mat tB = B.t();

        //  // (1 x p2) * (p2 x k) = (1 x k)
        //   arma::mat trans_LV_design = trans(LV_design);

        // Sigma_zeta^-1
        // sampling Sigma_zeta_inv for speed
        //    arma::mat Sigma_zeta_inv = inv(Sigma_zeta);

        // Calulate sigma_itheta ^2 (1 x 1). Note inv() = ()^(-1)
        double sigma_sq_itheta = arma::conv_to<double>::from(
            1.0 /
            (trans(a0) * a0 + tB * Sigma_zeta_inv * B + 1.0 / sigma2_theta));
        double mu_itheta = arma::conv_to<double>::from(
            sigma_sq_itheta *
            (Z_tilda_yi * a0 + tB * Sigma_zeta_inv * zeta_i.t() +
             mu_theta / sigma2_theta));

        // End eq 35 part 1

        // Perform the theta draw
        // We need to take the square root of the variance since rnorm requires
        // mu and std
        theta1(i) = R::rnorm(mu_itheta, sqrt(sigma_sq_itheta));
    }

    // arma::mat Sigma_xi_inv = inv(Sig_xi0);

    /* Albert's way to obtain A and B */
    arma::vec m1 = arma::ones<arma::vec>(N) * -1;
    arma::mat X_theta = join_rows(theta1, m1);
    arma::mat XX_xi = X_theta.t() * X_theta; // X^T*X
    arma::mat XZ = X_theta.t() * Z0;         // X^T * Z_y (augmented)

    arma::mat Sig_xi_star = inv(Sig_xi0 + XX_xi);
    // arma::vec mu_xi_star(2); // this temporary

    arma::vec as(J);
    arma::vec bs(J);

    double pa, mb_a, vb_a;

    for (unsigned j = 0; j < J; j++) {

        arma::vec mu_xi_star = Sig_xi_star * (Sig_xi0 * mu_xi0 + XZ.col(j));

        double sd_a = sqrt(Sig_xi_star(0, 0));
        pa = R::pnorm(0.0, mu_xi_star(0), sd_a, 1, 0);
        as(j) = R::qnorm(pa + R::runif(0.0, 1.0) * (1.0 - pa), mu_xi_star(0),
                         sd_a, 1, 0);

        mb_a = mu_xi_star(1) +
               Sig_xi_star(0, 1) / Sig_xi_star(1, 1) * (as(j) - mu_xi_star(0));
        vb_a = Sig_xi_star(1, 1) -
               Sig_xi_star(0, 1) * Sig_xi_star(0, 1) / Sig_xi_star(0, 0);
        bs(j) = R::qnorm(R::runif(0.0, 1.0), mb_a, sqrt(vb_a), 1, 0);
    }

    return Rcpp::List::create(Rcpp::Named("ai1", as), Rcpp::Named("bi1", bs),
                              Rcpp::Named("theta1", theta1),
                              Rcpp::Named("Z_c", Z_c),
                              Rcpp::Named("Wzeta_0", Wzeta_0));
}

//' Probit Hierarchical Level Model
//'
//' Performs modeling procedure for a Probit Hierarchical Level Model.
//'
//' @param unique_subject_ids   A `vector` with length N x 1 containing
//'                             unique subject IDs.
//' @param subject_ids          A `vector` with length N*K x 1 containing
//'                             subject IDs.
//' @param choices_nk           A `vector` with length N*K x 1 containing
//'                             subject choices.
//' @param fixed_effects_design A `matrix` with dimensions N*K x P containing
//'                             fixed effect variables.
//' @param rv_effects_design    A `matrix` with dimensions N*K x V containing
//'                             random effect variables.
//' @param B_elem_plus1         A `V[[1]]` dimensional column `vector`
//'                             indicating which zeta_i relate to theta_i.
//' @param gamma                A `vector` with dimensions P_1 x 1 containing
//'                             fixed parameter estimates.
//' @param beta                 A `vector` with dimensions P_2 x 1 containing
//'                             random parameter estimates.
//' @param theta                A `vector` with dimensions N x 1 containing
//'                             subject understanding estimates.
//' @param zeta_rv              A `matrix` with dimensions N x V containing
//'                             random parameter estimates.
//' @param WtW                  A `field<matrix>` P x P x N contains the
//'                             caching for direct sum.
//' @param Z_c                  A `vector` with dimensions N*K x 1
//' @param Wzeta_0              A `vector` with dimensions N*K x 1
//' @param inv_Sigma_gamma      A `matrix` with dimensions P x P that is the
//'                             prior inverse sigma matrix for gamma.
//' @param mu_gamma             A `vector` with length P x 1 that is the prior
//'                             mean vector for gamma.
//' @param Sigma_zeta_inv       A `matrix` with dimensions V x V that is the
//'                             prior inverse sigma matrix for zeta.
//' @param S0                   A `matrix` with dimensions V x V that is the
//'                             prior sigma matrix for zeta.
//' @param mu_beta              A `vector` with dimensions P_2 x 1, that is
//'                             the mean of beta.
//' @param sigma_beta_inv       A `matrix` with dimensions P_2 x P_2, that is
//'                             the inverse sigma matrix of beta.
//' @return
//' A `list` that contains:
//' \describe{
//'   \item{\code{zeta_1}}{A `vector` of length N}
//'   \item{\code{sigma_zeta_inv_1}}{A `matrix` of dimensions V x V}
//'   \item{\code{gamma_1}}{A `vector` of length P}
//'   \item{\code{beta_1}}{A `vector` of length V}
//'   \item{\code{B}}{A `matrix` of length V}
//' }
//' 
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @details
//' The function is implemented to decrease the amount of vectorizations
//' necessary.
//'
//' @seealso
//' [rwishart()] and [TwoPLChoicemcmc()]
//'
//' @export
// [[Rcpp::export]]
Rcpp::List
probitHLM(const arma::vec &unique_subject_ids, const arma::vec &subject_ids,
          const arma::vec &choices_nk, const arma::mat &fixed_effects_design,
          const arma::mat &rv_effects_design, const arma::uvec &B_elem_plus1,
          const arma::mat &gamma, const arma::vec &beta, const arma::vec &theta,
          const arma::mat &zeta_rv,
          /* Trading Memory for Speed */
          const arma::field<arma::mat> &WtW, const arma::vec Z_c,
          const arma::vec Wzeta_0,
          /* Bunch of priors */
          const arma::mat &inv_Sigma_gamma, const arma::mat &mu_gamma,
          const arma::mat &Sigma_zeta_inv, const arma::mat &S0,
          const arma::vec &mu_beta, const arma::mat &sigma_beta_inv)
{

    /* Number of random effects */
    unsigned int V = rv_effects_design.n_cols;

    /* Total number of choices to process (N*K) */
    // unsigned int NK = choices_nk.n_elem;

    // Z_ci = (Z_ci1, ... Z_ciK), which is K * 1 */
    // Z_c = (Z_c1, ... Z_cN), which is N*K x 1 */
    // Number of Subjects
    unsigned int N = unique_subject_ids.n_elem;

    // Zeta is N x V matrix
    arma::mat zeta_1(N, V);

    // Save time
    // arma::mat inv_Sigma_gamma = inv(Sigma_gamma);

    /* START 1.3 Full Conditional for Fixed Effects */

    // \Sigma^star_gamma update (P x P)
    arma::mat Sigma_star = inv(
        trans(fixed_effects_design) * fixed_effects_design + inv_Sigma_gamma);

    // \tilda{Z}_gamma update (NK x 1)
    arma::mat tilda_Z_gamma = Z_c - Wzeta_0;

    // Mu update (P x 1)
    arma::mat mu_star =
        Sigma_star * (trans(fixed_effects_design) * tilda_Z_gamma +
                      inv_Sigma_gamma * mu_gamma);

    // Gamma update (P x 1)
    arma::mat gamma_1 = trans(rmvnorm(1, mu_star, Sigma_star));

    /* END 1.3 Full Conditional for Fixed Effects */

    /* START 1.3.1 Full Conditional for Random Effect */

    // Creating B vector
    arma::vec B = arma::zeros<arma::vec>(V);
    B(B_elem_plus1 - 1) = beta;

    // Cache Sigma_zeta_inv
    //    arma::mat Sigma_zeta_inv = inv(Sigma_zeta);//maybe sample
    //    Sigma_zeta_inv from wishart instead
    // Build with new data
    for (unsigned int i = 0; i < N; i++) {
        // Subsetting for the current subject being considered
        arma::uvec current_subject = find(subject_ids == unique_subject_ids(i));

        // extracting theta_i
        double theta_i = theta(i);

        // Sigma_zeta_star_i = inv((V x K) * (K x V) + (V x V)) = inv(V x V)
        arma::mat Sigma_zeta_star_i = inv(WtW(i) + Sigma_zeta_inv);

        // Z_tilda_zeta_i = (K x 1) - (K x P)*(P X 1) = (K x 1)
        arma::vec Z_tilda_zeta_i =
            Z_c.elem(current_subject) -
            fixed_effects_design.rows(current_subject) * gamma_1;

        // Mu_star_zeta_i = (V x V) * (V x K) * (K x 1) = (V x 1)
        arma::vec Mu_star_zeta_i =
            Sigma_zeta_star_i *
            (trans(rv_effects_design.rows(current_subject)) * Z_tilda_zeta_i +
             Sigma_zeta_inv * B * theta_i);

        // Update zeta_1 by i*V (start) to V*i+V-1 (end) by (1 x V)
        zeta_1.row(i) = rmvnorm(1, Mu_star_zeta_i, Sigma_zeta_star_i);
    }

    /* END 1.3.1 Full Conditional for Random Effect */

    /* START 1.3.2 Full Conditional for beta in random-effects model*/
    // pull elements from Sigma_zeta_inv and zeta
    arma::mat Sigma_zeta_inv_ast =
        Sigma_zeta_inv(B_elem_plus1 - 1, B_elem_plus1 - 1);
    arma::mat zeta_ast = zeta_1.cols(B_elem_plus1 - 1);

    arma::mat sigma_beta_star = inv(
        Sigma_zeta_inv_ast * arma::conv_to<double>::from(theta.t() * theta) +
        sigma_beta_inv);
    arma::mat mu_beta_star =
        sigma_beta_star *
        (Sigma_zeta_inv_ast * zeta_ast.t() * theta + sigma_beta_inv * mu_beta);

    // beta update (V x 1)
    arma::mat beta_1 = trans(rmvnorm(1, mu_beta_star, sigma_beta_star));

    /* END 1.3.2 Full Conditional for beta in random-effects model*/

    /* START 1.3.2 Full Conditional for the Variance-Covariance*/

    /* V_o = ncol(rand_effects) + 2 // we added one here.
     * N + V_o = N + V + 2
     */
    //  arma::mat sigma_zeta_1 = riwishart(N + V + 2, (zeta_1-theta*B.t()).t() *
    //  (zeta_1-theta*B.t()) + S0);
    arma::mat sigma_zeta_inv_1 = rwishart(
        N + V + 2,
        inv((zeta_1 - theta * B.t()).t() * (zeta_1 - theta * B.t()) + S0));

    /* END 1.3.2 Full Conditional for the Variance-Covariance*/

    return Rcpp::List::create(Rcpp::Named("zeta_1", zeta_1),
                              Rcpp::Named("sigma_zeta_inv_1", sigma_zeta_inv_1),
                              Rcpp::Named("gamma_1", gamma_1),
                              Rcpp::Named("beta_1", beta_1),
                              Rcpp::Named("B", B));
}

//' Generic Implementation of Choice IRT MCMC
//'
//' Builds a model using MCMC
//'
//' @param subject_ids    A `vector` that contains subject IDs for each line of
//'                       data in the choice vector (e.g. For 1 subject that
//'                       made 5 choices, we would have the number 1 appear
//'                       five times consecutively.)
//' @param fixed_effects  A `matrix` with NK x P1 dimensions that acts as the
//'                       design matrix for terms WITHOUT theta.
//' @param B_elem_plus1   A `V[[1]]` dimensional column `vector` indicating
//'                       which zeta_i relate to theta_i.
//' @param rv_effects     A `matrix` with NK x V dimensions for random effects
//'                       design matrix.
//' @param trial_matrix   A `matrix` with N x J dimensions, where J denotes the
//'                       number of items presented. The matrix MUST contain
//'                       only 1's and 0's.
//' @param choices_nk     A `vector` with NK length that contains the choice
//'                       value e.g. 0 or 1.
//' @param chain_length   An `int` that controls how many MCMC draws there are.
//'                       (> 0)
//' @param burnit         An `int` that describes how many MCMC draws should be
//'                       discarded.
//'
//' @return
//' A `list` that contains:
//'
//' \describe{
//'   \item{\code{as}}{A `matrix` of dimension chain_length x J}
//'   \item{\code{bs}}{A `matrix` of dimension chain_length x J}
//'   \item{\code{gs}}{A `matrix` of dimension chain_length x P_1}
//'   \item{\code{Sigma_zeta_inv}}{An `array` of dimension V x V x chain_length}
//'   \item{\code{betas}}{A `matrix` of dimension chain_length x P_2}
//' }
//'
//' @author
//' Steven Andrew Culpepper and James Joseph Balamuta
//'
//' @seealso
//' [TwoPLChoicemcmc()], [probitHLM()], [center_matrix()],
//' [rmvnorm()], [rwishart()], and [riwishart()]
//'
//' @export
//' @examples
//' \dontrun{
//' # Variables
//' # Y = trial matix
//' # C = KN vector of binary choices
//' # N = #of subjects
//' # J = # of items
//' # K= # of choices
//' # atrue = true item discriminations
//' # btrue = true item locations
//' # thetatrue = true thetas/latent performance
//' # gamma = fixed effects coefficients
//' # Sig = random-effects variance-covariance
//' # subid = id variable for subjects
//' 
//' # Load the Package
//' library(cIRT)
//' 
//' # Load the Data
//' data(trial_matrix)
//' data(choice_matrix)
//'
//' # Thurstone design matrices
//' all_nopractice = subset(all_data_trials, experiment_loop.thisN > -1)
//' hard_items = choice_matrix$hard_q_id
//' easy_items = choice_matrix$easy_q_id
//'
//' D_easy = model.matrix( ~ -1 + factor(easy_items))
//' D_hard = -1 * model.matrix( ~ -1 + factor(hard_items))[, -c(5, 10, 15)]
//'
//' # Defining effect-coded contrasts
//' high_contrasts = rbind(-1, diag(4))
//' rownames(high_contrasts) = 12:16
//' low_contrasts = rbind(-1, diag(2))
//' rownames(low_contrasts) = 4:6
//'
//' # Creating high & low factors
//' high = factor(choice_matrix[, 'high_value'])
//' low = factor(choice_matrix[, 'low_value'])
//' contrasts(high) = high_contrasts
//' contrasts(low) = low_contrasts
//'
//' fixed_effects = model.matrix( ~ high + low)
//' fixed_effects_base = fixed_effects[, 1]
//' fixed_effects_int = model.matrix( ~ high * low)
//'
//'
//' # Model with Thurstone D Matrix
//' system.time({
//'  out_model_thurstone = cIRT(
//'    choice_matrix[, 'subject_id'],
//'    cbind(fixed_effects[, -1], D_easy, D_hard),
//'    c(1:ncol(fixed_effects)),
//'    as.matrix(fixed_effects),
//'    as.matrix(trial_matrix),
//'    choice_matrix[, 'choose_hard_q'],
//'    20000,
//'    25000
//'  )
//' })
//'
//'
//' vlabels_thurstone = colnames(cbind(fixed_effects[, -1], D_easy, D_hard))
//' G_thurstone = t(apply(
//'  out_model_thurstone$gs0,
//'  2,
//'  FUN = quantile,
//'  probs = c(.5, .025, .975)
//' ))
//'
//' rownames(G_thurstone) = vlabels_thurstone
//' B_thurstone = t(apply(
//'  out_model_thurstone$beta,
//'  2,
//'  FUN = quantile,
//'  probs = c(.5, 0.025, .975)
//' ))
//'
//' rownames(B_thurstone) = colnames(fixed_effects)
//'
//' S_thurstone = solve(
//'   apply(out_model_thurstone$Sigma_zeta_inv, c(1, 2), FUN = mean)
//' )
//'
//' inv_sd = diag(1 / sqrt(diag(solve(
//'  apply(out_model_thurstone$Sigma_zeta_inv, c(1, 2),
//'        FUN = mean)
//' ))))
//'
//' inv_sd %*% S_thurstone %*% inv_sd
//' apply(out_model_thurstone$as, 2, FUN = mean)
//' apply(out_model_thurstone$bs, 2, FUN = mean)
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List cIRT(arma::vec subject_ids, arma::mat fixed_effects,
                arma::uvec B_elem_plus1, arma::mat rv_effects,
                arma::mat trial_matrix, arma::vec choices_nk,
                unsigned int burnit, unsigned int chain_length = 10000)
{
    // Dimensions for cognitive items
    // N is subject and J is items
    unsigned int N = trial_matrix.n_rows;
    unsigned int J = trial_matrix.n_cols;
    // Dimensions for fixed design matrix
    // NK is N subjects * K choices
    unsigned int NK = fixed_effects.n_rows;

    // Number of choices.
    unsigned int K = NK / N;

    // P1 gives the number of fixed gamma parameters
    //  unsigned int P1 = fixed_effects.n_cols;
    unsigned int P = fixed_effects.n_cols;
    // V1 gives the number of slopes relating theta to zeta
    unsigned int V1 = B_elem_plus1.n_elem;
    // P gives the total number of parameters
    //  unsigned int P = P1+P2;
    // V gives the number of variables for the random effects matrix W
    unsigned int V = rv_effects.n_cols;
    arma::vec B = arma::zeros<arma::vec>(V);
    arma::vec beta = arma::zeros<arma::vec>(V1);

    // Create an index of unique values so that we can avoid having "find"
    // difficulties.
    arma::vec unique_subject_ids = unique(subject_ids);
    // Returns (N x 1) down from (NK x 1)

    /* The next snip of code does two things:
     * 1. Center the W_i matrix
     * 2. Create the WtW field
     */

    // There are N people, so we will have N W_i matrices.
    arma::field<arma::mat> WtW(N);
    for (unsigned int i = 0; i < N; i++) {
        // Random effects for subject i (K x P)
        arma::uvec current_subject = find(subject_ids == unique_subject_ids(i));
        arma::mat Wi = rv_effects.rows(current_subject);
        // arma::mat center_matrix_Wi =
        // center_matrix(rv_effects.rows(current_subject));
        // center_matrix_Wi.col(0) =
        // arma::ones<arma::vec>(current_subject.n_elem); // hack job
        // rv_effects.rows(current_subject) = center_matrix_Wi;
        WtW(i) = trans(Wi) * Wi;
    }

    // Center the fixed effects matrix
    // fixed_effects_n = center_matrix(fixed_effects_n);
    // fixed_effects_theta = center_matrix(fixed_effects_theta);

    // Force it to be one after centering.
    // fixed_effects_theta.col(0) = arma::ones<arma::vec>(NK);

    // 2PL Ogive Model
    // Set up theta0
    arma::vec theta0 = arma::randn<arma::vec>(N);

    // Set up xi0
    arma::vec ai0 = arma::ones<arma::vec>(J);  // a is 1
    arma::vec bi0 = arma::zeros<arma::vec>(J); // b is 0
    arma::mat Sig_xi0 = arma::eye<arma::mat>(2, 2);
    arma::vec mu_xi0 = arma::vec("0 0");

    // Starting values for choice model
    // Random effect
    arma::mat zeta_rv = arma::zeros<arma::mat>(N, V);
    arma::mat Sigma_zeta_inv = arma::eye(V, V);

    // Prior for Sigb
    arma::mat S0 = arma::eye(V, V);

    // Prior for beta
    arma::mat sigma_beta_inv = arma::eye(V1, V1);
    arma::mat mu_beta = arma::zeros<arma::vec>(V1);

    // Fixed effect
    arma::vec gamma = arma::zeros<arma::vec>(P);

    // Prior for gamma
    arma::vec mu_gamma = arma::zeros<arma::vec>(P);
    arma::mat inv_Sigma_gamma = arma::eye(P, P);

    // Storage
    arma::mat all_as(chain_length, J);
    arma::mat all_bs(chain_length, J);
    arma::mat all_gamma(chain_length, P);
    arma::mat all_beta(chain_length, V1);
    arma::mat nC(K + 1, chain_length - burnit);
    arma::mat nY(J + 1, chain_length - burnit);
    arma::mat Cs(NK, chain_length - burnit);
    /* Dimension for the elements of a upper triangle matrix */
    //  arma::mat all_Sigma_zeta(chain_length,V*(V+1)/2);
    arma::cube all_Sigma_zeta_inv(V, V, chain_length);

    // fixed_design = cbind(fixed_effects_n, fixed_effects_theta)
    arma::mat fixed_design = arma::zeros(NK, P);

    for (unsigned int i = 0; i < chain_length; i++) {

        // Rcpp::Rcout << "Calling TwoPLChoicemcmc! Iter:" << i << std::endl;
        Rcpp::List outit =
            TwoPLChoicemcmc(/* Start new parameters */
                            unique_subject_ids, subject_ids, choices_nk,
                            fixed_effects, B, rv_effects, gamma, beta, zeta_rv,
                            Sigma_zeta_inv,
                            /* End new parameters */
                            trial_matrix, theta0, ai0, bi0, mu_xi0, Sig_xi0);
        // Y, theta0, a0, b0, mu_xi0, Sig_xi0

        ai0 = Rcpp::as<arma::mat>(outit[0]);    // ai1
        bi0 = Rcpp::as<arma::mat>(outit[1]);    // bi1
        theta0 = Rcpp::as<arma::vec>(outit[2]); // theta1

        //  Rcpp::Rcout << "Theta Mean:" << mean(theta0) << " - var: " <<
        //  var(theta0) << std::endl;

        //    mu_xi0	= as<arma::vec>(outit[3]); // mu_xi1
        //    Sig_xi0	= as<arma::mat>(outit[4]); // Sig_xi1
        // New Add
        arma::vec Z_c = Rcpp::as<arma::vec>(outit[3]);
        arma::vec Wzeta_0 = Rcpp::as<arma::vec>(outit[4]);

        // Rcpp::Rcout << "Building Design Matrix: " << i << std::endl;

        /*
        // The loop builds the fixed design matrix that includes thetas.
        for(unsigned int w = 0; w < N; w++){
        // Subsetting for the current subject being considered
        arma::uvec current_subject = find(subject_ids == unique_subject_ids(w));

        // Fixed effects without theta for subject i (k x p_1)
        arma::mat Xi_ntheta = fixed_effects.rows(current_subject);

        // Fixed effects with theta for subject i (k x p_2)
        double theta_i = theta0(w);
        arma::mat Xi_theta = theta_i*fixed_effects_theta.rows(current_subject);

        // Fixed effects combined
        fixed_design.rows(current_subject) = arma::join_rows(Xi_ntheta,
        Xi_theta);
        }
        */

        // Rcpp::Rcout << "Calling probitHLM! Iter:" << i << std::endl;
        // DTA, w, Bs0, gs0, Sigb0, g0_prior, Sigg0_prior, S0
        Rcpp::List outchit = probitHLM(
            unique_subject_ids, subject_ids, /* Subsetting */
            choices_nk,    /* Type of question choice (n*k x 1) */
            fixed_effects, /* Fixed Design w/ Theta appended */
            rv_effects,    /* Random Effects design matrix */
            B_elem_plus1,  /* vector indicating nonzero beta elements */
            gamma,         /* Estimated Gamma */
            beta, theta0,  /* Estimated vector of ability parameters */
            zeta_rv,       /* Estimated Random Effect */
            /* Trading Memory for Speed */
            WtW, Z_c, Wzeta_0,
            /* Priors */
            inv_Sigma_gamma, mu_gamma, Sigma_zeta_inv, S0, mu_beta,
            sigma_beta_inv);

        // Bs0
        zeta_rv = Rcpp::as<arma::mat>(outchit[0]);
        Sigma_zeta_inv = Rcpp::as<arma::mat>(outchit[1]);

        // gs0
        gamma = Rcpp::as<arma::vec>(outchit[2]);

        // beta
        beta = Rcpp::as<arma::vec>(outchit[3]);
        B = Rcpp::as<arma::vec>(outchit[4]);

        // Rcpp::Rcout << "probitHLM Success!" << std::endl;

        /*indicator of where we are:
        if( i% 10==0){   Rcpp::Rcout << "Iteration " << i << " completed" <<
        std::endl;}*/

        all_as.row(i) = (ai0).t();
        all_bs.row(i) = (bi0).t();
        all_gamma.row(i) = gamma.t();
        all_beta.row(i) = beta.t();
        // Rcpp::Rcout<< "Sigma_zeta" << Sigma_zeta<< std::endl;

        //    all_Sigma_zeta.row(i) =
        //    (Sigma_zeta.elem(find(trimatu(Sigma_zeta)))).t();
        all_Sigma_zeta_inv.slice(i) = Sigma_zeta_inv;

        // PPPs for Y and C total scores
        if (i > burnit - 1) {
            Rcpp::List pppdat = Generate_Choice(
                N, J, K, theta0, ai0, bi0, zeta_rv, gamma, fixed_effects,
                rv_effects, subject_ids, unique_subject_ids);

            nC.col(i - burnit) = Rcpp::as<arma::vec>(pppdat[0]);
            nY.col(i - burnit) = Rcpp::as<arma::vec>(pppdat[1]);
            Cs.col(i - burnit) = Rcpp::as<arma::vec>(pppdat[2]);
        }
    }

    //  all_as.save("as.csv", arma::csv_ascii);
    //  all_bs.save("bs.csv", arma::csv_ascii);
    //  all_gamma.save("gs.csv", arma::csv_ascii);
    //  all_beta.save("betas.csv", arma::csv_ascii);
    //  all_Sigma_zeta.save("Sigma_zeta.csv", arma::csv_ascii);

    return Rcpp::List::create(
        Rcpp::Named("as", all_as), Rcpp::Named("bs", all_bs),
        Rcpp::Named("gs0", all_gamma),
        Rcpp::Named("Sigma_zeta_inv", all_Sigma_zeta_inv),
        Rcpp::Named("betas", all_beta), Rcpp::Named("nC", nC),
        Rcpp::Named("nY", nY), Rcpp::Named("Cs", Cs));
}
