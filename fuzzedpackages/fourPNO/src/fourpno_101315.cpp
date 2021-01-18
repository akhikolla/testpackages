#include <RcppArmadillo.h>

//' Generate Random Multivariate Normal Distribution
//'
//' Creates a random Multivariate Normal when given number of
//' obs, mean, and sigma.
//'
//' @param n     An `int`, which gives the number of observations.  (> 0)
//' @param mu    A `vector` length m that represents the means of
//'              the normals.
//' @param sigma A `matrix` with dimensions m x m that provides the
//'              covariance matrix.
//'
//' @return
//' A `matrix` that is a Multivariate Normal distribution
//'
//' @author
//' James J Balamuta
//'
//' @export
//' @examples
//' # Call with the following data:
//' rmvnorm(2, c(0,0), diag(2))
// [[Rcpp::export]]
arma::mat rmvnorm(unsigned int n, const arma::vec &mu, const arma::mat &sigma)
{
    unsigned int ncols = sigma.n_cols;
    Rcpp::RNGScope scope;
    arma::mat Y(n, ncols);
    Y.imbue(norm_rand);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

//' Initialize Thresholds
//'
//' Internal function for initializing item thresholds.
//'
//' @param Ms A `vector` with the number of scale values.
//'
//' @return
//' A `matrix` that is a Multivariate Normal distribution
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @seealso
//' [Gibbs_4PNO()]
//'
//' @noRd
arma::mat kappa_initialize(const arma::vec &Ms)
{
    unsigned int J = Ms.n_elem;
    unsigned int M = max(Ms);
    arma::mat KAP0(M + 1, J);

    (KAP0.row(0)).fill(-arma::datum::inf);
    KAP0.row(1).fill(0.0);

    for (unsigned int j = 0; j < J; j++) {
        for (unsigned int m = 2; m < M + 1; m++) {

            if (m < Ms(j)) {
                KAP0(m, j) = (m - 1) * 2;
            }

            if (m == Ms(j)) {
                KAP0(m, j) = arma::datum::inf;
            }

            if (m > Ms(j)) {
                KAP0(m, j) = arma::datum::nan;
            }
        }
    }

    return KAP0;
}

//' Internal Function for Updating Theta in Gibbs Sampler
//'
//' Update theta in Gibbs sampler
//'
//' @param N               An `int`, which gives the number of observations.
//'                        (> 0)
//' @param Z               A `matrix` N by J of continuous augmented data.
//' @param as              A `vector` of item discrimination parameters.
//' @param bs              A `vector` of item threshold parameters.
//' @param theta           A `vector` of prior thetas.
//' @param mu_theta        The prior mean for theta.
//' @param Sigma_theta_inv The prior inverse variance for theta.
//'
//' @return
//' A `vector` of thetas.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @seealso
//' [Gibbs_4PNO()]
//'
//' @noRd
void update_theta(unsigned int N, const arma::mat &Z, const arma::vec &as,
                        const arma::vec &bs, arma::vec &theta,
                        const double &mu_theta, const double &Sigma_theta_inv)
{

    double apb = dot(as, bs);
    arma::vec oneN = arma::ones<arma::vec>(N);
    double vartheta = 1.0 / (dot(as, as) + Sigma_theta_inv);
    //  apb = dot(as,bs);
    arma::vec zetas(N, arma::fill::randn);
    //  zetas.imbue( norm_rand );

    theta = zetas * sqrt(vartheta) +
            vartheta * (Z * as + oneN * (apb + mu_theta * Sigma_theta_inv));

}

//' Update a and b Parameters of 2PNO, 3PNO, 4PNO
//'
//' Update item slope and threshold
//'
//' @param N            An `int`, which gives the number of observations. (> 0)
//' @param J            An `int`, which gives the number of items. (> 0)
//' @param Z            A `matrix` N by J of continuous augmented data.
//' @param as           A `vector` of item discrimination parameters.
//' @param bs           A `vector` of item threshold parameters.
//' @param theta        A `vector` of prior thetas.
//' @param mu_xi        A two dimensional `vector` of prior item parameter
//'                     means.
//' @param Sigma_xi_inv A two dimensional identity `matrix` of prior item
//'                     parameter VC matrix.
//'
//' @return
//' A `list` of item parameters.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @seealso
//' [Gibbs_4PNO()]
//'
//' @noRd
void update_ab_NA(unsigned int N, unsigned int J, const arma::mat &Z,
                        arma::vec &as, arma::vec &bs, const arma::vec &theta,
                        const arma::vec &mu_xi, const arma::mat &Sigma_xi_inv)
{

    arma::vec m1(N);
    m1.fill(-1.0);
    arma::mat X_theta = join_rows(theta, m1);
    arma::mat XX_xi = X_theta.t() * X_theta;
    arma::mat XZ = X_theta.t() * Z;

    arma::mat Sig_xi_star = inv(Sigma_xi_inv + XX_xi);
    arma::vec mu_xi_star(2);
    double pa, mb_a, vb_a, ua, ub;

    for (unsigned j = 0; j < J; j++) {

        mu_xi_star = Sig_xi_star * (Sigma_xi_inv * mu_xi + XZ.col(j));

        ua = R::runif(0.0, 1.0);
        pa = R::pnorm(0.0, mu_xi_star(0), sqrt(Sig_xi_star(0, 0)), 1, 0);
        as(j) = R::qnorm(pa + ua * (1.0 - pa), mu_xi_star(0),
                         sqrt(Sig_xi_star(0, 0)), 1, 0);

        ub = R::runif(0.0, 1.0);
        mb_a = mu_xi_star(1) +
               Sig_xi_star(0, 1) / Sig_xi_star(1, 1) * (as(j) - mu_xi_star(0));
        vb_a = Sig_xi_star(1, 1) -
               Sig_xi_star(0, 1) * Sig_xi_star(0, 1) / Sig_xi_star(0, 0);
        bs(j) = R::qnorm(ub, mb_a, sqrt(vb_a), 1, 0);
    }

}

//' Update a and b Parameters of 4pno without alpha > 0 Restriction
//'
//' Update item slope and threshold
//'
//' @param N            An `int`, which gives the number of observations. (> 0)
//' @param J            An `int`, which gives the number of items. (> 0)
//' @param Z            A `matrix` N by J of continuous augmented data.
//' @param as           A `vector` of item discrimination parameters.
//' @param bs           A `vector` of item threshold parameters.
//' @param theta        A `vector` of prior thetas.
//' @param mu_xi        A two dimensional `vector` of prior item parameter
//'                     means.
//' @param Sigma_xi_inv A two dimensional identity `matrix` of prior item
//'                     parameter VC matrix.
//'
//' @return
//' A `list` of item parameters.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @seealso
//' [Gibbs_4PNO()]
//'
//' @noRd
void update_ab_norestriction(unsigned int N, unsigned int J,
                                   const arma::mat &Z, arma::vec &as,
                                   arma::vec &bs, const arma::vec &theta,
                                   const arma::vec &mu_xi,
                                   const arma::mat &Sigma_xi_inv)
{

    arma::vec m1(N);
    m1.fill(-1.0);
    arma::mat X_theta = join_rows(theta, m1);
    arma::mat XX_xi = X_theta.t() * X_theta;
    arma::mat XZ = X_theta.t() * Z;

    arma::mat Sig_xi_star = inv(Sigma_xi_inv + XX_xi);
    arma::vec mu_xi_star(2);
    double pa, mb_a, vb_a, ua, ub;

    for (unsigned j = 0; j < J; j++) {
        mu_xi_star = Sig_xi_star * (Sigma_xi_inv * mu_xi + XZ.col(j));

        ua = R::runif(0.0, 1.0);
        pa = 0;
        as(j) = R::qnorm(pa + ua * (1.0 - pa), mu_xi_star(0),
                         sqrt(Sig_xi_star(0, 0)), 1, 0);

        ub = R::runif(0.0, 1.0);
        mb_a = mu_xi_star(1) +
               Sig_xi_star(0, 1) / Sig_xi_star(1, 1) * (as(j) - mu_xi_star(0));
        vb_a = Sig_xi_star(1, 1) -
               Sig_xi_star(0, 1) * Sig_xi_star(0, 1) / Sig_xi_star(0, 0);
        bs(j) = R::qnorm(ub, mb_a, sqrt(vb_a), 1, 0);
    }

}

//' Update Lower and Upper Asymptote Parameters of 4PNO
//'
//' Internal function to update item lower and upper asymptote
//'
//' @param Y        A N by J `matrix` of item responses.
//' @param Ysum     A `vector` of item total scores.
//' @param Z        A `matrix` N by J of continuous augmented data.
//' @param as       A `vector` of item discrimination parameters.
//' @param bs       A `vector` of item threshold parameters.
//' @param gs       A `vector` of item lower asymptote parameters.
//' @param ss       A `vector` of item upper asymptote parameters.
//' @param theta    A `vector` of prior thetas.
//' @param Kaps     A `matrix` for item thresholds
//'                 (used for internal computations).
//' @param alpha_c  The lower asymptote prior 'a' parameter.
//' @param beta_c   The lower asymptote prior 'b' parameter.
//' @param alpha_s  The upper asymptote prior 'a' parameter.
//' @param beta_s   The upper asymptote prior 'b' parameter.
//' @param gwg_reps The number of Gibbs within Gibbs MCMC samples for
//'                 marginal distribution of gamma.
//'
//' @return
//' A `list` of item threshold parameters.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @seealso
//' [Gibbs_4PNO()]
//'
//' @noRd
Rcpp::List update_WKappaZ_NA(const arma::mat &Y, const arma::vec &Ysum,
                             arma::mat &Z, const arma::vec &as,
                             const arma::vec &bs, const arma::vec &gs,
                             const arma::vec &ss, const arma::vec &theta,
                             const arma::mat &Kaps, double alpha_c,
                             double beta_c, double alpha_s, double beta_s,
                             unsigned int gwg_reps)
{

    unsigned int N = Y.n_rows;
    unsigned int J = Y.n_cols;

    arma::vec eta(N);

    arma::vec p4pno(2);
    arma::mat W = arma::zeros<arma::mat>(N, J);
    arma::vec oneN = arma::ones<arma::vec>(N);
    double us, ug, uZ, u4pno;
    double Phi_eta;
    double pg_temp, ps_temp, s_temp, g_temp;
    arma::vec gs_new = arma::zeros<arma::vec>(J);
    arma::vec ss_new = arma::zeros<arma::vec>(J);
    double n1dot, n0dot, n01, n10;
    double ps, pg;

    for (unsigned int j = 0; j < J; j++) {

        // Calculate eta_ij for all of the models and individuals with data
        eta = as(j) * theta - oneN * bs(j);

        ////////////////////////////////////////////////
        // 4PNO model
        ////////////////////////////////////////////////
        //    if(model(j)=="4PNO"){

        for (unsigned int i = 0; i < N; i++) {

            //        if(MISS(i,j)==1.0){
            uZ = R::runif(0.0, 1.0);
            u4pno = R::runif(0.0, 1.0);
            Phi_eta = R::pnorm(eta(i), 0.0, 1.0, 1, 0);

            if ((1.0 - ss(j)) * Phi_eta /
                        (gs(j) + (1.0 - ss(j) - gs(j)) * Phi_eta) >
                    u4pno and
                Y(i, j) == 1.0) {
                W(i, j) = 1.0;
            }

            if (ss(j) * Phi_eta /
                        (1.0 - gs(j) - (1.0 - ss(j) - gs(j)) * Phi_eta) >
                    u4pno and
                Y(i, j) == 0.0) {
                W(i, j) = 1.0;
            }
            //          sj = (1-W(i))*Y(i,j) + sj;

            p4pno(0) = R::pnorm(Kaps(W(i, j), j), eta(i), 1.0, 1, 0);
            p4pno(1) = R::pnorm(Kaps(W(i, j) + 1.0, j), eta(i), 1.0, 1, 0);
            Z(i, j) = R::qnorm(p4pno(0) + uZ * (p4pno(1) - p4pno(0)), eta(i),
                               1.0, 1, 0);
            //        }
        }

        // update gs and ss
        us = R::runif(0, 1);
        ug = R::runif(0, 1);

        n1dot = sum(W.col(j));
        n0dot = N - n1dot;
        n01 = Ysum(j) - dot(W.col(j), Y.col(j));
        n10 = n1dot - dot(W.col(j), Y.col(j));

        // draw g from marginal distribution w Gibbs within Gibbs
        s_temp = R::rbeta(n10 + alpha_s, n1dot - n10 + beta_s);
        for (unsigned int r = 0; r < gwg_reps; r++) {
            pg_temp = R::pbeta(1.0 - s_temp, n01 + alpha_c,
                               n0dot - n01 + beta_c, 1, 0);
            g_temp = R::qbeta(R::runif(0, 1) * pg_temp, n01 + alpha_c,
                              n0dot - n01 + beta_c, 1, 0);
            ps_temp = R::pbeta(1.0 - g_temp, n10 + alpha_s,
                               n1dot - n10 + beta_s, 1, 0);
            s_temp = R::qbeta(R::runif(0, 1) * ps_temp, n10 + alpha_s,
                              n1dot - n10 + beta_s, 1, 0);
        }
        pg = R::pbeta(1.0 - s_temp, n01 + alpha_c, n0dot - n01 + beta_c, 1, 0);
        gs_new(j) =
            R::qbeta(ug * pg, n01 + alpha_c, n0dot - n01 + beta_c, 1, 0);

        // draw s conditoned upon g
        ps = R::pbeta(1.0 - gs_new(j), n10 + alpha_s, n1dot - n10 + beta_s, 1,
                      0);
        ss_new(j) =
            R::qbeta(us * ps, n10 + alpha_s, n1dot - n10 + beta_s, 1, 0);

        //      }
    }
    return Rcpp::List::create(Rcpp::Named("ss_new", ss_new),
                              Rcpp::Named("gs_new", gs_new));
}

//' Compute 4PNO Deviance
//'
//' Internal function to -2LL
//'
//' @param N     An `int`, which gives the number of observations.  (> 0)
//' @param J     An `int`, which gives the number of items.  (> 0)
//' @param Y     A N by J `matrix` of item responses.
//' @param as    A `vector` of item discrimination parameters.
//' @param bs    A `vector` of item threshold parameters.
//' @param gs    A `vector` of item lower asymptote parameters.
//' @param ss    A `vector` of item upper asymptote parameters.
//' @param theta A `vector` of prior thetas.
//'
//' @return
//' -2LL.
//' @author
//' Steven Andrew Culpepper
//'
//' @seealso
//' [Gibbs_4PNO()]
//'
//' @export
// [[Rcpp::export]]
double min2LL_4pno(unsigned int N, unsigned int J, const arma::mat &Y,
                   const arma::vec &as, const arma::vec &bs,
                   const arma::vec &gs, const arma::vec &ss,
                   const arma::vec &theta)
{

    double m2ll = 0.0;
    double eta, p_eta;

    for (unsigned int j = 0; j < J; j++) {
        for (unsigned int i = 0; i < N; i++) {

            eta = as(j) * theta(i) - bs(j);
            p_eta = R::pnorm(eta, 0.0, 1.0, 1, 0);

            if (Y(i, j) == 0.0) {
                m2ll = -2.0 * log(1.0 - gs(j) - (1.0 - ss(j) - gs(j)) * p_eta) +
                       m2ll;
            }

            if (Y(i, j) == 1.0) {
                m2ll = -2.0 * log(gs(j) + (1.0 - ss(j) - gs(j)) * p_eta) + m2ll;
            }
        }
    }
    return m2ll;
}

//' Simulate from 4PNO Model
//'
//' Generate item responses under the 4PNO
//'
//' @param N     An `int`, which gives the number of observations. (> 0)
//' @param J     An `int`, which gives the number of items. (> 0)
//' @param as    A `vector` of item discrimination parameters.
//' @param bs    A `vector` of item threshold parameters.
//' @param gs    A `vector` of item lower asymptote parameters.
//' @param ss    A `vector` of item upper asymptote parameters.
//' @param theta A `vector` of prior thetas.
//'
//' @return
//' A N by J `matrix` of dichotomous item responses.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @seealso
//' [Gibbs_4PNO()]
//'
//' @export
// [[Rcpp::export]]
arma::mat Y_4pno_simulate(unsigned int N, unsigned int J, const arma::vec &as,
                          const arma::vec &bs, const arma::vec &gs,
                          const arma::vec &ss, const arma::vec &theta)
{

    double u, eta, p_eta;
    arma::mat Ysim(N, J);

    for (unsigned int j = 0; j < J; j++) {
        for (unsigned int i = 0; i < N; i++) {

            u = R::runif(0.0, 1.0);
            eta = as(j) * theta(i) - bs(j);
            p_eta = R::pnorm(eta, 0.0, 1.0, 1, 0);

            if (gs(j) + (1.0 - ss(j) - gs(j)) * p_eta > u) {
                Ysim(i, j) = 1.0;
            }

            else {
                Ysim(i, j) = .0;
            }
        }
    }
    return Ysim;
}

//' Calculate Tabulated Total Scores
//'
//' Internal function to -2LL
//'
//' @param N  An `int`, which gives the number of observations. (> 0)
//' @param J  An `int`, which gives the number of items. (> 0)
//' @param Y  A N by J `matrix` of item responses.
//'
//' @return
//' A vector of tabulated total scores.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @seealso
//' [Gibbs_4PNO()]
//'
//' @export
// [[Rcpp::export]]
arma::uvec Total_Tabulate(unsigned int N, unsigned int J, const arma::mat Y)
{

    arma::vec T(N);
    //  arma::vec H(J);
    T = sum(Y, 1);
    arma::uvec H = arma::hist(T, arma::linspace<arma::vec>(0, J, J + 1));

    return H;
}

//' Gibbs Implementation of 4PNO
//'
//' Internal function to -2LL
//'
//' @param Y               A N by J `matrix` of item responses.
//' @param mu_xi           A two dimensional `vector` of prior item parameter
//'                        means.
//' @param Sigma_xi_inv    A two dimensional identity `matrix` of prior item
//'                        parameter VC matrix.
//' @param mu_theta        The prior mean for theta.
//' @param Sigma_theta_inv The prior inverse variance for theta.
//' @param alpha_c         The lower asymptote prior 'a' parameter.
//' @param beta_c          The lower asymptote prior 'b' parameter.
//' @param alpha_s         The upper asymptote prior 'a' parameter.
//' @param beta_s          The upper asymptote prior 'b' parameter.
//' @param burnin          The number of MCMC samples to discard.
//' @param cTF             A J dimensional `vector` indicating which
//'                        lower asymptotes to estimate.
//'                        0 = exclude lower asymptote and
//'                        1 = include lower asymptote.
//' @param sTF             A J dimensional `vector` indicating which
//'                        upper asymptotes to estimate.
//'                        0 = exclude upper asymptote and
//'                        1 = include upper asymptote.
//' @param gwg_reps        The number of Gibbs within Gibbs MCMC samples for
//'                        marginal distribution of gamma. Values between
//'                        5 to 10 are adequate.
//' @param chain_length    The number of MCMC samples.
//'
//' @return
//' Samples from posterior.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @export
//' @examples
//' # Simulate small 4PNO dataset to demonstrate function
//' J = 5
//' N = 100
//'
//' # Population item parameters
//' as_t = rnorm(J,mean=2,sd=.5)
//' bs_t = rnorm(J,mean=0,sd=.5)
//'
//' # Sampling gs and ss with truncation
//' gs_t = rbeta(J,1,8)
//' ps_g = pbeta(1-gs_t,1,8)
//' ss_t = qbeta(runif(J)*ps_g,1,8)
//' theta_t <- rnorm(N)
//' Y_t = Y_4pno_simulate(N,J,as=as_t,bs=bs_t,gs=gs_t,ss=ss_t,theta=theta_t)
//'
//' # Setting prior parameters
//' mu_theta=0
//' Sigma_theta_inv=1
//' mu_xi = c(0,0)
//' alpha_c=alpha_s=beta_c=beta_s=1
//' Sigma_xi_inv = solve(2*matrix(c(1,0,0,1),2,2))
//' burnin = 1000
//'
//' # Execute Gibbs sampler
//' out_t = Gibbs_4PNO(Y_t,mu_xi,Sigma_xi_inv,mu_theta,
//'                    Sigma_theta_inv,alpha_c,beta_c,alpha_s,
//'                    beta_s,burnin,rep(1,J),rep(1,J),
//'                    gwg_reps=5,chain_length=burnin*2)
//'
//' # Summarizing posterior distribution
//' OUT = cbind(apply(out_t$AS[,-c(1:burnin)],1,mean),
//'             apply(out_t$BS[,-c(1:burnin)],1,mean),
//'             apply(out_t$GS[,-c(1:burnin)],1,mean),
//'             apply(out_t$SS[,-c(1:burnin)],1,mean),
//'             apply(out_t$AS[,-c(1:burnin)],1,sd),
//'             apply(out_t$BS[,-c(1:burnin)],1,sd),
//'             apply(out_t$GS[,-c(1:burnin)],1,sd),
//'             apply(out_t$SS[,-c(1:burnin)],1,sd) )
//'
//' OUT = cbind(1:J,OUT)
//' colnames(OUT) = c('Item', 'as', 'bs', 'gs', 'ss', 'as_sd', 'bs_sd',
//'                   'gs_sd', 'ss_sd')
//' print(OUT, digits = 3)
// [[Rcpp::export]]
Rcpp::List Gibbs_4PNO(const arma::mat &Y, const arma::vec &mu_xi,
                      const arma::mat &Sigma_xi_inv, const double &mu_theta,
                      const double &Sigma_theta_inv, double alpha_c,
                      double beta_c, double alpha_s, double beta_s,
                      unsigned int burnin, const arma::vec &cTF,
                      const arma::vec &sTF, unsigned int gwg_reps,
                      unsigned int chain_length = 10000)
{

    // arma::vec theta,arma::vec as,arma::vec bs,
    unsigned int N = Y.n_rows;
    unsigned int J = Y.n_cols;
    //  unsigned int K = max(area);//assuming area is coded from 1 to # of areas

    // save model deviances over iteractions
    arma::vec deviance(chain_length - burnin);

    // sum of Y's columns & vc matrix
    arma::vec oneN = arma::ones<arma::vec>(N);
    arma::vec Ysum = Y.t() * oneN;

    double tmburn;
    arma::mat Ysimt(N, J);
    arma::umat HistS(J + 1, chain_length - burnin);

    // compute average kappas, thetas, as, and bs
    arma::vec EAPtheta = arma::zeros<arma::vec>(N);

    // Setting 2 categories per item
    arma::vec Ms(J);
    Ms.fill(2.0);

    // Savinging as,bs,gs,kappas, means and vcs of theta
    arma::mat AS(J, chain_length);
    arma::mat BS(J, chain_length);
    arma::mat GSS(J, chain_length);
    arma::mat SSS(J, chain_length);
    arma::vec ms_thetas(chain_length);
    arma::vec SD_thetas(chain_length);

    // need to initialize, theta, as, bs, kappas, Z; eventually gs
    arma::vec theta = arma::randn<arma::vec>(N);
    arma::vec oneJ = arma::ones<arma::vec>(J);
    arma::vec as = arma::ones<arma::vec>(J) + 0.5 * arma::randu<arma::vec>(J);
    arma::vec bs = arma::randn<arma::vec>(J);
    arma::vec uc = arma::randu<arma::vec>(J);
    arma::vec gs = (oneJ - arma::sqrt(uc)) % cTF;
    arma::vec ss = (arma::randu<arma::vec>(J) % (oneJ - gs)) % sTF;
    arma::mat KAPS = kappa_initialize(Ms);
    arma::mat Z = arma::zeros<arma::mat>(N, J);

    // Start chain
    for (unsigned int t = 0; t < chain_length; t++) {
        // Generate Z matrix. Note that Z will be overwritten to disk here.
        Rcpp::List step1Z =
            update_WKappaZ_NA(Y, Ysum, Z, as, bs, gs, ss, theta, KAPS, alpha_c,
                              beta_c, alpha_s, beta_s, gwg_reps);

        // update value for ss and gs
        ss = Rcpp::as<arma::vec>(step1Z[0]) % sTF;
        gs = Rcpp::as<arma::vec>(step1Z[1]) % cTF;

        // Update a, b; as and bs should be automatically updated here by
        // writing to disk
        update_ab_NA(N, J, Z, as, bs, theta, mu_xi, Sigma_xi_inv);
        
        // Update theta. Theta is stored in memory and function should update
        // each iteration too.
        update_theta(N, Z, as, bs, theta, mu_theta, Sigma_theta_inv);
        
        // saving means and vcs of thetas
        SD_thetas(t) = stddev(theta);
        ms_thetas(t) = mean(theta);

        // Storing output for as, bs, and kappas
        AS.col(t) = as;
        BS.col(t) = bs;
        GSS.col(t) = gs;
        SSS.col(t) = ss;

        // /*
        if (t > burnin - 1) {
            tmburn = t - burnin;
            // Compute EAPs here
            EAPtheta = (tmburn * EAPtheta + theta) / (tmburn + 1.0);

            // Compute Simulated VC matrix
            Ysimt = Y_4pno_simulate(N, J, as, bs, gs, ss, theta);
            HistS.col(tmburn) = Total_Tabulate(N, J, Ysimt);
            // Compute -2LL for DIC or Bayes Factors
            deviance(tmburn) = min2LL_4pno(N, J, Y, as, bs, gs, ss, theta);
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("EAPtheta", EAPtheta), Rcpp::Named("HistS", HistS),
        Rcpp::Named("AS", AS), Rcpp::Named("BS", BS), Rcpp::Named("GS", GSS),
        Rcpp::Named("SS", SSS), Rcpp::Named("ms_thetas", ms_thetas),
        Rcpp::Named("SD_thetas", SD_thetas), Rcpp::Named("Ds", deviance));
}

//' Update 2PNO Model Parameters
//'
//' Internal function to update 2PNO parameters
//'
//' @param N               The number of observations.
//' @param J               The number of items.
//' @param Y               A N by J `matrix` of item responses.
//' @param Z               A `matrix` N by J of continuous augmented data.
//' @param as              A `vector` of item discrimination parameters.
//' @param bs              A `vector` of item threshold parameters.
//' @param theta           A `vector` of prior thetas.
//' @param Kaps            A `matrix` for item thresholds
//'                        (used for internal computations).
//' @param mu_xi           Prior mean for item parameters.
//' @param Sigma_xi_inv    Prior item parameter inverse variance-covariance
//'                        matrix.
//' @param mu_theta        Prior mean for theta.
//' @param Sigma_theta_inv Prior inverse variance for theta.
//'
//' @return
//' A `list` of item parameters.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @seealso
//' [Gibbs_2PNO()]
//'
//' @noRd
Rcpp::List update_2pno(unsigned int N, unsigned int J, const arma::mat &Y,
                       arma::mat &Z, const arma::vec &as, const arma::vec &bs,
                       const arma::vec &theta, const arma::mat &Kaps,
                       const arma::vec &mu_xi, const arma::mat &Sigma_xi_inv,
                       const double &mu_theta, const double &Sigma_theta_inv)
{

    arma::vec eta(N);
    arma::vec p4pno(2);
    arma::vec oneN = arma::ones<arma::vec>(N);
    double uZ;

    // parms for alpha and beta
    arma::vec m1(N);
    m1.fill(-1.0);
    arma::mat X_theta = join_rows(theta, m1);
    arma::mat XX_xi = X_theta.t() * X_theta;
    arma::mat Sig_xi_star = inv(Sigma_xi_inv + XX_xi);
    arma::vec mu_xi_star(2);
    double pa, mb_a, vb_a, ua, ub;
    arma::vec as1(J), bs1(J);

    for (unsigned int j = 0; j < J; j++) {
        // Calculate eta_ij for all of the models and individuals with data
        eta = as(j) * theta - oneN * bs(j);
        // Update Z, W
        for (unsigned int i = 0; i < N; i++) {
            uZ = R::runif(0.0, 1.0);
            p4pno(0) = R::pnorm(Kaps(Y(i, j), j), eta(i), 1.0, 1, 0);
            p4pno(1) = R::pnorm(Kaps(Y(i, j) + 1.0, j), eta(i), 1.0, 1, 0);
            Z(i, j) = R::qnorm(p4pno(0) + uZ * (p4pno(1) - p4pno(0)), eta(i),
                               1.0, 1, 0);
        }

        // update alpha,beta
        arma::mat XZ_col = X_theta.t() * Z.col(j);
        mu_xi_star = Sig_xi_star * (Sigma_xi_inv * mu_xi + XZ_col);
        ua = R::runif(0.0, 1.0);
        pa = R::pnorm(0.0, mu_xi_star(0), sqrt(Sig_xi_star(0, 0)), 1, 0);
        as1(j) = R::qnorm(pa + ua * (1.0 - pa), mu_xi_star(0),
                          sqrt(Sig_xi_star(0, 0)), 1, 0);
        ub = R::runif(0.0, 1.0);
        mb_a = mu_xi_star(1) +
               Sig_xi_star(0, 1) / Sig_xi_star(1, 1) * (as(j) - mu_xi_star(0));
        vb_a = Sig_xi_star(1, 1) -
               Sig_xi_star(0, 1) * Sig_xi_star(0, 1) / Sig_xi_star(0, 0);
        bs1(j) = R::qnorm(ub, mb_a, sqrt(vb_a), 1, 0);
    }

    // update theta
    double apb = dot(as1, bs1);
    double vartheta = 1.0 / (dot(as1, as1) + Sigma_theta_inv);
    arma::vec zetas(N, arma::fill::randn);
    arma::vec theta1 =
        zetas * sqrt(vartheta) +
        vartheta * (Z * as1 + oneN * (apb + mu_theta * Sigma_theta_inv));

    return Rcpp::List::create(Rcpp::Named("as1", as1), Rcpp::Named("bs1", bs1),
                              Rcpp::Named("theta1", theta1));
}

//' Gibbs Implementation of 2PNO
//'
//' Implement Gibbs 2PNO Sampler
//'
//' @param Y                A N by J `matrix` of item responses.
//' @param mu_xi            A two dimensional `vector` of prior item parameter
//'                         means.
//' @param Sigma_xi_inv     A two dimensional identity `matrix` of prior item
//'                         parameter VC matrix.
//' @param mu_theta         The prior mean for theta.
//' @param Sigma_theta_inv  The prior inverse variance for theta.
//' @param burnin           The number of MCMC samples to discard.
//' @param chain_length     The number of MCMC samples.
//'
//' @return
//' Samples from posterior.
//'
//' @author
//' Steven Andrew Culpepper
//'
//' @export
//' @examples
//' # simulate small 2PNO dataset to demonstrate function
//' J = 5
//' N = 100
//'
//' # Population item parameters
//' as_t = rnorm(J,mean=2,sd=.5)
//' bs_t = rnorm(J,mean=0,sd=.5)
//'
//' # Sampling gs and ss with truncation
//' gs_t = rbeta(J,1,8)
//' ps_g = pbeta(1-gs_t,1,8)
//' ss_t = qbeta(runif(J)*ps_g,1,8)
//' theta_t = rnorm(N)
//' Y_t = Y_4pno_simulate(N,J,as=as_t,bs=bs_t,gs=gs_t,ss=ss_t,theta=theta_t)
//'
//' # Setting prior parameters
//' mu_theta = 0
//' Sigma_theta_inv = 1
//' mu_xi = c(0,0)
//' alpha_c = alpha_s = beta_c = beta_s = 1
//' Sigma_xi_inv = solve(2*matrix(c(1,0,0,1), 2, 2))
//' burnin = 1000
//'
//' # Execute Gibbs sampler. This should take about 15.5 minutes
//' out_t = Gibbs_4PNO(Y_t,mu_xi,Sigma_xi_inv,mu_theta,Sigma_theta_inv,
//'                     alpha_c,beta_c,alpha_s, beta_s,burnin,
//'                     rep(1,J),rep(1,J),gwg_reps=5,chain_length=burnin*2)
//'
//' # Summarizing posterior distribution
//' OUT = cbind(
//'     apply(out_t$AS[, -c(1:burnin)], 1, mean),
//'     apply(out_t$BS[, -c(1:burnin)], 1, mean),
//'     apply(out_t$GS[, -c(1:burnin)], 1, mean),
//'     apply(out_t$SS[, -c(1:burnin)], 1, mean),
//'     apply(out_t$AS[, -c(1:burnin)], 1, sd),
//'     apply(out_t$BS[, -c(1:burnin)], 1, sd),
//'     apply(out_t$GS[, -c(1:burnin)], 1, sd),
//'     apply(out_t$SS[, -c(1:burnin)], 1, sd)
//' )
//' OUT = cbind(1:J, OUT)
//' colnames(OUT) = c('Item','as','bs','gs','ss','as_sd','bs_sd',
//'                   'gs_sd','ss_sd')
//' print(OUT, digits = 3)
// [[Rcpp::export]]
Rcpp::List Gibbs_2PNO(const arma::mat &Y, const arma::vec &mu_xi,
                      const arma::mat &Sigma_xi_inv, const double &mu_theta,
                      const double &Sigma_theta_inv, unsigned int burnin,
                      unsigned int chain_length = 10000)
{

    // arma::vec theta,arma::vec as,arma::vec bs,
    unsigned int N = Y.n_rows;
    unsigned int J = Y.n_cols;
    //  unsigned int K = max(area);//assuming area is coded from 1 to # of areas

    // save model deviances over iteractions
    arma::vec deviance(chain_length - burnin);

    // sum of Y's columns & vc matrix
    arma::vec oneN = arma::ones<arma::vec>(N);
    arma::vec Ysum = Y.t() * oneN;

    double tmburn;
    arma::mat Ysimt(N, J);
    arma::umat HistS(J + 1, chain_length - burnin);

    // compute average kappas, thetas, as, and bs
    arma::vec EAPtheta = arma::zeros<arma::vec>(N);

    // Setting 2 categories per item
    arma::vec Ms(J);
    Ms.fill(2.0);

    // Savinging as,bs,gs,kappas, means and vcs of theta
    arma::mat AS(J, chain_length);
    arma::mat BS(J, chain_length);
    arma::vec ms_thetas(chain_length);
    arma::vec SD_thetas(chain_length);

    // need to initialize, theta, as, bs, kappas, Z; eventually gs
    arma::vec theta = arma::randn<arma::vec>(N); // arma::zeros<arma::vec>(N);
    arma::vec oneJ = arma::ones<arma::vec>(J);
    arma::vec as = 0.5 * arma::randu<arma::vec>(J) + arma::ones<arma::vec>(J);
    arma::vec bs = arma::randn<arma::vec>(J);
    arma::vec gs = arma::zeros<arma::vec>(J);
    arma::vec ss = arma::zeros<arma::vec>(J);
    arma::mat KAPS = kappa_initialize(Ms);
    arma::mat Z = arma::zeros<arma::mat>(N, J);

    // Start chain
    for (unsigned int t = 0; t < chain_length; t++) {
        // Generate Z matrix. Note that Z will be overwritten to disk here.
        Rcpp::List step1Z =
            update_2pno(N, J, Y, Z, as, bs, theta, KAPS, mu_xi, Sigma_xi_inv,
                        mu_theta, Sigma_theta_inv);
        // update value for as, bs, theta
        as = Rcpp::as<arma::vec>(step1Z[0]);
        bs = Rcpp::as<arma::vec>(step1Z[1]);
        theta = Rcpp::as<arma::vec>(step1Z[2]);

        SD_thetas(t) = stddev(theta);
        ms_thetas(t) = mean(theta);

        // Storing output for as, bs, and kappas
        AS.col(t) = as;
        BS.col(t) = bs;

        if (t > burnin - 1) {
            tmburn = t - burnin;
            // Compute EAPs here
            EAPtheta = (tmburn * EAPtheta + theta) / (tmburn + 1.0);

            // Compute Simulated VC matrix
            Ysimt = Y_4pno_simulate(N, J, as, bs, gs, ss, theta);
            HistS.col(tmburn) = Total_Tabulate(N, J, Ysimt);
            // Compute -2LL for DIC or Bayes Factors
            deviance(tmburn) = min2LL_4pno(N, J, Y, as, bs, gs, ss, theta);
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("EAPtheta", EAPtheta), Rcpp::Named("HistS", HistS),
        Rcpp::Named("AS", AS), Rcpp::Named("BS", BS),
        Rcpp::Named("ms_thetas", ms_thetas),
        Rcpp::Named("SD_thetas", SD_thetas), Rcpp::Named("Ds", deviance));
}
