#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

#define EPSILON 1e-10
#define HALFMAX 25

//' @useDynLib higlasso
// [[Rcpp::export]]
arma::field <arma::mat> generate_Xi(arma::field <arma::mat> Xm)
{
    uword s = Xm.n_elem;
    uword n = Xm(0).n_rows;
    field <mat> Xi = field <mat> (s, s);
    for (uword j$ = 0; j$ < s; ++j$) {
        for (uword j = 0; j < j$; ++j) {
            mat Xmk1 = Xm(j);
            mat Xmk2 = Xm(j$);

            uword n_Xmk1_col = Xmk1.n_cols;
            uword n_Xmk2_col = Xmk2.n_cols;
            cube C = cube(n, n_Xmk2_col, n_Xmk1_col);
            for (uword k = 0; k < n_Xmk1_col; ++k)
                C.slice(k) = Xmk2.each_col() % Xmk1.col(k);

            // Collapse C into n x (pj * p$) matrix */
            C.reshape(C.n_rows, C.n_cols * C.n_slices, 1);
            // Make Xi(j,j$) have unit column variance
            Xi(j, j$) = C.slice(0);
        }
    }
    return Xi;
}

field <vec> initalize_eta(field <vec> eta_init, uword s)
{
    uword i = 0;
    field <vec> eta = field <vec> (s, s);
    for (uword j$ = 0; j$ < s; ++j$)
        for (uword j = 0; j < j$; ++j)
            eta(j, j$) = eta_init(i++);

    return eta;
}

std::function <mat (uword, uword)> initalize_Xi(field <mat> &Xi_init, uword s)
{
    std::function <mat (uword, uword)> f = [s, &Xi_init] (uword j, uword k) ->
        mat {
            return Xi_init(j + s * k);
        };
    return f;
}

mat calculate_Xt(std::function <mat (uword, uword)> Xi, field <vec> beta)
{
    uword p = 0;
    for (uword k = 0; k < beta.n_elem; ++k)
        for (uword j = 0; j < k; ++j)
            p += Xi(j, k).n_cols;

    mat Xt = mat(Xi(0,1).n_rows, p);

    p = 0;
    for (uword k = 0; k < beta.n_elem; ++k)
        for (uword j = 0; j < k; ++j) {
            uword p$ = p + Xi(j, k).n_cols;
            Xt.cols(p,  p$ - 1) = Xi(j, k) * diagmat(kron(beta(j), beta(k)));
            p = p$;
        }
    return Xt;
}

vec calculate_Yt(vec residuals, std::function <mat (uword, uword)> Xi,
                     field <vec> beta, field <vec> eta)
{
    vec Yt = residuals;
    for (uword j$ = 0; j$ < beta.n_elem; ++j$)
        for (uword j = 0; j < j$; ++j)
                Yt += Xi(j, j$) * (eta(j, j$) % kron(beta(j), beta(j$)));

    return Yt;
}

vec calculate_D(vec v, double sigma)
{
    vec D = abs(v);

    if (norm(D, "inf") < EPSILON)
        D.fill(EPSILON);

    double l2_norm = norm(D);
    // calculate D from d_k
    double inf_norm = norm(D, "inf");
    double w        = exp(-inf_norm / sigma);
    for (uword k = 0; k < D.n_elem; ++k) {
        if (inf_norm - D(k) < 1e-12)
            D(k) =  w * (1.0 / l2_norm - l2_norm / (D(k) * sigma));
        else
            D(k) = w / l2_norm;
    }

    return D;
}

void update_beta(vec &residuals, field <vec> &beta, field <vec> eta, field <mat>
                     Xm, std::function <mat (uword, uword)> Xi, double l1,
                     double sigma)
{
    //update beta
    int n = residuals.n_elem;
    for (uword j$ = 0; j$ < beta.n_elem; ++j$) {
        mat Xtj  = Xm(j$);
        vec Ytj = residuals + Xm(j$) * beta(j$);

        mat I_pj = eye <mat> (beta(j$).n_elem, beta(j$).n_elem);
        for (uword k = 0; k < j$; ++k) {
            Xtj += Xi(k, j$) * diagmat(eta(k, j$)) * kron(beta(k), I_pj);
            Ytj += Xi(k, j$) * (eta(k, j$) % kron(beta(k), beta(j$)));
        }
        for (uword l = j$ + 1; l < beta.n_elem; ++l) {
            Xtj += Xi(j$, l) * diagmat(eta(j$, l)) * kron(I_pj, beta(l));
            Ytj += Xi(j$, l) * (eta(j$, l) % kron(beta(j$), beta(l)));
        }

        vec D = calculate_D(beta(j$), sigma);
        vec C = (abs(D) - D) % beta(j$);

        vec new_beta;
        mat M = Xtj.t() * Xtj +  n * l1 * diagmat(D);
        if (!solve(new_beta, M, Xtj.t() * Ytj + l1 * C))
            Rcpp::warning("Armadillo solve() failed.\n");

        vec Ytj_orig = Ytj;
        auto ppen_lik = [&](double omega)
        {
            vec b = omega * new_beta + (1.0 - omega) * beta(j$);
            double penlik = l1 * exp(-norm(b, "inf") / sigma) * norm(b);
            Ytj = Ytj_orig - Xm(j$) * b;
            for (uword k = 0; k < j$; ++k)
                Ytj -= Xi(k, j$) * (eta(k, j$) % kron(beta(k), b));
            for (uword l = j$ + 1; l < beta.n_elem; ++l)
                Ytj -= Xi(j$, l) * (eta(j$, l) % kron(b, beta(l)));
            return penlik += 0.5 * dot(Ytj, Ytj) / n;
        };

        double pen_zero = 0.5 * dot(Ytj, Ytj) / n;
        double pen0     = ppen_lik(0.0);
        double pen1     = ppen_lik(1.0);

        int i = HALFMAX;
        double omega = 1.0;
        while (i-- && pen1 > pen0)
            pen1 = ppen_lik(omega /= 2.0);

        switch (index_min(vec({pen_zero, pen0, pen1}))) {
        case 0:
            beta(j$).fill(0.0);
            residuals = Ytj_orig;
            continue;
        case 1:
            continue;
        case 2:
            beta(j$) = omega * new_beta + (1.0 - omega) * beta(j$);
            residuals = Ytj;
        }
    }
}

field <vec> update_eta(mat Xt, vec Yt, field <vec> eta, double l2, double sigma)
{
    int n = Yt.n_elem;

    vec e = vec(Xt.n_cols);
    vec D = vec(Xt.n_cols);
    uword p = 0;
    for (uword k = 0; k < eta.n_rows; ++k) {
        for (uword j = 0; j < k; ++j) {
            uword p$ = p + eta(j, k).n_elem;
            e.subvec(p, p$ - 1) = eta(j, k);
            D.subvec(p, p$ - 1) = calculate_D(eta(j, k), sigma);
            p = p$;
        }
    }

    vec C = (abs(D) - D) % e;
    if (!solve(e, Xt.t() * Xt + n * l2 * diagmat(D), Xt.t() * Yt + l2 * C)) {
        Rcpp::warning("Armadillo solve() failed.\n");
        return eta;
    }

    p = 0;
    field <vec> new_eta = eta;
    for (uword k = 0; k < eta.n_rows; ++k) {
        for (uword j = 0; j < k; ++j) {
            uword p$ = p + eta(j, k).n_elem;
            new_eta(j, k) = e.subvec(p, p$ - 1);
            p = p$;
        }
    }
    return new_eta;
}

// performs backtracking line search for eta. Yt and eta are modified.
void bls_eta(vec &Yt, field <vec> new_eta, field <vec> &eta, field <vec> beta,
                 std::function <mat (uword, uword)> Xi, double l2, double sigma)
{
    int n = Yt.n_elem;

    vec Yt_org = Yt;
    double eta_reg = 0.0;
    auto ppen_lik = [&](double omega)
    {
        Yt = Yt_org;
        eta_reg = 0.0;
        for (uword k = 0; k < beta.n_elem; ++k) {
            for (uword j = 0; j < k; ++j) {
                vec e = omega * new_eta(j, k) + (1.0 - omega) * eta(j, k);
                Yt -= Xi(j, k) * (e % kron(beta(j), beta(k)));
                eta_reg += exp(-norm(e, "inf") / sigma) * norm(e);
            }
        }
        return 0.5 * dot(Yt, Yt) / n + l2 * eta_reg;
    };

    double pen0 = ppen_lik(0.0);
    double pen1 = ppen_lik(1.0);
    double omega = 1.0;

    int i = HALFMAX;
    while (i-- && pen1 > pen0)
        pen1 = ppen_lik(omega /= 2.0);

    if (pen1 > pen0)
        Yt = Yt_org;
    else {
        for (uword k = 0; k < beta.n_elem; ++k)
            for (uword j = 0; j < k; ++j) {
                eta(j, k) = omega * new_eta(j, k) + (1.0 - omega) * eta(j, k);
                vec Ytz = Yt + Xi(j, k) * (eta(j, k) % kron(beta(j), beta(k)));
                double pen_zero = dot(Ytz, Ytz) / (2.0 * n) + l2 * (eta_reg
                                  - exp(-norm(eta(j, k), "inf") / sigma)
                                  * norm(eta(j, k)));
                if (pen_zero < pen1)
                    eta(j, k).fill(0.0);
            }
    }
    for (uword k = 0; k < beta.n_elem; ++k)
        for (uword j = 0; j < k; ++j) {
            if (norm(beta(j)) < DBL_EPSILON || norm(beta(k)) < DBL_EPSILON
                                            || norm(eta(j, k)) < DBL_EPSILON)
                eta(j, k).fill(0.0);
        }
}

double penalized_likelihood(vec residuals, field <vec> beta, field <vec> eta,
                                double sigma, double l1, double l2)
{
    int n = residuals.n_elem;
    double beta_reg = 0.0;
    for (uword j = 0; j < beta.n_elem; ++j)
        beta_reg += exp(-norm(beta(j), "inf") / sigma) * norm(beta(j));

    double eta_reg = 0.0;
    for (uword k = 0; k < beta.n_elem; ++k)
        for (uword j = 0; j < k; ++j)
            eta_reg += exp(-norm(eta(j, k), "inf") / sigma) * norm(eta(j, k));

    return 0.5 * dot(residuals, residuals) / n + l1 * beta_reg + l2 * eta_reg;
}

//' @useDynLib higlasso
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
Rcpp::List higlasso_internal(arma::vec Y, arma::field <arma::mat> Xm,
                                 arma::field <arma::mat> Xi_init, arma::mat Z,
                                 arma::field <arma::vec> beta, arma::field
                                 <arma::vec> eta_init, double l1, double l2,
                                 double sigma, int maxit, double delta)
{
    field <vec> eta = initalize_eta(eta_init, beta.n_elem);
    auto Xi = initalize_Xi(Xi_init, beta.n_elem);

    field <vec> new_eta = eta;

    // initialize residuals
    vec residuals = Y;
    for (uword j = 0; j < Xm.n_elem; ++j)
        residuals -= Xm(j) * beta(j);

    for (uword k = 0; k < beta.n_elem; ++k)
        for (uword j = 0; j < k; ++j)
            residuals -= Xi(j, k) * (eta(j, k) % kron(beta(j), beta(k)));

    // psuedo inverse of Z
    mat Zinv  = pinv(Z);
    vec alpha = vec(Z.n_cols, fill::zeros);

    double pen_lik0 = 0.0;
    double pen_lik1 = penalized_likelihood(residuals, beta, eta, sigma, l1, l2);
    int it = 0;
    do {
        pen_lik0 = pen_lik1;

        // update non regularized coefficients
        vec new_alpha = Zinv * (residuals + Z * alpha);
        residuals += Z * (alpha - new_alpha);
        alpha = new_alpha;

        update_beta(residuals, beta, eta, Xm, Xi, l1, sigma);
        // update eta
        vec Yt = calculate_Yt(residuals, Xi, beta, eta);
        mat Xt = calculate_Xt(Xi, beta);
        new_eta = update_eta(Xt, Yt, eta, l2, sigma);

        bls_eta(Yt, new_eta, eta, beta, Xi, l2, sigma);
        residuals = Yt;

        pen_lik1 = penalized_likelihood(residuals, beta, eta, sigma, l1, l2);
        // check penalized likelihood
    } while (it++ < maxit && (pen_lik0 - pen_lik1) >= delta);

    if (it >= maxit && (pen_lik0 - pen_lik1) / pen_lik0 >= delta)
        Rcpp::warning("'maxit' reached without convergence: %f > %f\n",
        (pen_lik0 - pen_lik1) / pen_lik0, delta);

    // clean up betas
    for (uword j = 0; j < beta.n_elem; ++j)
        if (norm(beta(j), "inf") < 10 * EPSILON)
            beta(j).fill(0);

    // clean up etas and convert them to gammas
    for (uword k = 0; k < beta.n_elem; ++k) {
        for (uword j = 0; j < k; ++j) {
            if (norm(eta(j, k), "inf") < 10 * EPSILON)
                eta(j, k).fill(0);
            else
                eta(j, k) = eta(j, k) % kron(beta(j), beta(k));
        }

    }

    double mse = dot(residuals, residuals) / residuals.n_elem;
    return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
        Rcpp::Named("beta") = beta, Rcpp::Named("gamma") = eta,
        Rcpp::Named("mse") = mse);
}
