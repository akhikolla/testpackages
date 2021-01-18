#ifndef EDINA_EDINA_MEAT_H
#define EDINA_EDINA_MEAT_H

namespace edina {

inline arma::mat OddsRatio(unsigned int N, unsigned int J, const arma::mat &Yt)
{
    arma::mat M2_temp = arma::zeros<arma::mat>(J, J);
    for (unsigned int j1 = 0; j1 < J - 1; ++j1) {
        for (unsigned int j2 = j1 + 1; j2 < J; ++j2) {
            double n11 = arma::accu(Yt.col(j1) % Yt.col(j2));
            double n00 = arma::accu((1. - Yt.col(j1)) % (1. - Yt.col(j2)));
            double n10 = arma::accu(Yt.col(j1) % (1. - Yt.col(j2)));
            double n01 = N - n11 - n00 - n10;
            M2_temp(j1, j2) = (n11 * n00) / (n10 * n01);
        }
    }
    return M2_temp;
}


inline arma::vec bijectionvector(unsigned int K)
{
    arma::vec vv(K);

    for (unsigned int k = 0; k < K; ++k) {
        vv(k) = std::pow(2.0, static_cast<double>(K - k) - 1.0);
    }

    return vv;
}

inline arma::vec inv_bijectionvector(unsigned int K, double CL)
{
    arma::vec alpha(K);

    for (unsigned int k = 0; k < K; ++k) {
        double twopow = std::pow(2.0, static_cast<double>(K - k) - 1.0);
        alpha(k) = (twopow <= CL);
        CL = CL - twopow * alpha(k);
    }

    return alpha;
}

inline arma::mat random_Q(unsigned int J, unsigned int K)
{

    // Generate two identity matrices
    arma::vec one_K = arma::ones<arma::vec>(K);
    arma::mat I_K = arma::diagmat(one_K);
    arma::mat Two_I_K = arma::join_cols(I_K, I_K);

    // Generate Q1
    unsigned int Jm2K = J - 2 * K;
    unsigned int J1max = K;

    if (K > Jm2K) {
        J1max = Jm2K;
    }

    unsigned int J1 =
        arma::as_scalar(arma::randi<arma::vec>(1, arma::distr_param(1, J1max)));

    arma::mat U1 = arma::randu<arma::mat>(J1, K);
    arma::mat Q1 = arma::zeros<arma::mat>(J1, K);

    // Fix elements so rows are nonzero
    arma::vec col_ks = arma::randi<arma::vec>(J1, arma::distr_param(0, K - 1));

    for (unsigned int j = 0; j < J1; ++j) {
        Q1(j, col_ks(j)) = 1;
    }

    // Fix elements so columns are nonzero
    arma::vec row_ks = arma::randi<arma::vec>(K, arma::distr_param(0, J1 - 1));

    for (unsigned int k = 0; k < K; ++k) {
        Q1(row_ks(k), k) = 1;
    }

    Q1.elem(arma::find(Q1 > .5)).fill(1.0);

    arma::mat Q = arma::join_cols(Two_I_K, Q1);

    // Generating the remaining elements of Q in Q2
    unsigned int Jm2KmJ1 = Jm2K - J1;
    arma::mat Q2 = arma::zeros<arma::mat>(Jm2KmJ1, K);
    if (Jm2KmJ1 > 0) {
        arma::mat U2 = arma::randu<arma::mat>(Jm2KmJ1, K);

        // Use sample()
        arma::vec jks_w_1s =
            arma::randi<arma::vec>(Jm2KmJ1, arma::distr_param(0, K - 1));
        for (unsigned int j = 0; j < Jm2KmJ1; ++j) {
            U2(j, jks_w_1s(j)) = 1;
        }

        Q2.elem(arma::find(U2 > .5)).fill(1.0);
        Q = arma::join_cols(Q, Q2);
    }

    // Q
    arma::uvec P = arma::uvec(J);
    for (unsigned int j = 0; j < J; ++j) {
        P(j) = j;
    }

    P = arma::shuffle(P);
    return Q.rows(P);
}


inline bool check_identifiability(const arma::mat& Q)
{
    unsigned int K = Q.n_cols;
    unsigned int J = Q.n_rows;

    arma::mat ones_zero_on_diag = -1 * arma::ones<arma::mat>(K, K);
    arma::vec zeros_K = arma::zeros<arma::vec>(K);
    ones_zero_on_diag.diag() = zeros_K;

    arma::vec c_sum = (arma::sum(Q, 0)).t();
    arma::vec r_sum = arma::sum(Q, 1);
    arma::mat I_check = Q * ones_zero_on_diag;
    arma::mat I_count = arma::zeros<arma::mat>(J, K);
    I_count.elem(arma::find(I_check > -1)).fill(1.0);
    arma::vec n_ek = (arma::sum(I_count, 0)).t();

    double min_c = (arma::min(c_sum) > 2);
    double min_r = (arma::min(r_sum) > 0);
    double min_ek = (arma::min(n_ek) > 1);

    return (min_c + min_r + min_ek > 2);
}

inline arma::mat alpha_matrix(const arma::mat& Q) {

    int K = Q.n_cols;

    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    arma::mat alpha_mat(nClass, K);

    for(unsigned int cc = 0; cc < nClass; cc++){
        alpha_mat.row(cc) = inv_bijectionvector(K, cc).t();
    }

    return alpha_mat;
}


inline arma::mat ETAmat(unsigned int K, unsigned int J, const arma::mat &Q)
{
    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    arma::mat ETA(J, nClass);

    for (unsigned int cc = 0; cc < nClass; ++cc) {
        arma::vec alpha_c = inv_bijectionvector(K, cc);

        for (unsigned int j = 0; j < J; ++j) {
            arma::rowvec qj = Q.row(j);
            // Switch to as_scalar
            double compare = arma::as_scalar(qj * alpha_c - qj * qj.t());
            ETA(j, cc) = (compare >= 0);
        }
    }

    return ETA;
}

inline arma::cube ETAmat_nok(unsigned int K)
{
    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    unsigned int nqjs =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K) - 1.0));

    arma::vec zero_to_Km1 = arma::linspace(0, K - 1, K);
    arma::cube ETA_nok(K, nqjs, nClass);
    for (unsigned int k = 0; k < K; ++k) {
        arma::uvec ks = find(zero_to_Km1 != k);
        for (unsigned int cc = 0; cc < nClass; ++cc) {
            arma::vec alpha_c = inv_bijectionvector(K, cc);
            arma::vec alpha_c_nok = alpha_c(ks);
            for (unsigned int cq = 0; cq < nqjs; ++cq) {
                arma::vec qj = inv_bijectionvector(K - 1, cq);
                // switched to as_scalar
                double compare =
                    arma::as_scalar(qj.t() * alpha_c_nok - qj.t() * qj);
                ETA_nok(k, cq, cc) = (compare >= 0);
            }
        }
    }
    return ETA_nok;
}

inline arma::cube ETAmat_nok_one_m_ac(unsigned int K)
{
    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    unsigned int nqjs =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K) - 1.0));

    arma::vec zero_to_Km1 = arma::linspace(0, K - 1, K);
    arma::cube ETA_nok(K, nqjs, nClass);
    for (unsigned int k = 0; k < K; ++k) {
        arma::uvec ks = find(zero_to_Km1 != k);
        for (unsigned int cc = 0; cc < nClass; ++cc) {
            arma::vec alpha_c = inv_bijectionvector(K, cc);
            double alpha_c_k = alpha_c(k);
            arma::vec alpha_c_nok = alpha_c(ks);
            for (unsigned int cq = 0; cq < nqjs; ++cq) {
                arma::vec qj = inv_bijectionvector(K - 1, cq);
                // as_scalar
                double compare =
                    arma::as_scalar(qj.t() * alpha_c_nok - qj.t() * qj);
                ETA_nok(k, cq, cc) = (1. - alpha_c_k) * (compare >= 0);
            }
        }
    }
    return ETA_nok;
}

inline arma::vec abcount_old(unsigned int K, unsigned int k, const arma::vec &qj,
                             const arma::vec &Yj, const arma::mat &alpha)
{

    arma::vec abcounts(2);
    arma::vec zero_to_Km1 = arma::linspace(0, K - 1, K);

    arma::uvec ks = find(zero_to_Km1 != k);
    arma::mat alpha_no_k = alpha.cols(ks);
    arma::vec qj_no_k = qj(ks);
    arma::vec alpha_qj_no_k = alpha_no_k * qj_no_k;
    arma::vec eta_no_k = arma::zeros<arma::vec>(alpha_qj_no_k.n_elem);

    // Converted to as_scalar
    double cut = arma::as_scalar((qj_no_k.t() * qj_no_k));
    eta_no_k.elem(arma::find(alpha_qj_no_k >= cut)).fill(1.0);
    arma::vec one_m_alpha_k = 1 - alpha.col(k);
    arma::vec one_m_alpha_k_times_eta_no_k = eta_no_k % one_m_alpha_k;

    arma::uvec ys_eq0 = find(Yj == 0);
    abcounts(0) = arma::accu(one_m_alpha_k_times_eta_no_k(ys_eq0));
    arma::uvec ys_eq1 = find(Yj == 1);
    abcounts(1) = arma::accu(one_m_alpha_k_times_eta_no_k(ys_eq1));
    return abcounts;
}

inline arma::vec abcounts(unsigned int N, const arma::vec &Yj, const arma::vec &CLASS,
                          const arma::vec &ETAtnokimes1ma)
{

    arma::vec ab = arma::zeros<arma::vec>(2);
    for (unsigned int i = 0; i < N; ++i) {
        ab(Yj(i)) += ETAtnokimes1ma(CLASS(i));
    }
    return ab;
}

inline arma::mat ClassbyQmat(unsigned int K)
{
    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    arma::mat ETAbyQ(nClass, nClass - 1);
    for (unsigned int cc = 0; cc < nClass; ++cc) {
        arma::vec alpha_c = inv_bijectionvector(K, cc);
        for (unsigned int r = 0; r < nClass - 1; ++r) {
            arma::vec qj = inv_bijectionvector(K, r + 1);
            // converted to as_scalar
            double compare = arma::as_scalar(qj.t() * alpha_c - qj.t() * qj);
            ETAbyQ(cc, r) = (compare >= 0);
        }
    }
    return ETAbyQ;
}

inline double llj(unsigned int N, const arma::vec &Yj, const arma::vec &ETAj,
                  const arma::vec &CLASS, double gj, double sj)
{
    arma::mat abmat = arma::zeros<arma::mat>(2, 2);
    for (unsigned int i = 0; i < N; ++i) {
        abmat(Yj(i), ETAj(CLASS(i))) += 1.;
    }
    double log_lik = abmat(0, 1) * log(sj) + abmat(1, 1) * log(1. - sj) +
        abmat(1, 0) * log(gj) + abmat(0, 0) * log(1. - gj);
    return log_lik;
}

inline double lnlik_dina_condclass(unsigned int N, unsigned int J, const arma::mat &Y,
                                   const arma::mat &ETA, const arma::vec &CLASS,
                                   const arma::vec &pis, const arma::vec &gs,
                                   const arma::vec &ss)
{
    double log_lik = 0.;
    for (unsigned int i = 0; i < N; ++i) {
        arma::rowvec Yi = Y.row(i);
        arma::vec eta_i = ETA.col(CLASS(i));
        for (unsigned int j = 0; j < J; ++j) {
            double sj = ss(j);
            double gj = gs(j);
            double Yij = Y(j);
            double eta_ij = eta_i(j);
            double ometaij = 1. - eta_ij;
            log_lik += Yij * (eta_ij * log(1 - sj) + ometaij * log(gj)) +
                (1. - Yij) * (eta_ij * log(sj) + ometaij * log(1. - gj));
        }
    }
    return log_lik;
}

inline double pYjeq1(const arma::vec &ETAbyQ, const arma::vec &pis, double nClass,
                     double sj, double gj)
{
    double one_m_s = 1. - sj;
    double pj1 = 0;

    for (unsigned int cc = 0; cc < nClass; ++cc) {
        double eta_c = ETAbyQ(cc);
        pj1 += pis(cc) * (eta_c * one_m_s + (1. - eta_c) * gj);
    }

    return pj1;
}

inline double pYit(const arma::vec &ETA_it, const arma::vec &Y_it, const arma::vec &ss,
                   const arma::vec &gs)
{

    arma::vec one_m_ss = 1. - ss;
    arma::vec one_m_gs = 1. - gs;
    arma::vec one_m_ETA_it = 1. - ETA_it;
    arma::vec one_m_Y_it = 1. - Y_it;

    arma::vec ps = Y_it % (one_m_ss % ETA_it + gs % one_m_ETA_it) +
        one_m_Y_it % (ss % ETA_it + one_m_gs % one_m_ETA_it);

    return arma::prod(ps);
}

inline double lnlik_dina(unsigned int N, unsigned int J, unsigned int nClass,
                         const arma::mat &Y, const arma::mat &ETA,
                         const arma::vec &pis, const arma::vec &gs,
                         const arma::vec &ss)
{

    double log_lik = 0.;

    for (unsigned int i = 0; i < N; ++i) {
        arma::rowvec Yi = Y.row(i);
        arma::rowvec pyi_c(nClass);

        for (unsigned int cc = 0; cc < nClass; ++cc) {
            pyi_c(cc) = pYit(ETA.col(cc), Yi.t(), ss, gs);
        }

        // as_scalar
        double lik_i = arma::as_scalar(pyi_c * pis);
        log_lik += log(lik_i);
    }

    return log_lik;
}

inline void updateQ_DINA(arma::mat &Q, const arma::mat &Y, const arma::mat &alpha,
                         const arma::vec &ss, const arma::vec &gs)
{

    unsigned int K = Q.n_cols;
    unsigned int J = Q.n_rows;
    double qjk, u, flag1, flag0;

    arma::vec qj(K);
    arma::vec flag(K);
    arma::vec zero_to_Km1 = arma::linspace(0, K - 1, K);

    for (unsigned int j = 0; j < J; ++j) {

        for (unsigned int k = 0; k < K; ++k) {
            // qjk = Q(j,k);
            // checking whether 1 is possible
            arma::mat Q1 = Q;
            Q1(j, k) = 1;
            flag1 = check_identifiability(Q1);

            if (flag1 == true) {
                // checking whether 0 is possible
                arma::mat Q0 = Q;
                Q0(j, k) = 0;
                flag0 = check_identifiability(Q0);

                // update based upon posterior for qjk
                if (flag0 == true) {
                    double s_d_1mg = ss(j) / (1.0 - gs(j));
                    double Onems_d_g = (1.0 - ss(j)) / gs(j);
                    arma::vec qj = (Q.row(j)).t();
                    arma::uvec ks = find(zero_to_Km1 != k);
                    arma::mat alpha_no_k = alpha.cols(ks);
                    arma::vec qj_no_k = qj(ks);
                    arma::vec alpha_qj_no_k = alpha_no_k * qj_no_k;
                    arma::vec eta_no_k =
                        arma::zeros<arma::vec>(alpha_qj_no_k.n_elem);
                    // switched to as_scalar
                    double cut = arma::as_scalar((qj_no_k.t() * qj_no_k));
                    //*arma::ones<arma::vec>(alpha_qj_no_k.n_elem);
                    eta_no_k.elem(arma::find(alpha_qj_no_k >= cut)).fill(1.0);
                    arma::vec one_m_alpha_k = 1 - alpha.col(k);
                    arma::vec one_m_alpha_k_times_eta_no_k =
                        eta_no_k % one_m_alpha_k;

                    arma::uvec ys_eq0 = find(Y.col(j) == 0);
                    double as0_m_as1 =
                        arma::accu(one_m_alpha_k_times_eta_no_k(ys_eq0));
                    arma::uvec ys_eq1 = find(Y.col(j) == 1);
                    double bs0_m_bs1 =
                        arma::accu(one_m_alpha_k_times_eta_no_k(ys_eq1));
                    //       Rcpp::Rcout << as0_m_as1<<bs0_m_bs1 << std::endl;
                    u = R::runif(0.0, 1.0);
                    qjk = 1.0 * (log(1. - u) - log(u) >
                                     as0_m_as1 * log(s_d_1mg) +
                                     bs0_m_bs1 * log(Onems_d_g));
                    // update Q
                    Q(j, k) = qjk;
                }
            }
        } // k loop
    }     // j loop
}

inline void updateQ_DINA_new(unsigned int N, unsigned int K, unsigned int J,
                             arma::mat &Q, const arma::mat &Y, const arma::vec &CLASS,
                             const arma::vec &ss, const arma::vec &gs,
                             const arma::vec &vj, const arma::cube &ETAmatnok,
                             const arma::mat &a_by_q, const arma::vec &vv)
{
    double qjk, flag1;
    arma::vec zero_to_Km1 = arma::linspace(0, K - 1, K);
    arma::vec abn;

    for (unsigned int j = 0; j < J; ++j) {

        arma::vec Yj = Y.col(j);
        double s_d_1mg = ss(j) / (1.0 - gs(j));
        double Onems_d_g = (1.0 - ss(j)) / gs(j);

        for (unsigned int k = 0; k < K; ++k) {
            qjk = Q(j, k);
            // checking whether 1 is possible
            arma::mat Q0 = Q;
            Q0(j, k) = 1. - qjk;
            flag1 = check_identifiability(Q0);

            if (flag1 == true) {

                /*
                 // Old
                 double sj = ss(j);
                 double gj = gs(j);
                 arma::rowvec qj1 = Q.row(j);
                 qj1(k) = 1.;
                 unsigned int qj1_biject = arma::conv_to< unsigned int >::from(
                 qj1*vv); arma::rowvec qj0 = qj1; qj0(k) = 0.; unsigned int
                 qj0_biject = arma::conv_to< unsigned int >::from( qj0*vv);
                 double llj1=llj(N,Yj,a_by_q.col(qj1_biject-1),CLASS,gj,sj);
                 double llj0=llj(N,Yj,a_by_q.col(qj0_biject-1),CLASS,gj,sj);
                 double u = R::runif(0.0,1.0);
                 qjk = 1.0*(log(1.-u) -log(u) > llj0 - llj1);
                 // end Old
                 */

                arma::uvec ks = find(zero_to_Km1 != k);
                arma::rowvec qj = Q.row(j);
                arma::vec qjnok = qj(ks);
                // switched to as_scalar
                double qjv = arma::as_scalar(qjnok.t() * vj);
                arma::vec ETAtnokimes1ma = ETAmatnok.tube(k, qjv);
                abn = abcounts(N, Yj, CLASS, ETAtnokimes1ma);
                double u = R::runif(0.0, 1.0);
                qjk = 1.0 * (log(1. - u) - log(u) >
                                 abn(0) * log(s_d_1mg) + abn(1) * log(Onems_d_g));
                // update Q
                Q(j, k) = qjk;
            }
        } // k loop
    }     // j loop
}

inline double cond_threshold(unsigned int k, unsigned int j, unsigned int n_noks,
                             unsigned int N, unsigned int K, const arma::vec &Yj,
                             const arma::vec &CLASS, const arma::mat &Q, double gj,
                             double sj, const arma::cube &ETAmatnok)
{
    double s_d_1mg = sj / (1.0 - gj);
    double Onems_d_g = (1.0 - sj) / gj;
    double mean_p1 = 0.;
    arma::vec I0pI1 = arma::zeros<arma::vec>(n_noks);
    arma::vec I0(n_noks);
    arma::vec I1(n_noks);
    arma::vec qj(K);
    arma::vec zero_to_Km1 = arma::linspace(0, K - 1, K);
    arma::uvec ks = find(zero_to_Km1 != k);
    arma::mat Q0 = Q;

    arma::vec abn(2);

    for (unsigned int ck = 0; ck < n_noks; ++ck) {
        qj(ks) = inv_bijectionvector(K - 1, ck);
        qj(k) = 0;
        Q0.row(j) = qj.t();
        I0(ck) = check_identifiability(Q0);
        Q0(j, k) = 1.;
        I1(ck) = check_identifiability(Q0);
        I0pI1(ck) = I0(ck) + I1(ck);
    }

    if (arma::accu(I0pI1) > 0) {
        arma::uvec valid_indices = find(I0pI1 > 0);
        unsigned int n_valid = valid_indices.n_elem;

        arma::vec p1(n_valid);
        for (unsigned int h = 0; h < n_valid; ++h) {
            unsigned int ck = valid_indices(h);
            arma::vec ETAtnokimes1ma = ETAmatnok.tube(k, ck);
            abn = abcounts(N, Yj, CLASS, ETAtnokimes1ma);

            double thres = abn(0) * log(s_d_1mg) + abn(1) * log(Onems_d_g);

            p1(h) = I0(ck) * I1(ck) * (1. / (1. + exp(thres))) +
                (1 - I0(ck)) * I1(ck);
        }

        mean_p1 = arma::mean(p1);
    }

    return mean_p1;
}

inline arma::mat sim_Y_dina(unsigned int N, unsigned int J, const arma::vec &CLASS,
                            const arma::mat &ETA, const arma::vec &gs,
                            const arma::vec &ss)
{
    arma::mat Y(N, J);
    for (unsigned int i = 0; i < N; ++i) {
        double class_i = CLASS(i);
        arma::vec ETA_i = ETA.col(class_i);
        for (unsigned int j = 0; j < J; ++j) {
            double u = R::runif(0, 1);
            Y(i, j) =
                1. * (gs(j) * (1. - ETA_i(j)) + (1. - ss(j)) * ETA_i(j) > u);
        }
    }
    return Y;
}

inline void parm_update_nomiss(unsigned int N, unsigned int J, unsigned int K,
                               unsigned int nClass, const arma::mat &Y,
                               const arma::mat &ETA, arma::vec &gs, arma::vec &ss,
                               arma::vec &CLASS, arma::vec &pis)
{
    arma::vec pY(nClass);
    arma::cube ab_tilde = arma::zeros<arma::cube>(J, 2, 2);
    // update alpha
    for (unsigned int i = 0; i < N; ++i) {
        arma::vec Yi = (Y.row(i)).t();

        for (unsigned int cc = 0; cc < nClass; ++cc) {
            arma::vec ETA_it = ETA.col(cc);
            pY(cc) = pYit(ETA_it, Yi, ss, gs);
        }

        arma::vec numerator = pY % pis;
        arma::vec PS = numerator / arma::sum(numerator);
        double class_i = rgen::rmultinomial(PS);
        CLASS(i) = class_i;

        // update guess and slip full conditional beta parms
        arma::vec ETA_c = ETA.col(class_i);
        for (unsigned int j = 0; j < J; ++j) {
            ab_tilde(j, ETA_c(j), Yi(j)) += 1.;
        }
    }
    // update pis
    arma::uvec class_sum =
        arma::hist(CLASS, arma::linspace<arma::vec>(0, nClass - 1, nClass));
    arma::vec deltatilde = arma::conv_to<arma::vec>::from(class_sum) + 1.;
    pis = rgen::rdirichlet(deltatilde);

    // update guess and slip probabilities
    for (unsigned int j = 0; j < J; ++j) {
        double us = R::runif(0, 1);
        double ug = R::runif(0, 1);
        double sold = ss(j);
        // draw g conditoned upon s_t-1
        double ab_g1 = ab_tilde(j, 0, 1);
        double ab_g0 = ab_tilde(j, 0, 0);
        double pg = R::pbeta(1.0 - sold, ab_g1 + 1., ab_g0 + 1., 1, 0);
        double gnew = R::qbeta(ug * pg, ab_g1 + 1., ab_g0 + 1., 1, 0);
        // draw s conditoned upon g
        double ab_s1 = ab_tilde(j, 1, 1);
        double ab_s0 = ab_tilde(j, 1, 0);
        double ps = R::pbeta(1.0 - gnew, ab_s0 + 1., ab_s1 + 1., 1, 0);
        double snew = R::qbeta(us * ps, ab_s0 + 1., ab_s1 + 1., 1, 0);
        gs(j) = gnew;
        ss(j) = snew;
    }
    /*return gs;*/
}

inline Rcpp::List edina_Gibbs_Q(const arma::mat &Y, unsigned int K,
                                unsigned int burnin = 1000,
                                unsigned int chain_length = 10000)
{

    // --- Initialize configuration

    // Number of Observations
    unsigned int N = Y.n_rows;

    // Number of Assessment items
    unsigned int J = Y.n_cols;

    // Number of Attributes
    // unsigned int K = Q.n_cols;

    // Number of Classes
    unsigned int nClass =
        static_cast<unsigned int>(std::pow(2.0, static_cast<double>(K)));

    // Total Chain length with burn
    // Prevents overflow
    unsigned int iter_total = chain_length + burnin;

    // Log-likelihood sum over iteration
    double loglike_summed = 0;

    // --- Saving output
    // Q Matrices

    arma::mat Q_summed = arma::zeros<arma::mat>(J, K);

    // arma::cube QS(J, K, chain_m_burn);

    // Latent probabilities
    arma::mat PIs(nClass, chain_length);

    // Slipping
    arma::mat SS(J, chain_length);

    // Guessing
    arma::mat GS(J, chain_length);

    // Compute the sample Odds Ratio
    arma::mat Sample_OR = OddsRatio(N, J, Y);

    arma::mat OR_tested_summed = arma::zeros<arma::mat>(J, J);

    // --- Setup

    // alphas, theta, pis
    arma::vec CLASS =
        arma::randi<arma::vec>(N, arma::distr_param(0, nClass - 1));

    // Slipping
    arma::vec ss = arma::randu<arma::vec>(J);

    // Guessing
    arma::vec gs = (arma::ones<arma::vec>(J) - ss) % arma::randu<arma::vec>(J);

    // 2^k ewk!
    // Conjugate prior for latent probabilities
    arma::vec delta0 = arma::ones<arma::vec>(nClass);

    // Latent Probabilities
    arma::vec pis = rgen::rdirichlet(delta0);

    // Q matrix
    arma::mat Q = random_Q(J, K);

    // ETA matrix
    arma::mat ETA = ETAmat(K, J, Q);

    // ETA Cube
    arma::cube ETAmatnokonemac = ETAmat_nok_one_m_ac(K);

    // Bijection vectors
    arma::vec vj = bijectionvector(K - 1);
    arma::vec vv = bijectionvector(K);

    // arma::mat alpha(N,K);

    // Classification by Q Matrix
    arma::mat a_by_q = ClassbyQmat(K);

    // --- Start Markov chain

    for (unsigned int t = 0; t < iter_total; ++t) {

        parm_update_nomiss(N, J, K, nClass, Y, ETA, gs, ss, CLASS, pis);

        /*
         // Old code
         // update alpha matrix w inv_bijection formula
         for(unsigned int i = 0 ;i < N; ++i) {
         arma::vec alpha_i = inv_bijectionvector(K,CLASS(i));
         alpha.row(i) = alpha_i.t();
         }
         updateQ_DINA(Q,Y,alpha,ss,gs);
         */

        updateQ_DINA_new(N, K, J, Q, Y, CLASS, ss, gs, vj, ETAmatnokonemac,
                         a_by_q, vv);

        ETA = ETAmat(K, J, Q);

        if (t > burnin - 1) {
            int tmburn = t - burnin;
            // update parameter value via pointer. save classes and PIs
            SS.col(tmburn) = ss;
            GS.col(tmburn) = gs;
            PIs.col(tmburn) = pis;
            Q_summed += Q / chain_length;

            // Simulate new data
            arma::mat Yt = sim_Y_dina(N, J, CLASS, ETA, gs, ss);

            // Compute the loglikelihood of estimated iteration and add it to
            // total
            loglike_summed += lnlik_dina(N, J, nClass, Y, ETA, pis, gs, ss);

            // Sum up the positive OR tests (upper diagonal)
            OR_tested_summed += arma::conv_to<arma::mat>::from(
                OddsRatio(N, J, Yt) > Sample_OR) /
                    chain_length;
        }
    }

    // Take means by row (e.g. 1)
    arma::mat coefs(J, 4);

    // Guessing Parameter Estimates
    coefs.col(0) = mean(GS, 1);
    // norm_type = 0: Using N-1 in STD denominator
    coefs.col(1) = stddev(GS, 0, 1);

    // Slipping parameter estimates
    coefs.col(2) = mean(SS, 1);
    // norm_type = 0: Using N-1 in STD denominator
    coefs.col(3) = stddev(SS, 0, 1);

    arma::vec PI_summed = mean(PIs, 1);

    // Estimated Q value
    arma::mat Qest = arma::conv_to<arma::mat>::from(Q_summed > .5);

    // Error here?
    double loglike_pmean = lnlik_dina(N, J, nClass, Y, ETAmat(K, J, Qest),
                                      PI_summed, coefs.col(0), coefs.col(2));

    // Release
    return Rcpp::List::create(
        Rcpp::Named("coefficients", coefs),
        Rcpp::Named("loglike_summed", loglike_summed),
        Rcpp::Named("loglike_pmean", loglike_pmean),
        Rcpp::Named("pis", PI_summed), Rcpp::Named("avg_q", Q_summed),
        Rcpp::Named("est_q", Qest), Rcpp::Named("or_tested", OR_tested_summed),
        Rcpp::Named("sample_or", Sample_OR));
}

}

#endif
