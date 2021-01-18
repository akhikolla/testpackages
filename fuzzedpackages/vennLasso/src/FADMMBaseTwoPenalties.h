#ifndef FADMMBASETWOPENALTIES_H
#define FADMMBASETWOPENALTIES_H

#include <RcppEigen.h>
#include "Linalg/BlasWrapper.h"


// General problem setting
//   minimize f(x) + g(z)
//   s.t. Ax + Bz = c
//
// x(n, 1), z(m, 1), A(p, n), B(p, m), c(p, 1)
//
template<typename VecTypeBeta, typename VecTypeGamma, typename VecTypeNu>
class FADMMBaseTwo
{
protected:
    typedef typename VecTypeNu::RealScalar Yscalar;

    const int dim_main;         // dimension of x
    const int dim_aux_gamma;    // dimension of z
    const int dim_aux_eta;
    const int dim_dual;         // dimension of Ax + Bz - c

    VecTypeBeta main_beta;      // parameters to be optimized
    VecTypeGamma aux_gamma;     // auxiliary parameters
    VecTypeGamma aux_eta;       // auxiliary parameters
    VecTypeNu dual_nu;          // Lagrangian multiplier

    VecTypeGamma adj_gamma;     // adjusted z vector, used for acceleration
    VecTypeGamma adj_eta;       // adjusted z vector, used for acceleration
    VecTypeNu adj_nu;           // adjusted y vector, used for acceleration
    VecTypeGamma old_gamma;     // z vector in the previous iteration, used for acceleration
    VecTypeGamma old_eta;       // z vector in the previous iteration, used for acceleration
    VecTypeNu old_nu;           // y vector in the previous iteration, used for acceleration
    double adj_a;               // coefficient used for acceleration
    double adj_c;               // coefficient used for acceleration

    double rho;           // augmented Lagrangian parameter
    const double eps_abs; // absolute tolerance
    const double eps_rel; // relative tolerance

    double eps_primal;    // tolerance for primal residual
    double eps_dual;      // tolerance for dual residual

    double resid_primal;  // primal residual
    double resid_dual;    // dual residual

    virtual void A_mult (VecTypeNu &res, VecTypeBeta &x) = 0;   // operation res -> Ax, x can be overwritten
    virtual void At_mult(VecTypeNu &res, VecTypeNu &y) = 0;   // operation res -> A'y, y can be overwritten
    virtual void B_mult (VecTypeNu &res, VecTypeGamma &z) = 0;   // operation res -> Bz, z can be overwritten
    virtual double c_norm() = 0;                            // L2 norm of c

    // res = Ax + Bz - c
    virtual void next_residual(VecTypeNu &res) = 0;
    // res = x in next iteration
    virtual void next_beta(VecTypeBeta &res) = 0;
    // res = z in next iteration
    virtual void next_gamma(VecTypeGamma &res) = 0;
    virtual void next_eta(VecTypeGamma &res) = 0;
    // action when rho is changed, e.g. re-factorize matrices
    virtual void rho_changed_action() {}

    // calculating eps_primal
    // eps_primal = sqrt(p) * eps_abs + eps_rel * max(||Ax||, ||Bz||, ||c||)
    virtual double compute_eps_primal()
    {
        VecTypeNu betares, gammares, etares;
        VecTypeBeta betacopy   = main_beta;
        VecTypeGamma gammacopy = aux_gamma;
        VecTypeGamma etacopy   = aux_eta;
        A_mult(betares, betacopy);
        B_mult(gammares, gammacopy);
        B_mult(etares, etacopy);
        double r = std::max(etares.norm(), std::max(betares.norm(), gammares.norm()) );
        r = std::max(r, c_norm());
        return r * eps_rel + std::sqrt(double(dim_dual)) * eps_abs;
    }
    // calculating eps_dual
    // eps_dual = sqrt(n) * eps_abs + eps_rel * ||A'y||
    virtual double compute_eps_dual()
    {
        VecTypeNu nures, nucopy = dual_nu;

        At_mult(nures, nucopy);

        return nures.norm() * eps_rel + std::sqrt(double(dim_main)) * eps_abs;
    }
    // calculating dual residual
    // resid_dual = rho * A'B(auxz - oldz)
    virtual double compute_resid_dual()
    {
        VecTypeGamma gammadiff = aux_gamma - old_gamma;
        VecTypeGamma etadiff   = aux_eta   - old_eta;
        VecTypeNu tmp;
        B_mult(tmp, gammadiff);

        VecTypeNu tmp2;
        B_mult(tmp2, etadiff);

        VecTypeNu dual;
        At_mult(dual, tmp);

        VecTypeNu dual2;
        At_mult(dual2, tmp2);

        return rho * (dual.norm() + dual2.norm());
    }
    // calculating combined residual
    // resid_combined = rho * ||resid_primal||^2 + rho * ||auxz - adjz||^2
    virtual double compute_resid_combined()
    {
        VecTypeGamma tmp = aux_gamma - adj_gamma;
        VecTypeNu tmp2;
        B_mult(tmp2, tmp);

        VecTypeGamma tmp3 = aux_eta - adj_eta;
        VecTypeNu tmp4;
        B_mult(tmp4, tmp3);

        return rho * resid_primal * resid_primal + rho * (tmp2.squaredNorm() + tmp4.squaredNorm());
    }
    // increase or decrease rho in iterations
    virtual void update_rho()
    {
        if(resid_primal / eps_primal > 10 * resid_dual / eps_dual)
        {
            rho *= 2;
            rho_changed_action();
        }
        else if(resid_dual / eps_dual > 10 * resid_primal / eps_primal)
        {
            rho /= 2;
            rho_changed_action();
        }

        if(resid_primal < eps_primal)
        {
            rho /= 1.2;
            rho_changed_action();
        }

        if(resid_dual < eps_dual)
        {
            rho *= 1.2;
            rho_changed_action();
        }
    }
    // Debugging residual information
    void print_header(std::string title)
    {
        const int width = 80;
        const char sep = ' ';

        Rcpp::Rcout << std::endl << std::string(width, '=') << std::endl;
        Rcpp::Rcout << std::string((width - title.length()) / 2, ' ') << title << std::endl;
        Rcpp::Rcout << std::string(width, '-') << std::endl;

        Rcpp::Rcout << std::left << std::setw(7)  << std::setfill(sep) << "iter";
        Rcpp::Rcout << std::left << std::setw(13) << std::setfill(sep) << "eps_primal";
        Rcpp::Rcout << std::left << std::setw(13) << std::setfill(sep) << "resid_primal";
        Rcpp::Rcout << std::left << std::setw(13) << std::setfill(sep) << "eps_dual";
        Rcpp::Rcout << std::left << std::setw(13) << std::setfill(sep) << "resid_dual";
        Rcpp::Rcout << std::left << std::setw(13) << std::setfill(sep) << "rho";
        Rcpp::Rcout << std::endl;

        Rcpp::Rcout << std::string(width, '-') << std::endl;
    }
    void print_row(int iter)
    {
        const char sep = ' ';

        Rcpp::Rcout << std::left << std::setw(7)  << std::setfill(sep) << iter;
        Rcpp::Rcout << std::left << std::setw(13) << std::setfill(sep) << eps_primal;
        Rcpp::Rcout << std::left << std::setw(13) << std::setfill(sep) << resid_primal;
        Rcpp::Rcout << std::left << std::setw(13) << std::setfill(sep) << eps_dual;
        Rcpp::Rcout << std::left << std::setw(13) << std::setfill(sep) << resid_dual;
        Rcpp::Rcout << std::left << std::setw(13) << std::setfill(sep) << rho;
        Rcpp::Rcout << std::endl;
    }
    void print_footer()
    {
        const int width = 80;
        Rcpp::Rcout << std::string(width, '=') << std::endl << std::endl;
    }

public:
    FADMMBaseTwo(int n_, int m_, int f_, int p_,
                 double eps_abs_ = 1e-6, double eps_rel_ = 1e-6) :
        dim_main(n_),  dim_aux_gamma(m_), dim_aux_eta(f_), dim_dual(p_),
        main_beta(n_),     aux_gamma(m_),     aux_eta(f_), dual_nu(p_),  // allocate space but do not set values
        adj_gamma(m_), adj_eta(f_), adj_nu(p_),
        old_gamma(m_), old_eta(f_), old_nu(p_),
        adj_a(1.0), adj_c(9999999),
        eps_abs(eps_abs_), eps_rel(eps_rel_)
    {}

    virtual ~FADMMBaseTwo() {}

    void update_beta()
    {
        eps_primal = compute_eps_primal();
        eps_dual   = compute_eps_dual();

        VecTypeBeta newbeta(dim_main);
        next_beta(newbeta);

        main_beta.swap(newbeta);
    }
    void update_gamma()
    {
        VecTypeGamma newgamma(dim_aux_gamma);
        next_gamma(newgamma);
        aux_gamma.swap(newgamma);

        resid_dual = compute_resid_dual();
    }
    void update_eta()
    {
        VecTypeGamma neweta(dim_aux_eta);
        next_eta(neweta);
        aux_eta.swap(neweta);

        //resid_dual = compute_resid_dual();
    }
    void update_nu()
    {
        VecTypeNu newr(dim_dual);
        next_residual(newr);

        resid_primal = newr.norm();

        // dual_nu.noalias() = adj_nu + rho * newr;
        std::copy(adj_nu.data(), adj_nu.data() + dim_dual, dual_nu.data());
        dual_nu += Yscalar(rho) * newr;
        // Linalg::vec_add(dual_nu.data(), Yscalar(rho), newr.data(), dim_dual);
    }

    bool converged()
    {
        return (resid_primal < eps_primal) &&
               (resid_dual < eps_dual);
    }

    virtual int solve(int maxit)
    {
        int i;

        for(i = 0; i < maxit; i++)
        {
            old_gamma = aux_gamma;
            old_eta   = aux_eta;
            // old_nu = dual_nu;
            std::copy(dual_nu.data(), dual_nu.data() + dim_dual, old_nu.data());

            update_beta();
            update_eta();
            update_gamma();
            update_nu();

            // print_row(i);

            if(converged())
                break;

            double old_c = adj_c;
            adj_c = compute_resid_combined();

            if(adj_c < 0.999 * old_c)
            {
                double old_a = adj_a;
                adj_a            = 0.5 + 0.5 * std::sqrt(1 + 4.0 * old_a * old_a);
                double ratio     = (old_a - 1.0) / adj_a;
                // ratio = 0;
                adj_gamma        = (1 + ratio) * aux_gamma - ratio * old_gamma;
                adj_eta          = (1 + ratio) * aux_eta   - ratio * old_eta;
                adj_nu.noalias() = (1 + ratio) * dual_nu   - ratio * old_nu;
            } else {
                adj_a = 1.0;
                adj_gamma = old_gamma;
                adj_eta   = old_eta;
                // adj_nu = old_nu;
                std::copy(old_nu.data(), old_nu.data() + dim_dual, adj_nu.data());
                adj_c = old_c / 0.999;
            }
            // only update rho after a few iterations and after every 40 iterations.
            // too many updates makes it slow.
            if(i > 5 && i % 2500 == 0)
                update_rho();
        }

        // print_footer();
        return std::min(i + 1, maxit);
    }


    virtual Eigen::MatrixXd get_hessian() { return Eigen::MatrixXd(1,1); }
    virtual Eigen::VectorXd get_xty() { return Eigen::VectorXd(1); }
    virtual VecTypeBeta get_beta() { return main_beta; }
    virtual VecTypeGamma get_gamma() { return aux_gamma; }
    virtual VecTypeGamma get_aux_gamma() { return aux_gamma; }
    virtual VecTypeGamma get_eta() { return aux_eta; }
    virtual VecTypeGamma get_aux_eta() { return aux_eta; }
    virtual VecTypeNu get_nu() { return dual_nu; }

    virtual double get_loss() { return 1e99; }
    virtual double get_lambda_zero() const { return 0; }

    virtual void update_adaptive_group_weights(Eigen::VectorXd &weights_) {}
    virtual void update_adaptive_fused_weights(Eigen::VectorXd &weights_) {}
    virtual void init(double lambda_, double rho_, double lambdaF_) {}
    virtual void init_warm(double lambda_, double lambdaF_) {}

};



#endif // FADMMBASETWOPENALTIES_H
