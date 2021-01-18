#ifndef ADMMBASE_H
#define ADMMBASE_H

#include <RcppEigen.h>
#include "Linalg/BlasWrapper.h"

// General problem setting
//   minimize f(x) + g(z)
//   s.t. Ax + Bz = c
//
// x(n, 1), z(m, 1), A(p, n), B(p, m), c(p, 1)
//
template<typename VecTypeBeta, typename VecTypeGamma, typename VecTypeNu>
class ADMMBase
{
protected:
    const int dim_main;   // dimension of x
    const int dim_aux;    // dimension of z
    const int dim_dual;   // dimension of Ax + Bz - c

    VecTypeBeta main_beta;      // parameters to be optimized
    VecTypeGamma aux_gamma;       // auxiliary parameters
    VecTypeNu dual_nu;      // Lagrangian multiplier

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
    // action when rho is changed, e.g. re-factorize matrices
    virtual void rho_changed_action() {}

    // calculating eps_primal
    // eps_primal = sqrt(p) * eps_abs + eps_rel * max(||Ax||, ||Bz||, ||c||)
    virtual double compute_eps_primal()
    {
        VecTypeNu betares, gammares;
        VecTypeBeta betacopy = main_beta;
        VecTypeGamma gammacopy = aux_gamma;
        A_mult(betares, betacopy);
        B_mult(gammares, gammacopy);
        double r = std::max(betares.norm(), gammares.norm());
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
    virtual double compute_resid_dual(const VecTypeGamma &new_gamma)
    {
        VecTypeGamma zdiff = new_gamma - aux_gamma;
        VecTypeNu tmp;
        B_mult(tmp, zdiff);

        VecTypeNu dual;
        At_mult(dual, tmp);

        return rho * dual.norm();
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
    ADMMBase(int n_, int m_, int p_,
             double eps_abs_ = 1e-6, double eps_rel_ = 1e-6) :
        dim_main(n_), dim_aux(m_), dim_dual(p_),
        main_beta(n_), aux_gamma(m_), dual_nu(p_),  // allocate space but do not set values
        eps_abs(eps_abs_), eps_rel(eps_rel_)
    {}

    virtual ~ADMMBase() {}

    void update_beta()
    {
        eps_primal = compute_eps_primal();
        eps_dual = compute_eps_dual();

        VecTypeBeta newx(dim_main);
        next_beta(newx);
        main_beta.swap(newx);
    }
    void update_gamma()
    {
        VecTypeGamma newgamma(dim_aux);
        next_gamma(newgamma);

        resid_dual = compute_resid_dual(newgamma);

        aux_gamma.swap(newgamma);
    }
    void update_nu()
    {
        VecTypeNu newr(dim_dual);
        next_residual(newr);

        resid_primal = newr.norm();

        // dual_nu.noalias() += rho * newr;
        Linalg::vec_add(dual_nu.data(), typename VecTypeNu::RealScalar(rho), newr.data(), dim_dual);
    }

    bool converged()
    {
        return (resid_primal < eps_primal) &&
               (resid_dual < eps_dual);
    }

    int solve(int maxit)
    {
        int i;

        // print_header("ADMM iterations");

        for(i = 0; i < maxit; i++)
        {
            update_beta();
            update_gamma();
            update_nu();

            // print_row(i);

            if(converged())
                break;

            if(i > 3)
                update_rho();
        }

        // print_footer();

        return i + 1;
    }

    virtual VecTypeBeta get_beta() { return main_beta; }
    virtual VecTypeGamma get_gamma() { return aux_gamma; }
    virtual VecTypeGamma get_aux_gamma() { return aux_gamma; }
    virtual VecTypeNu get_nu() { return dual_nu; }

    virtual double get_loss() { return 1e99; }
    virtual double get_lambda_zero() const { return 0; }

    virtual void init(double lambda_, double rho_) {}
    virtual void init_warm(double lambda_, int iternum) {}
};



#endif // ADMMBASE_H
