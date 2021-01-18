#ifndef FADMMBASE_H
#define FADMMBASE_H

#include <RcppEigen.h>
#include "Linalg/BlasWrapper.h"

// [[Rcpp::plugins(cpp17)]]

// General problem setting
//   minimize f(x) + g(z)
//   s.t. Ax + Bz = c
//
// x(n, 1), z(m, 1), A(p, n), B(p, m), c(p, 1)
//
template<typename VecTypeBeta, typename VecTypeGamma, typename VecTypeNu>
class FADMMBase
{
protected:
    typedef typename VecTypeNu::RealScalar Yscalar;
    
    double eps_primal;    // tolerance for primal residual
    double eps_dual;      // tolerance for dual residual
    
    double resid_primal;  // primal residual
    double resid_dual;    // dual residual

    const int dim_main;   // dimension of x
    const int dim_aux;    // dimension of z
    const int dim_dual;   // dimension of Ax + Bz - c

    VecTypeBeta main_beta;      // parameters to be optimized
    VecTypeGamma aux_gamma;       // auxiliary parameters
    VecTypeNu dual_nu;      // Lagrangian multiplier

    VecTypeGamma adj_gamma;       // adjusted z vector, used for acceleration
    VecTypeNu adj_nu;       // adjusted y vector, used for acceleration
    VecTypeGamma old_gamma;       // z vector in the previous iteration, used for acceleration
    VecTypeNu old_nu;       // y vector in the previous iteration, used for acceleration
    double adj_a;         // coefficient used for acceleration
    double adj_c;         // coefficient used for acceleration

    double rho;           // augmented Lagrangian parameter
    const double eps_abs; // absolute tolerance
    const double eps_rel; // relative tolerance

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
    
    
    virtual bool converged()
    {
        // std::cout << "resid_primal:" << resid_primal << "eps_primal:" << eps_primal <<
        //     "resid_dual:" << resid_dual << "eps_dual:" << eps_dual << std::endl;
        // 
        // if (std::isinf(resid_primal) || std::isinf(eps_primal) || std::isinf(resid_dual) ||
        //     std::isinf(eps_dual) || resid_primal != resid_primal || eps_primal != eps_primal ||
        //     resid_dual != resid_dual || eps_dual != eps_dual)
        // {
        //     bool is_conv = false;
        //     return is_conv;
        // }
            
        bool is_conv = (resid_primal < eps_primal) && (resid_dual < eps_dual);
        return is_conv;
    }

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
    virtual double compute_resid_dual()
    {
        VecTypeGamma gammadiff = aux_gamma - old_gamma;
        VecTypeNu tmp;
        B_mult(tmp, gammadiff);

        VecTypeNu dual;
        At_mult(dual, tmp);

        return rho * dual.norm();
    }
    // calculating combined residual
    // resid_combined = rho * ||resid_primal||^2 + rho * ||auxz - adjz||^2
    virtual double compute_resid_combined()
    {
        VecTypeGamma tmp = aux_gamma - adj_gamma;
        VecTypeNu tmp2;
        B_mult(tmp2, tmp);

        return rho * resid_primal * resid_primal + rho * tmp2.squaredNorm();
    }
    // increase or decrease rho in iterations
    virtual void update_rho()
    {
        if(resid_primal / eps_primal > 10.0 * resid_dual / eps_dual)
        {
            rho *= 2.0;
            rho_changed_action();
        }
        else if(resid_dual / eps_dual > 10.0 * resid_primal / eps_primal)
        {
            rho /= 2.0;
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
    FADMMBase(int n_, int m_, int p_,
              double eps_abs_ = 1e-6, double eps_rel_ = 1e-6,
              double smallval_ = 1e-15, double bigval_ = 1e99) :
        eps_primal(smallval_), eps_dual(smallval_),
        resid_primal(bigval_), resid_dual(bigval_),
        dim_main(n_), dim_aux(m_), dim_dual(p_),
        main_beta(n_), aux_gamma(m_), dual_nu(p_),  // allocate space but do not set values
        adj_gamma(m_), adj_nu(p_),
        old_gamma(m_), old_nu(p_),
        adj_a(1.0), adj_c(bigval_),
        eps_abs(eps_abs_), eps_rel(eps_rel_)
        
    {}

    virtual ~FADMMBase() {}

    void update_beta()
    {
        eps_primal = compute_eps_primal();
        eps_dual = compute_eps_dual();

        VecTypeBeta newbeta(dim_main);
        next_beta(newbeta);

        main_beta.swap(newbeta);
    }
    void update_gamma()
    {
        VecTypeGamma newgamma(dim_aux);
        next_gamma(newgamma);
        aux_gamma.swap(newgamma);

        resid_dual = compute_resid_dual();
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


    virtual int solve(int maxit)
    {
        int i = 0;
        
        for(i = 0; i < maxit; i++)
        {
            old_gamma = aux_gamma;
            // old_nu = dual_nu;
            std::copy(dual_nu.data(), dual_nu.data() + dim_dual, old_nu.data());
            
            update_beta();
            update_gamma();
            update_nu();
            
            // print_row(i);
            
            if (i > 0)
            {
                if(converged())
                    break;
            }
            
            double old_c = adj_c;
            adj_c = compute_resid_combined();
            
            if(adj_c < 0.999 * old_c)
            {
                double old_a = adj_a;
                adj_a = 0.5 + 0.5 * std::sqrt(1 + 4.0 * old_a * old_a);
                double ratio = (old_a - 1.0) / adj_a;
                adj_gamma = (1 + ratio) * aux_gamma - ratio * old_gamma;
                adj_nu.noalias() = (1 + ratio) * dual_nu - ratio * old_nu;
            } else {
                adj_a = 1.0;
                adj_gamma = old_gamma;
                // adj_nu = old_nu;
                std::copy(old_nu.data(), old_nu.data() + dim_dual, adj_nu.data());
                adj_c = old_c / 0.999;
            }
            // only update rho after a few iterations and after every 40 iterations.
            // too many updates makes it slow.
            if(i > 5 && i % 500 == 0)
                update_rho();
        }
        
        // print_footer();
        
        return i + 1;
    }


    virtual Eigen::MatrixXd get_hessian() { return Eigen::MatrixXd(1,1); }
    virtual Eigen::VectorXd get_xty() { return Eigen::VectorXd(1); }
    virtual VecTypeBeta get_beta() { return main_beta; }
    virtual VecTypeGamma get_gamma() { return aux_gamma; }
    virtual VecTypeGamma get_aux_gamma() { return aux_gamma; }
    virtual VecTypeNu get_nu() { return dual_nu; }

    virtual double get_loss() { return 1e99; }
    virtual double get_lambda_zero() const { return 0; }
    virtual int get_nselected(VecTypeBeta &beta_vector) { return 0; }

    virtual void update_adaptive_group_weights(Eigen::VectorXd &weights_) {}
    virtual void init(double lambda_, double rho_) {}
    virtual void init_warm(double lambda_) {}

};



#endif // FADMMBASE_H
