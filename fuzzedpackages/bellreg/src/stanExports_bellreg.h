// Generated by rstantools.  Do not edit by hand.

/*
    bellreg is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    bellreg is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with bellreg.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_bellreg_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_bellreg");
    reader.add_event(1, 1, "include", "chunks/mylib.stan");
    reader.add_event(1, 0, "start", "chunks/mylib.stan");
    reader.add_event(67, 66, "end", "chunks/mylib.stan");
    reader.add_event(67, 2, "restart", "model_bellreg");
    reader.add_event(112, 45, "end", "model_bellreg");
    return reader;
}
template <typename T0__>
typename boost::math::tools::promote_args<T0__>::type
lambertW(const T0__& x, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 6;
        local_scalar_t__ y(DUMMY_VAR__);
        (void) y;  // dummy to suppress unused var warning
        stan::math::initialize(y, DUMMY_VAR__);
        stan::math::fill(y, DUMMY_VAR__);
        current_statement_begin__ = 7;
        local_scalar_t__ w(DUMMY_VAR__);
        (void) w;  // dummy to suppress unused var warning
        stan::math::initialize(w, DUMMY_VAR__);
        stan::math::fill(w, DUMMY_VAR__);
        current_statement_begin__ = 8;
        stan::math::assign(y, stan::math::sqrt((1 + (stan::math::exp(1) * x))));
        current_statement_begin__ = 9;
        stan::math::assign(w, (-(1) + (2.036 * stan::math::log(((1 + (1.14956131 * y)) / (1 + (0.45495740 * stan::math::log((1 + y)))))))));
        current_statement_begin__ = 10;
        stan::math::assign(w, ((w / (1 + w)) * (1 + stan::math::log((x / w)))));
        current_statement_begin__ = 11;
        stan::math::assign(w, ((w / (1 + w)) * (1 + stan::math::log((x / w)))));
        current_statement_begin__ = 12;
        stan::math::assign(w, ((w / (1 + w)) * (1 + stan::math::log((x / w)))));
        current_statement_begin__ = 13;
        return stan::math::promote_scalar<fun_return_scalar_t__>(w);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct lambertW_functor__ {
    template <typename T0__>
        typename boost::math::tools::promote_args<T0__>::type
    operator()(const T0__& x, std::ostream* pstream__) const {
        return lambertW(x, pstream__);
    }
};
double
bellnumber(const int& n, std::ostream* pstream__) {
    typedef double local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        current_statement_begin__ = 17;
        if (as_bool(logical_lt(n, 2))) {
            current_statement_begin__ = 18;
            return stan::math::promote_scalar<fun_return_scalar_t__>(1);
        } else {
            {
            current_statement_begin__ = 20;
            int k(0);
            (void) k;  // dummy to suppress unused var warning
            stan::math::fill(k, std::numeric_limits<int>::min());
            current_statement_begin__ = 21;
            validate_non_negative_index("B", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> B(n);
            stan::math::initialize(B, DUMMY_VAR__);
            stan::math::fill(B, DUMMY_VAR__);
            current_statement_begin__ = 22;
            validate_non_negative_index("Bneu", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> Bneu(n);
            stan::math::initialize(Bneu, DUMMY_VAR__);
            stan::math::fill(Bneu, DUMMY_VAR__);
            current_statement_begin__ = 23;
            stan::model::assign(B, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        1, 
                        "assigning variable B");
            current_statement_begin__ = 24;
            for (int i = 1; i <= (n - 1); ++i) {
                current_statement_begin__ = 25;
                stan::math::assign(k, i);
                current_statement_begin__ = 26;
                stan::model::assign(Bneu, 
                            stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                            get_base1(B, i, "B", 1), 
                            "assigning variable Bneu");
                current_statement_begin__ = 27;
                for (int j = 2; j <= (i + 1); ++j) {
                    current_statement_begin__ = 28;
                    stan::model::assign(Bneu, 
                                stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), 
                                (get_base1(B, (j - 1), "B", 1) + get_base1(Bneu, (j - 1), "Bneu", 1)), 
                                "assigning variable Bneu");
                }
                current_statement_begin__ = 30;
                for (int j = 1; j <= n; ++j) {
                    current_statement_begin__ = 31;
                    stan::model::assign(B, 
                                stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), 
                                get_base1(Bneu, j, "Bneu", 1), 
                                "assigning variable B");
                }
            }
            current_statement_begin__ = 34;
            return stan::math::promote_scalar<fun_return_scalar_t__>(get_base1(Bneu, (k + 1), "Bneu", 1));
            }
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct bellnumber_functor__ {
            double
    operator()(const int& n, std::ostream* pstream__) const {
        return bellnumber(n, pstream__);
    }
};
template <bool propto, typename T1__>
typename boost::math::tools::promote_args<T1__>::type
bell_lpmf(const int& x,
              const T1__& theta, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T1__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 40;
        local_scalar_t__ Bx(DUMMY_VAR__);
        (void) Bx;  // dummy to suppress unused var warning
        stan::math::initialize(Bx, DUMMY_VAR__);
        stan::math::fill(Bx, DUMMY_VAR__);
        current_statement_begin__ = 41;
        local_scalar_t__ lprob(DUMMY_VAR__);
        (void) lprob;  // dummy to suppress unused var warning
        stan::math::initialize(lprob, DUMMY_VAR__);
        stan::math::fill(lprob, DUMMY_VAR__);
        current_statement_begin__ = 42;
        stan::math::assign(Bx, bellnumber(x, pstream__));
        current_statement_begin__ = 43;
        stan::math::assign(lprob, (((((x * stan::math::log(theta)) - stan::math::exp(theta)) + 1) + stan::math::log(Bx)) - stan::math::lgamma((x + 1))));
        current_statement_begin__ = 44;
        return stan::math::promote_scalar<fun_return_scalar_t__>(lprob);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
template <typename T1__>
typename boost::math::tools::promote_args<T1__>::type
bell_lpmf(const int& x,
              const T1__& theta, std::ostream* pstream__) {
    return bell_lpmf<false>(x,theta, pstream__);
}
struct bell_lpmf_functor__ {
    template <bool propto, typename T1__>
        typename boost::math::tools::promote_args<T1__>::type
    operator()(const int& x,
              const T1__& theta, std::ostream* pstream__) const {
        return bell_lpmf(x, theta, pstream__);
    }
};
template <typename T1__>
typename boost::math::tools::promote_args<T1__>::type
loglik_bell(const std::vector<int>& x,
                const std::vector<T1__>& theta, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T1__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 59;
        validate_non_negative_index("lprob", "num_elements(x)", num_elements(x));
        std::vector<local_scalar_t__  > lprob(num_elements(x), local_scalar_t__(DUMMY_VAR__));
        stan::math::initialize(lprob, DUMMY_VAR__);
        stan::math::fill(lprob, DUMMY_VAR__);
        current_statement_begin__ = 60;
        for (int i = 1; i <= num_elements(x); ++i) {
            current_statement_begin__ = 61;
            stan::model::assign(lprob, 
                        stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                        ((get_base1(x, i, "x", 1) * stan::math::log(get_base1(theta, i, "theta", 1))) - stan::math::exp(get_base1(theta, i, "theta", 1))), 
                        "assigning variable lprob");
        }
        current_statement_begin__ = 63;
        return stan::math::promote_scalar<fun_return_scalar_t__>(sum(lprob));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct loglik_bell_functor__ {
    template <typename T1__>
        typename boost::math::tools::promote_args<T1__>::type
    operator()(const std::vector<int>& x,
                const std::vector<T1__>& theta, std::ostream* pstream__) const {
        return loglik_bell(x, theta, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_bellreg
  : public stan::model::model_base_crtp<model_bellreg> {
private:
        int n;
        int p;
        std::vector<int> y;
        matrix_d X;
        row_vector_d x_mean;
        vector_d x_sd;
        int approach;
        double mu_beta;
        double sigma_beta;
public:
    model_bellreg(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_bellreg(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_bellreg_namespace::model_bellreg";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 70;
            context__.validate_dims("data initialization", "n", "int", context__.to_vec());
            n = int(0);
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            n = vals_i__[pos__++];
            check_greater_or_equal(function__, "n", n, 1);
            current_statement_begin__ = 71;
            context__.validate_dims("data initialization", "p", "int", context__.to_vec());
            p = int(0);
            vals_i__ = context__.vals_i("p");
            pos__ = 0;
            p = vals_i__[pos__++];
            check_greater_or_equal(function__, "p", p, 1);
            current_statement_begin__ = 72;
            validate_non_negative_index("y", "n", n);
            context__.validate_dims("data initialization", "y", "int", context__.to_vec(n));
            y = std::vector<int>(n, int(0));
            vals_i__ = context__.vals_i("y");
            pos__ = 0;
            size_t y_k_0_max__ = n;
            for (size_t k_0__ = 0; k_0__ < y_k_0_max__; ++k_0__) {
                y[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 73;
            validate_non_negative_index("X", "n", n);
            validate_non_negative_index("X", "p", p);
            context__.validate_dims("data initialization", "X", "matrix_d", context__.to_vec(n,p));
            X = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n, p);
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_j_2_max__ = p;
            size_t X_j_1_max__ = n;
            for (size_t j_2__ = 0; j_2__ < X_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_j_1_max__; ++j_1__) {
                    X(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 74;
            validate_non_negative_index("x_mean", "p", p);
            context__.validate_dims("data initialization", "x_mean", "row_vector_d", context__.to_vec(p));
            x_mean = Eigen::Matrix<double, 1, Eigen::Dynamic>(p);
            vals_r__ = context__.vals_r("x_mean");
            pos__ = 0;
            size_t x_mean_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < x_mean_j_1_max__; ++j_1__) {
                x_mean(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 75;
            validate_non_negative_index("x_sd", "p", p);
            context__.validate_dims("data initialization", "x_sd", "vector_d", context__.to_vec(p));
            x_sd = Eigen::Matrix<double, Eigen::Dynamic, 1>(p);
            vals_r__ = context__.vals_r("x_sd");
            pos__ = 0;
            size_t x_sd_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < x_sd_j_1_max__; ++j_1__) {
                x_sd(j_1__) = vals_r__[pos__++];
            }
            check_greater_or_equal(function__, "x_sd", x_sd, 0);
            current_statement_begin__ = 76;
            context__.validate_dims("data initialization", "approach", "int", context__.to_vec());
            approach = int(0);
            vals_i__ = context__.vals_i("approach");
            pos__ = 0;
            approach = vals_i__[pos__++];
            check_greater_or_equal(function__, "approach", approach, 0);
            check_less_or_equal(function__, "approach", approach, 1);
            current_statement_begin__ = 77;
            context__.validate_dims("data initialization", "mu_beta", "double", context__.to_vec());
            mu_beta = double(0);
            vals_r__ = context__.vals_r("mu_beta");
            pos__ = 0;
            mu_beta = vals_r__[pos__++];
            current_statement_begin__ = 78;
            context__.validate_dims("data initialization", "sigma_beta", "double", context__.to_vec());
            sigma_beta = double(0);
            vals_r__ = context__.vals_r("sigma_beta");
            pos__ = 0;
            sigma_beta = vals_r__[pos__++];
            check_greater_or_equal(function__, "sigma_beta", sigma_beta, 0);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 82;
            validate_non_negative_index("beta_std", "p", p);
            num_params_r__ += p;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_bellreg() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 82;
        if (!(context__.contains_r("beta_std")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta_std missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta_std");
        pos__ = 0U;
        validate_non_negative_index("beta_std", "p", p);
        context__.validate_dims("parameter initialization", "beta_std", "vector_d", context__.to_vec(p));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_std(p);
        size_t beta_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_std_j_1_max__; ++j_1__) {
            beta_std(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta_std);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta_std: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 82;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_std;
            (void) beta_std;  // dummy to suppress unused var warning
            if (jacobian__)
                beta_std = in__.vector_constrain(p, lp__);
            else
                beta_std = in__.vector_constrain(p);
            // transformed parameters
            current_statement_begin__ = 87;
            validate_non_negative_index("beta", "p", p);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta(p);
            stan::math::initialize(beta, DUMMY_VAR__);
            stan::math::fill(beta, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 88;
            if (as_bool(logical_eq(p, 1))) {
                current_statement_begin__ = 89;
                stan::model::assign(beta, 
                            stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                            (get_base1(beta_std, 1, "beta_std", 1) / get_base1(x_sd, 1, "x_sd", 1)), 
                            "assigning variable beta");
            } else {
                current_statement_begin__ = 91;
                stan::model::assign(beta, 
                            stan::model::cons_list(stan::model::index_min_max(2, p), stan::model::nil_index_list()), 
                            elt_divide(stan::model::rvalue(beta_std, stan::model::cons_list(stan::model::index_min_max(2, p), stan::model::nil_index_list()), "beta_std"), stan::model::rvalue(x_sd, stan::model::cons_list(stan::model::index_min_max(2, p), stan::model::nil_index_list()), "x_sd")), 
                            "assigning variable beta");
                current_statement_begin__ = 92;
                stan::model::assign(beta, 
                            stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                            ((get_base1(beta_std, 1, "beta_std", 1) / get_base1(x_sd, 1, "x_sd", 1)) - multiply(stan::model::rvalue(x_mean, stan::model::cons_list(stan::model::index_min_max(2, p), stan::model::nil_index_list()), "x_mean"), stan::model::rvalue(beta, stan::model::cons_list(stan::model::index_min_max(2, p), stan::model::nil_index_list()), "beta"))), 
                            "assigning variable beta");
            }
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 87;
            size_t beta_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(beta(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: beta" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable beta: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            {
            current_statement_begin__ = 97;
            validate_non_negative_index("theta", "n", n);
            std::vector<local_scalar_t__  > theta(n, local_scalar_t__(DUMMY_VAR__));
            stan::math::initialize(theta, DUMMY_VAR__);
            stan::math::fill(theta, DUMMY_VAR__);
            current_statement_begin__ = 98;
            validate_non_negative_index("eta", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> eta(n);
            stan::math::initialize(eta, DUMMY_VAR__);
            stan::math::fill(eta, DUMMY_VAR__);
            stan::math::assign(eta,multiply(X, beta_std));
            current_statement_begin__ = 99;
            validate_non_negative_index("mu", "n", n);
            std::vector<local_scalar_t__  > mu(n, local_scalar_t__(DUMMY_VAR__));
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 100;
            for (int i = 1; i <= n; ++i) {
                current_statement_begin__ = 101;
                stan::model::assign(mu, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            stan::math::exp(get_base1(eta, i, "eta", 1)), 
                            "assigning variable mu");
                current_statement_begin__ = 102;
                stan::model::assign(theta, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            lambertW(get_base1(mu, i, "mu", 1), pstream__), 
                            "assigning variable theta");
            }
            current_statement_begin__ = 105;
            lp_accum__.add(loglik_bell(y, theta, pstream__));
            current_statement_begin__ = 106;
            if (as_bool(logical_eq(approach, 1))) {
                current_statement_begin__ = 107;
                lp_accum__.add(normal_log<propto__>(beta_std, mu_beta, sigma_beta));
            }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("beta_std");
        names__.push_back("beta");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(p);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(p);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_bellreg_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_std = in__.vector_constrain(p);
        size_t beta_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_std_j_1_max__; ++j_1__) {
            vars__.push_back(beta_std(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 87;
            validate_non_negative_index("beta", "p", p);
            Eigen::Matrix<double, Eigen::Dynamic, 1> beta(p);
            stan::math::initialize(beta, DUMMY_VAR__);
            stan::math::fill(beta, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 88;
            if (as_bool(logical_eq(p, 1))) {
                current_statement_begin__ = 89;
                stan::model::assign(beta, 
                            stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                            (get_base1(beta_std, 1, "beta_std", 1) / get_base1(x_sd, 1, "x_sd", 1)), 
                            "assigning variable beta");
            } else {
                current_statement_begin__ = 91;
                stan::model::assign(beta, 
                            stan::model::cons_list(stan::model::index_min_max(2, p), stan::model::nil_index_list()), 
                            elt_divide(stan::model::rvalue(beta_std, stan::model::cons_list(stan::model::index_min_max(2, p), stan::model::nil_index_list()), "beta_std"), stan::model::rvalue(x_sd, stan::model::cons_list(stan::model::index_min_max(2, p), stan::model::nil_index_list()), "x_sd")), 
                            "assigning variable beta");
                current_statement_begin__ = 92;
                stan::model::assign(beta, 
                            stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                            ((get_base1(beta_std, 1, "beta_std", 1) / get_base1(x_sd, 1, "x_sd", 1)) - multiply(stan::model::rvalue(x_mean, stan::model::cons_list(stan::model::index_min_max(2, p), stan::model::nil_index_list()), "x_mean"), stan::model::rvalue(beta, stan::model::cons_list(stan::model::index_min_max(2, p), stan::model::nil_index_list()), "beta"))), 
                            "assigning variable beta");
            }
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t beta_j_1_max__ = p;
                for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
                    vars__.push_back(beta(j_1__));
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_bellreg";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_std_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_std" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t beta_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_std_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_std_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_std" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t beta_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_bellreg_namespace::model_bellreg stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
