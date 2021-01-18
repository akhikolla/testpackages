// Generated by rstantools.  Do not edit by hand.

/*
    bayesbr is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    bayesbr is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with bayesbr.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_bayesbr_namespace {
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
    reader.add_event(0, 0, "start", "model_bayesbr");
    reader.add_event(70, 68, "end", "model_bayesbr");
    return reader;
}
#include <stan_meta_header.hpp>
class model_bayesbr
  : public stan::model::model_base_crtp<model_bayesbr> {
private:
        int n;
        int p;
        int q;
        double a;
        double b;
        vector_d Y;
        matrix_d X;
        matrix_d W;
        vector_d mean_betas;
        vector_d variance_betas;
        vector_d mean_gammas;
        vector_d variance_gammas;
public:
    model_bayesbr(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_bayesbr(stan::io::var_context& context__,
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
        static const char* function__ = "model_bayesbr_namespace::model_bayesbr";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "n", "int", context__.to_vec());
            n = int(0);
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            n = vals_i__[pos__++];
            check_greater_or_equal(function__, "n", n, 0);
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "p", "int", context__.to_vec());
            p = int(0);
            vals_i__ = context__.vals_i("p");
            pos__ = 0;
            p = vals_i__[pos__++];
            check_greater_or_equal(function__, "p", p, 0);
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "q", "int", context__.to_vec());
            q = int(0);
            vals_i__ = context__.vals_i("q");
            pos__ = 0;
            q = vals_i__[pos__++];
            check_greater_or_equal(function__, "q", q, 0);
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "a", "double", context__.to_vec());
            a = double(0);
            vals_r__ = context__.vals_r("a");
            pos__ = 0;
            a = vals_r__[pos__++];
            check_greater_or_equal(function__, "a", a, 0);
            current_statement_begin__ = 7;
            context__.validate_dims("data initialization", "b", "double", context__.to_vec());
            b = double(0);
            vals_r__ = context__.vals_r("b");
            pos__ = 0;
            b = vals_r__[pos__++];
            check_greater_or_equal(function__, "b", b, 0);
            current_statement_begin__ = 8;
            validate_non_negative_index("Y", "n", n);
            context__.validate_dims("data initialization", "Y", "vector_d", context__.to_vec(n));
            Y = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
            vals_r__ = context__.vals_r("Y");
            pos__ = 0;
            size_t Y_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < Y_j_1_max__; ++j_1__) {
                Y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 9;
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
            current_statement_begin__ = 10;
            validate_non_negative_index("W", "n", n);
            validate_non_negative_index("W", "q", q);
            context__.validate_dims("data initialization", "W", "matrix_d", context__.to_vec(n,q));
            W = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n, q);
            vals_r__ = context__.vals_r("W");
            pos__ = 0;
            size_t W_j_2_max__ = q;
            size_t W_j_1_max__ = n;
            for (size_t j_2__ = 0; j_2__ < W_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < W_j_1_max__; ++j_1__) {
                    W(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 11;
            validate_non_negative_index("mean_betas", "p", p);
            context__.validate_dims("data initialization", "mean_betas", "vector_d", context__.to_vec(p));
            mean_betas = Eigen::Matrix<double, Eigen::Dynamic, 1>(p);
            vals_r__ = context__.vals_r("mean_betas");
            pos__ = 0;
            size_t mean_betas_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < mean_betas_j_1_max__; ++j_1__) {
                mean_betas(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 12;
            validate_non_negative_index("variance_betas", "p", p);
            context__.validate_dims("data initialization", "variance_betas", "vector_d", context__.to_vec(p));
            variance_betas = Eigen::Matrix<double, Eigen::Dynamic, 1>(p);
            vals_r__ = context__.vals_r("variance_betas");
            pos__ = 0;
            size_t variance_betas_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < variance_betas_j_1_max__; ++j_1__) {
                variance_betas(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 13;
            validate_non_negative_index("mean_gammas", "q", q);
            context__.validate_dims("data initialization", "mean_gammas", "vector_d", context__.to_vec(q));
            mean_gammas = Eigen::Matrix<double, Eigen::Dynamic, 1>(q);
            vals_r__ = context__.vals_r("mean_gammas");
            pos__ = 0;
            size_t mean_gammas_j_1_max__ = q;
            for (size_t j_1__ = 0; j_1__ < mean_gammas_j_1_max__; ++j_1__) {
                mean_gammas(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 14;
            validate_non_negative_index("variance_gammas", "q", q);
            context__.validate_dims("data initialization", "variance_gammas", "vector_d", context__.to_vec(q));
            variance_gammas = Eigen::Matrix<double, Eigen::Dynamic, 1>(q);
            vals_r__ = context__.vals_r("variance_gammas");
            pos__ = 0;
            size_t variance_gammas_j_1_max__ = q;
            for (size_t j_1__ = 0; j_1__ < variance_gammas_j_1_max__; ++j_1__) {
                variance_gammas(j_1__) = vals_r__[pos__++];
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 18;
            validate_non_negative_index("betas", "p", p);
            num_params_r__ += p;
            current_statement_begin__ = 19;
            validate_non_negative_index("gammas", "q", q);
            num_params_r__ += q;
            current_statement_begin__ = 20;
            num_params_r__ += 1;
            current_statement_begin__ = 21;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_bayesbr() { }
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
        current_statement_begin__ = 18;
        if (!(context__.contains_r("betas")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable betas missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("betas");
        pos__ = 0U;
        validate_non_negative_index("betas", "p", p);
        context__.validate_dims("parameter initialization", "betas", "vector_d", context__.to_vec(p));
        Eigen::Matrix<double, Eigen::Dynamic, 1> betas(p);
        size_t betas_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < betas_j_1_max__; ++j_1__) {
            betas(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(betas);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable betas: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 19;
        if (!(context__.contains_r("gammas")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable gammas missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("gammas");
        pos__ = 0U;
        validate_non_negative_index("gammas", "q", q);
        context__.validate_dims("parameter initialization", "gammas", "vector_d", context__.to_vec(q));
        Eigen::Matrix<double, Eigen::Dynamic, 1> gammas(q);
        size_t gammas_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gammas_j_1_max__; ++j_1__) {
            gammas(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(gammas);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable gammas: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 20;
        if (!(context__.contains_r("zeta_e")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable zeta_e missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("zeta_e");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "zeta_e", "double", context__.to_vec());
        double zeta_e(0);
        zeta_e = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, zeta_e);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable zeta_e: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 21;
        if (!(context__.contains_r("theta_e")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable theta_e missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("theta_e");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "theta_e", "double", context__.to_vec());
        double theta_e(0);
        theta_e = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0, 1, theta_e);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable theta_e: ") + e.what()), current_statement_begin__, prog_reader__());
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
            current_statement_begin__ = 18;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> betas;
            (void) betas;  // dummy to suppress unused var warning
            if (jacobian__)
                betas = in__.vector_constrain(p, lp__);
            else
                betas = in__.vector_constrain(p);
            current_statement_begin__ = 19;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> gammas;
            (void) gammas;  // dummy to suppress unused var warning
            if (jacobian__)
                gammas = in__.vector_constrain(q, lp__);
            else
                gammas = in__.vector_constrain(q);
            current_statement_begin__ = 20;
            local_scalar_t__ zeta_e;
            (void) zeta_e;  // dummy to suppress unused var warning
            if (jacobian__)
                zeta_e = in__.scalar_lb_constrain(0, lp__);
            else
                zeta_e = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 21;
            local_scalar_t__ theta_e;
            (void) theta_e;  // dummy to suppress unused var warning
            if (jacobian__)
                theta_e = in__.scalar_lub_constrain(0, 1, lp__);
            else
                theta_e = in__.scalar_lub_constrain(0, 1);
            // transformed parameters
            current_statement_begin__ = 27;
            validate_non_negative_index("theta", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> theta(n);
            stan::math::initialize(theta, DUMMY_VAR__);
            stan::math::fill(theta, DUMMY_VAR__);
            current_statement_begin__ = 28;
            validate_non_negative_index("zeta", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> zeta(n);
            stan::math::initialize(zeta, DUMMY_VAR__);
            stan::math::fill(zeta, DUMMY_VAR__);
            current_statement_begin__ = 29;
            validate_non_negative_index("lpredt", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lpredt(n);
            stan::math::initialize(lpredt, DUMMY_VAR__);
            stan::math::fill(lpredt, DUMMY_VAR__);
            current_statement_begin__ = 30;
            validate_non_negative_index("lpredz", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lpredz(n);
            stan::math::initialize(lpredz, DUMMY_VAR__);
            stan::math::fill(lpredz, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 32;
            if (as_bool(logical_neq(p, 0))) {
                current_statement_begin__ = 33;
                stan::math::assign(lpredt, multiply(X, betas));
                current_statement_begin__ = 34;
                stan::math::assign(theta, elt_divide(stan::math::exp(lpredt), add(1, stan::math::exp(lpredt))));
            }
            current_statement_begin__ = 36;
            if (as_bool(logical_neq(q, 0))) {
                current_statement_begin__ = 37;
                stan::math::assign(lpredz, multiply(W, gammas));
                current_statement_begin__ = 38;
                stan::math::assign(zeta, stan::math::exp(lpredz));
            }
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 27;
            size_t theta_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(theta(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: theta" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable theta: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 28;
            size_t zeta_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < zeta_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(zeta(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: zeta" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable zeta: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 29;
            size_t lpredt_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < lpredt_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(lpredt(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: lpredt" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable lpredt: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 30;
            size_t lpredz_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < lpredz_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(lpredz(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: lpredz" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable lpredz: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            current_statement_begin__ = 45;
            if (as_bool(logical_neq(p, 0))) {
                current_statement_begin__ = 46;
                if (as_bool(logical_neq(q, 0))) {
                    current_statement_begin__ = 47;
                    lp_accum__.add(beta_log<propto__>(Y, elt_multiply(zeta, theta), elt_multiply(zeta, subtract(1, theta))));
                } else {
                    current_statement_begin__ = 50;
                    for (int i = 1; i <= n; ++i) {
                        current_statement_begin__ = 51;
                        lp_accum__.add(beta_log<propto__>(get_base1(Y, i, "Y", 1), (zeta_e * get_base1(theta, i, "theta", 1)), (zeta_e * (1 - get_base1(theta, i, "theta", 1)))));
                    }
                }
            } else {
                current_statement_begin__ = 56;
                if (as_bool(logical_neq(p, 0))) {
                    current_statement_begin__ = 57;
                    for (int i = 1; i <= n; ++i) {
                        current_statement_begin__ = 58;
                        lp_accum__.add(beta_log<propto__>(get_base1(Y, i, "Y", 1), (get_base1(zeta, i, "zeta", 1) * theta_e), (get_base1(zeta, i, "zeta", 1) * (1 - theta_e))));
                    }
                }
            }
            current_statement_begin__ = 64;
            lp_accum__.add(normal_log<propto__>(betas, mean_betas, variance_betas));
            current_statement_begin__ = 65;
            lp_accum__.add(beta_log<propto__>(theta_e, a, b));
            current_statement_begin__ = 66;
            lp_accum__.add(normal_log<propto__>(gammas, mean_gammas, variance_gammas));
            current_statement_begin__ = 67;
            lp_accum__.add(gamma_log<propto__>(zeta_e, a, b));
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
        names__.push_back("betas");
        names__.push_back("gammas");
        names__.push_back("zeta_e");
        names__.push_back("theta_e");
        names__.push_back("theta");
        names__.push_back("zeta");
        names__.push_back("lpredt");
        names__.push_back("lpredz");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(p);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(q);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
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
        static const char* function__ = "model_bayesbr_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> betas = in__.vector_constrain(p);
        size_t betas_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < betas_j_1_max__; ++j_1__) {
            vars__.push_back(betas(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> gammas = in__.vector_constrain(q);
        size_t gammas_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gammas_j_1_max__; ++j_1__) {
            vars__.push_back(gammas(j_1__));
        }
        double zeta_e = in__.scalar_lb_constrain(0);
        vars__.push_back(zeta_e);
        double theta_e = in__.scalar_lub_constrain(0, 1);
        vars__.push_back(theta_e);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 27;
            validate_non_negative_index("theta", "n", n);
            Eigen::Matrix<double, Eigen::Dynamic, 1> theta(n);
            stan::math::initialize(theta, DUMMY_VAR__);
            stan::math::fill(theta, DUMMY_VAR__);
            current_statement_begin__ = 28;
            validate_non_negative_index("zeta", "n", n);
            Eigen::Matrix<double, Eigen::Dynamic, 1> zeta(n);
            stan::math::initialize(zeta, DUMMY_VAR__);
            stan::math::fill(zeta, DUMMY_VAR__);
            current_statement_begin__ = 29;
            validate_non_negative_index("lpredt", "n", n);
            Eigen::Matrix<double, Eigen::Dynamic, 1> lpredt(n);
            stan::math::initialize(lpredt, DUMMY_VAR__);
            stan::math::fill(lpredt, DUMMY_VAR__);
            current_statement_begin__ = 30;
            validate_non_negative_index("lpredz", "n", n);
            Eigen::Matrix<double, Eigen::Dynamic, 1> lpredz(n);
            stan::math::initialize(lpredz, DUMMY_VAR__);
            stan::math::fill(lpredz, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 32;
            if (as_bool(logical_neq(p, 0))) {
                current_statement_begin__ = 33;
                stan::math::assign(lpredt, multiply(X, betas));
                current_statement_begin__ = 34;
                stan::math::assign(theta, elt_divide(stan::math::exp(lpredt), add(1, stan::math::exp(lpredt))));
            }
            current_statement_begin__ = 36;
            if (as_bool(logical_neq(q, 0))) {
                current_statement_begin__ = 37;
                stan::math::assign(lpredz, multiply(W, gammas));
                current_statement_begin__ = 38;
                stan::math::assign(zeta, stan::math::exp(lpredz));
            }
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t theta_j_1_max__ = n;
                for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                    vars__.push_back(theta(j_1__));
                }
                size_t zeta_j_1_max__ = n;
                for (size_t j_1__ = 0; j_1__ < zeta_j_1_max__; ++j_1__) {
                    vars__.push_back(zeta(j_1__));
                }
                size_t lpredt_j_1_max__ = n;
                for (size_t j_1__ = 0; j_1__ < lpredt_j_1_max__; ++j_1__) {
                    vars__.push_back(lpredt(j_1__));
                }
                size_t lpredz_j_1_max__ = n;
                for (size_t j_1__ = 0; j_1__ < lpredz_j_1_max__; ++j_1__) {
                    vars__.push_back(lpredz(j_1__));
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
        return "model_bayesbr";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t betas_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < betas_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "betas" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t gammas_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gammas_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "gammas" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "zeta_e";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "theta_e";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t theta_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "theta" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t zeta_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < zeta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "zeta" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t lpredt_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < lpredt_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "lpredt" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t lpredz_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < lpredz_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "lpredz" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t betas_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < betas_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "betas" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t gammas_j_1_max__ = q;
        for (size_t j_1__ = 0; j_1__ < gammas_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "gammas" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "zeta_e";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "theta_e";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t theta_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "theta" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t zeta_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < zeta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "zeta" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t lpredt_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < lpredt_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "lpredt" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t lpredz_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < lpredz_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "lpredz" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_bayesbr_namespace::model_bayesbr stan_model;
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
