// Generated by rstantools.  Do not edit by hand.

/*
    qmix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    qmix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with qmix.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_qmixbr_namespace {
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
    reader.add_event(0, 0, "start", "model_qmixbr");
    reader.add_event(78, 76, "end", "model_qmixbr");
    return reader;
}
template <typename T0__, typename T1__>
typename boost::math::tools::promote_args<T0__, T1__>::type
pald2(const T0__& mu,
          const T1__& p, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 4;
        local_scalar_t__ prob(DUMMY_VAR__);
        (void) prob;  // dummy to suppress unused var warning
        stan::math::initialize(prob, DUMMY_VAR__);
        stan::math::fill(prob, DUMMY_VAR__);
        current_statement_begin__ = 5;
        if (as_bool(logical_lt(mu, 0))) {
            current_statement_begin__ = 6;
            stan::math::assign(prob, (p * stan::math::exp((mu * (1 - p)))));
        } else {
            current_statement_begin__ = 8;
            stan::math::assign(prob, (1 - ((1 - p) * stan::math::exp((-(mu) * p)))));
        }
        current_statement_begin__ = 10;
        return stan::math::promote_scalar<fun_return_scalar_t__>(prob);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct pald2_functor__ {
    template <typename T0__, typename T1__>
        typename boost::math::tools::promote_args<T0__, T1__>::type
    operator()(const T0__& mu,
          const T1__& p, std::ostream* pstream__) const {
        return pald2(mu, p, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_qmixbr
  : public stan::model::model_base_crtp<model_qmixbr> {
private:
        int N;
        int D;
        vector_d Y;
        matrix_d X;
        int k;
        double offset;
public:
    model_qmixbr(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_qmixbr(stan::io::var_context& context__,
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
        static const char* function__ = "model_qmixbr_namespace::model_qmixbr";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 18;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 19;
            context__.validate_dims("data initialization", "D", "int", context__.to_vec());
            D = int(0);
            vals_i__ = context__.vals_i("D");
            pos__ = 0;
            D = vals_i__[pos__++];
            current_statement_begin__ = 20;
            validate_non_negative_index("Y", "N", N);
            context__.validate_dims("data initialization", "Y", "vector_d", context__.to_vec(N));
            Y = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("Y");
            pos__ = 0;
            size_t Y_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < Y_j_1_max__; ++j_1__) {
                Y(j_1__) = vals_r__[pos__++];
            }
            check_greater_or_equal(function__, "Y", Y, 0);
            check_less_or_equal(function__, "Y", Y, 1);
            current_statement_begin__ = 21;
            validate_non_negative_index("X", "N", N);
            validate_non_negative_index("X", "D", D);
            context__.validate_dims("data initialization", "X", "matrix_d", context__.to_vec(N,D));
            X = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, D);
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_j_2_max__ = D;
            size_t X_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < X_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_j_1_max__; ++j_1__) {
                    X(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 22;
            context__.validate_dims("data initialization", "k", "int", context__.to_vec());
            k = int(0);
            vals_i__ = context__.vals_i("k");
            pos__ = 0;
            k = vals_i__[pos__++];
            current_statement_begin__ = 23;
            context__.validate_dims("data initialization", "offset", "double", context__.to_vec());
            offset = double(0);
            vals_r__ = context__.vals_r("offset");
            pos__ = 0;
            offset = vals_r__[pos__++];
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 33;
            validate_non_negative_index("beta", "D", D);
            validate_non_negative_index("beta", "k", k);
            num_params_r__ += (D * k);
            current_statement_begin__ = 34;
            validate_non_negative_index("p_tmp", "k", k);
            num_params_r__ += k;
            current_statement_begin__ = 35;
            validate_non_negative_index("theta", "k", k);
            num_params_r__ += (k - 1);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_qmixbr() { }
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
        current_statement_begin__ = 33;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "D", D);
        validate_non_negative_index("beta", "k", k);
        context__.validate_dims("parameter initialization", "beta", "vector_d", context__.to_vec(k,D));
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > beta(k, Eigen::Matrix<double, Eigen::Dynamic, 1>(D));
        size_t beta_j_1_max__ = D;
        size_t beta_k_0_max__ = k;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < beta_k_0_max__; ++k_0__) {
                beta[k_0__](j_1__) = vals_r__[pos__++];
            }
        }
        size_t beta_i_0_max__ = k;
        for (size_t i_0__ = 0; i_0__ < beta_i_0_max__; ++i_0__) {
            try {
                writer__.vector_unconstrain(beta[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        current_statement_begin__ = 34;
        if (!(context__.contains_r("p_tmp")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable p_tmp missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("p_tmp");
        pos__ = 0U;
        validate_non_negative_index("p_tmp", "k", k);
        context__.validate_dims("parameter initialization", "p_tmp", "vector_d", context__.to_vec(k));
        Eigen::Matrix<double, Eigen::Dynamic, 1> p_tmp(k);
        size_t p_tmp_j_1_max__ = k;
        for (size_t j_1__ = 0; j_1__ < p_tmp_j_1_max__; ++j_1__) {
            p_tmp(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_lub_unconstrain(0, 1, p_tmp);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable p_tmp: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 35;
        if (!(context__.contains_r("theta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable theta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("theta");
        pos__ = 0U;
        validate_non_negative_index("theta", "k", k);
        context__.validate_dims("parameter initialization", "theta", "vector_d", context__.to_vec(k));
        Eigen::Matrix<double, Eigen::Dynamic, 1> theta(k);
        size_t theta_j_1_max__ = k;
        for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
            theta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.simplex_unconstrain(theta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable theta: ") + e.what()), current_statement_begin__, prog_reader__());
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
            current_statement_begin__ = 33;
            std::vector<Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> > beta;
            size_t beta_d_0_max__ = k;
            beta.reserve(beta_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < beta_d_0_max__; ++d_0__) {
                if (jacobian__)
                    beta.push_back(in__.vector_constrain(D, lp__));
                else
                    beta.push_back(in__.vector_constrain(D));
            }
            current_statement_begin__ = 34;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> p_tmp;
            (void) p_tmp;  // dummy to suppress unused var warning
            if (jacobian__)
                p_tmp = in__.vector_lub_constrain(0, 1, k, lp__);
            else
                p_tmp = in__.vector_lub_constrain(0, 1, k);
            current_statement_begin__ = 35;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> theta;
            (void) theta;  // dummy to suppress unused var warning
            if (jacobian__)
                theta = in__.simplex_constrain(k, lp__);
            else
                theta = in__.simplex_constrain(k);
            // transformed parameters
            current_statement_begin__ = 40;
            validate_non_negative_index("p", "k", k);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> p(k);
            stan::math::initialize(p, DUMMY_VAR__);
            stan::math::fill(p, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 41;
            stan::math::assign(p, sort_asc(p_tmp));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 40;
            size_t p_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < p_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(p(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: p" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable p: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            {
            current_statement_begin__ = 46;
            local_scalar_t__ lik(DUMMY_VAR__);
            (void) lik;  // dummy to suppress unused var warning
            stan::math::initialize(lik, DUMMY_VAR__);
            stan::math::fill(lik, DUMMY_VAR__);
            current_statement_begin__ = 47;
            for (int i = 1; i <= k; ++i) {
                current_statement_begin__ = 48;
                lp_accum__.add(normal_log<propto__>(get_base1(beta, i, "beta", 1), 0, 10));
            }
            current_statement_begin__ = 51;
            lp_accum__.add(dirichlet_log<propto__>(theta, rep_vector(1.0, k)));
            current_statement_begin__ = 52;
            lp_accum__.add(uniform_log<propto__>(p_tmp, 0, 1));
            current_statement_begin__ = 56;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 57;
                if (as_bool(logical_eq(get_base1(Y, i, "Y", 1), 1))) {
                    current_statement_begin__ = 58;
                    stan::math::assign(lik, 0);
                    current_statement_begin__ = 59;
                    for (int j = 1; j <= k; ++j) {
                        current_statement_begin__ = 60;
                        stan::math::assign(lik, (lik + (get_base1(theta, j, "theta", 1) * pald2(dot_product(stan::model::rvalue(X, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "X"), get_base1(beta, j, "beta", 1)), get_base1(p, j, "p", 1), pstream__))));
                    }
                    current_statement_begin__ = 62;
                    stan::math::assign(lik, (lik + offset));
                }
                current_statement_begin__ = 64;
                if (as_bool(logical_eq(get_base1(Y, i, "Y", 1), 0))) {
                    current_statement_begin__ = 65;
                    stan::math::assign(lik, 0);
                    current_statement_begin__ = 66;
                    for (int j = 1; j <= k; ++j) {
                        current_statement_begin__ = 67;
                        stan::math::assign(lik, (lik + (get_base1(theta, j, "theta", 1) * (1 - pald2(dot_product(stan::model::rvalue(X, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "X"), get_base1(beta, j, "beta", 1)), get_base1(p, j, "p", 1), pstream__)))));
                    }
                    current_statement_begin__ = 69;
                    stan::math::assign(lik, (lik + offset));
                }
                current_statement_begin__ = 72;
                lp_accum__.add(stan::math::log(lik));
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
        names__.push_back("beta");
        names__.push_back("p_tmp");
        names__.push_back("theta");
        names__.push_back("p");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(k);
        dims__.push_back(D);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(k);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(k);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(k);
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
        static const char* function__ = "model_qmixbr_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > beta;
        size_t beta_d_0_max__ = k;
        beta.reserve(beta_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < beta_d_0_max__; ++d_0__) {
            beta.push_back(in__.vector_constrain(D));
        }
        size_t beta_j_1_max__ = D;
        size_t beta_k_0_max__ = k;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < beta_k_0_max__; ++k_0__) {
                vars__.push_back(beta[k_0__](j_1__));
            }
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> p_tmp = in__.vector_lub_constrain(0, 1, k);
        size_t p_tmp_j_1_max__ = k;
        for (size_t j_1__ = 0; j_1__ < p_tmp_j_1_max__; ++j_1__) {
            vars__.push_back(p_tmp(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> theta = in__.simplex_constrain(k);
        size_t theta_j_1_max__ = k;
        for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
            vars__.push_back(theta(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 40;
            validate_non_negative_index("p", "k", k);
            Eigen::Matrix<double, Eigen::Dynamic, 1> p(k);
            stan::math::initialize(p, DUMMY_VAR__);
            stan::math::fill(p, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 41;
            stan::math::assign(p, sort_asc(p_tmp));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t p_j_1_max__ = k;
                for (size_t j_1__ = 0; j_1__ < p_j_1_max__; ++j_1__) {
                    vars__.push_back(p(j_1__));
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
        return "model_qmixbr";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = D;
        size_t beta_k_0_max__ = k;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < beta_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta" << '.' << k_0__ + 1 << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        size_t p_tmp_j_1_max__ = k;
        for (size_t j_1__ = 0; j_1__ < p_tmp_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "p_tmp" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t theta_j_1_max__ = k;
        for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t p_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < p_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "p" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = D;
        size_t beta_k_0_max__ = k;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < beta_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta" << '.' << k_0__ + 1 << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        size_t p_tmp_j_1_max__ = k;
        for (size_t j_1__ = 0; j_1__ < p_tmp_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "p_tmp" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t theta_j_1_max__ = (k - 1);
        for (size_t j_1__ = 0; j_1__ < theta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t p_j_1_max__ = k;
            for (size_t j_1__ = 0; j_1__ < p_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "p" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_qmixbr_namespace::model_qmixbr stan_model;
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
