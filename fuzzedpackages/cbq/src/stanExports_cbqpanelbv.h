// Generated by rstantools.  Do not edit by hand.

/*
    cbq is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cbq is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cbq.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_cbqpanelbv_namespace {
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
    reader.add_event(0, 0, "start", "model_cbqpanelbv");
    reader.add_event(56, 54, "end", "model_cbqpanelbv");
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
        current_statement_begin__ = 3;
        local_scalar_t__ prob(DUMMY_VAR__);
        (void) prob;  // dummy to suppress unused var warning
        stan::math::initialize(prob, DUMMY_VAR__);
        stan::math::fill(prob, DUMMY_VAR__);
        current_statement_begin__ = 4;
        if (as_bool(logical_lt(mu, 0))) {
            current_statement_begin__ = 5;
            stan::math::assign(prob, (p * stan::math::exp((mu * (1 - p)))));
        } else {
            current_statement_begin__ = 7;
            stan::math::assign(prob, (1 - ((1 - p) * stan::math::exp((-(mu) * p)))));
        }
        current_statement_begin__ = 9;
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
class model_cbqpanelbv
  : public stan::model::model_base_crtp<model_cbqpanelbv> {
private:
        int N;
        int D;
        vector_d Y;
        matrix_d X;
        double offset;
        double q;
        int N_person;
        std::vector<int> person;
        int N_wave;
        std::vector<int> wave;
public:
    model_cbqpanelbv(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_cbqpanelbv(stan::io::var_context& context__,
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
        static const char* function__ = "model_cbqpanelbv_namespace::model_cbqpanelbv";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 14;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 15;
            context__.validate_dims("data initialization", "D", "int", context__.to_vec());
            D = int(0);
            vals_i__ = context__.vals_i("D");
            pos__ = 0;
            D = vals_i__[pos__++];
            current_statement_begin__ = 16;
            validate_non_negative_index("Y", "N", N);
            context__.validate_dims("data initialization", "Y", "vector_d", context__.to_vec(N));
            Y = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("Y");
            pos__ = 0;
            size_t Y_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < Y_j_1_max__; ++j_1__) {
                Y(j_1__) = vals_r__[pos__++];
            }
            check_greater_or_equal(function__, "Y", Y, -(1));
            check_less_or_equal(function__, "Y", Y, 1);
            current_statement_begin__ = 17;
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
            current_statement_begin__ = 18;
            context__.validate_dims("data initialization", "offset", "double", context__.to_vec());
            offset = double(0);
            vals_r__ = context__.vals_r("offset");
            pos__ = 0;
            offset = vals_r__[pos__++];
            current_statement_begin__ = 19;
            context__.validate_dims("data initialization", "q", "double", context__.to_vec());
            q = double(0);
            vals_r__ = context__.vals_r("q");
            pos__ = 0;
            q = vals_r__[pos__++];
            check_greater_or_equal(function__, "q", q, 0);
            check_less_or_equal(function__, "q", q, 1);
            current_statement_begin__ = 20;
            context__.validate_dims("data initialization", "N_person", "int", context__.to_vec());
            N_person = int(0);
            vals_i__ = context__.vals_i("N_person");
            pos__ = 0;
            N_person = vals_i__[pos__++];
            current_statement_begin__ = 21;
            validate_non_negative_index("person", "N", N);
            context__.validate_dims("data initialization", "person", "int", context__.to_vec(N));
            person = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("person");
            pos__ = 0;
            size_t person_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < person_k_0_max__; ++k_0__) {
                person[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 22;
            context__.validate_dims("data initialization", "N_wave", "int", context__.to_vec());
            N_wave = int(0);
            vals_i__ = context__.vals_i("N_wave");
            pos__ = 0;
            N_wave = vals_i__[pos__++];
            current_statement_begin__ = 23;
            validate_non_negative_index("wave", "N", N);
            context__.validate_dims("data initialization", "wave", "int", context__.to_vec(N));
            wave = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("wave");
            pos__ = 0;
            size_t wave_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < wave_k_0_max__; ++k_0__) {
                wave[k_0__] = vals_i__[pos__++];
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 28;
            validate_non_negative_index("beta", "D", D);
            num_params_r__ += D;
            current_statement_begin__ = 29;
            validate_non_negative_index("beta_ind", "N_person", N_person);
            num_params_r__ += N_person;
            current_statement_begin__ = 30;
            validate_non_negative_index("beta_wave", "N_wave", N_wave);
            num_params_r__ += N_wave;
            current_statement_begin__ = 31;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_cbqpanelbv() { }
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
        current_statement_begin__ = 28;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "D", D);
        context__.validate_dims("parameter initialization", "beta", "vector_d", context__.to_vec(D));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta(D);
        size_t beta_j_1_max__ = D;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            beta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 29;
        if (!(context__.contains_r("beta_ind")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta_ind missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta_ind");
        pos__ = 0U;
        validate_non_negative_index("beta_ind", "N_person", N_person);
        context__.validate_dims("parameter initialization", "beta_ind", "vector_d", context__.to_vec(N_person));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_ind(N_person);
        size_t beta_ind_j_1_max__ = N_person;
        for (size_t j_1__ = 0; j_1__ < beta_ind_j_1_max__; ++j_1__) {
            beta_ind(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta_ind);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta_ind: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 30;
        if (!(context__.contains_r("beta_wave")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta_wave missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta_wave");
        pos__ = 0U;
        validate_non_negative_index("beta_wave", "N_wave", N_wave);
        context__.validate_dims("parameter initialization", "beta_wave", "vector_d", context__.to_vec(N_wave));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_wave(N_wave);
        size_t beta_wave_j_1_max__ = N_wave;
        for (size_t j_1__ = 0; j_1__ < beta_wave_j_1_max__; ++j_1__) {
            beta_wave(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta_wave);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta_wave: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 31;
        if (!(context__.contains_r("sigma_beta_ind")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma_beta_ind missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma_beta_ind");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "sigma_beta_ind", "double", context__.to_vec());
        double sigma_beta_ind(0);
        sigma_beta_ind = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, sigma_beta_ind);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma_beta_ind: ") + e.what()), current_statement_begin__, prog_reader__());
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
            current_statement_begin__ = 28;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.vector_constrain(D, lp__);
            else
                beta = in__.vector_constrain(D);
            current_statement_begin__ = 29;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_ind;
            (void) beta_ind;  // dummy to suppress unused var warning
            if (jacobian__)
                beta_ind = in__.vector_constrain(N_person, lp__);
            else
                beta_ind = in__.vector_constrain(N_person);
            current_statement_begin__ = 30;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_wave;
            (void) beta_wave;  // dummy to suppress unused var warning
            if (jacobian__)
                beta_wave = in__.vector_constrain(N_wave, lp__);
            else
                beta_wave = in__.vector_constrain(N_wave);
            current_statement_begin__ = 31;
            local_scalar_t__ sigma_beta_ind;
            (void) sigma_beta_ind;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma_beta_ind = in__.scalar_lb_constrain(0, lp__);
            else
                sigma_beta_ind = in__.scalar_lb_constrain(0);
            // model body
            {
            current_statement_begin__ = 36;
            local_scalar_t__ lik(DUMMY_VAR__);
            (void) lik;  // dummy to suppress unused var warning
            stan::math::initialize(lik, DUMMY_VAR__);
            stan::math::fill(lik, DUMMY_VAR__);
            current_statement_begin__ = 37;
            lp_accum__.add(normal_log<propto__>(beta, 0, 10));
            current_statement_begin__ = 39;
            lp_accum__.add(cauchy_log<propto__>(sigma_beta_ind, 0, 1));
            current_statement_begin__ = 40;
            lp_accum__.add(normal_log<propto__>(beta_ind, 0, sigma_beta_ind));
            current_statement_begin__ = 41;
            lp_accum__.add(normal_log<propto__>(beta_wave, 0, 10));
            current_statement_begin__ = 43;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 44;
                if (as_bool(logical_eq(get_base1(Y, i, "Y", 1), 1))) {
                    current_statement_begin__ = 45;
                    stan::math::assign(lik, ((1 - pald2(((dot_product(stan::model::rvalue(X, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "X"), beta) + get_base1(beta_ind, get_base1(person, i, "person", 1), "beta_ind", 1)) + get_base1(beta_wave, get_base1(wave, i, "wave", 1), "beta_wave", 1)), q, pstream__)) + offset));
                }
                current_statement_begin__ = 47;
                if (as_bool(logical_eq(get_base1(Y, i, "Y", 1), 0))) {
                    current_statement_begin__ = 48;
                    stan::math::assign(lik, (pald2(-(((dot_product(stan::model::rvalue(X, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "X"), beta) + get_base1(beta_ind, get_base1(person, i, "person", 1), "beta_ind", 1)) + get_base1(beta_wave, get_base1(wave, i, "wave", 1), "beta_wave", 1))), q, pstream__) + offset));
                }
                current_statement_begin__ = 50;
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
        names__.push_back("beta_ind");
        names__.push_back("beta_wave");
        names__.push_back("sigma_beta_ind");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(D);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N_person);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N_wave);
        dimss__.push_back(dims__);
        dims__.resize(0);
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
        static const char* function__ = "model_cbqpanelbv_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta = in__.vector_constrain(D);
        size_t beta_j_1_max__ = D;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            vars__.push_back(beta(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_ind = in__.vector_constrain(N_person);
        size_t beta_ind_j_1_max__ = N_person;
        for (size_t j_1__ = 0; j_1__ < beta_ind_j_1_max__; ++j_1__) {
            vars__.push_back(beta_ind(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_wave = in__.vector_constrain(N_wave);
        size_t beta_wave_j_1_max__ = N_wave;
        for (size_t j_1__ = 0; j_1__ < beta_wave_j_1_max__; ++j_1__) {
            vars__.push_back(beta_wave(j_1__));
        }
        double sigma_beta_ind = in__.scalar_lb_constrain(0);
        vars__.push_back(sigma_beta_ind);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
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
        return "model_cbqpanelbv";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = D;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_ind_j_1_max__ = N_person;
        for (size_t j_1__ = 0; j_1__ < beta_ind_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_ind" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_wave_j_1_max__ = N_wave;
        for (size_t j_1__ = 0; j_1__ < beta_wave_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_wave" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_beta_ind";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = D;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_ind_j_1_max__ = N_person;
        for (size_t j_1__ = 0; j_1__ < beta_ind_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_ind" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_wave_j_1_max__ = N_wave;
        for (size_t j_1__ = 0; j_1__ < beta_wave_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_wave" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_beta_ind";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_cbqpanelbv_namespace::model_cbqpanelbv stan_model;
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
