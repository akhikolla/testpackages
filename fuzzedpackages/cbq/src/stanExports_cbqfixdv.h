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
namespace model_cbqfixdv_namespace {
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
    reader.add_event(0, 0, "start", "model_cbqfixdv");
    reader.add_event(137, 135, "end", "model_cbqfixdv");
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
int
group_size(const std::vector<int>& ref,
               const int& value, std::ostream* pstream__) {
    typedef double local_scalar_t__;
    typedef int fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 15;
        int count(0);
        (void) count;  // dummy to suppress unused var warning
        stan::math::fill(count, std::numeric_limits<int>::min());
        current_statement_begin__ = 16;
        stan::math::assign(count, 0);
        current_statement_begin__ = 17;
        for (int ii = 1; ii <= size(ref); ++ii) {
            current_statement_begin__ = 18;
            if (as_bool(logical_eq(get_base1(ref, ii, "ref", 1), value))) {
                current_statement_begin__ = 19;
                stan::math::assign(count, (count + 1));
            }
        }
        current_statement_begin__ = 20;
        return stan::math::promote_scalar<fun_return_scalar_t__>(count);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct group_size_functor__ {
            int
    operator()(const std::vector<int>& ref,
               const int& value, std::ostream* pstream__) const {
        return group_size(ref, value, pstream__);
    }
};
std::vector<int>
subset_intarray(const std::vector<int>& y,
                    const std::vector<int>& ref,
                    const int& value, std::ostream* pstream__) {
    typedef double local_scalar_t__;
    typedef int fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 25;
        int jj(0);
        (void) jj;  // dummy to suppress unused var warning
        stan::math::fill(jj, std::numeric_limits<int>::min());
        current_statement_begin__ = 26;
        validate_non_negative_index("res", "group_size(ref, value, pstream__)", group_size(ref, value, pstream__));
        std::vector<int  > res(group_size(ref, value, pstream__), int(0));
        stan::math::fill(res, std::numeric_limits<int>::min());
        current_statement_begin__ = 27;
        if (as_bool(logical_neq(size(ref), size(y)))) {
            current_statement_begin__ = 28;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "illegal input: non-matching dimensions";
            throw std::domain_error(errmsg_stream__.str());
        }
        current_statement_begin__ = 29;
        stan::math::assign(jj, 1);
        current_statement_begin__ = 30;
        for (int ii = 1; ii <= size(ref); ++ii) {
            current_statement_begin__ = 31;
            if (as_bool(logical_eq(get_base1(ref, ii, "ref", 1), value))) {
                current_statement_begin__ = 32;
                stan::model::assign(res, 
                            stan::model::cons_list(stan::model::index_uni(jj), stan::model::nil_index_list()), 
                            get_base1(y, ii, "y", 1), 
                            "assigning variable res");
                current_statement_begin__ = 33;
                stan::math::assign(jj, (jj + 1));
            }
        }
        current_statement_begin__ = 36;
        return stan::math::promote_scalar<fun_return_scalar_t__>(res);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct subset_intarray_functor__ {
            std::vector<int>
    operator()(const std::vector<int>& y,
                    const std::vector<int>& ref,
                    const int& value, std::ostream* pstream__) const {
        return subset_intarray(y, ref, value, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_cbqfixdv
  : public stan::model::model_base_crtp<model_cbqfixdv> {
private:
        int N;
        int D_common;
        vector_d Y;
        matrix_d X_common;
        int N_indx;
        std::vector<int> ind;
        int N_wave;
        std::vector<int> wave;
        double q;
        double offset;
        std::vector<int> n_group;
public:
    model_cbqfixdv(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_cbqfixdv(stan::io::var_context& context__,
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
        static const char* function__ = "model_cbqfixdv_namespace::model_cbqfixdv";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 46;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 47;
            context__.validate_dims("data initialization", "D_common", "int", context__.to_vec());
            D_common = int(0);
            vals_i__ = context__.vals_i("D_common");
            pos__ = 0;
            D_common = vals_i__[pos__++];
            current_statement_begin__ = 48;
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
            current_statement_begin__ = 49;
            validate_non_negative_index("X_common", "N", N);
            validate_non_negative_index("X_common", "D_common", D_common);
            context__.validate_dims("data initialization", "X_common", "matrix_d", context__.to_vec(N,D_common));
            X_common = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, D_common);
            vals_r__ = context__.vals_r("X_common");
            pos__ = 0;
            size_t X_common_j_2_max__ = D_common;
            size_t X_common_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < X_common_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_common_j_1_max__; ++j_1__) {
                    X_common(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 50;
            context__.validate_dims("data initialization", "N_indx", "int", context__.to_vec());
            N_indx = int(0);
            vals_i__ = context__.vals_i("N_indx");
            pos__ = 0;
            N_indx = vals_i__[pos__++];
            current_statement_begin__ = 51;
            validate_non_negative_index("ind", "N", N);
            context__.validate_dims("data initialization", "ind", "int", context__.to_vec(N));
            ind = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("ind");
            pos__ = 0;
            size_t ind_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < ind_k_0_max__; ++k_0__) {
                ind[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 54;
            context__.validate_dims("data initialization", "N_wave", "int", context__.to_vec());
            N_wave = int(0);
            vals_i__ = context__.vals_i("N_wave");
            pos__ = 0;
            N_wave = vals_i__[pos__++];
            current_statement_begin__ = 55;
            validate_non_negative_index("wave", "N", N);
            context__.validate_dims("data initialization", "wave", "int", context__.to_vec(N));
            wave = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("wave");
            pos__ = 0;
            size_t wave_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < wave_k_0_max__; ++k_0__) {
                wave[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 56;
            context__.validate_dims("data initialization", "q", "double", context__.to_vec());
            q = double(0);
            vals_r__ = context__.vals_r("q");
            pos__ = 0;
            q = vals_r__[pos__++];
            current_statement_begin__ = 57;
            context__.validate_dims("data initialization", "offset", "double", context__.to_vec());
            offset = double(0);
            vals_r__ = context__.vals_r("offset");
            pos__ = 0;
            offset = vals_r__[pos__++];
            // initialize transformed data variables
            current_statement_begin__ = 61;
            validate_non_negative_index("n_group", "N_indx", N_indx);
            n_group = std::vector<int>(N_indx, int(0));
            stan::math::fill(n_group, std::numeric_limits<int>::min());
            // execute transformed data statements
            current_statement_begin__ = 62;
            for (int ii = 1; ii <= N_indx; ++ii) {
                current_statement_begin__ = 63;
                stan::model::assign(n_group, 
                            stan::model::cons_list(stan::model::index_uni(ii), stan::model::nil_index_list()), 
                            group_size(ind, ii, pstream__), 
                            "assigning variable n_group");
            }
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 68;
            validate_non_negative_index("beta", "D_common", D_common);
            num_params_r__ += D_common;
            current_statement_begin__ = 70;
            validate_non_negative_index("beta_wave", "N_wave", N_wave);
            num_params_r__ += N_wave;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_cbqfixdv() { }
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
        current_statement_begin__ = 68;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "D_common", D_common);
        context__.validate_dims("parameter initialization", "beta", "vector_d", context__.to_vec(D_common));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta(D_common);
        size_t beta_j_1_max__ = D_common;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            beta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 70;
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
            current_statement_begin__ = 68;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.vector_constrain(D_common, lp__);
            else
                beta = in__.vector_constrain(D_common);
            current_statement_begin__ = 70;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_wave;
            (void) beta_wave;  // dummy to suppress unused var warning
            if (jacobian__)
                beta_wave = in__.vector_constrain(N_wave, lp__);
            else
                beta_wave = in__.vector_constrain(N_wave);
            // model body
            {
            current_statement_begin__ = 80;
            int pos(0);
            (void) pos;  // dummy to suppress unused var warning
            stan::math::fill(pos, std::numeric_limits<int>::min());
            current_statement_begin__ = 81;
            validate_non_negative_index("xb_common", "N", N);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> xb_common(N);
            stan::math::initialize(xb_common, DUMMY_VAR__);
            stan::math::fill(xb_common, DUMMY_VAR__);
            current_statement_begin__ = 85;
            lp_accum__.add(normal_log<propto__>(beta, 0, 10));
            current_statement_begin__ = 87;
            lp_accum__.add(normal_log<propto__>(beta_wave, 0, 10));
            current_statement_begin__ = 91;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 92;
                stan::model::assign(xb_common, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            (multiply(stan::model::rvalue(X_common, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "X_common"), beta) + get_base1(beta_wave, get_base1(wave, i, "wave", 1), "beta_wave", 1)), 
                            "assigning variable xb_common");
            }
            current_statement_begin__ = 95;
            stan::math::assign(pos, 1);
            current_statement_begin__ = 97;
            for (int i = 1; i <= N_indx; ++i) {
                {
                current_statement_begin__ = 98;
                local_scalar_t__ lik0(DUMMY_VAR__);
                (void) lik0;  // dummy to suppress unused var warning
                stan::math::initialize(lik0, DUMMY_VAR__);
                stan::math::fill(lik0, DUMMY_VAR__);
                current_statement_begin__ = 99;
                local_scalar_t__ lik(DUMMY_VAR__);
                (void) lik;  // dummy to suppress unused var warning
                stan::math::initialize(lik, DUMMY_VAR__);
                stan::math::fill(lik, DUMMY_VAR__);
                current_statement_begin__ = 100;
                local_scalar_t__ lik2(DUMMY_VAR__);
                (void) lik2;  // dummy to suppress unused var warning
                stan::math::initialize(lik2, DUMMY_VAR__);
                stan::math::fill(lik2, DUMMY_VAR__);
                current_statement_begin__ = 101;
                local_scalar_t__ lik3(DUMMY_VAR__);
                (void) lik3;  // dummy to suppress unused var warning
                stan::math::initialize(lik3, DUMMY_VAR__);
                stan::math::fill(lik3, DUMMY_VAR__);
                current_statement_begin__ = 102;
                local_scalar_t__ lik4(DUMMY_VAR__);
                (void) lik4;  // dummy to suppress unused var warning
                stan::math::initialize(lik4, DUMMY_VAR__);
                stan::math::fill(lik4, DUMMY_VAR__);
                current_statement_begin__ = 103;
                validate_non_negative_index("y_g", "get_base1(n_group, i, \"n_group\", 1)", get_base1(n_group, i, "n_group", 1));
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> y_g(get_base1(n_group, i, "n_group", 1));
                stan::math::initialize(y_g, DUMMY_VAR__);
                stan::math::fill(y_g, DUMMY_VAR__);
                current_statement_begin__ = 104;
                validate_non_negative_index("xb_common_g", "get_base1(n_group, i, \"n_group\", 1)", get_base1(n_group, i, "n_group", 1));
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> xb_common_g(get_base1(n_group, i, "n_group", 1));
                stan::math::initialize(xb_common_g, DUMMY_VAR__);
                stan::math::fill(xb_common_g, DUMMY_VAR__);
                current_statement_begin__ = 105;
                validate_non_negative_index("ystar_g", "get_base1(n_group, i, \"n_group\", 1)", get_base1(n_group, i, "n_group", 1));
                Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> ystar_g(get_base1(n_group, i, "n_group", 1));
                stan::math::initialize(ystar_g, DUMMY_VAR__);
                stan::math::fill(ystar_g, DUMMY_VAR__);
                current_statement_begin__ = 107;
                stan::math::assign(y_g, segment(Y, pos, get_base1(n_group, i, "n_group", 1)));
                current_statement_begin__ = 108;
                stan::math::assign(xb_common_g, segment(xb_common, pos, get_base1(n_group, i, "n_group", 1)));
                current_statement_begin__ = 110;
                stan::math::assign(lik, 1);
                current_statement_begin__ = 111;
                for (int j = 1; j <= (get_base1(n_group, i, "n_group", 1) - 1); ++j) {
                    current_statement_begin__ = 113;
                    stan::math::assign(lik, (lik * (1 - pald2(-((get_base1(xb_common_g, get_base1(n_group, i, "n_group", 1), "xb_common_g", 1) - get_base1(xb_common_g, j, "xb_common_g", 1))), q, pstream__))));
                }
                current_statement_begin__ = 116;
                stan::math::assign(lik3, 1);
                current_statement_begin__ = 117;
                for (int j = 1; j <= (get_base1(n_group, i, "n_group", 1) - 1); ++j) {
                    current_statement_begin__ = 118;
                    stan::math::assign(lik0, 1);
                    current_statement_begin__ = 119;
                    for (int k = 1; k <= get_base1(n_group, i, "n_group", 1); ++k) {
                        current_statement_begin__ = 120;
                        if (as_bool(logical_neq(j, k))) {
                            current_statement_begin__ = 121;
                            stan::math::assign(lik0, (lik0 * (1 - pald2(-((get_base1(xb_common_g, j, "xb_common_g", 1) - get_base1(xb_common_g, k, "xb_common_g", 1))), q, pstream__))));
                        }
                    }
                    current_statement_begin__ = 124;
                    stan::math::assign(lik3, (lik3 * (1 - lik0)));
                }
                current_statement_begin__ = 127;
                stan::math::assign(lik2, ((lik + offset) * lik3));
                current_statement_begin__ = 129;
                lp_accum__.add(stan::math::log(lik2));
                current_statement_begin__ = 130;
                stan::math::assign(pos, (pos + get_base1(n_group, i, "n_group", 1)));
                }
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
        names__.push_back("beta_wave");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(D_common);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N_wave);
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
        static const char* function__ = "model_cbqfixdv_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta = in__.vector_constrain(D_common);
        size_t beta_j_1_max__ = D_common;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            vars__.push_back(beta(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta_wave = in__.vector_constrain(N_wave);
        size_t beta_wave_j_1_max__ = N_wave;
        for (size_t j_1__ = 0; j_1__ < beta_wave_j_1_max__; ++j_1__) {
            vars__.push_back(beta_wave(j_1__));
        }
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
        return "model_cbqfixdv";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = D_common;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_wave_j_1_max__ = N_wave;
        for (size_t j_1__ = 0; j_1__ < beta_wave_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_wave" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = D_common;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_wave_j_1_max__ = N_wave;
        for (size_t j_1__ = 0; j_1__ < beta_wave_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta_wave" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_cbqfixdv_namespace::model_cbqfixdv stan_model;
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
