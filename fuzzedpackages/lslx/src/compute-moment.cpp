// utility functions for package lslx
// written by Po-Hsien Huang psyphh@gmail.com

#include "utility-function.h"

// [[Rcpp::depends(RcppEigen)]]
// compute saturated moment via EM algorithm
// [[Rcpp::export]]
void compute_saturated_moment_cpp(
    Rcpp::List y_obs,
    Rcpp::List w,
    Rcpp::List m_idx,
    Rcpp::List saturated_mean,
    Rcpp::List saturated_cov,
    int iter_other_max,
    double tol_other) {
  Rcpp::List y_obs_i;
  Rcpp::List w_i;
  Rcpp::List m_idx_i;
  Eigen::MatrixXd saturated_mean_i, saturated_cov_i;
  Rcpp::IntegerVector m_idx_ik;
  Eigen::MatrixXd saturated_cov_ik, saturated_cov_ik_inv;
  Eigen::MatrixXd e_sum_i, c_sum_i;
  Eigen::MatrixXd a_ik, b_ik;
  Eigen::MatrixXd y_ik, y_ik_w;
  double moment_diff_max;
  
  int i, j, k;
  int n_group = y_obs.size();
  for (i = 0; i < n_group; i++) {
    y_obs_i = Rcpp::as<List>(y_obs[i]);
    m_idx_i = Rcpp::as<List>(m_idx[i]);
    w_i = Rcpp::as<List>(w[i]);
    for (j = 0; j < iter_other_max; j++) {
      saturated_cov_i = Rcpp::as< Eigen::MatrixXd >(saturated_cov[i]);
      saturated_mean_i = Rcpp::as< Eigen::MatrixXd >(saturated_mean[i]);
      e_sum_i = Eigen::MatrixXd::Zero(saturated_mean_i.rows(), 1);
      c_sum_i = Eigen::MatrixXd::Zero(saturated_cov_i.rows(), saturated_cov_i.cols());
      for (k = 0; k < y_obs_i.size(); k++) {
        Eigen::Map<Eigen::MatrixXd> y_obs_ik(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(y_obs_i[k]));
        Eigen::Map<Eigen::VectorXd> w_ik(Rcpp::as< Eigen::Map <Eigen::VectorXd> >(w_i[k]));
        m_idx_ik = Rcpp::as< Rcpp::IntegerVector >(m_idx_i[k]);
        saturated_cov_ik = slice_both(saturated_cov_i, m_idx_ik, m_idx_ik);
        saturated_cov_ik_inv = saturated_cov_ik.inverse();
        
        a_ik = saturated_mean_i.transpose() - 
          slice_row(saturated_mean_i, m_idx_ik).transpose() * 
          saturated_cov_ik_inv * slice_row(saturated_cov_i, m_idx_ik);
        b_ik = saturated_cov_ik_inv * slice_row(saturated_cov_i, m_idx_ik);
        y_ik = y_obs_ik * b_ik + Eigen::MatrixXd::Ones(y_obs_ik.rows(), 1)  * a_ik;
        y_ik_w = (y_ik.array().colwise() * w_ik.array()).matrix();
        e_sum_i += y_ik_w.colwise().sum().transpose();
        c_sum_i += w_ik.sum() * saturated_cov_i -
          w_ik.sum() * slice_col(saturated_cov_i, m_idx_ik) * b_ik +
          y_ik_w.transpose() * y_ik;
      }
      saturated_mean[i] = e_sum_i;
      saturated_cov[i] = c_sum_i - e_sum_i * e_sum_i.transpose();
      moment_diff_max = 
        std::max((Rcpp::as< Eigen::MatrixXd >(saturated_mean[i]) - saturated_mean_i).array().abs().maxCoeff(),
                 (Rcpp::as< Eigen::MatrixXd >(saturated_cov[i]) - saturated_cov_i).array().abs().maxCoeff());
      if (moment_diff_max < tol_other) {
        break;
      }
    }
  }
}

// compute asymptotic covariance of saturated moment via response (raw data)
// [[Rcpp::export]]
void compute_saturated_moment_acov_response_cpp(
    Rcpp::List y_obs,
    Rcpp::List w,
    Rcpp::List m_idx,
    Rcpp::List m2_idx,
    Rcpp::List saturated_mean,
    Rcpp::List saturated_cov,
    Rcpp::List saturated_moment_acov) {
  int n_response_i, n_moment_i, sample_size_i;
  Rcpp::List y_obs_i;
  Rcpp::List w_i;
  Rcpp::List m_idx_i;
  Rcpp::List m2_idx_i;
  Eigen::MatrixXd saturated_mean_i, saturated_cov_i;
  Eigen::MatrixXd duplication_i;
  Eigen::MatrixXd score2_sum_i;
  Eigen::MatrixXd hessian_sum_i, hessian_sum_i_inv;
  Eigen::MatrixXd saturated_mean_ij, saturated_cov_ij, saturated_cov_ij_vech, saturated_cov_ij_inv;
  Eigen::MatrixXd y_obs_ij;
  Eigen::VectorXd w_ij;
  Eigen::MatrixXd yc_obs_ij, yc_obs_ij_w;
  Eigen::MatrixXd yc2c_obs_ij, yc2c_obs_ijk;
  Eigen::MatrixXd score_ij, score_ij_w;
  Eigen::MatrixXd saturated_moment_acov_i;
  Rcpp::IntegerVector m_idx_ij, m2_idx_ij;
  Rcpp::IntegerVector idx_vech_i, idx_tvech_i, idx_vech_match_i;
  Rcpp::IntegerVector idx_vech_ij, idx_tvech_ij, idx_vech_match_ij;
  int n_response_ij, n_moment_ij, sample_size_ij;
  Eigen::MatrixXd duplication_ij;
  
  int i, j, k;
  int n_group = y_obs.size();
  for (i = 0; i < n_group; i++) {
    sample_size_i = 0;
    y_obs_i = Rcpp::as<List>(y_obs[i]);
    m_idx_i = Rcpp::as<List>(m_idx[i]);
    m2_idx_i = Rcpp::as<List>(m2_idx[i]);
    w_i = Rcpp::as<List>(w[i]);
    saturated_cov_i = Rcpp::as< Eigen::MatrixXd >(saturated_cov[i]);
    saturated_mean_i = Rcpp::as< Eigen::MatrixXd >(saturated_mean[i]);
    n_response_i = saturated_mean_i.rows();
    n_moment_i = (n_response_i * (n_response_i + 3)) / 2;
    score2_sum_i = Eigen::MatrixXd::Zero(n_moment_i, n_moment_i);
    hessian_sum_i = Eigen::MatrixXd::Zero(n_moment_i, n_moment_i);
    
    duplication_i = create_duplication(n_response_i);
    idx_vech_i = create_idx_vech(n_response_i, true);
    idx_tvech_i = create_idx_tvech(n_response_i, true);
    idx_vech_match_i = find_idx_match(idx_vech_i, idx_tvech_i);
    
    if ((y_obs_i.size() == 1) & ((Rcpp::as< Rcpp::IntegerVector >(m_idx_i[0])).size() == n_response_i)) {
      y_obs_ij = Rcpp::as<Eigen::MatrixXd>(y_obs_i[0]);
      w_ij = Rcpp::as<Eigen::VectorXd>(w_i[0]);
      sample_size_ij = y_obs_ij.rows();
      sample_size_i += sample_size_ij;
      saturated_mean_ij = saturated_mean_i;
      saturated_cov_ij = saturated_cov_i;
      saturated_cov_ij_vech = vech(saturated_cov_ij); 
      n_response_ij = saturated_mean_ij.rows();
      n_moment_ij = (n_response_ij * (n_response_ij + 3)) / 2;
      
      yc_obs_ij = y_obs_ij - 
        Eigen::MatrixXd::Ones(sample_size_ij, 1) * saturated_mean_ij.transpose();
      yc2c_obs_ij.resize(sample_size_ij, n_moment_ij - n_response_ij);
      
      for (k = 0; k < sample_size_ij; k++) {
        yc2c_obs_ijk = vech(yc_obs_ij.row(k).transpose() * yc_obs_ij.row(k));
        yc2c_obs_ij.row(k) = (yc2c_obs_ijk - saturated_cov_ij_vech).transpose();
      }
      saturated_cov_ij_inv = saturated_cov_ij.inverse();
      duplication_ij = create_duplication(n_response_ij);
      idx_vech_ij = create_idx_vech(n_response_ij, true);
      idx_tvech_ij = create_idx_tvech(n_response_ij, true);
      idx_vech_match_ij = find_idx_match(idx_vech_ij, idx_tvech_ij);
      
      score_ij.resize(sample_size_ij, n_moment_i);
      score_ij.leftCols(n_response_i) = yc_obs_ij;
      score_ij.rightCols((n_moment_i - n_response_i)) = yc2c_obs_ij;
      
      score_ij_w = (score_ij.array().colwise() * w_ij.array()).matrix();
      score2_sum_i += score_ij_w.transpose() * score_ij;
      saturated_moment_acov_i = score2_sum_i / double(sample_size_i);
    } else {
      for (j = 0; j < y_obs_i.size(); j++) {
        y_obs_ij = Rcpp::as<Eigen::MatrixXd>(y_obs_i[j]);
        w_ij = Rcpp::as<Eigen::VectorXd>(w_i[j]);
        m_idx_ij = Rcpp::as< Rcpp::IntegerVector >(m_idx_i[j]);
        m2_idx_ij = Rcpp::as< Rcpp::IntegerVector >(m2_idx_i[j]);
        sample_size_ij = y_obs_ij.rows();
        sample_size_i += sample_size_ij;
        saturated_mean_ij = slice_row(saturated_mean_i, m_idx_ij);
        saturated_cov_ij = slice_both(saturated_cov_i, m_idx_ij, m_idx_ij);
        saturated_cov_ij_vech = vech(saturated_cov_ij); 
        n_response_ij = saturated_mean_ij.rows();
        n_moment_ij = (n_response_ij * (n_response_ij + 3)) / 2;
        
        yc_obs_ij = y_obs_ij - 
          Eigen::MatrixXd::Ones(sample_size_ij, 1) * saturated_mean_ij.transpose();
        yc2c_obs_ij.resize(sample_size_ij, n_moment_ij - n_response_ij);
        
        for (k = 0; k < sample_size_ij; k++) {
          yc2c_obs_ijk = vech(yc_obs_ij.row(k).transpose() * yc_obs_ij.row(k));
          yc2c_obs_ij.row(k) = (yc2c_obs_ijk - saturated_cov_ij_vech).transpose();
        }
        saturated_cov_ij_inv = saturated_cov_ij.inverse();
        
        duplication_ij = create_duplication(n_response_ij);
        idx_vech_ij = create_idx_vech(n_response_ij, true);
        idx_tvech_ij = create_idx_tvech(n_response_ij, true);
        idx_vech_match_ij = find_idx_match(idx_vech_ij, idx_tvech_ij);
        
        score_ij.resize(sample_size_ij, n_moment_i);
        score_ij.leftCols(n_response_i) = 
          expand_col((yc_obs_ij * saturated_cov_ij_inv), m_idx_ij, n_response_i);
        score_ij.rightCols((n_moment_i - n_response_i)) = 
          expand_col((0.5 * yc2c_obs_ij * deduplify_both(
              Eigen::kroneckerProduct(saturated_cov_ij_inv, saturated_cov_ij_inv), 
              idx_vech_ij, idx_tvech_ij, idx_vech_match_ij)), 
          m2_idx_ij, (n_moment_i - n_response_i));
        
        score_ij_w = (score_ij.array().colwise() * w_ij.array()).matrix();
        score2_sum_i += score_ij_w.transpose() * score_ij;
        
        yc_obs_ij = expand_col(yc_obs_ij, m_idx_ij, n_response_i);
        yc_obs_ij_w = (yc_obs_ij.array().colwise() * w_ij.array()).matrix();
        yc2c_obs_ij = expand_col(yc2c_obs_ij, m2_idx_ij, (n_moment_i - n_response_i));
        saturated_cov_ij_inv = expand_both(saturated_cov_ij_inv, m_idx_ij, m_idx_ij,
                                           n_response_i, n_response_i);
        
        hessian_sum_i.block(0, 0, n_response_i, n_response_i) +=
          w_ij.sum() * saturated_cov_ij_inv;
        hessian_sum_i.block(n_response_i, n_response_i, 
                            (n_moment_i - n_response_i), 
                            (n_moment_i - n_response_i)) += 
                              deduplify_both(
                                (Eigen::kroneckerProduct(saturated_cov_ij_inv,
                                                         (saturated_cov_ij_inv * 
                                                           (yc_obs_ij_w.transpose() * yc_obs_ij) * 
                                                           saturated_cov_ij_inv - 0.5 * w_ij.sum() * 
                                                           saturated_cov_ij_inv))), 
                                                           idx_vech_i, idx_tvech_i, idx_vech_match_i);
        hessian_sum_i.block(0, n_response_i, 
                            n_response_i, 
                            (n_moment_i - n_response_i)) += 
                              deduplify_right((Eigen::kroneckerProduct(saturated_cov_ij_inv, 
                                                                       yc_obs_ij_w.colwise().sum() * saturated_cov_ij_inv)),
                                                                       idx_vech_i, idx_tvech_i, idx_vech_match_i);
        hessian_sum_i.block(n_response_i, 0, 
                            (n_moment_i - n_response_i), 
                            n_response_i) = 
                              hessian_sum_i.block(0, n_response_i, 
                                                  n_response_i, 
                                                  (n_moment_i - n_response_i)).transpose(); 
      }
      hessian_sum_i_inv = hessian_sum_i.inverse();
      saturated_moment_acov_i = (hessian_sum_i_inv * score2_sum_i * hessian_sum_i_inv) / double(sample_size_i);
    }
    saturated_moment_acov[i] = saturated_moment_acov_i;
  }
}

// compute asymptotic covariance of saturated moment via moment (moment data)
// [[Rcpp::export]]
void compute_saturated_moment_acov_moment_cpp(
    int n_observation,
    Rcpp::List sample_proportion,
    Rcpp::List saturated_cov,
    Rcpp::List saturated_moment_acov) {
  Eigen::MatrixXd saturated_cov_i_inv, saturated_moment_acov_i;
  Eigen::MatrixXd duplication_i;
  Rcpp::IntegerVector idx_vech_i, idx_tvech_i, idx_vech_match_i;
  double sample_proportion_i;
  int n_response_i, n_moment_i;
  int i;
  for (i = 0; i < saturated_cov.size(); i++) {
    Eigen::Map<MatrixXd> saturated_cov_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(saturated_cov[i]));
    sample_proportion_i = Rcpp::as<double>(sample_proportion[i]);
    saturated_cov_i_inv = saturated_cov_i.inverse();
    n_response_i = saturated_cov_i.cols();
    n_moment_i = n_response_i * (n_response_i + 3) / 2;
    duplication_i = create_duplication(n_response_i);
    idx_vech_i = create_idx_vech(n_response_i, true);
    idx_tvech_i = create_idx_tvech(n_response_i, true);
    idx_vech_match_i = find_idx_match(idx_vech_i, idx_tvech_i);
    saturated_moment_acov_i = Eigen::MatrixXd::Zero(n_moment_i, n_moment_i);
    saturated_moment_acov_i.block(0, 0, n_response_i, n_response_i) = saturated_cov_i_inv;
    saturated_moment_acov_i.block(n_response_i, n_response_i,
                                  n_moment_i - n_response_i, n_moment_i - n_response_i) =
                                    0.5 * deduplify_both(Eigen::kroneckerProduct(saturated_cov_i_inv, saturated_cov_i_inv), 
                                                         idx_vech_i, idx_tvech_i, idx_vech_match_i);
    saturated_moment_acov_i = saturated_moment_acov_i.inverse() / 
      (sample_proportion_i *double(n_observation));
    saturated_moment_acov[i] = saturated_moment_acov_i;
  }
}
