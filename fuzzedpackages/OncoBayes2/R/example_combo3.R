#' @title Three-drug combination example
#'
#' @name example-combo3
#'
#' @description
#' Example using a combination of two experimental drugs, with EXNEX
#' and differential discounting.
#'
#' @details
#'
#' This dataset involves a hypothetical dose-escalation study of combination
#' therapy with three treatment components. From two previous studies
#' \code{HistAgent1} and \code{HistAgent2}, historical data is available on each
#' of the treatments as single-agents, as well as two of the two-way
#' combinations. However, due to a difference in treatment schedule between the
#' \code{Combo} study and the historical studies, a stratification (through \code{stratum})
#' is made between the groups to allow differential discounting of the
#' alternate-schedule data. The association is as below.
#'
#' \tabular{cc}{
#'  group_id (j):  \tab stratum (s_j): \cr
#'  Combo (1)      \tab BID (1)         \cr
#'  HistAgent1 (2) \tab QD (2)          \cr
#'  HistAgent2 (3) \tab QD (2)
#' }
#'
#' For additional robustness, EXNEX priors are used for all group-level
#' treatment component and interaction parameters, to limit the amount of
#' borrowing in case of significant heterogeneity across groups.
#'
#' The complete model is as follows. As a function of doses \eqn{d_1,d_2,d_3}, the
#' DLT rate in group \eqn{j} is, for \eqn{j = 1,\ldots,3},
#' \deqn{\mbox{logit}\, \pi_j(d_1,d_2,d_3) = \mbox{logit}\Bigl( 1 - \prod_{i=1}^3 (1-\pi_{ij}(d_i))\Bigr) + \eta_{j}^{(12)}\frac{d_1}{d_1^*}\frac{d_2}{d_2^*} + \eta_{j}^{(13)}\frac{d_1}{d_1^*}\frac{d_3}{d_3^*} + \eta_{j}^{(23)}\frac{d_2}{d_2^*}\frac{d_3}{d_3^*} + \eta_{j}^{(123)}\frac{d_1}{d_1^*}\frac{d_2}{d_2^*}\frac{d_3}{d_3^*}.}
#'
#' In group \eqn{j} each treatment component \eqn{i} toxicity is modeled with
#' logistic regression,
#' \deqn{\mbox{logit}\, \pi_{ij}(d_i) = \log\, \alpha_{ij} + \beta_{ij} \, \log\, \Bigl(\frac{d_i}{d_i^*}\Bigr).}
#' The intercept and log-slope parameters \eqn{\boldsymbol\theta_{ij} = (\log\, \alpha_{ij}, \log\, \beta_{ij})}
#' are are given an EXNEX prior
#' \deqn{\boldsymbol\theta_{ij} \sim p_{ij} \mbox{BVN}(\boldsymbol\mu_i, \boldsymbol\Sigma_{ij}) + (1-p_{ij}) \mbox{BVN}(\boldsymbol m_{ij}, \boldsymbol S_{ij}),}
#' where the exchangeability weights are all \eqn{p_{ij} = 0.9}.
#' The NEX parameters are set to \eqn{\boldsymbol m_{ij} = (\mbox{logit}(1/3), \log\, 1)},
#' \eqn{\boldsymbol S_{ij} = \mbox{diag}(2^2, 1^2)} for all components \eqn{i=1,2,3} and
#' groups \eqn{j = 1,2,3}, and the EX parameters are modeled hierarchically. The
#' mean of the exchangeable part has the distribution
#' \deqn{\boldsymbol\mu_i = (\mu_{\alpha i}, \mu_{\beta i}) \sim \mbox{BVN}(\boldsymbol m_i, \boldsymbol S_i),}
#' with \eqn{\boldsymbol m_i = (\mbox{logit}(1/3), \log 1)} and
#' \eqn{\boldsymbol S_i  = \mbox{diag}(2^2, 1^2)} for each component \eqn{i = 1,2,3}.
#' For differentially discounting data from each schedule (QD and BID), the
#' covariance parameters for the exchangeable part
#' \deqn{\Sigma_{ij} = \left( \begin{array}{cc}
#' \tau^2_{\alpha s_j i} & \rho_i \tau_{\alpha s_j i} \tau_{\beta s_j i}\\
#' \rho_i \tau_{\alpha s_j i} \tau_{\beta s_j i} & \tau^2_{\beta s_j i}
#' \end{array} \right).}
#' are allowed to vary across groups \eqn{j} depending on their mapping
#' to strata \eqn{s(j)} as described above. For stratum \eqn{s=1} (\code{BID},
#' which contains only the group \eqn{j = 1} (\code{Combo})), the standard
#' deviations are modeled as
#' \deqn{\tau_{\alpha 1 i} \sim \mbox{Log-Normal}(\log\,0.25, (\log 4 / 1.96)^2)}
#' \deqn{\tau_{\beta 1 i} \sim \mbox{Log-Normal}(\log\,0.125, (\log 4 / 1.96)^2).}
#' Whereas in stratum \eqn{s=2} (\code{QD}, which contains the historical groups
#' \eqn{j=2,3} (\code{HistData1}, \code{HistData2})), the standard deviations are
#' \deqn{\tau_{\alpha 2 i} \sim \mbox{Log-Normal}(\log\,0.5, (\log 4 / 1.96)^2)}
#' \deqn{\tau_{\beta 2 i} \sim \mbox{Log-Normal}(\log\,0.25, (\log 4 / 1.96)^2).}
#'
#' For all interaction parameters \eqn{\eta_{j}^{(12)}}, \eqn{\eta_{j}^{(13)}},
#' \eqn{\eta_{j}^{(23)}}, and \eqn{\eta_{j}^{(123)}} (\eqn{j = 1,2,3}), the following
#' prior is assumed:
#' \deqn{\eta_{j}^{(\cdot)} \sim p_{\eta j}^{(\cdot)} \mbox{N}(\mu_{\eta}^{(\cdot)},{\tau_{\eta s_j}^{(\cdot)}}^2) + (1-p_{\eta j}^{(\cdot)}) \mbox{N}(m_{\eta j}^{(\cdot)}, {s_{\eta j}^{(\cdot)}}^2).}
#' The exchangeability weights are \eqn{p_{\eta j}^{(\cdot)} = 0.9} for all interaction
#' parameters and all groups. Here, for each \eqn{\mu_{\eta}^{(12)}}, \eqn{\mu_{\eta}^{(13)}},
#' \eqn{\mu_{\eta}^{(23)}}, and \eqn{\mu_{\eta}^{(123)}}, we take
#' \deqn{\mu_{\eta}^{(\cdot)} \sim \mbox{N}(0, 1/2),}
#' and for each \eqn{\tau_{\eta s}^{(12)}}, \eqn{\tau_{\eta s}^{(13)}},
#' \eqn{\tau_{\eta s}^{(23)}}, and \eqn{\tau_{\eta s}^{(123)}},
#' \deqn{\tau_{\eta s}^{(\cdot)} \sim \mbox{Log-Normal}(\log(0.25), (\log 2 / 1.96)^2),}
#' for both strata \eqn{s = 1,2}. Furthermore, \eqn{m_{\eta j}^{(\cdot)} = 0} and
#' \eqn{{s_{\eta j}^{(\cdot)}}^2 = 1/2}, uniformly across all indices.
#'
#' Below is the syntax for specifying this model in \code{blrm_exnex}.
#'
#' @template ref-mac
#'
#' @template start-example
#' @template example-combo3
#' @template stop-example
NULL
