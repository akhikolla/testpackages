#' @title Two-drug combination example
#'
#' @name example-combo2
#'
#' @description
#' Example using a combination of two experimental drugs.
#'
#' @details
#'
#' The following example is described in the reference
#' Neuenschwander, B. et al (2016). The data are described
#' in the help page for \code{codata_combo2}. In the study
#' \code{trial_AB}, the risk of DLT was studied as a function of
#' dose for two drugs, drug A and drug B. Historical information
#' on the toxicity profiles of these two drugs was available from
#' single agent trials \code{trial_A} and \code{trial_B}. Another
#' study \code{IIT} was run concurrently to \code{trial_AB}, and
#' studies the same combination.
#'
#' The model described in Neuenschwander, et al (2016) is adapted as follows.
#' For groups \eqn{j = 1,\ldots, 4} representing each of the four sources
#' of data mentioned above,
#' \deqn{\mbox{logit}\, \pi_{1j}(d_1) = \log\, \alpha_{1j} + \beta_{1j} \, \log\, \Bigl(\frac{d_1}{d_1^*}\Bigr),}
#' and
#' \deqn{\mbox{logit}\, \pi_{2j}(d_2) = \log\, \alpha_{2j} + \beta_{2j} \, \log\, \Bigl(\frac{d_2}{d_2^*}\Bigr),}
#' are logistic regressions for the single-agent toxicity of drugs A and B,
#' respectively, when administered in group \eqn{j}. Conditional on the
#' regression parameters
#' \eqn{\boldsymbol\theta_{1j} = (\log \, \alpha_{1j}, \log \, \beta_{1j})} and
#' \eqn{\boldsymbol\theta_{2j} = (\log \, \alpha_{2j}, \log \, \beta_{2j})},
#' the toxicity \eqn{\pi_{j}(d_1, d_2)} for
#' the combination is modeled as the "no-interaction" DLT rate,
#' \deqn{\tilde\pi_{j}(d_1, d_2) = 1 - (1-\pi_{1j}(d_1) )(1- \pi_{2j}(d_2))}
#' with a single interaction term added on the log odds scale,
#' \deqn{\mbox{logit} \, \pi_{j}(d_1, d_2) = \mbox{logit} \, \tilde\pi_{j}(d_1, d_2) + \eta_j \frac{d_1}{d_1^*}\frac{d_2}{d_2^*}.}
#' A hierarchical model across the four groups \eqn{j} allows
#' dose-toxicity information to be shared through common hyperparameters.
#'
#' For the component parameters \eqn{\boldsymbol\theta_{ij}},
#' \deqn{\boldsymbol\theta_{ij} \sim \mbox{BVN}(\boldsymbol \mu_i, \boldsymbol\Sigma_i).}
#' For the mean, a further prior is specified as
#' \deqn{\boldsymbol\mu_i = (\mu_{\alpha i}, \mu_{\beta i}) \sim \mbox{BVN}(\boldsymbol m_i, \boldsymbol S_i),}
#' with \eqn{\boldsymbol m_i = (\mbox{logit}\, 0.1, \log 1)} and
#' \eqn{\boldsymbol S_i  = \mbox{diag}(3.33^2, 1^2)} for each \eqn{i = 1,2}.
#' For the standard deviations and correlation parameters in the covariance matrix,
#' \deqn{\boldsymbol\Sigma_i = \left( \begin{array}{cc}
#' \tau^2_{\alpha i} & \rho_i \tau_{\alpha i} \tau_{\beta i}\\
#' \rho_i \tau_{\alpha i} \tau_{\beta i} & \tau^2_{\beta i}
#' \end{array} \right), }
#' the specified priors are
#' \eqn{\tau_{\alpha i} \sim \mbox{Log-Normal}(\log\, 0.25, ((\log 4) / 1.96)^2)},
#'
#' \eqn{\tau_{\beta i} \sim \mbox{Log-Normal}(\log\, 0.125, ((\log 4) / 1.96)^2)},
#' and \eqn{\rho_i \sim \mbox{U}(-1,1)} for \eqn{i = 1,2}.
#'
#' For the interaction parameters \eqn{\eta_j} in each group, the hierarchical
#' model has
#' \deqn{\eta_j \sim \mbox{N}(\mu_\eta, \tau^2_\eta),}
#' for \eqn{j = 1,\ldots, 4}, with \eqn{\mu_\eta \sim \mbox{N}(0, 1.121^2)}
#' and \eqn{\tau_\eta \sim \mbox{Log-Normal}(\log\, 0.125, ((\log 4) / 1.96)^2).}
#'
#' Below is the syntax for specifying this fully exchangeable model in
#' \code{blrm_exnex}.
#'
#' @template ref-mac
#'
#' @template start-example
#' @template example-combo2
#' @template stop-example
NULL
