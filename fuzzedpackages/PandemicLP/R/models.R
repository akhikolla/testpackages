#' Models used in the PandemicLP package
#'
#' This document explains the models used in the PandemicLP package in detail.
#'
#' @name models
#' @section Growth curve for the mean cases:
#' The count data for number of cases or deaths is modeled according to an epidemiological model of growth.
#' In particular, the average counts are \eqn{mu(t)} modeled with a generalized logistic curve:
#' \deqn{\mu(t) = a c f \frac{e^{-c t}}{(b+e^{-c t})^{f+1}}.}
#' All parameters, that is \eqn{a,b,c} and \eqn{f} are positive.
#'
#' Parameter \eqn{c} is interpreted as the infection rate. Parameter \eqn{f} controls the asymmetry, so
#' if it is equal to 1, then the curve is symmetric. If it is lesser than 1, then the cases grow slower before
#' the peak than they decrease after. The behavior is inverted when \eqn{f} is greater than 1.
#'
#' The counts for the Covid-19 pandemic typically had a behavior with positive asymmetry, and so the default for the
#' package functions is to use a greater than 1 truncation for \eqn{f}.
#'
#' It was common in the early stages of the Covid-19 pandemic that the predictions would result in very high
#' and absurd values for the total number of cases (TNC). It is straightforward to show that
#' \deqn{TNC = \frac{a}{b^f}.}
#' Since all locations displayed a total number of cases that never exceeded 5\% of that location's population,
#' another truncation is applied, so that \eqn{a\le b^f 0.08 Pop}, where \eqn{Pop} is the location's population.
#' @section Probabilistic model:
#' The simplest probabilistic model for the counts is the Poisson model.
#' If \eqn{y_t} is the count at time \eqn{t}, then
#' \deqn{y_t | \theta \sim Poisson(\mu(t)),} where \eqn{\theta} represents the model parameters.
#' @section Advanced modeling:
#' Here we present some other forms for the growth curve in the mean.
#' \subsection{Seasonal effects}{
#' A weekly seasonal effect can be added. This is done by multiplying \eqn{\mu(t)} by a positive effect \eqn{d} when \eqn{t}
#' is the desired weekday. If \eqn{d<1} then that weekday represents underreporting. It is overreporting if \eqn{d>1}.}
#' \subsection{Multiple curves}{
#' Additionally, two curves can be fitted, as happened in the Covid-19 pandemic in many locations.
#' In this case the model is slightly different. In this case,
#' \deqn{\mu(t) = \mu_1(t)+\mu_2(t)}\deqn{\mu_j(t) = a_j c_j \frac{e^{-c_j t}}{(b_j+e^{-c_j t})^2}\Phi(\alpha_j (t-\delta_j)), j = 1, 2,}
#' where \eqn{\Phi(.)} is the probit function. The probit function induces asymmetry in the curve, similarly to
#' parameter \eqn{f}.}
#' @section Prior distribution:
#' Apart from the truncation mentioned above, the prior is defined as below. Note that every available model used in
#' the \code{\link{pandemic_model}} function uses only a subset of these parameters. The priors are described here
#' for reference.
#' \deqn{a\sim Gamma(0.1,0.1)} \deqn{b\sim LogNormal(0,20)} \deqn{c\sim Gamma(2,9)} \deqn{f\sim Gamma(0.01, 0.01).}
#' \deqn{d_k\sim Gamma(2,1),k=1,...,3} \deqn{\delta_j\sim Normal(0,100),j=1,2} \deqn{\alpha_j\sim Gamma(0.01,0.01),j=1,2}
#' @section Options for the \code{\link{pandemic_model}} function:
#' Two arguments in the function change the fitted model, as described below:
#' \itemize{
#'   \item 'seasonal': By leaving this argument \code{NULL}, the standard model is fitted. By supplying it with a vector
#'   of up to three weekdays, the desired seasonal effects are added to the model.
#'   \item 'n_waves': By leaving this argument equal to 1, the standard model is fitted. By changing it to
#'   2 implies a two waves model. Future versions of the package will allow for more waves.}
#' Seasonal and multiple waves cannot currently be used in the same model at the same time. Future versions of the
#' package will allow such mixtures and also other distributions for the data.
#' @seealso \code{\link{pandemic_model}} and \code{\link{posterior_predict.pandemicEstimated}}.
#' @references
#' CovidLP Team, 2020. CovidLP: Short and Long-term Prediction for COVID-19. Departamento de Estatistica. UFMG,
#' Brazil. URL: \url{http://est.ufmg.br/covidlp/home/en/}
NULL
