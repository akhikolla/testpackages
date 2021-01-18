# -----------------------------------------------------------------------------
#                            Functions to Simulate Data
# -----------------------------------------------------------------------------

#' Simulate data from the circle model.
#'
#' Simulate draws from a bernoulli distribution over \code{c(-1,1)}.  First, the
#' predictors \eqn{x} are drawn i.i.d. uniformly over the square in the two dimensional
#' plane centered at the origin with side length \code{2*outer_r}, and then the
#' response is drawn according to \eqn{p(y=1|x)}, which depends
#' on \eqn{r(x)}, the euclidean norm of \eqn{x}.  If
#' \eqn{r(x) \le inner_r}, then \eqn{p(y=1|x) = 1}, if \eqn{r(x) \ge outer_r}
#' then \eqn{p(y=1|x) = 1}, and \eqn{p(y=1|x) = (outer_r - r(x))/(outer_r - inner_r)}
#' when \eqn{inner_r <= r(x) <= outer_r}.  See Mease (2008).
#'
#'
#' @param n Number of points to simulate.
#' @param inner_r Inner radius of annulus.
#' @param outer_r Outer radius of annulus.
#'
#' @return Returns a list with the following components:
#'
#' \item{y}{Vector of simulated response in \code{c(-1,1)}.}
#' \item{X}{An \code{n}x\code{2} matrix of simulated predictors.}
#' \item{p}{The true conditional probability \eqn{p(y=1|x)}.}
#'
#' @references Mease, D., Wyner, A. and Buha, A. (2007). Costweighted
#' boosting with jittering and over/under-sampling:
#' JOUS-boost. J. Machine Learning Research 8 409-439.
#'
#' @examples
#' # Generate data from the circle model
#' set.seed(111)
#' dat = circle_data(n = 500, inner_r = 1, outer_r = 5)
#'
#' \dontrun{
#' # Visualization of conditional probability p(y=1|x)
#' inner_r = 0.5
#' outer_r = 1.5
#' x = seq(-outer_r, outer_r, by=0.02)
#' radius = sqrt(outer(x^2, x^2, "+"))
#' prob = ifelse(radius >= outer_r, 0, ifelse(radius <= inner_r, 1,
#'              (outer_r-radius)/(outer_r-inner_r)))
#' image(x, x, prob, main='Probability Density: Circle Example')
#' }
#'
#' @export
circle_data = function(n = 500, inner_r = 8, outer_r = 28){
  if(outer_r <= inner_r)
    stop('outer_r must be strictly larger than inner_r')

  X = matrix(stats::runif(2*n, -outer_r, outer_r), nrow=n, ncol=2)
  r = apply(X, 1, function(x) sqrt(sum(x^2)))
  p = 1*(r < inner_r) +
    (outer_r-r)/(outer_r-inner_r)*((inner_r < r) & (r < outer_r))
  y = 2*stats::rbinom(n, 1, p) - 1
  list(X=X, y=y, p=p)
}

#' Simulate data from the Friedman model
#'
#' Simulate draws from a bernoulli distribution over \code{c(-1,1)}, where the
#' log-odds is defined according to:
#' \deqn{log{p(y=1|x)/p(y=-1|x)} = gamma*(1 - x_1 + x_2 - ... + x_6)*(x_1 + x_2 + ... + x_6)}
#' and \eqn{x} is distributed as N(0, I_\code{d}x\code{d}).  See Friedman (2000).
#'
#' @param n Number of points to simulate.
#' @param d The dimension of the predictor variable \eqn{x}.
#' @param gamma A parameter controlling the Bayes error, with higher values of
#'        \code{gamma} corresponding to lower error rates.
#'
#' @return Returns a list with the following components:
#' \item{y}{Vector of simulated response in \code{c(-1,1)}.}
#' \item{X}{An \code{n}x\code{d} matrix of simulated predictors.}
#' \item{p}{The true conditional probability \eqn{p(y=1|x)}.}
#'
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2000). Additive logistic
#' regression: a statistical view of boosting (with discussion), Annals of
#' Statistics 28: 337-307.
#'
#' @examples
#' set.seed(111)
#' dat = friedman_data(n = 500, gamma = 0.5)
#'
#' @export
friedman_data = function(n = 500, d = 10, gamma = 10){
  if(d < 6)
    stop('d must be greater than 6')
  X = matrix(stats::rnorm(d*n), nrow=n, ncol=d)
  log_odds = gamma*(1 - X[,1] + X[,2] - X[,3] + X[,4] - X[,5] + X[,6])*
    rowSums(X[,1:6])
  p = exp(log_odds)/(1 + exp(log_odds))
  y = 2*stats::rbinom(n, 1, p) - 1
  list(X=X, y=y, p=p)

}
