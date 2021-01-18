\name{pred_rgasp}
\alias{pred_rgasp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Prediction for robust GaSP model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function to make prediction on robust GaSP models after the robust GaSP model has been constructed.
}
\usage{
pred_rgasp(beta, nu, input, X, zero_mean,output, testing_input,
           X_testing, L, LX, theta_hat, sigma2_hat, 
           q_025, q_975, r0, kernel_type, alpha,method,interval_data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
%%     ~~Describe \code{beta} here~~
inverse-range parameters.
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
noise-variance ratio parameter.
}
  \item{input}{
%%     ~~Describe \code{input} here~~
input matrix.
}
  \item{X}{
%%     ~~Describe \code{X} here~~
the mean basis function i.e. the trend function.
}
  \item{zero_mean}{
%%     ~~Describe \code{zero_mean} here~~
The mean basis function is zero or not.
}

  \item{output}{
%%     ~~Describe \code{output} here~~
output matrix.
}
  \item{testing_input}{
%%     ~~Describe \code{testing_input} here~~
testing input matrix.
}
  \item{X_testing}{
%%     ~~Describe \code{X_testing} here~~
mean/trend matrix of testing inputs.
}
  \item{L}{
%%     ~~Describe \code{L} here~~
a lower triangular matrix for the cholesky decomposition of \code{R}, the correlation matrix.
}
  \item{LX}{
%%     ~~Describe \code{LX} here~~
a lower triangular matrix for the cholesky decomposition of $X^tR^{-1}X$.
\ifelse{html}{\out{X^tR^{-1}X}}{\eqn{4(X^tR^{-1}X}{X^tR^{-1}X}}
}
  \item{theta_hat}{
%%     ~~Describe \code{theta_hat} here~~
estimated mean/trend parameters.
}
  \item{sigma2_hat}{
%%     ~~Describe \code{sigma2_hat} here~~
estimated variance parameter.
}
  \item{q_025}{
%%     ~~Describe \code{qt_025} here~~
0.025 quantile of \code{t} distribution.
}
  \item{q_975}{
%%     ~~Describe \code{qt_975} here~~
0.975 quantile of \code{t} distribution.
}
  \item{r0}{
%%     ~~Describe \code{r0} here~~
a matrix of absolute difference between inputs and testing inputs.
}
  \item{kernel_type}{
%%     ~~Describe \code{kernel_type} here~~
type of kernel. \code{matern_3_2} and \code{matern_5_2} are \code{Matern kernel} with roughness parameter 3/2 and 5/2 respectively. \code{pow_exp} is power exponential kernel with roughness parameter alpha. If \code{pow_exp} is to be used, one needs to specify its roughness parameter alpha.
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
Roughness parameters in the kernel functions.
}
\item{method}{  
method of parameter estimation. \code{post_mode} means the marginal posterior mode is used for estimation. \code{mle} means the maximum likelihood estimation is used. \code{mmle} means the maximum marginal likelihood estimation is used. The \code{post_mode}  is the default method. 
}
\item{interval_data}{
a boolean value. If \code{T}, the interval of the data will be calculated. If \code{F}, the  interval of the mean of the data will be calculated. 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of 4 elements. The first is a vector for predictive mean for testing inputs. The second is a vector for lower quantile for 95\% posterior credible interval and the  third is the upper quantile for 95\% posterior credible interval for these testing inputs. The last is a vector of standard deviation of each testing inputs.
}
\references{
Mengyang Gu. (2016). Robust Uncertainty Quantification and Scalable Computation for Computer Models with Massive Output. Ph.D. thesis. Duke University.
}
\author{
%%  ~~who you are~~
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{predict.rgasp}}
}




% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")

\keyword{internal}