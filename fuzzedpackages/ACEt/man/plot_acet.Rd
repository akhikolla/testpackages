\name{plot_acet}
\alias{plot_acet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot variance curves for various ACE(t) or ACE(t)-p models with 95\% CIs 
}
\description{
Plot variance curves or a heritability curve (with 95\% CIs) of the A, C and E components with respect to age modelled by B-splines or P-splines.
}
\usage{
plot_acet(acet, boot = FALSE, heri = FALSE, xlab, ylab, main, col, legend = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{acet}{
An object obtained from the B-splines or P-splines functions. For the P-splines functions, an object from the MCMC method must be used.
}
  \item{boot}{
An logical  indicator of whether the confidence bands estimated from the bootstrap method are plotted. The default is FALSE. Only available for the AtCtEt model.
}
  \item{heri}{
A logical indicator of whether to plot the dynamic heritability curve. The default is FALSE.
}
  \item{xlab}{
The 'xlab' argument in the plot funciton. The default is 'Age'.
}
  \item{ylab}{
The 'ylab' argument in the plot funciton. The default is 'Variance'.
}
  \item{main}{
The 'main' argument in the plot funciton. The default is 'Variance curves of the A, C, and E components' for variances and 'Dynamic heritability' for heritability.
}
  \item{col}{
The 'col' argument in the plot funciton.
}
  \item{legend}{
An logical indicator of whether the default legend is plotted.
}
}


\references{
He, L., Sillanpää, M.J., Silventoinen, K., Kaprio, J. and Pitkäniemi, J., 2016. Estimating Modifying Effect of Age on Genetic and Environmental Variance Components in Twin Models. Genetics, 202(4), pp.1313-1328.
}
\author{
Liang He
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

# data(data_ace)
# result <- AtCtEtp(data_ace$mz, data_ace$dz, knot_a = 5, knot_c = 4)
# result_mc <- acetp_mcmc(result, iter_num=10000, burnin = 500)
# plot_acet(result_mc)

# result <- AtCtEt(data_ace$mz, data_ace$dz, mod=c('d','c','c'), knot_a = 9)
# plot_acet(result)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

