\name{pqlseq} 
\alias{pqlseq} 
\title{Fit Generalized Linear Mixed Model with Known Kinship Matrices Through Penalized-quasi Likelihood} 
\description{
  Fit a generalized linear mixed model with a random intercept. The covariance matrix of the random intercept is proportional to a known kinship matrix. 
}
\usage{ 
  pqlseq(RawCountDataSet, Phenotypes, Covariates=NULL,
  RelatednessMatrix=NULL, LibSize=NULL, fit.model="PMM",
  fit.method = "AI.REML", fit.maxiter=500, fit.tol=1e-5, 
  numCore=1, filtering=TRUE, verbose=FALSE,...)  
} 
\arguments{   
\item{RawCountDataSet}{a data frame containing the read count.}   

\item{Phenotypes}{a vector containing the predictor of interest.}

\item{Covariates}{a data frame containing the covariates subject to adjustment (Default = NULL).}

\item{RelatednessMatrix}{a known relationship matrix (e.g. kinship matrix in genetic studies). When supplied with a matrix, this matrix should be a positive semi-definite matrix with dimensions equal to the
sample size in count data, and the order of subjects in this matrix should also match the order of subjects in count data. Currently there is no ID checking feature implemented, and it is the user's responsibility to match the orders. } 

\item{LibSize}{a data frame containing the total read count. For
possion mixed model, it will be calculated automatically if users do not provide. For binomial mixed model, it is required.  }

\item{fit.model}{a description of the error distribution and link function to be used in the model. Either "PMM" for possion model, or "BMM" for binomial model (default = "PMM").}

\item{fit.method}{method of fitting the generalized linear mixed model, currently only "REML" version is available.}

\item{fit.maxiter}{a positive integer specifying the maximum number of iterations when fitting the generalized linear mixed model (default = 500).}

\item{fit.tol}{a positive number specifying tolerance, the difference threshold for parameter estimates below which iterations should be stopped (default = 1e-5).}

\item{numCore}{a positive integer specifying the number of cores for parallel computing (default = 1).}

\item{filtering}{a logical switch for RNAseq data. By default, for each gene, at least two individuals should have read counts greater than 5. Otherwise, the gene is filtered (default = TRUE).}

\item{verbose}{a logical switch for printing detailed information (parameter estimates in each iteration) for testing and debugging purpose (default = FALSE).}

\item{\dots}{additional arguments that could be passed to glm.}
}

\details{
Generalized linear mixed models (GLMM) are fitted using the penalized quasi-likelihood (PQL) method proposed by Breslow and Clayton (1993). Statistical inference in GLMM is notoriously difficult because of an intractable high-dimensional integral in the likelihood (Chen, 2016 and Lea, 2015), and by default we use the Average Information REML algorithm (Gilmour, Thompson and Cullis, 1995; Yang et al., 2011) to fit the model. An eigen-decomposition is performed in each outer iteration and the estimate of the variance component parameter \eqn{\tau} is obtained by maximizing the profiled log restricted likelihood. When the Average Information REML algorithm fails to converge, a warning message is given and the algorithm is default to INLA approaches (Rue, 2009).
}
\value{
  \item{numIDV}{number of individuals with data being analyzed}
  \item{beta}{the fixed effect parameter estimate for the predictor of interest.}
  \item{se_beta}{the standard deviation of fixed effect.}
  \item{pvalue}{P value for the fixed effect, based on the wald test.}
  \item{h2}{heritability of the transformed rate.}
  \item{sigma2}{total variance component.}
  \item{overdisp}{dispersion parameter estimate}
  \item{converged}{a logical indicator for convergence.}
}

\references{
Breslow, N.E. and Clayton, D.G. (1993) Approximate Inference in Generalized Linear Mixed Models. Journal of the American Statistical Association 88, 9-25.

Chen, H., Wang, C., Conomos, M.P., Stilp, A.M., Li, Z., Sofer, T., Szpiro, A.A., Chen, W., Brehm, J.M., Celedon, J.C., Redline, S., Papanicolaou, G.J., Thornton, T.A., Laurie, C.C., Rice, K. and Lin, X. Control for population structure and relatedness for binary traits in genetic association studies using logistic mixed models. The American Journal of Human Genetics, 98, 653-666.

Gilmour, A.R., Thompson, R. and Cullis, B.R. (1995) Average Information REML: An Efficient Algorithm for Variance Parameter Estimation in Linear Mixed Models. Biometrics 51, 1440-1450.

Lea,A., Tung, J. and Zhou,X. (2015) A flexible, effcient binomial mixed model for identifying differential DNA methylation in bisulfite sequencing data. PLoS Genetics. 11: e1005650.

Rue, H., Martino, S. , and Chopin, N.(2009) Approximate Bayesian inference for latent Gaussian models using integrated nested Laplace approximations (with discussion). Journal of the Royal Statistical Society, Series B, 71(2):319-392

Yang, J., Lee, S.H., Goddard, M.E. and Visscher, P.M. (2011) GCTA: A Tool for Genome-wide Complex Trait Analysis. The American Journal of Human Genetics 88, 76-82.

Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis for association studies. Nature Genetics 44, 821-824.
}
\author{
Shiquan Sun, Jiaqiang Zhu, Xiang Zhou
}

\examples{
data(ExampleRNAseq)
attach(ExampleRNAseq)
model_RNA=pqlseq(RawCountDataSet=count, Phenotypes=predictor, 
  RelatednessMatrix=relatednessmatrix, LibSize=totalcount,
  fit.model="PMM",numCore=1)
head(model_RNA)
detach(ExampleRNAseq)
}

\keyword{function}
\concept{GLMMs}
\concept{pqlseq}
