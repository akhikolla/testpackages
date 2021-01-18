\name{MAT}
\alias{MAT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Multidimensional Adaptive Testing (MAT)
}
\description{
  MAT is a package to simulate multidimensional adaptive testing for the Multidimensional 3-Parameter Logistic (M3PL) model.
}
\usage{
  MAT(ipar, resp, cors,
      target.content.dist = NULL, content.cat = NULL, ncc = 1, 
      content.order = NULL, p = stop("p is required"), 
      selectionMethod = c("D", "A", "C", "R"), 
      selectionType = c("FISHER", "BAYESIAN"), c.weights = NA, 
      stoppingCriterion = c("CONJUNCTIVE", "COMPENSATORY"), 
      topN = 1, minNI = 10, maxNI = 30, minSE = 0.3, D = 1, 
      maxIter = 30, conv = 0.001, minTheta = -4, maxTheta = 4, 
      plot.audit.trail = T, theta.labels = NULL, easiness=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ipar}{a data frame containing M3PL item parameters, specifically a1, a2, ... , d, and c}
  \item{resp}{a data frame (that will be converted to a numeric matrix) of item responses, e.g., R1, R2, ..., R180}
  \item{cors}{a square matrix of the lower diagonal elements of a variance-covariance (VCV) matrix, 
   including 1's in the main diagonal}
  \item{target.content.dist}{an optional vector of target content distributions summed to 1.0, e.g., c(0.25,0.5,0.25)}
  \item{content.cat}{an optional vector specifying content designations}
  \item{ncc}{the number of content categories (default=1, i.e., no content balancing)}
  \item{content.order}{an optional vector specifying administration order of content categories, e.g., c(3,1,2)}
  \item{p}{the number of latent dimensions}
  \item{selectionMethod}{item selection criterion: "D"=D-optimality, "A"=A-optimality, "C"=C-optimality, "R"=Random (default="D")}
  \item{selectionType}{item selection method type: "FISHER"=Fisher information, "BAYESIAN"=adds inverse prior VCV}
  \item{c.weights}{an optional vector of weights of length p when selectionMethod="C"}
  \item{stoppingCriterion}{stopping criterion: "CONJUNCTIVE"=SEs for all dimensions must be met, "COMPENSATORY"=the generalized variance or SEs weighted by c-weights must be met}
  \item{topN}{Randomesque exposure control: selects an item randomly from the top N most informative items (default=1, no exposure control)}
  \item{minNI}{minimum number of items to administer (default=10)}
  \item{maxNI}{maximum number of items to administer (default=30)}
  \item{minSE}{minimum SE for stopping (default=0.3)}
  \item{D}{scaling constant: 1.7 or 1.0 (default=1.0)}
  \item{maxIter}{maximum number of Fisher scoring (default=30)}
  \item{conv}{convergence criterion for Fisher scoring (default=0.001)}
  \item{minTheta}{minimum theta value for plotting (default=-4)}
  \item{maxTheta}{maximum theta value for plotting (default=4)}
  \item{plot.audit.trail}{show CAT audit trail: T or F (default=T)}
  \item{theta.labels}{theta labels for plotting (default=c("Theta 1","Theta 2",...))}
  \item{easiness}{logical, T if d is related to the \emph{easiness} of items per Reckase, F otherwise}
}

\details{
  The purpose of this function is to simulate multidimensional adaptive testing based on the Multidimensional 3-Parameter Logistic (M3PL) model (Reckase, 2009):
  \deqn{P_i(\theta) \equiv P(U_i = 1|\boldsymbol{\theta}, \mathbf{a}_i, d_i, c_i) \equiv c_i + \frac{1-c_i}{1 + exp[-D(\mathbf{a}_i\cdot\boldsymbol{\theta} + d_i)]}}
  where \eqn{\mathbf{a}_i} is a vector of discrimination parameters of item \emph{i}, \eqn{\boldsymbol{\theta}} is a vector of abilities, \eqn{c_i} is a scalar representing the 
  guessing parameter of item \emph{i}, \eqn{d_i} is a scalar representing the easiness of item \emph{i}.
  Thetas are estimated using the Bayesian maximum a posteriori (MAP) estimator and the Fisher scoring method. Three item
  selection criteria are available: D-optimality, A-optimality, and C-optimality (Segall, 1996; van der Linden, 1999;
  Mulder & van der Linden, 2009). An option is provided to add the inverse of a prior variance-covariance
  matrix to the multivariate information matrix (selectionType="BAYESIAN"). The stopping condition can be specified as a conjunctive criterion or a compensatory
  criterion. Content balancing can be imposed by specifying target content distributions. An exposure control option is
  provided via the randomesque technique.
}
\value{
  Returns a list of class "MAT" with the following components:
%%  If it is a LIST, use
    \item{call}{function call stack}
    \item{items.used}{a matrix of items administered}
	\item{selected.item.resp}{a matrix containing item responses for selected items}
	\item{ni.administered}{a vector of the number of items administered}
	\item{theta.CAT}{a matrix of theta estimates from CAT}
	\item{se.CAT}{a matrix of SE estimates from CAT}
	\item{theta.history}{a matrix of theta history from CAT}
	\item{se.history}{a matrix of SE history from CAT}
	\item{theta.Full}{a matrix of theta estimates based on the full bank}
	\item{se.Full}{a matrix of SE estimates based on the full bank}
	\item{ipar}{a matrix of item parameters}
	\item{p}{the number of latent dimensions}
}
\note{
\enumerate{
  \item The MAT function performs a number of checks to determine if the arguments for content balancing and content ordering have been specified correctly. If the arguments have not been specified correctly, content balancing and/or content ordering will not be used for the simulation. Additionally, a warning message will be printed to the console detailing the misspecification.  
  \item Content ordering is only available for fixed-length CAT. Namely, to invoke a particular content order, the user must set the minimum number of items equal to the maximum number of items (e.g., minNI=30 & maxNI=30).
}
}

\references{
  \enumerate{
    \item Segall, D. O. (1996). Multidimensional adaptive testing, Psychometrika, 61(2), 331-354
    \item van der Linden, W. J. (1999). Multidimensional adaptive testing with a minimum error-variance criterion,
          Journal of Educational and Behavioral Statistics, 24(4), 398-412.
    \item Mulder, J., & van der Linden, W. J. (2009). Multidimensional adaptive testing with optimal design criteria
          for item selection, Psychometrika, 74(2), 273-296.
    \item Reckase, M. D. (2009). Multidimensional Item Response Theory. New York: Springer.
  }
}
\author{
  Seung W. Choi and David R. King
}
\note{
  requires \pkg{MASS}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
  \dontrun{MCAT.1<-MAT(ipar1,
              resp1,
              vcv1,
              target.content.dist=target.content.dist1,
              content.cat=content.cat1,
              ncc=3,
              p=3,
              selectionMethod="A",
              topN=1,
              selectionType="FISHER",
              stoppingCriterion="CONJUNCTIVE",
              minNI=10,
              maxNI=30)
	}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Psychometrics }
\keyword{ Computerized Adaptive Testing }% __ONLY ONE__ keyword per line
