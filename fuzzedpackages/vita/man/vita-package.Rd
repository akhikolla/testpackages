\name{vita-package}
\alias{vita-package}

\docType{package}
\title{
Variable importance testing approaches (vita)
}
\description{
Implements the novel testing approach by Janitza et al.(2015) for the permutation variable importance measure
in a random forest and the PIMP-algorithm by Altmann et al.(2010). Janitza  et al.(2015) do not use
the "standard" permutation variable importance but the cross-validated permutation variable
importance for the novel test approach. The cross-validated permutation variable importance
is not based on the out-of-bag observations but uses a similar strategy which is inspired by
the cross-validation procedure. The novel test approach can be applied for classification trees
as well as for regression trees.
However, the use of the novel testing approach has not been tested for regression trees so far, so
this routine is meant for the expert user only and its current state is rather experimental.
}
\details{


\strong{The novel test approach (\code{\link[vita]{NTA}})}:

The observed non-positive permutation variable importance values are used to approximate the distribution of
variable importance for non-relevant variables. The null distribution Fn0 is computed by mirroring the
non-positive variable importance values on the y-axis. Given the approximated null importance distribution,
the p-value is the probability of observing the \code{original PerVarImp} or a larger value. This testing
approach is suitable for data with large number of variables without any effect.

\code{PerVarImp} should be computed based on the hold-out permutation variable importance measures. If using
standard variable importance measures the results may be biased.

This function has not been tested for regression tasks so far, so this routine is meant for the expert user
only and its current state is rather experimental.

\strong{Cross-validated permutation variable importance (\code{\link[vita]{CVPVI}})}:

This method randomly splits the dataset into k sets of equal size. The method constructs k random forests, where the l-th forest is constructed based on observations that are not part of the l-th set. For each forest the fold-specific permutation variable importance measure is computed using all observations in the l-th data set: For each tree, the prediction error on the l-th data set is recorded. Then the same is done after permuting the values of each predictor variable.

The differences between the two prediction errors are then averaged over all trees. The cross-validated permutation variable importance is the average of all k-fold-specific permutation variable importances. For classification the mean decrease in accuracy over all classes is used and for regression the mean decrease in MSE.

\strong{PIMP testing approach (\code{\link[vita]{PIMP}})}:

The PIMP-algorithm by Altmann et al.(2010) permutes \eqn{S} times the response variable \eqn{y}.
For each permutation of the response vector \eqn{y^{*s}}, a new forest is grown and the permutation
variable importance measure (\eqn{VarImp^{*s}}) for all predictor variables \eqn{X} is computed.
The vector \code{perVarImp^{s}} for every predictor variables are used to approximate the null importance distributions.

Given the fitted null importance distribution, the p-value is the probability of observing the \emph{original VarImp} or a larger value.


}
\author{
Ender Celik
}

\references{
Breiman L. (2001), \emph{Random Forests}, Machine Learning 45(1),5-32, <doi:10.1023/A:1010933404324>

Altmann A.,Tolosi L., Sander O. and  Lengauer T. (2010),\emph{Permutation importance: a corrected feature importance measure}, Bioinformatics Volume 26 (10), 1340-1347, <doi:10.1093/bioinformatics/btq134>

Janitza S, Celik E, Boulesteix A-L, (2015),\emph{ A computationally fast variable importance test for random forest for high dimensional data},Technical Report 185, University of Munich <\url{http://nbn-resolving.de/urn/resolver.pl?urn=nbn:de:bvb:19-epub-25587-4}>
}

\keyword{ package }
\seealso{
\code{\link[vita]{PIMP}}, \code{\link[vita]{NTA}}, \code{\link[vita]{CVPVI}}, \code{\link[randomForest]{importance}}, \code{\link[randomForest]{randomForest}}
}
