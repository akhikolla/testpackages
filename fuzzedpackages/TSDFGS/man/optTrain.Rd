\name{optTrain}
\alias{optTrain}

\title{
Algorithm for optimal training set determination
}
\description{
It uses a genetic algorithm or simple exchange algorithm with three different criteria (r-score (J.H. Ou et al., (2019) <DOI:10.6342/NTU201802290>), PEV-score (Akdemir D. et al., (2015) <DOI:10.1186/s12711-015-0116-6>), CD-score (Laloe D. (1993) <DOI:10.1186/1297-9686-25-6-557>)) to determine an optimal training set.
}
\usage{
optTrain(geno, cand, n.train, subpop = NULL, 
        test = NULL, method = "rScore", min.iter = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{
A numeric matrix of principal components (rows: individuals; columns: PCs). To reduce computing time, one may use first k PCs by geno[,1:k].
}
  \item{cand}{
An integer vector of which rows of individuals are candidates of the training set in the geno matrix.
}
  \item{n.train}{
The size of the target training set.
}
  \item{subpop}{
A character vector of subpopulation's group name. The algorithm will ignore the population structure if it remains NULL.
}
  \item{test}{
An integer vector of which rows of individuals are in the test set in the geno matrix. The algorithm will use an un-target method if it remains NULL.
}
  \item{method}{
Choices are rScore, PEV and CD. rScore will be used by default.
}
  \item{min.iter}{
Minimum iteration of all methods can be appointed. One should always check if the algorithm is converged or not. A minimum iteration will set by considering the candidate and test set size if it remains NULL.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{OPTtrain}{An integer vector of the chosen optimal training set.}
  \item{TOPscore}{Score of each iteration. (Given by one of three criterions)}
  \item{ITERscore}{Score of the best solution in by far. (Given by one of three criterions.)}
}
\references{
Akdemir D., Sanchez JI., Jannink JL. (2015), Optimization of genomic selection training populations with a genetic algorithm. GenSetic Selection Evolution 47:38.\\

Laloe D. (1993), Precision and information in linear models of genetic evolution. Genetics Selection Evolution 25:557.\\

Holland J. H. (1975) Adaptation in Natural and Artificial Systems. University of Michigan Press.
}
\author{
\packageAuthor{TSDFGS}

Maintainer: \packageMaintainer{TSDFGS}
}
\note{
Both genetic algorithm and simple exchange algorithms do not assure convergence to global optimal, and it is highly recommended to draw the convergence plot to check it converges to the local optimal.
}
%\seealso{}
\examples{
  \donttest{
    ## LOAD EXAMPLE DATA ##
    data("rice44kPCA")
    out.RNN = optTrain(geno, cand = 1:404, n.train = 100)
  } 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.S
%\keyword{Optimal training set}
