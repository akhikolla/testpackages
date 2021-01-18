\name{IMaGES}
\alias{IMaGES}
\encoding{UTF-8}
\concept{Independent Multi-sample Greedy Equivalence Search}
\concept{essential graph}
\title{Estimate the DAG and Markov equivalence class of multiple datasets using
the IMaGES algorithm}
\description{
  Estimate the DAG using the IMaGES algorithm of Ramsey, Hanson, Hanson, Halchenko, Poldrack, Glymour (2010). Adapted from 'GES' as implemented in the 'pcalg' package.
}
\usage{
  IMaGES(matrices, penalty=3, num.markovs=5,
          use.verbose=FALSE)
}
\arguments{
  \item{matrices}{A list of (at least one) matrix object(s) to run in the algorithm}
  \item{penalty}{Node labels; by default, they are determined from the scoring
    object.}
  \item{num.markovs}{Number of graphs to store in the Markov equivalence class.}
  \item{use.verbose}{boolean value indicating whether or not the run will output things such as step direction, current IMScore, step number, and whether or not early stopping is kicking in.}
}
\details{

  IMaGES (Independent Multi-sample Greedy Equivalence Search) is a score-based 
  algorithm that greedily maximizes a score function similar to the one used
  in the pcalg implementation of GES. It modifies the scoring by creating a
  global score across all datasets and uses this score to determine which 
  step from the individual datasets best represents all of the datasets. It 
  accomplishes this by using forward, backward, and turning steps as described
  below:
  \describe{
    \item{Forward phase}{In the forward phase, IMaGES moves through the space of
    essential graphs in steps that correspond to the addition of a single edge
    in the space of DAGs; the phase is aborted as soon as the score cannot be
    augmented any more.}
    \item{Backward phase}{In the backward phase, the algorithm performs moves
    that correspond to the removal of a single edge in the space of DAGs until
    the score cannot be augmented any more.}
    \item{Turning phase}{In the turning phase, the algorithm performs moves
    that correspond to the reversal of a single arrow in the space of DAGs
    until the score cannot be augmented any more.}
    
    
  }
    During each step, the IMaGES algorithm simulates the addition, deletion, or turning
    of an edge for each individual dataset. The step that most augments the score (each
    edge is assigned a SEM BIC score) for the individual dataset is selected, and the
    global step across all datasets is selected by finding the most commonly recommended
    step. The algorithm then executes that step and updated the IMScore accordingly. This
    repeats for (number of nodes)**2 or until the algorithm detects that no steps augment 
    the score for five consecutive iterations (also known as 'early stopping').
    
    These operations are carried out and result in a global representative graph
    and a Markov Equivalence Class.

}
\value{
  \code{IMaGES} returns a list with the following components:
  \item{results}{A named list containing: 
       \describe{
       \item{\code{global}}{a named list containing \code{.graph}, the global graphNEL
       object, and \code{.params}, the Structural Equation Modeling data for the global graph}
       }
       \describe{
       \item{\code{.single.graphs}}{a list containing named lists of the same structure as above that        corresponds to each individual dataset passsed into IMaGES}
       }
       \describe{
       \item{\code{.markovs}}{a list containing named lists of the same structure as above that              corresponds to each Markov Equivalence Class (the size of which is specified by the user)}
       }
       \describe{
       \item{\code{.means}}{a list containing the mean values for the SEM data for each edge in
       the graph}
       }
       \describe{
       \item{\code{.std.errs}}{a list containing the standard errors for the SEM data for each edge in        the graph}
       }}
  \item{.graphs}{A list containing the individual graph information, including raw graph data, score information, and intermediary variables}
}
\references{
  
  J.D. Ramsey, S.J. Hanson, et al. (2010) Six problems for causal inference
  from fMRI. \emph{NeuroImage} \bold{49}, 1545-1558.

  D.M. Chickering (2002).  Optimal structure identification with greedy search.
  \emph{Journal of Machine Learning Research} \bold{3}, 507--554

  A. Hauser and P. Bühlmann (2012).  Characterization and greedy learning of
  interventional Markov equivalence classes of directed acyclic graphs.
  \emph{Journal of Machine Learning Research} \bold{13}, 2409--2464.

  P. Nandy, A. Hauser and M. Maathuis (2015).  Understanding consistency in 
  hybrid causal structure learning.  \emph{arXiv preprint} 1507.02608
  
  P. Spirtes, C.N. Glymour, and R. Scheines (2000).  \emph{Causation,
  Prediction, and Search}, MIT Press, Cambridge (MA).
}
\author{
  Noah Frazier-Logue (\email{n.frazier.logue@nyu.edu})
  
  Stephen Jose Hanson (\email{jose@rubic.rutgers.edu})
}
\seealso{
  \code{\link{plotIMGraph}}, \code{\link{plotAll}}, \code{\link{plotMarkovs}}
}
\examples{
library(graph)
## Load predefined data
data(IMData)
data(IMTrue)

#run IMaGES
im.fits <- IMaGES(matrices=IMData, penalty=3, num.markovs=5)

#plot global graph and all individual graphs with own SEM data
plotAll(im.fits)

#plot true individual graphs as a comparison
for (i in 1:length(IMTrue)) {
  plot(IMTrue[[1]])
}

#plot graph with SEM data from first dataset
plotIMGraph(im.fits$.single.graph[[1]])

#plot MEC
plotMarkovs(im.fits)

#print mean SEM values for each edge
print(im.fits$.means)

#print standard errors for each edge
print(im.fits$.std.errs)

}

\keyword{IMaGES}
\keyword{graphs}
