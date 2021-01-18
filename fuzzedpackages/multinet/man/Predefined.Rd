\name{multinet.predefined}
\alias{multinet.predefined}
\alias{ml_empty}
\alias{ml_aucs}
\alias{ml_bankwiring}
\alias{ml_florentine}
\alias{ml_monastery}
\alias{ml_tailorshop}
\alias{ml_toy}
\title{
Loading predefined multilayer networks
}
\description{
Creates predefined multilayer networks from the literature.
\itemize{
\item \code{ml_empty} returns an empty multilayer network, not containing any actor, layer, node or edge.
\item \code{ml_aucs} returns the AUCS multiplex network described in \cite{Rossi and Magnani, 2015. "Towards effective visual analytics on multiplex networks". Chaos, Solitons and Fractals. Elsevier.}
\item \code{ml_bankwiring} returns Padgett's Florentine Families multiplex network.
\item \code{ml_florentine} returns Padgett's Florentine Families multiplex network.
\item \code{ml_monastery} returns Sampson's monastery multiplex network.
\item \code{ml_tailorshop} returns Kapferer's' tailorshop multiplex network.
\item \code{ml_toy} returns the toy network used as a running example in \cite{Dickison, Magnani and Rossi. "Multilayer Social Networks". Cambridge University Press.}
}
}
\usage{
ml_empty(name="")
ml_aucs()
ml_bankwiring()
ml_florentine()
ml_monastery()
ml_tailorshop()
ml_toy()
}
\arguments{
\item{name}{The name of the new multilayer network.}
}
\value{
All these functions return a multilayer network.
}
\references{
\itemize{
\item Rossi, Luca, and Magnani, Matteo (2015). Towards effective visual analytics on multiplex and multilayer networks. Chaos, Solitons and Fractals, 72, 68-76. (for ml_aucs()).
\item Padgett, John F., and McLean, Paul D. (2006). Organizational Invention and Elite Transformation: The Birth of Partnership Systems in Renaissance Florence. American Journal of Sociology, 111(5), 1463-1568. (for ml_florentine()).
\item Breiger, R. and Boorman, S. and Arabic, P. (1975). An algorithm for clustering relational data with applications to social network analysis and comparison with multidimensional scaling. Journal of Mathematical Psychology, 12 (for ml_monastery() and ml_bankwiring() - these authors prepared the datasets, see multilayer.it.uu.se/datasets.html for references to the data collectors).
\item Kapferer, Bruce (1972). Strategy and Transaction in an African Factory: African Workers and Indian Management in a Zambian Town. Manchester University Press (for ml_tailorshop()).
}
}
\seealso{
\link{multinet.IO}, \link{multinet.generation}
}
\examples{
empty <- ml_empty("new network")
aucs <- ml_aucs()
bankwiring <- ml_bankwiring()
florentine <- ml_florentine()
monastery <- ml_monastery()
tailorshop <- ml_tailorshop()
}
