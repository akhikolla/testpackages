\name{IrishDirectoratesData}
\alias{IrishDirectoratesData}
\docType{data}
\title{
Board Composition For Companies Quoted On The Irish Stock Exchange From 2003 To 2013
}
\description{
Board composition for companies quoted on Irish Stock Exchange from 2003 to 2013.
Board compositions are only observed at the end of each year.
}
\usage{data("IrishDirectoratesData")}
\format{
  \code{IrishDirectoratesData} is a list containing:
  \describe{
    \item{\code{edgelist}}{the edgelist for a bipartite dynamic network. Each row of this dataframe corresponds to an undirected edge in the network. For each row, the first entry identifies the time frame where the edge occurs, the second entry represents the director, whereas the third identifies the company. The presence of an edge at a certain time frame between a director and a company means that the director was part of the company's board at the end of the corresponding year. }
    \item{\code{years}}{lookup table for the time frame labels.}
    \item{\code{directors_names}}{lookup table for directors' names.}
    \item{\code{companies_names}}{lookup table for companies' names.}
  }
}
\details{
The adjacency cube can be constructed from the edgelist. Please see example for sample code.
}
\source{
Irish Stock Exchange (http://www.ise.ie/).
}
\references{
Friel, N., Rastelli, R., Wyse, J. and Raftery, A.E. (2016) <DOI:10.1073/pnas.1606295113>.
}
\examples{
data(IrishDirectoratesData)
attach(IrishDirectoratesData)

N <- length(directors_names)
M <- length(companies_names)
tframes <- length(years)

# construct the binary adjacency cube
adj <- array(0,c(N,M,tframes))
for (l in 1:nrow(edgelist)) adj[edgelist[l,2],edgelist[l,3],edgelist[l,1]] = 1
dimnames(adj) = list(directors_names, companies_names, years)

# calculate the degrees of directors and boards
out_degrees <- apply(adj,c(1,3),sum)
in_degrees <- apply(adj,c(2,3),sum)

# create a binary matrix with ones corresponding to interlocked directors
interlocked_directors <- ifelse(out_degrees > 1, 1, 0)

# create a binary matrix with ones corresponding to interlocking companies
interlocking_companies <- matrix(0,M,tframes)
for (t in 1:tframes) for (i in 1:N) for (j in 1:M) if (adj[i,j,t] == 1) {
  if (interlocked_directors[i,t] > 0) interlocking_companies[j,t] = 1
}

# extract labels of interlocking companies
selected_companies <- which(rowSums(interlocking_companies) > 0)

# extract labels of remaining active directors
new_out_degrees <- apply(adj[,selected_companies,], c(1,3), sum)
selected_directors <- which(rowSums(new_out_degrees) > 0)

# create the new adjacency cube for the reduced data, as shown in the referenced paper
adj_reduced <- adj[selected_directors, selected_companies, ]
}
\keyword{datasets}
