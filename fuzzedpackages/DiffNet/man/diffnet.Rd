\name{diffnet}
\alias{diffnet}
\title{
Detection of Statistically Significant Changes in Paired Biological Networks
}
\description{
Performs differential network analysis for paired biological networks to identify statistically significant changes between two graphs. Currently, the approaches available for doing this include the "closed-form", "original" (dGHD) and the "fast-approx" techniques described in the paper of Mall et al. The methods works better for large-scale complex biological networks (in pairs).
}
\usage{
diffnet(g_A = sample_grg(6, 0.15, torus = TRUE, coords = TRUE), 
        g_B = permute(g_A, c(sample(5), 6)), p = 6, 
        threshold = 1e-50, approach = "closed-form")
}
\arguments{
  \item{g_A}{
  An igraph object representing graph g_A
}
  \item{g_B}{
  An igraph object representing the second graph B with same number of nodes.
}
  \item{p}{
  Represents the number of nodes in graph g_A which is the same as number of nodes in graph g_B 
}
  \item{threshold}{
  Threshold after which the "closed-form" and "fast-approx" techniques switch to use a model selection criterion similar to the "original" approach to identify statistically significant changes between two networks. By default its value is 1e-250 and a good range for this value is between 1e-50 to 1e-250.
}
  \item{approach}{
  Either "closed-form"/"original"/"fast-approx". By default its "closed-form"
}
}
\value{
  An ordered vector representing the p-value for each node. Nodes whose p-values are less than 0.01 form the differential sub-networks in paired graphs g_A and g_B.
}
\references{
\url{https://www.researchgate.net/publication/305355965_Detection_of_statistically_significant_network_changes_in_complex_biological_networks}
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}

\examples{
library("DiffNet")
# this step is optional, it helps speed up calculations, run in parallel on 2 processors
library(doParallel)
registerDoParallel(2)
# Run the differential network analysis tehcnique on sample data
p <- diffnet()
}
