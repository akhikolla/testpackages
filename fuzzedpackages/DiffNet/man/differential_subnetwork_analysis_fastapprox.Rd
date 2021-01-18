\name{differential_subnetwork_analysis_fastapprox}
\alias{differential_subnetwork_analysis_fastapprox}
\title{
Fast-Approximation Approach for Identifying Differential Sub-networks in Paired Graphs
}
\description{
This method identifies the differential sub-network between two graphs using the proposed Fast-Approximation approach of Mall et al paper.
}
\usage{
differential_subnetwork_analysis_fastapprox(ghd_val, mu_perm, p, 
                                            matrixA, matrixB, threshold)
}

\arguments{
  \item{ghd_val}{
  Generalized Hamming Distance value calculated using topological graphs of g_A and g_B.
}
  \item{mu_perm}{
  Asymptotic value of mean permutation for graph g_A. 
}
  \item{p}{
  Represents the number of nodes in graph g_A which is the same as number of nodes in graph g_B.
}
  \item{matrixA}{
  Topological matrix obtained from graph g_A.
}
  \item{matrixB}{
  Topological matrix obtained from graph g_B.
}
  \item{threshold}{
    Threshold after which the "fast-approx" technique switches to use a model selection criterion similar to the "original" approach to identify statistically significant changes between two networks. By default its value is 1e-250 and a good range for this value is between 1e-50 to 1e-250.
}
}

\value{
  A data frame comprising of:
  \item{actual_id}{
    Id of a node from the set of nodes in g_A
  }
  \item{dim_name}{
    Name associated with a node from the set of nodes in g_A.
  }
  \item{p_val}{
    P-value associated with that node.
  }
  \item{ghd_val}{
    Generalized Hamming Distance between the topological matrices after removal of that node.
  }
  \item{mu_perm}{
    Asymptotic first order moment: mean value.
  }
  \item{std_perm}{
    Asymptotic second order moment: standard deviation value.
  }
  \item{V7}{
    Adjusted p-value associated with that node.
  }
}
\references{
\url{https://www.researchgate.net/publication/305355965_Detection_of_statistically_significant_network_changes_in_complex_biological_networks}
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}

\seealso{
\code{\link{differential_subnetwork_analysis_original}}, \code{\link{differential_subnetwork_analysis_closedform}}
}