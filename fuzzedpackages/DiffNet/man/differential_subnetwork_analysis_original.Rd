\name{differential_subnetwork_analysis_original}
\alias{differential_subnetwork_analysis_original}
\title{
Original dGHD Approach for Identifying Differential Sub-networks in Paired Graphs
}
\description{
This method identifies the differential sub-network between two graphs using the original dGHD approach of Ruan et al paper.
}
\usage{
differential_subnetwork_analysis_original(ghd_val, mu_perm, p, 
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
    Not used in the original approach.
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
Ruan, D., Young, A. and Montana, G., 2015. 
Differential analysis of biological networks. \emph{BMC bioinformatics}, 
\bold{16(1)}, 327-334.
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}

\seealso{
\code{\link{differential_subnetwork_analysis_closedform}}, \code{\link{differential_subnetwork_analysis_fastapprox}}
}
