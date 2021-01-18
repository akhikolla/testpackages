\name{mvsvdl1}
\alias{mvsvdl1}
\title{Multi-view bi-clustering via SSVD}
\usage{mvsvdl1(datasets, lvs, lz, maxOuter=100000, thresOuter=.00001, maxInner=10000, 
thresInner=.00001, logLvl=0)}
\description{
	Identify consistent sample cluster among all views and simultaneously associated feature clusters per view. Clusters are obtained via multi-view sparse singular value decomposition (SSVD). One sample cluster and its associated feature clusters are identified and returned through each call of this function. If multiple clusters are desired, call this function repeatedly with samples left unclustered.
}
\details{
	This method is multi-view sparse singular value decomposition (SSVD) based. Consistent sample cluster among views and associated view specific feature clusters are simultaneously identified. Sample and feature clusters are associated in the sense that they help determine each other. The sparsity of the singular vectors in the decomposition is enforced using L1-norm. Please refer to the referenced paper for more details.
}
\arguments{
	\item{datasets}{
		List of input data, where each element is a matrix. For all the matrices, rows represent samples and columns represent features. These matrices are the characterization of a same set of samples from different perspectives(views), one matrix per view. So all the matrices should have the exact same rows, but can have different columns. The rows in the same position in all matrices represent a same sample.
	}
	\item{lvs}{
		A numerical vector with length that equals to the number of view, which controls the sparsity of the right singular vector in the SVD of each view. This in turn controls the size of the feature clusters. Larger value indicates stronger sparsity and thus leads to smaller feature clusters.
	}
	\item{lz}{
		A number, which controls the sparsity of the left singular vector in the SSVD of all views. It helps to tune the size of the identified sample cluster. A larger value means stronger sparsity and thus a smaller sample cluster.		
	}
	\item{maxOuter}{
		(Optional) Maximum number of outer loop iterations, which is one of the criteria that controls when to terminate the outer loop in the process of searching for an optimal solution. Please read the referenced paper for details. The default value is 100000.
	}
	\item{thresOuter}{
		(Optional) the other criteria (besides the above 'maxOuter' argument) for terminating the outer loop. When the sum of squares of the difference between two consecutive outer loop iterations of vector z passes (is smaller than) this threshold, the loop is terminated. The default value is 0.00001.
	}
	\item{maxInner}{
		(Optional) Maximum number of inner loop iterations. It works the same way as above 'maxOuter' argument, but controls the inner loop in the optimization process. The default value is 10000.
	}
	\item{thresInner}{
		(Optional) This works the same way as the above 'thresOuter' argument, but looks at the sum of squares of the difference between two consecutive iterations of vector u (the other multiplier (besides vector z) of left singular vector) and controls the inner loop in the optimization process. The default value is 0.00001.
	}
	\item{logLvl}{
		(Optional) Logging level, which can be set to 0, 1 and 2. This controls the amount of printing, where a larger value means more printing.  The default value is 0, which turns off the logging. 
	}
}
\value{
	A list with following named fields:
	\item{Cluster}{A binary vector (with 1 or 0 entries) of length equal to the sample size. It indicates whether a sample is in the identified cluster (with 1 in its corresponding position) or not (with 0).} 
	\item{FeatClusters}{A list of binary vectors that give identified feature cluster for each view. Value 1 indicates the corresponding feature belongs to the identified feature cluster.}
	\item{U}{A matrix with size of n x m, where n is the number of samples, m is the number of views. So each column of this matrix is one of the two multipliers of the left singular vector in the SSVD of corresponding view.}
	\item{V}{A list of vectors that give the right singular vector in the SSVD of each view.} 
	\item{z}{The common multiplier of the left singular vector in SSVD of all views.}
}
\references{
	Jiangwen Sun, Jinbo Bi, and Henry R. Kranzler \emph{Multi-view Singular Value Decomposition for Disease Subtyping and Genetic Associations}
	BMC Genetics, 15(73):1-12, 2014	
}
\examples{
		library(mvcluster)
		data(phe)
		data(gen)
        views <- list(phe,gen)
        result <- mvsvdl1(views,c(0.45,0.65),0.016)
}