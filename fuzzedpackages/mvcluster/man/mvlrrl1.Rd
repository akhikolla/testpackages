\name{mvlrrl1}
\alias{mvlrrl1}
\title{Multi-view bi-clustering via L1-norm enforced sparse LRR}
\usage{mvlrrl1(datasets, lus, lvs, lz, maxOuter=100000, thresOuter=.00001, maxInner=10000,
 thresInner=.00001, logLvl=0)}
\description{
	Identify consistent sample cluster among all views and simultaneously associated feature clusters per view. Clusters are obtained via finding sparse low rank representation (LRR) of input data matrices, where the sparsity is enforced using L1-norm. One sample cluster and its associated feature clusters are identified and returned each time this function is used. If multiple clusters are desired, call this function repeatedly with samples left unclustered.
}
\details{
	This method identifies clusters via finding multi-view sparse low rank representations (LRR) of input data matrices. The LRR is obtained via matrix decomposition in the form: (zu)v, where both z and u are column vectors and their length are equal to the number of samples, (zu) represents the element-wise product between z and u, and v is a column vector and its length is equal to the number of features in the view. Vector z is a multiplier shared across all views, while both u and v are view specific. Sample cluster is read off directly from z by assigning samples with non-zero component in z to the cluster. Because z is shared across all views, the obtained sample cluster is a common cluster among all views. Feature cluster is obtained by assigning features with non-zero component in v to the cluster. Because v is view specific, feature cluster is view specific as well. Sample cluster and feature cluster are associated in the sense that they help determine each other. The sparsity of vector z, u and v is enforced using L1-norm.
}
\arguments{
	\item{datasets}{
		List of input data, where each element is a matrix. For all the matrices, rows represent samples and columns represent features. These matrices are the characterization of a same set of samples from different perspectives (views), one matrix per view. So all the matrices should have the exact same rows, but can have different columns. The rows in the same position in all matrices represent a same sample.
	}
	\item{lus}{
		A numerical vector with length that equals to the number of views, which controls the sparsity of vector u in the decomposition of each view. This parameter helps control the size of the sample cluster. A larger value indicates stronger sparsity and thus leads to a smaller sample cluster.
	}
	\item{lvs}{
		A numerical vector with length that equals to the number of views, which controls the sparsity of vector v in the decomposition of each view. This parameter helps control the size of the feature clusters. A larger value indicates stronger sparsity and thus leads to smaller feature clusters.
	}
	\item{lz}{
		A number, which controls the sparsity of vector z in the decomposition of all views. It helps to tune the size of the identified sample cluster. A larger value means stronger sparsity and thus leads to a smaller sample cluster.		
	}
	\item{maxOuter}{
		(Optional) Maximum number of outer loop iterations, which is one of the criteria that controls when to terminate the outer loop in the process of searching for an optimal solution. The default value is 100000.
	}
	\item{thresOuter}{
		(Optional) The other criteria (besides the above 'maxOuter' argument) for terminating the outer loop. When the sum of squares of the difference between two consecutive outer loop iterations of vector z passes (is smaller than) this threshold, the loop is terminated. The default value is 0.00001.
	}
	\item{maxInner}{
		(Optional) Maximum number of inner loop iterations. It works the same way as above 'maxOuter' argument, but controls the inner loop in the optimization process. The default value is 10000.
	}
	\item{thresInner}{
		(Optional) This works the same way as above 'thresOuter' argument, but looks at the sum of squares of the difference between consecutive vector u (the other multiplier (besides vector z) of left singular vector) and controls the inner loop in the optimization process. The default value is 0.00001.
	}
	\item{logLvl}{
		(Optional) Logging level, which can be set to 0, 1 or 2.  It controls the amount of printing, where a larger value means more printing.  The default value is 0, which turns off the logging. 
	}
}
\value{
	A list with following named fields:
	\item{Cluster}{A binary vector (with 1 or 0 entries) of length equal to the sample size. It indicates whether a sample is in the identified cluster (with 1 in its corresponding position) or not (with 0).}
	\item{FeatClusters}{A list of binary vectors that give identified feature cluster for each view. Value 1 indicates the corresponding feature belongs to the identified feature cluster.}
	\item{U}{A matrix with size of n x m, where n is the number of samples, m is the number of views. So each column gives the u in the decomposition of the corresponding view.}
	\item{V}{A list of vectors that give the v in the decomposition of each view.}
	\item{z}{The common multiplier shared across all views in their decomposition.}
}
\examples{
		library(mvcluster)
		data(view1)
		data(view2)
        views <- list(view1,view2)
        result <- mvlrrl1(views,c(1.2,1.2),c(166.6667,133.3333),1.2)
}