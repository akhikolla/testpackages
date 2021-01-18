\name{mvlrrl0}
\alias{mvlrrl0}
\title{Multi-view bi-clustering via L0-norm enforced sparse LRR}
\usage{mvlrrl0(datasets, svs, sz, seed=0, maxIter=1000, thres=0.00001, logLvl=0)}
\description{
	Identify consistent sample cluster among all views and simultaneously associated feature clusters per view. Clusters are obtained via finding sparse low rank representation (LRR) of input data matrices, where the sparsity is enforced using L0-norm. One sample cluster and its associated feature clusters are identified and returned each time this function is used. If multiple clusters are desired, call this function repeatedly with samples left unclustered.
}
\details{
	This method identifies clusters via finding multi-view sparse low rank representations (LRR) of input data matrices. The LRR is obtained via matrix decomposition in the form: (zu)v, where both z and u are column vectors and their length is equal to the number of samples, (zu) represents the element-wise product between z and u, and v is a column vector and its length is equal to the number of features in the view. Vector z is a multiplier shared across all views, while both u and v are view specific. Sample cluster is read off directly from z by assigning samples with non-zero component in z to the cluster. Because z is shared across all views, the obtained sample cluster is a common cluster among all views. Feature cluster is obtained by assigning features with non-zero component in v to the cluster. Because v is view specific, feature cluster is view specific as well. Sample cluster and feature cluster are associated in the sense that they help determine each other. The sparsity of vector z and v is enforced using L0-norm.
}
\arguments{
	\item{datasets}{
		List of input data, where each element is a matrix. For all the matrices, rows represent samples and columns represent features. These matrices are the characterization of a same set of samples from different perspectives (views), one matrix per view. So all the matrices should have the exact same rows, but can have different columns. The rows in the same position in all matrices represent a same sample.
	}
	\item{svs}{
		A numerical vector with length that equals to the number of view, which controls the sparsity of vector v in the decomposition of each view. It specifies the maximum number of non-zero elements in vector v.
	}
	\item{sz}{
		A number, which controls the sparsity of vector z in the decomposition of all views. It specifies the maximum number of non-zero elements in vector z.
	}
	\item{seed}{
		(Optional) A number, which gives the feature (the index) in the first view that will be used as a seed during initialization of vector v in the decomposition. Vector v is initialized with zero in all positions except the one corresponding to the given feature, which is initialized with one. If a seed is not provided, v is initialized with one in all positions. Because the optimization problem this function solves is non-convex, the returned solution is a local optimizer and could vary while different initializations are used. Based on the observation in our experiments, the seed feature is included in the returned feature cluster in most cases. So if there is a feature needed to be included in the feature cluster, it can be set as the seed feature using this argument.
	}
	\item{maxIter}{
		(Optional) Maximum number of loop iterations, which is one of the criteria that controls when to terminate the iterative process for searching for the optimal solution. Please read the referenced paper for details. The default value is 1000.
	}
	\item{thres}{
		(Optional) The other criteria (besides above 'maxIter' argument) for terminating the optimizing loop. When the difference between the two z from the two most recent consecutive iterations passes (is smaller than) this threshold, the loop is terminated. The default value is 0.00001.
	}
	\item{logLvl}{
		(Optional) Logging level, which can be set to either 0 or 1, controls the amount of printing, where a larger value means more printing.  Default value is 0, which turns off the logging. 
	}
}
\value{
	A list with following named fields:
	\item{Cluster}{A binary vector (with 1 or 0 entries) of length equal to the sample size. It indicates whether a sample is in the identified cluster (with 1 in its corresponding position) or not (with 0).  }
	\item{FeatClusters}{A list of binary vectors that give identified feature cluster for each view. Value 1 indicates the corresponding feature belongs to the identified feature cluster.}
	\item{U}{A matrix of with size of n x m, where n is the number of samples, m is the number of views. So each column gives the u in the decomposition of the corresponding view.}
	\item{V}{A list of vectors that give the v in the decomposition of each view. }
	\item{z}{The common multiplier shared across all views in their decomposition.}
}
\references{
	Jiangwen Sun, Jin Lu, Tingyang Xu, and Jinbo Bi \emph{Multi-view Sparse Co-clustering via Proximal Alternating Linearized Minimization}
	Proceedings of the 32nd International Conference on Machine Learning (ICML), pp. 757-766
}
\examples{
		library(mvcluster)
		data(phe)
		data(gen)
        views <- list(phe,gen)
        result <- mvlrrl0(views,c(3,10),320)
}