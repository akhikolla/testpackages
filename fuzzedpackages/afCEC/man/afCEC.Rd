\name{afCEC}
\alias{afCEC}

\title{
afCEC
}

\description{
Performs the Active Function Cross-Entropy Clustering on the data set.
}

\usage{
afCEC (
    points,
    maxClusters,
    initialLabels="k-means++",
    cardMin=0.01,
    costThreshold=-0.000001,
    minIterations=10,
    maxIterations=100,
    numberOfStarts=10,
    method="Hartigan",
    values="quadratic",
    interactive=FALSE
)
}

\arguments{

\item{points}{
A (#points x n) matrix of data containing n-dimensional points stored row-by-row.
}

\item{maxClusters}{
Indicates the maximal number of clusters, that afCEC can partition the data into.
}

\item{initialLabels}{
Initial labelling determining the data membership to the particular clusters. There are 3 options allowed:
\itemize{
\item{"random" - causes the labelling to be randomly generated.}
\item{"k-means++" - causes the labelling to be generated using the k-means++ heuristics.}
\item{A (#points x numberOfStarts) matrix of type integer containing values in the range: 0...maxClusters - 1. The value x in the row i and column j indicates that in the j-th start, the i-th point will be initially assigned to the (x + 1)-th cluster.}
}
The defalt value is "k-means++".
}

\item{cardMin}{
Value so that if the number of points in particular cluster relatively to the size of the whole data drops bellow, then that cluster gets removed and it's data is assigned to the other clusters. The default value is 0.01.
}

\item{costThreshold}{
Negative value so that if the difference between the calculated cost in the subsequent iterations is greater than it then the afCEC terminates returning the solution from the terminated iteration as the final solution (in the present start), provided that at least minIterations already passed. The default value is -0.000001;
}

\item{minIterations}{
Value indicating the minimal number of iterations needed before afCEC can terminate, provided that no error occurred. The default value is 10.
}

\item{maxIterations}{
Value indicating the maximal number of iterations in one start, that afCEC cannot exceed. The default value is 100.
}

\item{numberOfStarts}{
Value indicating the number of runs of the algorithm to be performed. The best solution is chosen out of the solutions given in the particular steps, so that it has the lowest value of the cost function among them. The default value is 10.
}

\item{method}{
Heuristics used to perform the clustering. There are 2 options allowed:
\itemize{
\item{"Lloyd" - indicates that the Lloyd heuristics will be used to perform the clustering.}
\item{"Hartigan" - indicates that the Hartigan heuristics will be used to perform the clustering.}
}
The default value is "Hartigan".
}

\item{values}{
Definition of the active function family used to perform the clustering. There are 2 options allowed:
\itemize{
\item{"quadratic" - indicates that the function family of the form:
\cr
\cr
f(x_1,...,x_(n-1))=a_1*x_1^2+...+a_(n-1)*x_(n-1)^2+a_n*x_1+...+a_(2*n-2)*x_(n-1)+a_(2*n-1)*1,
\cr
\cr
where n is the dimensionality of the data and a_i are the optimal coefficients determined by afCEC using the linear least square fitting, will be used to perform the clustering.
}
\item{A ((m*n) x #points) matrix containing the values of the particular components of the active function on the data set, placed according to the following layout (row by row):
\cr
\cr
First m rows:
\cr
\cr
[f_1(x_1_2,...,x_1_n), f_1(x_2_2,...,x_2_n), ..., f_1(x_#points_2,...,x_#points_n)]\cr
[f_2(x_1_2,...,x_1_n), f_2(x_2_2,...,x_2_n), ..., f_2(x_#points_2,...,x_#points_n)]\cr
...\cr
[f_m(x_1_2,...,x_1_n), f_m(x_2_2,...,x_2_n), ..., f_m(x_#points_2,...,x_#points_n)]\cr
\cr
Second m rows:
\cr
\cr
[f_1(x_1_1,x_1_3,...,x_1_n), f_1(x_2_1,x_2_3,...,x_2_n), ...]\cr
[f_2(x_1_1,x_1_3,...,x_1_n), f_2(x_2_1,x_2_3,...,x_2_n), ...]\cr
...\cr
[f_m(x_1_1,x_1_3,...,x_1_n), f_m(x_2_1,x_2_3,...,x_2_n), ...]\cr
\cr
Last m rows:
\cr
\cr
[f_1(x_1_1,...,x_1_(n-1)), f_1(x_2_1,...,x_2_(n-1)), ...]\cr
[f_2(x_1_1,...,x_1_(n-1)), f_2(x_2_1,...,x_2_(n-1)), ...]\cr
...\cr
[f_m(x_1_1,...,x_1_(n-1)), f_m(x_2_1,...,x_2_(n-1)), ...],\cr
\cr
where: n - dimensionality of the data, x_i_j - j-th coordinate of the i-th point of data, m - number of components of the active function.
\cr
\cr
In the foregoing case, the active function family consists of the functions of the form:
\cr
\cr
f(x_1,...,x_(n-1))=a_1*f_1(x_1,...,x_(n-1))+...+a_m*f_m(x_1,...,x_(n-1)),
\cr
\cr
where n is the dimensionality of the data and a_i are the optimal coefficients determined by afCEC using the linear least square fitting.
}
}
The default value is "quadratic".
\cr
\cr
\strong{Remarks:}
\cr
\cr
In the case of the second way of defining the function, the returned object doesn't contain the information about the means coordinates corresponding to the active direction. Moreover, the plotting capabilities of the afCEC package are severely impaired in that case. See \code{\link{plot}}.
}

\item{interactive}{
Indicates if the algorithm runs in the "interactive" mode. The "interactive" mode allows to track the intermediate steps of the afCEC. Instead of one object of the afCEC class, the whole list of such objects is returned, where each item of the list corresponds to the intermediate step of the algorithm. See \code{value} section for more details. The default value is FALSE.
}

}
\value{
\itemize{
\item{Empty list - if clustering failed.}
\item{Object of class afCEC containing the parameters of the best fitted mixed generalized multivariate normal distribution - if argument interactive=FALSE and clustering succeded.}
\item{List of k objects of class afCEC containing the parameters of the best fitted mixed generalized multivariate normal distribution across the k subsequent steps of the algorithm - if argument interactive=TRUE and clustering succeeded.}
}
}

\seealso{
\code{\link{plot}}
}

\examples{
	\dontshow{
		library(afCEC);
		data(dog);
		res<-afCEC(dog,5,numberOfStarts=1);
	}
	\dontrun{
		# The following two examples demonstrate two equivalent ways of passing the same
		# 3D quadratic active function to the afCEC routine:
		
		# 1) Using "quadratic" value (default):
		
		library(afCEC);
		data(airplane);
		result <- afCEC(airplane, 17);
		
		# what is equivalent to:
		
		library(afCEC);
		data(airplane);
		result <- afCEC(airplane, 17, values="quadratic");
		
		# 2) Using the matrix containing the explicit values of the active function
		#    across the subsequent dimensions:
		
		library(afCEC);
		data(airplane);
		values = matrix(rep(0, 5*3*dim(airplane)[1]), 5*3, dim(airplane)[1]);
		for (i in 1:dim(airplane)[1]) {
			tmp <- airplane[i,2:3];
			for (j in 1:3) {
				values[((j - 1) * 5) + 1, i] <- tmp[1]^2;
				values[((j - 1) * 5) + 2, i] <- tmp[2]^2;
				values[((j - 1) * 5) + 3, i] <- tmp[1];
				values[((j - 1) * 5) + 4, i] <- tmp[2];
				values[((j - 1) * 5) + 5, i] <- 1;
				if (j < 3) tmp[j] <- airplane[i,j];
			}
		}
		result <- afCEC(airplane, 17, values=values);
	}
}