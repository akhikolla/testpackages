\encoding{UTF-8}
\name{RearrangeJoseph}
\alias{RearrangeJoseph}

\title{
Makes the initialization of auxiliary matrices and calls to sequence of functions to perform the reordering of the elements (leaves) of the hierarchical tree according to class labels of the data objects.
}
\description{
The function makes the preprocessing of the hierarchical merging matrix in order to simplify the hierarchical tree. For each node of the hierarchical tree forms two vectors, consisting of elements of its left and write subtrees. The function initialize the four auxiliary matrices, which are necessary for the main dynamic programming algorithm \code{\link{OrderingJoseph}} for reordering the leaves of the hierarchical tree according to available class labels. The function calls the main function code{\link{OrderingJoseph}} and than make use of the A, maxI and maxJ matrices to recive the optimal sequence of the leaves in the tree by calling the function \code{\link{funMerge}}. It outputs the vector of row indices of the merge matrix, which coressponds to the nodes with switched subtrees. The merge matrix with switched subtrees present the optimal reordering of the hierarchical tree.
}
\usage{
RearrangeJoseph(hc, dis,class,cpp)
}

\arguments{
 \item{hc}{
An object of class hclust which describes the tree produced by the clustering process.There are such components: merge, height, order, labels,call,method,dist.method.
}
\item{dis}{
an n by n distance matrix, each (i,j) element  corresponds to the Euclidean distance between ith and jth data object.
}
  \item{class}{
a vector of length n, which stores the class labels of the dataset objects. 
}
  \item{cpp}{a binary variable, which allows to switch between dynamic programming algorithm, realized in R language \code{\link{OrderingJoseph}} and the same algorithm as the C++ function \code{\link{OrderingJosephC}}
}
}

\details{
This function initializes four auxiliary matrices A, r , maxI and maxJ. Matrix A of dimension nxn stores the results of the calculation of the evaluation function values for each subtree of the hierarchical tree. Each cell (i,j) of the matrix A stores the best ordering (according to the evaluation function value) of the subtree starting with element i and ending with element j.
The algorithms use the dynamic programming concept and works, analyzing the hierarchical tree from the bottom to the top (root node). At the beginning all cells of the matrix A and r are equal 1, than the best ordering is calculated for the nodes, consisting of two elements. At each following iteration the best ordering of the more complex node is estimated on the basis of the best selected orderings received for its left and right subtrees. The process stops when the root node is reached. At each stage of this bottom up recursive process of computing the optimised function value for each inner node we updata the marices maxI,maxJ in order to store the intermediate leaves that provide the optimized value for each pair of (i,j) elements of this node. To evaluate the complex node not only the function values for its left and right subtrees are taken into account but also the coincidence of the class labels of the rightmost element of the left subtree and the leftmost element of the write subtree. When the class labels are equal the resulting value for the node is not simply the sum of the values of two its subtrees but is recalculated using the values of the auxiliary matrix r. Matrix r with dimensions nxn stores for each subtree the number of elements with the same class label starting from leftmost or the rightmost element of the subtree.

Each cell (i,j) of the matrix maxI stores the code of the intermediate element i1, that provides the best evaluation value A(i,j) for the subtree or node(i,j), which has element i as the leftmost and element j as the rightmost. Element i1 is the rightmost element of the left subtree of this node(i,j).
Each cell (i,j) of the matrix maxJ stores the code of the intermediate element j1, that provides the best evaluation value A(i,j) for the subtree or node(i,j), which has element i as the leftmost and element j as the rightmost. Element j1 is the leftmost element of the right subtree of this node(i,j).
After performing the reordering procedure \code{\link{OrderingJoseph}} using the dynamic programming concept the recovering of the sequence of leaves which gives the best evaluation function value is performed \code{\link{funMerge}}. At first we find out the leftmost and rightmost elements of the root node or the whole tree T. After that the sequence of leaves is recovered using a backtracking scheme with the help of matrices maxI and maxJ, which store the codes of the intermediate elements, providing the best value of the evaluation function for each subtree. The sequence corresponds to the optimized leaf ordering according to evaluation function.
As a result the optimized merge matrix is formed, with some nodes, having their subtrees switched.
}
\value{
\item{hcl}{A modified merging matrix of dimension nx2 which describes the tree produced optimal reordering of the tree leaves.
}
\item{A}{An auxiliary matrix A of dimension nxn, which stores the results of the calculation of the evaluation function values for each subtree of the hierarchical tree.
}
\item{max}{An maximal element of the matrix A.
}
}

\references{
Therese Biedl,Brona Brejova, Erik D,Demaine, Angele M.Hanmel and Tomas Vinar:Optimal Arrangement of Leaves in the Tree Representing Hierarchical Clustering of Gene Expression Data/Technical report 2001-14
}
\author{
Natalia Novoselova,Frank Klawonn,Junxi wang

Maintainer: Natalia Novoselova <novos65@mail.ru>
}

\seealso{
\code{\link{OrderingJoseph}}, \code{\link{OrderingJosephC}}, \code{\link{funMerge}}
}

\examples{
data(leukemia)
rownames(leukemia)=leukemia[,1]
leukemia=leukemia[,-1]
matr=leukemia[,-101]
class=leukemia[,101]
  
matr=as.matrix(matr)
dist=dist(matr)
hc <- hclust(dist)

res<-RearrangeJoseph(hc,as.matrix(dist),class,TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ reorder }
\keyword{ clustering }
\keyword{ backtracking }