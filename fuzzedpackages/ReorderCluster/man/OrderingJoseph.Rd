\name{OrderingJoseph}
\alias{OrderingJoseph}

\title{
Makes the calculation of the evaluation function for each subtree of the hierarchical tree using the dynamic programming approach
}
\description{
The function realizes the dynamic programming approach in order to find the optimal reordering of the leaves of the hierarchical tree. The optimal reordering is made according to the available class labels. As an output the four auxiliary matrices A, R, maxI, maxJ are returned to the \code{\link{RearrangeJoseph}} function. The cells of the A matrix store the values of the evaluation function for each subtree of the hierarchical tree.
}
\usage{
OrderingJoseph(ind, hc, node, A, r, maxI, maxJ, class, coef)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ind}{
the value of the last row of the merge matrix, which is equal to n-1, where n is the number of data objects.
}
  \item{hc}{
An object of class hclust which describes the tree produced by the clustering process.There are such components: merge, height, order, labels,call,method,dist.method.
}
  \item{node}{
a list of lists of length n-1 \code{\link{testBar}}, where each single list stores two vectors, consisting of elements of the left and right subtrees of the corresponding node in the merging matrix \code{hcmatr}
}
  \item{A}{
a numerical nxn matrix, which stores the values of the evaluation function fir the subtrees. Each cell (i,j) of the matrix A stores the best ordering (according to the evaluation function value) of the subtree starting with element i and ending with element j.
}
  \item{r}{
a numerical nxn matrix, each cell stores the number of the elements with the same class label, starting from the leftmost or the rightmost element of the subtree.
}
  \item{maxI}{
a nxn numeric matrix. Each cell (i,j) of the matrix maxI stores the code of the intermediate element i1, that provides the best value A(i,j) of the evaluation function for the subtree or node(i,j), which has element i as the leftmost and element j as the rightmost. Element i1 is the rightmost element of the left subtree of this node(i,j).
}
  \item{maxJ}{
a nxn numeric matrix. Each cell (i,j) of the matrix maxJ stores the code of the intermediate element j1, that provides the best value A(i,j) of the evaluation function for the subtree or node(i,j), which has element i as the leftmost and element j as the rightmost. Element j1 is the leftmost element of the right subtree of this node(i,j).
}
  \item{class}{
a numerical or character vector, containing the class labels of the dataset objects.
}
  \item{coef}{
a parameter used in the expression for the calculation of the evaluation function. It reinforces the influence of longer sequences. Its values are in the interval ]1,2].
}
}
\details{
The function performs the main part of the reordering process, forming the matrix A of values of the evaluation function for each subtree of the initially formed hierarchical tree.
At the beginning all cells of the matrix A and r are equal 1, than the best ordering is calculated for the nodes, consisting of two elements. At each following iteration the best ordering of the more complex node is estimated on the basis of the best selected orderings received for its left and right subtrees. The process stops when the root node is reached. At each stage of this bottom up recursive process of computing the optimised function value for each inner node we updata the marices maxI,maxJ in order to store the intermediate leaves that provide the optimized value for each pair of (i,j) elements of this node. To evaluate the complex node not only the function values for its left and right subtrees are taken into account but also the coincidence of the class labels of the rightmost element of the left subtree and the leftmost element of the write subtree. When the class labels are equal the resulting value for the node is not simply the sum of the values of two its subtrees but is recalculated using the values of the auxiliary matrix r. Matrix r with dimensions nxn stores for each subtree the number of elements with the same class label and are not interrupted by elements of other class labels starting from leftmost or the rightmost element of the subtree.
}
\value{
\item{A}{
a numerical nxn matrix, which stores the values of the evaluation function fir the subtrees. Each cell (i,j) of the matrix A stores the best ordering (according to the evaluation function value) of the subtree starting with element i and ending with element j.
}
\item{maxI}{
a nxn numeric matrix. Each cell (i,j) of the matrix maxI stores the code of the intermediate element i1, that provides the best value A(i,j) of the evaluation function for the subtree or node(i,j), which has element i as the leftmost and element j as the rightmost. Element i1 is the rightmost element of the left subtree of this node(i,j).
}
  \item{maxJ}{
a nxn numeric matrix. Each cell (i,j) of the matrix maxJ stores the code of the intermediate element j1, that provides the best value A(i,j) of the evaluation function for the subtree or node(i,j), which has element i as the leftmost and element j as the rightmost. Element j1 is the leftmost element of the right subtree of this node(i,j).
}
\item{r}{
a numerical nxn matrix, each cell stores the number of the elements with the same class label, starting from the leftmost or the rightmost element of the subtree.
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
\code{\link{RearrangeJoseph}}, \code{\link{funMerge}}
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

coef=1.5
num=dim(hc$merge)[1]
A=array(1,c(num+1,num+1))
r=array(1,c(num+1,num+1))

maxI=array(0,c(num+1,num+1))
maxJ=array(0,c(num+1,num+1))
ind=num
  
node=testBar(hc)
flag=CalcMerge(hc,node,class)
## change matrix hc$merge
hO<-hc
nodeO<-node
out=SubTree(ind,flag,node,hc,A,r,coef)
hc=out$hc
node=out$node
A=out$A
r=out$r

res=OrderingJoseph(ind, hc, node, A, r, maxI, maxJ, class, coef)
}

\keyword{ reordering }
\keyword{ clustering }
\keyword{tree}

