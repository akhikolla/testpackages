\name{funMerge}
\alias{funMerge}

\title{
Recover the optimal sequence of leaves in the hierarchical tree according to available class labels.
}
\description{
This is the auxiliary function of the package. It recovers the optimal sequence of leaves in the tree. Recovering is made using a backtracking scheme with the help of matrices maxI and maxJ, which store the codes of the intermediate elements, providing the best value of the evaluation function for each subtree. The sequence corresponds to the optimised leaf ordering according to evaluation function. During the recovering process for some nodes the coressponding subtrees are switched and the modified merge matrix is formed.
It returns a \dQuote{matrix} object, which presents the modified merge matrix with some nodes with switched subtrees. The function output is used in the function \code{\link{RearrangeJoseph}}. 
}
\usage{
funMerge(ind,row,col,hc,node,maxI,maxJ,cpp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ind}{
the value of the last row of the merge matrix, which is equal to n-1, where n is the number of data objects.
}
  \item{row}{
the leftmost element of the root tree T. It is defined as a row of the maximal element of the matrix A, holding the values of evaluation function for each subtree of the hierarchical tree.
}
  \item{col}{
the rightmost element of the root tree T. It is defined as a column of the maximal element of the matrix A, holding the values of evaluation function for each subtree of the hierarchical tree.
}
  \item{hc}{
An object of class hclust which describes the tree produced by the clustering process.There are such components: merge, height, order, labels,call,method,dist.method.
}
  \item{node}{
a list of lists of length n-1 \code{\link{testBar}}, where each single list stores two vectors, consisting of elements of the left and right subtrees of the corresponding node in the merging matrix \code{hc}
}
  \item{maxI}{
a nxn numeric matrix. Each cell (i,j) of the matrix maxI stores the code of the intermediate element i1, that provides the best value A(i,j) of the evaluation function for the subtree or node(i,j), which has element i as the leftmost and element j as the rightmost. Element i1 is the rightmost element of the left subtree of this node(i,j).
}
  \item{maxJ}{
a nxn numeric matrix. Each cell (i,j) of the matrix maxJ stores the code of the intermediate element j1, that provides the best value A(i,j) of the evaluation function for the subtree or node(i,j), which has element i as the leftmost and element j as the rightmost. Element j1 is the leftmost element of the right subtree of this node(i,j).
}
 \item{cpp}{a binary variable, which allows to switch between dynamic programming algorithm, realized in R language \code{\link{OrderingJoseph}} and the same algorithm as the C++ function \code{\link{OrderingJosephC}}
}
}
\details{
the function recovers the optimal sequence of leaves of the hierarchicla tree. The leftmost and rightmost elements of the root tree are used as the inputs function parameters. Than with the help of the matrices maxI, maxJ, the sequence of leaves which gives the best function value and corresponds to the optimised leaf ordering can be recovered using a backtracking scheme.
}
\value{
\item{hc}{modified (n-1)x2 merge matrix with some nodes switched in the process of backtracking}
}
\references{
Therese Biedl,Brona Brejova, Erik D,Demaine, Angele M.Hanmel and Tomas Vinar:Optimal Arrangement of Leaves in the Tree Representing Hierarchical Clustering of Gene Expression Data/Technical report 2001-14
}
\author{
Natalia Novoselova,Frank Klawonn,Junxi wang

Maintainer: Natalia Novoselova <novos65@mail.ru>
}

\seealso{
\code{\link{RearrangeJoseph}}, \code{\link{OrderingJoseph}}
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

res=OrderingJosephC(ind-1, hc$merge, node, A, r, maxI, maxJ, class, coef)
A=res$A
maxI=res$maxI
maxJ=res$maxJ
r=res$r
	
col=which.max(apply(A[node[[ind]]$left,node[[ind]]$right,drop=FALSE],2,max))
row=which.max(apply(A[node[[ind]]$left,node[[ind]]$right,drop=FALSE],1,max))
col=node[[ind]]$right[col]
row=node[[ind]]$left[row]

hcl=funMerge(ind,row,col,hO,node,maxI,maxJ,TRUE)
}
\keyword{ recovering }
\keyword{ optimal }
