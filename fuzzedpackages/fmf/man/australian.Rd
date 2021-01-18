\name{australian}
\alias{australian}
\docType{data}
\title{
Australian Credit Approval
}
\description{
This is the famous Australian Credit Approval dataset, originating from the StatLog project. It concerns credit card applications. 
All attribute names and values have been changed to meaningless symbols to protect the confidentiality of the data.

}
\usage{data(australian)}
\format{
  A data frame with 690 Instances and 15 attributes (including the class attribute, "Class")
}
\details{
There are 6 numerical and 8 categorical attributes, all normalized to [-1,1]. The original formatting was as follows:
A1: A,B class attribute (formerly: +,-)
A2: 0,1 CATEGORICAL (formerly: a,b)
A3: continuous.
A4: continuous.
A5: 1,2,3 CATEGORICAL (formerly: p,g,gg)
A6: 1, 2,3,4,5, 6,7,8,9,10,11,12,13,14 CATEGORICAL (formerly: ff,d,i,k,j,aa,m,c,w, e, q, r,cc, x)
A7: 1, 2,3, 4,5,6,7,8,9 CATEGORICAL (formerly: ff,dd,j,bb,v,n,o,h,z)
A8: continuous.
A9: 1, 0 CATEGORICAL (formerly: t, f)
A10: 1, 0 CATEGORICAL (formerly: t, f)
A11: continuous.
A12: 1, 0 CATEGORICAL (formerly t, f)
A13: 1, 2, 3 CATEGORICAL (formerly: s, g, p)
A14: continuous.
A15: continuous.

}
\source{
Confidential. Donated by Ross Quinlan
}
\references{
[LibSVM] (https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html), UCI - 1987 
}
\examples{

data(australian)

X = australian[, -1]

y = australian[, 1]
}
\keyword{datasets}
