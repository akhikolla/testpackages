\name{multiplication_speed}
\alias{multiplication_speed}
\title{Matrix multiplication acceleration algorithm.}
\description{
Matrix multiplication acceleration algorithm.}
\usage{
multiplication_speed(A,B)
}
\arguments{
      \item{A}{matrix A.}
	\item{B}{matrix B.}
	}
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, and Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
\dontrun{
A<-matrix(1:10,2,5)
B<-matrix(1:10,5:2)
result<-multiplication_speed(A,B)
}
}
