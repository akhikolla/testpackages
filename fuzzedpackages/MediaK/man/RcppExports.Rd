\name{RcppExports}

\alias{dis_value}
\alias{permute}

\title{Set of functions in example RcppEigen package}
\description{
  These two functions are created to calculate different values about MeDiA_k
}
\usage{

dis_value(iTest,jTest,select)
permute(iTest,jTest,times,selectvec)

}

\arguments{
  \item{iTest}{a numeric matrix}
  \item{jTest}{a numeric matrix which has the same rows as iTest}
  \item{select}{an integer  which means get k values you want to get }
  \item{selectvec}{an integer vector which contains a series of K values you want to get}
  \item{times}{an integer indicates permutation times}
}
\value{


  \code{dis_val(iTest,jTest,select)} returns a distance
  \code{permute(iTest,jTest,times,selectvec)} return a list contains mean and sd values based on the selectvec(K value you want to select)





}
\examples{

iTest <- matrix(rnorm(1000),100)
jTest<- matrix(rnorm(200),100)
select<- 5
times<- 100
selectvec <- c(1:10)
dis_value(iTest,jTest,select)
permute(iTest,jTest,times,selectvec)


}
