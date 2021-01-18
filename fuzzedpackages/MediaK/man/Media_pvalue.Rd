\name{Media_pvalue}

\alias{Media_pvalue}

\title{ functions to calculate pvalue}
\description{
  This function is created to calculate p values about MeDiA_k
}
\usage{

Media_pvalue(iTest,jTest,times,selectvec)

}

\arguments{
  \item{iTest}{a numeric matrix}
  \item{jTest}{a numeric matrix which has the same rows as iTest}
  \item{selectvec}{an integer vector which contains a series of K values you want to get}
  \item{times}{an integer indicates permutation times}
}
\value{

  \code{Media_pvalue(iTest,jTest,times,selectvec)} return a numeric vectors contains p values based on the selectvec(K value you want to select)


}
\examples{

iTest <- matrix(rnorm(1000),100)
jTest<- matrix(rnorm(200),100)
times<- 100
selectvec <- c(1:10)
Media_pvalue(iTest,jTest,times,selectvec)


}
