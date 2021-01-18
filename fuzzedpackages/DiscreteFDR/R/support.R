#'@title Computing discrete p-values and their support for binomial and Fisher's exact tests
#'
#'@description
#'Computes discrete raw p-values and their support
#'for binomial test or Fisher's exact test applied to 2 x 2 contingency tables
#'summarizing counts coming from two categorical measurements.
#'
#'@details
#'Assume that each contingency tables compares 2 variables and resumes the counts of association
#'or not with a condition. This can be resumed in the following table:
#'\tabular{lccc}{
#'\tab Association \tab No association \tab Total \cr
#'Variable 1 \tab X1 \tab Y1 \tab N1 \cr
#'Variable 2 \tab X2 \tab Y2 \tab N2 \cr
#'Total \tab X1 + X2 \tab Y1 + Y2 \tab N1 + N2
#'}
#'If \code{input="noassoc"}, \code{counts} has 4 columns which respectively contain X1, Y1, X2 and Y2.
#'If \code{input="marginal"}, \code{counts} has 4 columns which respectively contain X1, N1, X2 and N2.
#'
#'If \code{input="HG2011"}, we are in the situation of the \code{\link{amnesia}} data set as in
#'Heller & Gur (2011, see References). Each contingency table is obtained from one variable which is compared 
#'to all other variables of the study. That is, counts for "second variable" are replaced by the sum of the counts
#'of the other variables:
#'\tabular{lccc}{
#'\tab Association \tab No association \tab Total \cr
#'Variable j \tab Xj \tab Yj \tab Nj \cr
#'Variables !=j \tab SUM(Xi) - Xj  \tab SUM(Yi) - Yj \tab SUM(Ni) - Nj \cr
#'Total \tab SUM(Xi) \tab SUM(Yi)  \tab SUM(Ni)
#'}
#'Hence \code{counts} needs to have only 2 columns which respectively contain Xj and Yj.
#'
#'\code{binomial.pvalues.support} and \code{fisher.pvalues.support} are wrapper functions for \code{pvalues.support},
#'setting \code{test.type = "binomial"} and \code{test.type = "fisher"}, respectively.
#'
#'The code for the computation of the p-values of Fisher's
#'exact test is inspired by the example in the help page of
#'\code{\link[discreteMTP]{p.discrete.adjust}}.
#'
#'See the Wikipedia article about Fisher's
#'exact test, paragraph Example, for a good depiction
#'of what the code does for each possible value
#'of \code{alternative}.
#'
#'The binomial test simply tests for p = 0.5 by using X1
#'as the test statistic and N1 as the number of trials.
#'
#'This version: 2019-11-15.
#'
#'@seealso
#'\code{\link[discreteMTP]{p.discrete.adjust}}, \code{\link{fisher.test}}
#'
#'@param counts        a data frame of 2 or 4 columns and any number of lines,
#'                     each line representing a 2 x 2 contingency table to
#'                     test. The number of columns and what they must contain
#'                     depend on the value of the \code{input} argument, see
#'                     Details.
#'@param alternative   same argument as in \code{\link{fisher.test}}. The three
#'                     possible values are \code{"greater"} (default),
#'                     \code{"two.sided"} or \code{"less"} and you can specify
#'                     just the initial letter.
#'@param input         the format of the input data frame, see Details. The
#'                     three possible values are \code{"noassoc"} (default),
#'                     \code{"marginal"} or \code{"HG2011"} and you can specify
#'                     just the initial letter.
#'
#'@template example
#'@template exampleHG
#'
#'@return
#'A list of two elements:
#'\item{raw}{raw discrete p-values.}
#'\item{support}{a list of the supports of the CDFs of the p-values.
#'Each support is represented by a vector in increasing order.}
#'
#' @section References:
#' R. Heller and H. Gur (2011). False discovery rate controlling procedures for discrete tests. arXiv preprint arXiv:1112.4627v2 \href{https://arxiv.org/abs/1112.4627v2}{link}.
#' 
#' "Fisher's exact test", Wikipedia, The Free Encyclopedia,
#' accessed 2018-03-20,
#' \href{https://en.wikipedia.org/w/index.php?title=Fisher\%27s_exact_test&oldid=823327889}{link}.
#'
#'@importFrom stats dhyper phyper pbinom
#'@export
fisher.pvalues.support <- function(counts, alternative = "greater", input = "noassoc"){
  input <- match.arg(input, c("noassoc", "marginal", "HG2011"))
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  
  number.items <- nrow(counts)
  switch (input, 
          noassoc = {
            A11 <- counts[, 1]
            A12 <- counts[, 2]
            n <- A11 + A12
            #A21 <- counts[, 3]
            #A22 <- counts[, 4]
            A1. <- A11 + counts[, 3]
            A2. <- A12 + counts[, 4]
            #N <- A1. + A2.
            # Entry j in each of the above vectors is a count for the j-th test of no association
            # between two variables and a condition :  
            #                   Association    No association   Marginal counts
            #  Variable 1       A11[j]         A12[j]           n[j]
            #  Variable 2       A21[j]         A22[j]           N[j]-n[j]
            #  All variables    A1.[j]         A2.[j]           N[j]
          },
          marginal = {
            A11 <- counts[, 1]
            n <- counts[, 2]
            A12 <- n - A11
            #A21 <- counts[, 3]
            #A22 <- counts[, 4] - A21
            A1. <- A11 + counts[, 3]
            A2. <- counts[, 4] + n - A1.
            #N <- A1. + A2.
            # Entry j in each of the above vectors is a count for the j-th test of no association
            # between two variables and a condition :  
            #                   Association    No association   Marginal counts
            #  Variable 1       A11[j]         A12[j]           n[j]
            #  Variable 2       A21[j]         A22[j]           N[j]-n[j]
            #  All variables    A1.[j]         A2.[j]           N[j]
          },
          HG2011 = {
            A11 <- counts[, 1]
            #A21 <- sum(counts[, 1]) - A11
            A12 <- counts[, 2]
            #A22 <- sum(counts[, 2]) - A12
            A1. <- rep(sum(counts[, 1]), number.items)
            A2. <- rep(sum(counts[, 2]), number.items)
            n <- A11 + A12
            #N <- A1. + A2.
            # Entry j in each of the above vectors is a count for the test of no association
            # between variable j and the condition, compared to all other variables :  
            #                    Association    No association   Marginal counts
            #  Variable j        A11[j]         A12[j]           n[j]
            #  Other variables   A21[j]         A22[j]           N[j]-n[j]
            #  All variables     A1.[j]         A2.[j]           N[j]
          }
  )
  
  pCDFlist <- vector("list", number.items)
  raw.pvalues <- numeric(number.items)
  
  k <- pmin(n, A1.)
  l <- pmax(0, n - A2.)
  switch(alternative, greater =
           for (i in 1:number.items){
             x <- l[i]:k[i]
             # the "-1" below is because lower.tail = FALSE computes P[X > x],
             # and we want P[X >= x]=P[X > x-1]
			       # pmin/pmax below is to account for machine rounding issues
             pCDFlist[[i]] <- pmax(0, pmin(1, phyper(x-1, A1.[i], A2.[i], n[i], lower.tail = FALSE)))
             # the "+1" below is because vectors start with index 1 in R: x[A11[i]+1]=A11[i]
             raw.pvalues[i] <- pCDFlist[[i]][A11[i] + 1]
             # we want to have pCDFlist[[i]] in increasing order:
             pCDFlist[[i]] <- rev(pCDFlist[[i]])
           },
         less =
           for (i in 1:number.items){
             x <- l[i]:k[i]
			       # pmin/pmax below is to account for machine rounding issues
             pCDFlist[[i]] <- pmax(0, pmin(1, phyper(x, A1.[i], A2.[i], n[i], lower.tail = TRUE)))
             # the "+1" below is because vectors start with index 1 in R: x[A11[i]+1]=A11[i]
             raw.pvalues[i] <- pCDFlist[[i]][A11[i]+1]
           },
         two.sided =
           for (i in 1:number.items){
             x <- l[i]:k[i]
             atoms <- dhyper(x, A1.[i], A2.[i], n[i])
             # ensure that probabilities sum up to 1 (needs to be done multiple times sometimes)
             s <- sum(atoms)
             while(s != 1){
               atoms <- atoms/s
               s <- sum(atoms)
             }
			       # pmin/pmax below is to account for machine rounding issues
             pCDFlist[[i]] <- pmax(0, pmin(1, sapply(x, function(nu){sum(atoms[which(atoms <= atoms[nu + 1])])})))
             # the "+1" above and below is because vectors start with index 1 in R: x[A11[i]+1]=A11[i]
             raw.pvalues[i] <- pCDFlist[[i]][A11[i] + 1]
             # we want to have pCDFlist[[i]] in increasing order:
             pCDFlist[[i]] <- sort(pCDFlist[[i]])
           }
  )
  
  return(list(raw = raw.pvalues, support = pCDFlist))
}