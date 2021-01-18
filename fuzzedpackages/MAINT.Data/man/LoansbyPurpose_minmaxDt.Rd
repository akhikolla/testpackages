\name{LoansbyPurpose_minmaxDt}
\docType{data}
\alias{LoansbyPurpose_minmaxDt}
\title{Loans by purpose: minimum and maximum Data Set}

\description{This data set consist of the lower and upper bounds of the intervals for four interval characteristics of the loans aggregated by their purpose. The original microdata is available at the Kaggle Data Science platform and consists of 887 383 loan records characterized by 75 descriptors.
Among the large set of variables available, we focus on borrowers' income and account and loan information aggregated by the 14 loan purposes, wich are considered as the units of interest.}

\usage{data(LoansbyPurpose_minmaxDt)}

\format{A data frame containing 14 observations on the following 8 variables.
  \describe{
    \item{ln-inc_min}{The minimum, for the current loan purpose, of natural logarithm of the self-reported annual income provided by the borrower during registration.}
    \item{ln-inc_max}{The maximum, for the current loan purpose, of natural logarithm of the self-reported annual income provided by the borrower during registration.}
    \item{ln-revolbal_min}{The minimum, for the current loan purpose, of natural logarithm of the total credit revolving balance.}
    \item{ln-revolbal_max}{The maximum, for the current loan purpose, of natural logarithm of the total credit revolving balance.}
    \item{open-acc_min}{The minimum, for the current loan purpose, of the number of open credit lines in the borrower's credit file.}
    \item{open-acc_max}{The maximum, for the current loan purpose, of the number of open credit lines in the borrower's credit file.}
    \item{total-acc_min}{The minimum, for the current loan purpose, of the total number of credit lines currently in the borrower's credit file.}
    \item{total-acc_max}{The maximum, for the current loan purpose, of the total number of credit lines currently in the borrower's credit file.}
  }
}

\source{https:www.kaggle.com/wendykan/lending-club-loan-data}

\keyword{datasets}
\keyword{interval data}

