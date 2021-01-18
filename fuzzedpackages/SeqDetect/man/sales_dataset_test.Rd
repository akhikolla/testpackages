\name{sales_dataset_test}
\alias{sales_dataset_test}
\title{Sales time-series test}
\usage{
sales_dataset_test(learning_set = 1:20, testing_set = 21:40,
  th_increment = 1, max_th = NULL)
}
\arguments{
\item{learning_set}{(vector) - A set of products to learn ETTs in the Sequence Detector.}

\item{testing_set}{(vector) - A set of products to test previously learned sales numbers.}

\item{th_increment}{(integer) - A threshold increment between two tests.}

\item{max_th}{(integer) - Maximal thershold for testing. When reached, no further tests and no further threshold increment is done. If NULL, re-testing is done while there are some sequences detected.}
}
\value{
A list that comprises sequence statistics for all tests and thresholds.
}
\description{
Sales dataset taken from [2], which comprises 811 product one year sales quantities. We applied this dataset to test the Sequence Detector.
The results are available in [1]. The results of the test are various statistics on detected sequences.
The testing set of products is re-tested by simultaneously rising the projection threshold, until no more sequences could be detected or \emph{max_th} parameter is reached.
}
\references{
[1] D. Krleža, B. Vrdoljak, and M. Brčić, Latent Process Discovery using Evolving Tokenized Transducer, \emph{IEEE Access}, vol. 7, pp. 169657 - 169676, Dec. 2019 \cr
[2] S. C. Tan and J. P. San Lau, Time series clustering: A superior alternative for market basket analysis, in Proceedings of the First International Conference on Advanced Data and Information Engineering (DaEng-2013), Singapore, 2014, pp. 241–248.
}
