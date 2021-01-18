\name{reg_2d}
\Rdversion{1.1}
\alias{reg_2d}
\title{
    Isotonic Regression on 2D input.
}
\description{
    Isotonic regression on weighted or unweighted 2D input with L1, L2 metric and other options.
}
\usage{
    reg_2d(y_vec, w_vec, metric)
}
\arguments{
    \item{y_vec}{
        The 2D NumericMatrix of input data that we use to regression. It must 
        be the same size as the w_vec argument.
    }
    \item{w_vec}{
        The 2D NumericMatrix of the weight of the input data. The default 
        value is 1 for every entry. It must be the same size as y_vec.
    }
    \item{metric}{
        This is an integer input, metric = 1 stands for using L1 metric,
        metric = 2 stands for using L2 metric
    }
}
\section{Error Messages}{
    \itemize{
        \item The size of y_vec is 0: Empty data.
        \item The rows of w_vec doesn't match the rows of y_vec: Data and weight have different number of rows
        \item The columns of w_vec doesn't match the rows of y_vec: Data and weight have different number of columns
        \item The entry of w_vec has negative value: Negative weight detected
        \item Metric input is not in {1,2,3}: Metric does not exist
    }
}
\details{
    See the paper about 2D regression in the reference.
}
\value{
    A 2D NumericMatrix of the regression result which has the same size of y_vec.
}
\references{
    Q.F. Stout, Isotonic median regression via partitioning, Algorithmica 66 (2013), pp. 93-112 doi.org/10.1007/s00453-012-9628-4
}
\author{Zhipeng Xu, Chenkai Sun, Aman Karunakaran, Quentin Stout
    \email{xzhipeng@umich.edu}
    \url{https://github.com/xzp1995/UniIsoRegression}
}
\examples{
    library(UniIsoRegression)
    #===2d monotonic===
    y=matrix(c(2, 4, 3, 1, 5, 7,9,0), nrow=2, ncol=4, byrow = TRUE)
    weight=matrix(c(1, 10, 3, 9, 5, 7,9,10), nrow=2, ncol=4, byrow = TRUE)

    #l_1 metric
    temp=UniIsoRegression::reg_2d(y, weight, metric = 1)
    print(temp)

    #l_2 metric
    temp=UniIsoRegression::reg_2d(y, weight, metric = 2)
    print(temp)
}








