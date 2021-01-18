\name{reg_1d}
\Rdversion{1.1}
\alias{reg_1d}
\title{
    Isotonic and Unimodal Regression on 1D input.
}
\description{
    Isotonic and unimodal regression on weighted or unweighted 1D input with L1, L2 and Linf metric and other options.
}
\usage{
    reg_1d(y_vec, w_vec, metric, unimodal = FALSE, decreasing = FALSE)
}
\arguments{
    \item{y_vec}{
        The vector of input data that we use to regression. It must be the 
        same size as the w_vec argument.
    }
    \item{w_vec}{
        The vector of the weight of the input data. The default value is 1 for
        every entry. It must be the same size as y_vec. It's only avaliable 
        for L1 and L2.
    }
    \item{metric}{
        This is an integer input, metric = 1 stands for using L1 metric,
        metric = 2 stands for using L2 metric, metric = 3 stands for
        using Linf metric.
    }
    \item{unimodal}{
        This is a boolean input, unimodal = false or 0 stands for isotonic
        regression and unimodal = true or 1 stands for unimodal regression
    }
    \item{decreasing}{
        This is a boolean input, decreasing = false or 0 stands for increasing
        model and decreasing = true or 1 stands for decreasing model. 
    }
}
\section{Error Messages}{
    \itemize{
        \item The size of y_vec is 0: Empty data.
        \item The size of w_vec doesn't match the size of y_vec: Data and weight have different number of entries
        \item The entry of w_vec has negative value: Negative weight detected
        \item Metric input is not in {1,2,3}: Metric does not exist
    }
}
\details{
    See the paper about unimodal regression via prefix isotonic regression in the reference.
}
\value{
    A vector of the regression result which has the same size of y_vec.
}
\references{
    Quentin F.Stout; Unimodal Regression via Prefix Isotonic Regression Computational Statistics and Data Analysis 53 (2008), pp. 289-297;
    Spouge, J., Wan, H. & Wilbur, W. Journal of Optimization Theory and Applications (2003) 117: 585-605 doi.org/10.1023/A:1023901806339
}
\author{Zhipeng Xu, Chenkai Sun, Aman Karunakaran, Quentin Stout
    \email{xzhipeng@umich.edu}
    \url{https://github.com/xzp1995/UniIsoRegression}
}
\examples{
    library(UniIsoRegression)

    #===1d monotonic===
    y=c(1,3,6,7,-1)
    weight=c(1,3,4,9,10)

    #l_1 metric decreasing
    temp=UniIsoRegression::reg_1d(y, weight, metric = 1, decreasing = TRUE)
    print(temp)

    #l_2 metric unimodel
    temp=UniIsoRegression::reg_1d(y, weight, metric = 2, unimodal = TRUE)
    print(temp)

    #l_infinity metric increasing
    temp=UniIsoRegression::reg_1d(y, weight, metric = 3)
    print(temp)
}








