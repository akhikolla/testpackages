\name{BigQuic}
\alias{BigQuic}
\docType{package}
\title{
\packageTitle{BigQuic}
}
\description{
\packageDescription{BigQuic}
}

\usage{
BigQuic (X = NULL, inputFileName = NULL, outputFileName = NULL, lambda = 0.5, 
         numthreads = 4, maxit = 5, epsilon = 1e-3, k = 0, memory_size = 8000, 
         verbose = 0, isnormalized = 1, seed = NULL, use_ram = FALSE)
}
\arguments{
  \item{X}{
    An n rows by p columns matrix of the data without the response vector (Y).  
  }
  \item{inputFileName}{
    Full path to a file containing X in the format used for input to BigQUIC, 
    p n
    X
    An example is installed, you can get its path with: 
    paste(path.package("BigQuic"), "/extdata/testInput", sep = "")
  }
  \item{outputFileName}{
    Location and name of output file that will be extrapolated for their 
    naming, e.g. /home/username/test when 3 files are being output will 
    result in /home/username/test.1.output /home/username/test.2.output 
    and /home/username/test.3.output
  }
  \item{lambda}{
    The tuning parameter 0 <= lambda <= 1, but small values should not be used 
    for performance reasons, e.g. < .4 or so.  A vector of lambdas may also be 
    input, in which case BigQUIC will be run for each lambda.  Yes, the 
    examples shows lambda as small as 0.1, but that is only because the 
    testInput matrix is very small so the small lambdas can still finish in a 
    sensible amount of time.  
  }
  \item{numthreads}{
    Number of threads to use for this computation.  
  }
  \item{maxit}{
    Maximum number of Newton iterations.
  }
  \item{epsilon}{
    Convergence tolerance.
  }
  \item{k}{
    Number of memory blocks to use, ideally should be the smallest k such that 
    p/k columns fit in the memory_size.  
  }
  \item{memory_size}{
    The amount of memory this computation is constrained to.  
  }
  \item{verbose}{
    Controls how verbose messages should be printed during execution.
    Valid value range: 0--4. Higher numbers will give more messages for 
    debugging.  
  }
  \item{isnormalized}{
    Whether or not the input is already normalized.  
  }
  \item{seed}{
    A seed for the random number generation, useful for replicating results.
  }
  \item{use_ram}{
    By default the results are written into files, using this option will 
    load those files back to R and return them instead of their paths 
    (the default behavior).  When doing this there is a possibility that 
    R will crash if you don't have enough RAM, use with caution on larger 
    data sets or with many lambdas.  
  }
}
\details{
BigQUIC is finally here!  The original authors of QUIC and BigQUIC brought 
QUIC to Matlab (MEX), Standalone (C++), and R, but BigQUIC was delivered for 
Matlab and Standalone only with no R package.  There are also some other 
features to the package, including sample data generation, inverse selection, 
and plotting.  
IMPORTANT: Due to the practicalities of formatting and working with large data 
sets, files are written to disk at various times when using BigQuic.  The 
locations of the files BigQuic wrote to disk are kept in the object returned 
by BigQuic.  They can be deleted when you're finished with the BigQuic_object 
manually by using the cleanFiles() function as shown in the examples.
There are basiclly 8 cases for file creation, the following will give you an 
idea of where they are in case R crashes completely and loses the references 
to the files so you need to delete them manually.  Files created in tmp are 
deleted on reboot, so no worries if you're having trouble finding them.  

1. X, output file, use_ram = TRUE
length(lambda) output files created in output location
1 file created for X in tmp
Note: this is the same as 5, use_ram doesn't matter in this case

2. input file, no output file, use_ram = FALSE
length(lambda) output files in location of input file

3. input file, output file, use_ram = FALSE
length(lambda)a output files in location of output file
Same as 8, use_ram doesn't matter in this case

4. X, no output file, use_ram = FALSE
length(lambda) output files in tmp
1 file created for X in tmp
Also same as 1 and 5

5. X, output file, use_ram = FALSE
length(lambda) output files created in output file location
1 file created for X in tmp

6. X, no output file, use_ram = TRUE
1 file created for X in tmp

7. input file, no output file, use_ram = TRUE
no files created

8. input file, output file, use_ram = TRUE
length(lambda) output files created in output file location
}

\value{
An object with Reference Class \code{"BigQuic_object"}
  \item{X}{The X input for BigQuic, if given}
  \item{inputFileName}{The file name input for BigQuic, if given}
  \item{isnormalized}{Whether or not the input data was previously normalized}
  \item{k}{k used in BigQuic}
  \item{epsilon}{The epsilon that was used in this run of BigQuic}
  \item{lambda}{lambda used in BigQuic}
  \item{maxit}{maxit used in BigQuic}
  \item{memory_size}{memory_size used in BigQuic}
  \item{numthreads}{numthreads used in BigQuic}
  \item{seed}{seed used in BigQuic}
  \item{use_ram}{use_ram used in BigQuic}
  \item{verbose}{level of verbosity used in BigQuic}


  \item{opt.lambda}{The selected optimal lambda value, initially empty, it will be filled
                    in by running BigQuic.select on the object, see the use in the Examples below}
  \item{precision_matrices}{The precision matrix for each of the lambdas in a list, so to access
                            the one for the 1st lambda in the example: 
                            \code{exampleResult$precision_matrices[[1]]}}


  \item{output_file_names}{Lists files created by the class}
  \item{clean}{Indicates whether or not cleanFiles() has been called on this object before}
  \item{inFlag}{An internal indicator for the class}
  \item{outFlag}{An internal indicator for the class}
  \item{getClass}{Returns Class method definition}
  \item{cleanFiles}{Deletes files created by the class, except for those intentionally 
                    output by specifying an output file name}
  \item{setX}{Used internall to set X}
  \item{setOptLambda}{used internally to set opt.lambda}
  \item{setSeed}{used internally to set the seed}
  \item{.self}{returns the object itself again}
  \item{.refClassDef}{Lists fields and methods of the reference class}
}


\author{
\packageAuthor{BigQuic}

Maintainer: \packageMaintainer{BigQuic}
}
\references{
BIG & QUIC: Sparse Inverse Covariance Estimation for a Million Variables (pdf)
C. Hsieh, M. Sustik, I. Dhillon, P. Ravikumar, R. Poldrack.
In Neural Information Processing Systems (NIPS), December 2013. (Oral)
http://www.cs.utexas.edu/~cjhsieh/hugeQUIC.pdf
http://bigdata.ices.utexas.edu/publication/big-quic-sparse-inverse-covariance-estimation-for-a-million-variables-2/

QUIC: Quadratic Approximation for Sparse Inverse Covariance Matrix Estimation (pdf)
C. Hsieh, M. Sustik, I. Dhillon, P. Ravikumar.
Journal of Machine Learning Research (JMLR), October 2014.
http://jmlr.org/papers/volume15/hsieh14a/hsieh14a.pdf
http://bigdata.ices.utexas.edu/publication/quic-quadratic-approximation-for-sparse-inverse-covariance-matrix-estimation-2/

METIS:"A Fast and Highly Quality Multilevel Scheme for Partitioning Irregular Graphs".
George Karypis and Vipin Kumar. 
SIAM Journal on Scientific Computing, Vol. 20, No. 1, pp. 359-392, 1999.
http://glaros.dtc.umn.edu/gkhome/fetch/papers/mlSIAMSC99.pdf

PCG: A Family of Simple Fast Space-Efficient Statistically Good Algorithms for Random Number Generation. 
This paper is currently submitted to ACM Transactions on Mathematical Software, where it is currently under review.
http://www.pcg-random.org/pdf/toms-oneill-pcg-family-v1.02.pdf
}
\keyword{BigQUIC}
\keyword{package}
\keyword{covariance}
\keyword{regularization}
\keyword{graphical model}
\keyword{lasso}
\keyword{Newton}
\keyword{coordinate descent}
\seealso{
\code{\link[QUIC:QUIC]{QUIC}}
}
\examples{
lambda <- 0.91
exampleResult <- BigQuic(inputFileName = paste(path.package("BigQuic"), 
                         "/extdata/testInput", sep = ""), 
                         outputFileName = tempfile(pattern = 
                         "BigQuic_output_matrix", fileext = ".Bmat"),
                         lambda = lambda, numthreads = 1, memory_size = 512, 
                         seed = 1, use_ram = TRUE)
BigQuic.select(exampleResult)
plot(exampleResult)
exampleResult$cleanFiles()
\dontrun{
If you have the hdi package installed: 
library(hdi)
data(riboflavin)
lambda <- seq(from = 0.9, to = 0.99, by = 0.01)
exampleResult <- BigQuic(as.matrix(riboflavin), lambda = lambda, 
                         numthreads = 1, memory_size = 512, seed = 1, 
                         use_ram = TRUE)
BigQuic.select(exampleResult)
plot(exampleResult)
}
}
