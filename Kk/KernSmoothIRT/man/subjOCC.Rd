\name{subjOCC}
\alias{subjOCC}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ksIRT - kernel smoothing in Item Response Theory
}
\description{
Returns a list containing a matrix for each item. Each matrix in the list contains a row for each option with each column representing a subject with the probability of selecting that option for each subject. 
}
\usage{
subjOCC(x, stype = c("ObsScore","ExpectedScore","MLScore","Theta","MLTheta"))

}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{x}{
a \code{ksIRT} object to be analyzed.
}

\item{stype}{
the scale on which to evaluate each subject. \code{stype = "ObsScore"} uses the subject's observed test score. \code{stype = "ExpectedScore"} uses the subject's expected test score.\code{stype = "MLScore"} uses the maximum likelihood estimate for the subject's overall score.\code{stype = "Theta"} uses the subject's rank on the \code{thetadist} scale. \code{stype = "MLTheta"} uses the maximum likelihood estimate for the subject on the \code{thetadist} scale. 
}

}



\references{
Ramsay, J.O. (2000). 
\pkg{TestGraf}: A program for the graphical analysis of multiple choice test and questionnaire data. 
   
 Silverman, B.W. (1986). 
 \emph{Density Estimation for Statistics and Data Analysis}. 
 Chapman & Hall, London.
}




