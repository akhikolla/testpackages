\name{subjOCCDIF}
\alias{subjOCCDIF}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{ksIRT - kernel smoothing in Item Response Theory
}
\description{
Applies only when the \code{groups} argument is specified. It returns a list containing a matrix for each item for each of the different groups. Each matrix in the list contains a row for each option with each column representing a subject with the probability of selecting that option for each subject. 
}
\usage{
subjOCCDIF(x, stype = c("ObsScore","Theta","ThetaML","ScoreML"))

}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{x}{
a \code{ksIRT} object to be analyzed.
}

\item{stype}{
the scale on which to evaluate each subject. \code{stype = "ObsScore"} uses the subject's total score. \code{stype = "Theta"} uses the subject's rank on the \code{thetadist} scale. \code{stype = "ThetaML"} uses the maximum likelihood estimate for the subject on the \code{thetadist} scale. code{stype = "ScoreML"} uses the maximum likelihood estimate for the subject on the overall test score scale.
}


}


\references{
Ramsay, J.O. (2000). 
\pkg{TestGraf}: A program for the graphical analysis of multiple choice test and questionnaire data. 
   
 Silverman, B.W. (1986). 
 \emph{Density Estimation for Statistics and Data Analysis}. 
 Chapman & Hall, London.
}




