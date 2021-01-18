\name{Psych101}
\alias{Psych101}
\alias{Psychkey}
\alias{Psychresponses}
\docType{data}
\title{
The Introductory Psychology Data

}
\description{
A multiple-choice examination given to 379 students in an introductory course in psychology in the Christmas exam period of 1989 at McGill. 
The test consisted of 100 multiple-choice items, each having four response options. 

}
\usage{
data(Psych101)

}
\format{
	\describe{
	\item{\code{Psychkey}}{ 
	is a numeric vector with the correct responses for each of the 100 questions on the exam.
	}
	\item{\code{Psychresponses}}{ 
	is a numeric matrix with students' answers. 
	Each column represents a question and each row represents a student. }
		}

}
\details{

See \code{ksIRT()} help page for more information and examples on how to use these data with the \pkg{KernSmoothIRT} package. 
}
\source{
Ramsay, J.O. (2000). 
\pkg{TestGraf}: A program for the graphical analysis of multiple choice test and questionnaire data. 

}

\references{
Ramsay, J.O. and Abrahamowicz M. (1989). 
Binomial regression with monotone splines: A psychometric application. 
\emph{Journal of the American Statistical Association}, 84(408), 906-915.

Ramsay, J.O. (2000). 
\pkg{TestGraf}: A program for the graphical analysis of multiple choice test and questionnaire data. 

}


\keyword{datasets}