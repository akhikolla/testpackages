\name{BDI}
\alias{BDI}
\alias{BDIkey}
\alias{BDIresponses}
\docType{data}
\title{
The Beck Depression Inventory Data
}
\description{
The Beck Depression Inventory (Beck, Rush, Shaw & Emery, 1979) is a 21-question self-report inventory, one of the most widely used instruments for measuring the severity of depression. 
Each question has a set of four possible answer choices, ranging in intensity. 
In this dataset, the respondents are 242 McGill University students. 
The data were collected by Prof. David Zuroff.
}
\usage{
data(BDI)

}

\format{
	\describe{
	\item{\code{BDIkey}}{a numeric vector with the highest level ordinal response for each of the 21 questions.}
	\item{\code{BDIresponses}}{is a numeric matrix with examinees responses. Each column represents a question and each row represents a subject; thus, \code{BDIresponses} has 242 rows and 21 columns.}
}
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~

Ramsay JO (2000). 
\pkg{TestGraf}: A program for the graphical analysis of multiple choice test and questionnaire data. 
}



\references{
Beck, A. T., Rush, A., Shaw, B. and Emery, G. (1979).
\emph{Cognitive Therapy of Depression}, 
New York: Guilford Press.

Santor, D. A., Ramsay, J. O. and Zuroff, D. C. (1994). 
Nonparametric item analyses of the Beck Depression Inventory: Evaluating gender item bias and response option weights. 
\emph{Psychological Assessment}, \bold{6}, 255-270.

}

\keyword{datasets}
