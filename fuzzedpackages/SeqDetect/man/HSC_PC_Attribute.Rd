\name{HSC_PC_Attribute}
\alias{HSC_PC_Attribute}
\title{Attribute pre-classifier}
\usage{
HSC_PC_Attribute(field)
}
\arguments{
\item{field}{(character) - Field taken as the classification value from the input event stream.}
}
\description{
Extends the \link{HSC_PC} abstract class.
}
\details{
A pre-classifier takes classification from the predefined field in the input event stream and copies these values to the
\emph{.clazz} field. The rest of the input event stream remains unmodified.
}
\examples{
event_stream <- data.frame(product=c("P45","P134","P45","P134","P134","P45","P134"),
                           sales=c(2,12,18,16,18,24,8),
                           alert=c(NA,NA,NA,NA,NA,"Alert P45","Alert P134"))
pc <- HSC_PC_Attribute("sales")
cons_stream <- classify(pc,event_stream)
}
