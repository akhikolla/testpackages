\name{HSC_PC_Binning}
\alias{HSC_PC_Binning}
\title{Binning pre-classifier}
\usage{
HSC_PC_Binning(min_value, max_value, bins, value_field)
}
\arguments{
\item{min_value}{(numeric) - Minimal value.}

\item{max_value}{(numeric) - Maximal value:}

\item{bins}{(integer) - A number of bins that needs to be created.}

\item{value_field}{(character) - The name of the value field in the input event stream.}
}
\description{
Extends the \link{HSC_PC} abstract class.
}
\details{
A pre-classifier takes performs binning on a value field of the input event stream.
}
\examples{
event_stream <- data.frame(product=c("P45","P134","P45","P134","P134","P45","P134"),
                           sales=c(2,12,18,16,18,24,8),
                           alert=c(NA,NA,NA,NA,NA,"Alert P45","Alert P134"))
pc <- HSC_PC_Binning(0,100,40,"sales")
cons_stream <- classify(pc,event_stream)
# Minimal value = 0, Maximal value = 100, 40 bins, values taken from the field named *sales*
}
