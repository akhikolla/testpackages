\name{HSC_PC_None}
\alias{HSC_PC_None}
\title{Straight-through pre-classifier}
\usage{
HSC_PC_None()
}
\description{
Extends the \link{HSC_PC} abstract class.
}
\details{
A pre-classifier class that does not contain any classifier. It passes an input event stream straight through without any modifications.
The only thing is to check whether the input event stream contains \emph{.clazz} field, which should carry classification and input symbols for Sequence Detector ETTs.
}
\examples{
event_stream <- data.frame(product=c("P45","P134","P45","P134","P134","P45","P134"),
                           sales=c(2,12,18,16,18,24,8),
                           alert=c(NA,NA,NA,NA,NA,"Alert P45","Alert P134"),
                           .clazz=c(2,12,18,16,18,24,8))
pc <- HSC_PC_None()
cons_stream <- classify(pc,event_stream)
}
