\docType{methods}
\name{printMachines,HybridSequenceClassifier-method}
\alias{printMachines,HybridSequenceClassifier-method}
\title{Sequence Detector Method: printout machines (ETTs)}
\usage{
\S4method{printMachines}{HybridSequenceClassifier}(machine_id=NULL,state=NULL,print_cache=TRUE,print_keys=TRUE)
}
\arguments{
\item{machine_id}{(character) - If defined, printout only machine that has the supplied identifier. If NULL, printout all machines.}
\item{state}{(character) - If defined, printout only states that have the supplied identifier. If NULL, printout all states.}
\item{print_cache}{(logical) - Switch for printout of the cache. If FALSE, the cache printout is omitted. The cache can be quite big for each machine and state, and could potentially blur the printout.}
\item{print_keys}{(logical) - Switch for printout of the current token set. If FALSE, the token set printout is omitted. The number of tokens can be considerable, and could potentially blur the printout.}
}
\description{
Sequence Detector method for printing out the machines (ETTs) in the Sequence Detector object. See The SeqDetect vignette for proper usage and cases.
}
\seealso{
\code{\link{HybridSequenceClassifier}}
}
\examples{
ldf1 <- data.frame(product=c("P1","P1","P1","P1"),sequence_id=c(1,3,5,7),
                   sales=c(5,76,123,1),alert=c(NA,NA,NA,"Alert P1"))
ldf2 <- data.frame(product=c("P2","P2","P2","P2"),sequence_id=c(2,4,6,8),
                   sales=c(21,76,123,42),alert=c(NA,NA,NA,"Alert P2"))
input_streams <- list(stream1=ldf1,stream2=ldf2)
pp <- HSC_PP(c("product","sales","alert","sequence_id"),"sequence_id")
pc <- HSC_PC_Attribute("sales")
hsc <- HybridSequenceClassifier(c("sequence_id","product","sales","alert"),
                                         "sequence_id","sequence_id",context_field="product",
                                         preclassifier=pc,preprocessor=pp,reuse_states=TRUE,
                                         pattern_field="alert")
hsc$process(input_streams,learn=TRUE)
hsc$cleanKeys()
hsc$mergeMachines()
hsc$printMachines()
}
