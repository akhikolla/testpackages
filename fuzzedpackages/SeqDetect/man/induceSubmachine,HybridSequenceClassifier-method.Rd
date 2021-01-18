\docType{methods}
\name{induceSubmachine,HybridSequenceClassifier-method}
\alias{induceSubmachine,HybridSequenceClassifier-method}
\title{Sequence Detector Method: ETT projection}
\usage{
\S4method{induceSubmachine}{HybridSequenceClassifier}(threshold,isolate=FALSE)
}
\arguments{
\item{threshold}{(integer) - A threshold for the ETT projection. All transitions that have invocation statistic above the threshold are moved to a submachine.}
\item{isolate}{(logical) - After the regular sequences are moved the the submachine, the original parent can be removed, leaving only the most regular sequences. If TRUE, the parent ETT is removed and only the most regular sequences are left.}
}
\description{
Sequence Detector method for ETT projections. See the SeqDetect vignette for proper usage and cases. All projection changes are performed on the same Sequence Detector object.
}
\value{
Returns:\cr
TRUE - projection was performed successfully\cr
FALSE - no projection was performed.
}
\seealso{
\code{\link{HybridSequenceClassifier}}
}
\examples{
st <- data.frame(product=c("P1","P2"),sales=c(5,76),alert=c(NA,NA))
for(i in 1:400) {
  st <- rbind(st,data.frame(product=c("P1","P2"),sales=c(10,58),alert=c(NA,NA)))
  st <- rbind(st,data.frame(product=c("P1","P2"),sales=c(20,31),alert=c(NA,NA)))
}
st <- rbind(st,data.frame(product=c("P1","P2"),sales=c(30,11),
                          alert=c("Sequence 1","Sequence 2")))
input_streams <- list(stream=st)
pp <- HSC_PP(c("product","sales","alert"),"sequence_id",auto_id=TRUE)
pc <- HSC_PC_Attribute("sales")
hsc <- HybridSequenceClassifier(c("sequence_id","product","sales","alert"),"sequence_id",
                                         "sequence_id",context_field="product",preclassifier=pc,
                                         preprocessor=pp,reuse_states=TRUE,pattern_field="alert")
hsc$process(input_streams,learn=TRUE)
hsc$cleanKeys()
hsc$induceSubmachine(200,isolate=TRUE)
hsc$printMachines()
}
