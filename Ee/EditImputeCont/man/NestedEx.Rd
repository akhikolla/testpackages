\name{NestedEx}
 \docType{data}
 \alias{NestedEx}
 \title{Data Example With 11 Variables}
 \description{
   This data set gives an example consisting of observed data, ratio edits, range restrictions and nested balance edits.
 }
 \usage{data(NestedEx)}
 \format{A list containing the followings:
 \itemize{
 \item \bold{D.obs}: \code{data.frame} of observed data with \eqn{n=1000} records and \eqn{p=11} variables.
 \item \bold{Ratio.edit}: \code{data.frame} of ratio edits.
 \item \bold{Range.edit}: \code{data.frame} of range restrictions.
 \item \bold{Balance.edit}: \code{data.frame} of balance edits of nested forms, which represent V1=V2+V3+V4, V5=V6+0.4 V10+0.6 V11 and V7=0.4 V10+0.6 V11.
 }
 }
 \seealso{\code{\link{readData}}} 
 \keyword{datasets}
