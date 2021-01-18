
\name{dataPreprocess}

\alias{dataPreprocess}

\title{Preprocess Data to Be Fed into Joint Models}

\description{
\code{dataPreprocess} is a function to preprocess data to be used in fitting joint models. Suppose the situation is that 
the longitudinal measurements are recorded in a data frame with one row per measurment and the survival information 
are recorded in another data frame with one row per subject. This function merges the two data frames by subject identification 
and generate three new columns: \code{start}, \code{stop}, \code{event}. See \bold{Value}.
}

\usage{
dataPreprocess(long, surv, id.col, long.time.col, surv.time.col, surv.event.col, 
               surv.event.indicator = list(censored = 0, event = 1), suffix = ".join")
}

\arguments{
  \item{long}{a data frame for the longitudinal data, one row per measurment, with subject identification, time of measurement, and longitudinal measurements, etc.}
  \item{surv}{a data frame for the survival data, one row per subject, with subject identification (column name should match that in \code{long}), possibly censored time-to-event, and event indicator (normally 0=censored, 1=event), etc.}
  \item{id.col}{a \code{character} string specifying the subject identification column in both \code{long} and \code{surv}.}
  \item{long.time.col}{a \code{character} string specifying the time of measurement column in \code{long}.}
  \item{surv.time.col}{a \code{character} string specifying the possibly censored time-to-event column in \code{surv}.}
  \item{surv.event.col}{a \code{character} string specifying the event status column in \code{surv}.}
  \item{surv.event.indicator}{a \code{list} specifying the values in column \code{surv.event.col} corresponding to censored and event status.}
  \item{suffix}{a optional \code{character} string specifying the suffix to be added to the \code{start}, \code{stop}, \code{event} columns in case  \code{long} or \code{surv} already have columns with these names.}
}

\value{
A data frame merging \code{long} and \code{surv} by subject identification, with one row per longitudinal measurment, 
and generate three new columns: \code{start}, \code{stop}, \code{event} (column names are added with suffix specified by \code{suffix}:
  \describe{
    \item{\code{start}}{starting time of the interval which contains the time of the longitudinal measurements.}
    \item{\code{stop}}{ending time of the interval which contains the time of the longitudinal measurements.}
    \item{\code{event}}{event indicator suggesting whether the event-of-interest, e.g. death, happens in the interval given by \code{start} and \code{stop}.}
  }
}

\author{
Cong Xu \email{helenxu1112@gmail.com}
}

\note{
1. If \code{long} and \code{surv} have columns sharing the same column names, the columns from \code{long} and \code{surv} would be named with suffixes ".long" and ".surv", respectively, in the output data frame.
2. The time of measurement of the longitudinal measurements and possibly censored time-to-event should be recorded consistently for each subject, i.e. time 0 means the same time point for the longitudinal and survival measurements. 
}

\examples{
\dontrun{
liver.join <- dataPreprocess(liver.long, liver.surv, 'ID', 'obstime', 'Time', 'death')
}
}

\keyword{datagen}
